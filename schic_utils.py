# Author: Lissa Cruz-Saavedra
# Date: 24-09-2026
"""
schic_utils.py - shared helpers for the scHiC scripts.

Barcode linker (palindromic), found ANYWHERE in a read:

  [arm A] GATC+BC1tail  BC2  BC3tail | BC2  BC1(+GATC) [arm B]
          '--- bcr1 ---''bcr2''rc(bcr3)[4:]' '-- mirror copy --'

  BC1  bcr1.txt, 10 bp, starts with GATC (Sau3AI overhang)
  BC2  bcr2.txt, 14 bp, ends with AGAG
  BC3  bcr3.txt, 10 bp, ends with CTCT; anneals to AGAG and is filled in,
       so only revcomp(bcr3)[4:] (6 bp) appears in the read.

The forward module (GATC + BC1 tail + BC2 + BC3 tail = 30 bp) identifies
the cell. Each piece is matched with <= 1 mismatch.
"""
import gzip
import sys

ANCHOR = "GATC"
TRUSEQ_READTHROUGH = "AGATCGGAAGAGCACACGTCTGAACTCC"

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(_COMP)[::-1]


def cell_id(i1, i2, i3):
    """0-based round indices -> cell_BC1_BC2_BC3 (1-based, as in the whitelist)."""
    return f"cell_{i1 + 1:02d}_{i2 + 1:02d}_{i3 + 1:02d}"


# ----------------------------------------------------------------------------
# Barcode files
# ----------------------------------------------------------------------------
def load_round(path):
    with open(path) as fh:
        barcodes = [line.strip().upper() for line in fh if line.strip()]
    if len({len(b) for b in barcodes}) != 1:
        sys.exit(f"ERROR: {path} has barcodes of different lengths")
    if len(barcodes) != len(set(barcodes)):
        sys.exit(f"ERROR: {path} contains duplicate barcodes")
    return barcodes


def load_barcodes(bcr1_path, bcr2_path, bcr3_path):
    bcr1, bcr2, bcr3 = load_round(bcr1_path), load_round(bcr2_path), load_round(bcr3_path)
    if not all(b.startswith(ANCHOR) for b in bcr1):
        sys.exit("ERROR: every BC1 must start with GATC")
    if not all(b.endswith("AGAG") for b in bcr2):
        sys.exit("ERROR: every BC2 must end with AGAG")
    if not all(b.endswith("CTCT") for b in bcr3):
        sys.exit("ERROR: every BC3 must end with CTCT")
    bc3_tails = [revcomp(b)[4:] for b in bcr3]
    if len(set(bc3_tails)) != len(bc3_tails):
        sys.exit("ERROR: two BC3 barcodes give the same 6 bp tail")
    return {"bcr1": bcr1, "bcr2": bcr2, "bcr3": bcr3, "bc3_tails": bc3_tails}


def add_barcode_args(parser):
    parser.add_argument("--bcr1", required=True, help="round-1 barcodes (10 bp, GATC + 6)")
    parser.add_argument("--bcr2", required=True, help="round-2 barcodes (14 bp)")
    parser.add_argument("--bcr3", required=True, help="round-3 barcodes (10 bp, end CTCT)")


# ----------------------------------------------------------------------------
# FASTQ
# ----------------------------------------------------------------------------
def open_text(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode, compresslevel=4) if "w" in mode else gzip.open(path, mode)
    return open(path, mode)


def iter_fastq(path, n_reads=None):
    """Yield (read_id, seq, qual). read_id = first header token, '/1' '/2' removed."""
    with open_text(path) as fh:
        n = 0
        while True:
            header = fh.readline()
            if not header:
                return
            seq = fh.readline().rstrip()
            fh.readline()
            qual = fh.readline().rstrip()
            read_id = header[1:].split()[0]
            if read_id.endswith(("/1", "/2")):
                read_id = read_id[:-2]
            yield read_id, seq, qual
            n += 1
            if n_reads and n >= n_reads:
                return


# ----------------------------------------------------------------------------
# Matching
# ----------------------------------------------------------------------------
def one_mismatch_index(barcodes):
    """Map every sequence within 1 substitution of a barcode to its index.
    Exact matches always win; sequences 1 mismatch from 2+ barcodes are left
    out (ambiguous)."""
    index = {b: i for i, b in enumerate(barcodes)}
    near = {}
    for i, b in enumerate(barcodes):
        for p in range(len(b)):
            for c in "ACGTN":
                if c == b[p]:
                    continue
                v = b[:p] + c + b[p + 1:]
                if v in index:
                    continue
                prev = near.get(v)
                near[v] = i if prev is None or prev == i else -1
    index.update({v: i for v, i in near.items() if i >= 0})
    return index


def within_one(query, target):
    if len(query) != len(target):
        return False
    return sum(a != b for a, b in zip(query, target)) <= 1


def chance_rate(k, n_barcodes, read_len):
    """P(>= 1 of n random k-mers occurs in a random read of length L)."""
    positions = max(int(read_len) - k + 1, 0)
    return 1 - (1 - min(n_barcodes / 4 ** k, 1.0)) ** positions


class JunctionFinder:
    """Finds the first barcode junction in a read and the end of its mirror copy."""

    def __init__(self, bc):
        self.bc1_full = bc["bcr1"]
        self.bc2 = bc["bcr2"]
        bc1_tails = [b[4:] for b in bc["bcr1"]]
        self.l1, self.l2, self.l3 = len(bc1_tails[0]), len(self.bc2[0]), len(bc["bc3_tails"][0])
        self.idx1 = one_mismatch_index(bc1_tails)
        self.idx2 = one_mismatch_index(self.bc2)
        self.idx3 = one_mismatch_index(bc["bc3_tails"])
        self.module_len = 4 + self.l1 + self.l2 + self.l3  # 30 bp

    def find(self, seq):
        """Return (start, genomic_b_start, i1, i2, i3, mirror) or None.

        start           first base of GATC (arm A = seq[:start])
        genomic_b_start first base after the linker (arm B = seq[genomic_b_start:])
        mirror          'forward', 'revcomp', '<orientation>_BC2_only' or 'none'
        """
        p = seq.find(ANCHOR)
        while p != -1:
            end = p + self.module_len
            if end > len(seq):
                return None
            a = p + 4
            b = a + self.l1
            c = b + self.l2
            i1 = self.idx1.get(seq[a:b], -1)
            if i1 >= 0:
                i2 = self.idx2.get(seq[b:c], -1)
                if i2 >= 0:
                    i3 = self.idx3.get(seq[c:end], -1)
                    if i3 >= 0:
                        b_start, mirror = self._skip_mirror(seq, end, i1, i2)
                        return p, b_start, i1, i2, i3, mirror
            p = seq.find(ANCHOR, p + 1)
        return None

    def _skip_mirror(self, seq, pos, i1, i2):
        """Skip the mirror copy (BC2 then BC1, same indices) after the module.
        Both orientations are tried: forward (BC2 + BC1) and reverse
        complement (rc(BC2) + rc(BC1)). A read too short to hold the mirror,
        or a mirror that does not match, returns pos unchanged."""
        bc2, bc1 = self.bc2[i2], self.bc1_full[i1]
        for label, first, second in (("forward", bc2, bc1),
                                     ("revcomp", revcomp(bc2), revcomp(bc1))):
            if within_one(seq[pos:pos + len(first)], first):
                cur = pos + len(first)
                if within_one(seq[cur:cur + len(second)], second):
                    return cur + len(second), label
                return cur, f"{label}_BC2_only"
        return pos, "none"
