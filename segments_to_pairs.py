#!/usr/bin/env python3
# Author: Lissa Cruz-Saavedra
# Date: 24-09-2026
"""
segments_to_pairs.py  (step 5b, run by slurm/05_map_segments.sh)

Turns mapped segments into per-cell Hi-C contacts. Segments from the same
original read pair (same read_id, same cell) that map uniquely are paired
with each other:
  junction  armA x armB of the same mate (the two sides of one ligation)
  mate      any other combination (R1 segment x R2 segment)

Input: SAM text from bwa mem -C (the CB:Z: tag is carried from the FASTQ
comment), in the ORDER bwa wrote it (not coordinate-sorted), so that all
segments of a read pair are next to each other. Header (-h) is needed.

Output: 4DN .pairs file (upper triangle, unsorted) with extra columns
cell_id and contact_type, an optional folder of one .pairs file per cell,
and a per-cell contact table.

Usage:
    samtools view -h sample_segments.filtered.bam | \
        python3 segments_to_pairs.py --sam - --out-pairs sample.pairs.gz \
            --per-cell-dir per_cell_pairs --stats-out sample_contacts_per_cell.tsv
"""
import argparse
import os
import re
import sys
from collections import Counter, defaultdict
from itertools import combinations

from schic_utils import open_text

CIGAR_REF = re.compile(r"(\d+)([MDN=X])")
COLUMNS = "readID chrom1 pos1 chrom2 pos2 strand1 strand2 cell_id contact_type"


def ref_length(cigar):
    return sum(int(n) for n, _ in CIGAR_REF.findall(cigar))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sam", required=True, help="SAM file or '-' for stdin")
    ap.add_argument("--out-pairs", required=True)
    ap.add_argument("--per-cell-dir", default=None)
    ap.add_argument("--stats-out", default=None)
    ap.add_argument("--min-mapq", type=int, default=10)
    args = ap.parse_args()

    sam = sys.stdin if args.sam == "-" else open_text(args.sam)
    chroms = []  # (name, length) in header order
    out = open_text(args.out_pairs, "wt")
    stats = Counter()
    per_cell = Counter()
    buffers = defaultdict(list)
    header_written = False

    def header():
        lines = ["## pairs format v1.0", "#sorted: none", "#shape: upper triangle"]
        lines += [f"#chromsize: {c} {l}" for c, l in chroms]
        lines.append(f"#columns: {COLUMNS}")
        return "\n".join(lines) + "\n"

    def flush_cell(cell):
        path = os.path.join(args.per_cell_dir, f"{cell}.pairs")
        new = not os.path.exists(path)
        with open(path, "a") as fh:
            if new:
                fh.write(header())
            fh.writelines(buffers[cell])
        buffers[cell].clear()

    order = {}

    def flush_group(read_id, segs):
        stats["read_pairs_with_mapped_segment"] += 1
        if len(segs) < 2:
            stats["read_pairs_single_segment"] += 1
            return
        stats["read_pairs_with_contact"] += 1
        for a, b in combinations(segs, 2):
            if (order[a[2]], a[3]) > (order[b[2]], b[3]):
                a, b = b, a
            kind = "junction" if a[0] == b[0] and {a[1], b[1]} == {"armA", "armB"} else "mate"
            cell = a[5]
            line = f"{read_id}\t{a[2]}\t{a[3]}\t{b[2]}\t{b[3]}\t{a[4]}\t{b[4]}\t{cell}\t{kind}\n"
            out.write(line)
            stats[f"contacts_{kind}"] += 1
            stats["contacts_cis" if a[2] == b[2] else "contacts_trans"] += 1
            per_cell[cell] += 1
            if args.per_cell_dir:
                buffers[cell].append(line)
                if len(buffers[cell]) >= 10000:
                    flush_cell(cell)

    if args.per_cell_dir:
        os.makedirs(args.per_cell_dir, exist_ok=True)

    current, segs = None, []
    for line in sam:
        if line.startswith("@"):
            if line.startswith("@SQ"):
                tags = dict(t.split(":", 1) for t in line.rstrip("\n").split("\t")[1:])
                chroms.append((tags["SN"], tags["LN"]))
                order[tags["SN"]] = len(order)
            continue
        if not header_written:
            out.write(header())
            header_written = True
        f = line.rstrip("\n").split("\t")
        flag, mapq = int(f[1]), int(f[4])
        stats["alignments_read"] += 1
        if flag & 0x904 or mapq < args.min_mapq:  # unmapped, secondary, supplementary
            stats["alignments_skipped"] += 1
            continue
        read_id, mate, arm = f[0].rsplit("_", 2)
        cell = next((t[5:] for t in f[11:] if t.startswith("CB:Z:")), "unknown")
        strand = "-" if flag & 16 else "+"
        pos = int(f[3]) + (ref_length(f[5]) - 1 if strand == "-" else 0)  # 5' end
        if read_id != current:
            if segs:
                flush_group(current, segs)
            current, segs = read_id, []
        segs.append((mate, arm, f[2], pos, strand, cell))
    if segs:
        flush_group(current, segs)
    if not header_written:
        out.write(header())
    out.close()

    if args.per_cell_dir:
        for cell in list(buffers):
            if buffers[cell]:
                flush_cell(cell)

    if args.stats_out:
        with open(args.stats_out, "w") as fh:
            fh.write("cell_id\tcontacts\n")
            for cell, n in per_cell.most_common():
                fh.write(f"{cell}\t{n}\n")

    for k in sorted(stats):
        print(f"  {k:<34} {stats[k]:>12,}", file=sys.stderr)
    print(f"  {'cells_with_contacts':<34} {len(per_cell):>12,}", file=sys.stderr)


if __name__ == "__main__":
    main()
