#!/usr/bin/env python3
# Author: Lissa Cruz-Saavedra
# Date: 24-09-2026
"""
build_whitelist.py  (step 1)

Builds every possible cell barcode from the three round files and a table
that links each 30 bp barcode to its cell ID (cell_BC1_BC2_BC3).
The cell IDs match the CB:Z: tags written by extract_schic_reads.py.

  whitelist entry = bcr1 (10 bp) + bcr2 (14 bp) + revcomp(bcr3)[4:] (6 bp)

Usage:
    python3 build_whitelist.py --bcr1 bcr1.txt --bcr2 bcr2.txt --bcr3 bcr3.txt \
        --whitelist-out chromap_whitelist.txt --translate-out barcode_translate.tsv
"""
import argparse
import sys

from schic_utils import add_barcode_args, cell_id, load_barcodes


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    add_barcode_args(ap)
    ap.add_argument("--whitelist-out", default="chromap_whitelist.txt")
    ap.add_argument("--translate-out", default="barcode_translate.tsv")
    args = ap.parse_args()

    bc = load_barcodes(args.bcr1, args.bcr2, args.bcr3)
    n = len(bc["bcr1"]) * len(bc["bcr2"]) * len(bc["bc3_tails"])

    # Error correction allows 1 mismatch per piece. Barcodes fewer than 3
    # substitutions apart can make a read ambiguous (it is then not assigned).
    rounds = {"BC1": bc["bcr1"], "BC2": bc["bcr2"], "BC3 tail": bc["bc3_tails"]}
    for name, barcodes in rounds.items():
        close = [(a, b, d) for i, a in enumerate(barcodes) for b in barcodes[i + 1:]
                 if (d := sum(x != y for x, y in zip(a, b))) < 3]
        for a, b, d in close:
            print(f"WARNING: {name} {a} and {b} differ at only {d} position(s); "
                  f"reads with a mismatch between them stay unassigned", file=sys.stderr)

    with open(args.whitelist_out, "w") as wl, open(args.translate_out, "w") as tr:
        tr.write("barcode\tcell_id\tbc1_index\tbc2_index\tbc3_index\n")
        for i1, b1 in enumerate(bc["bcr1"]):
            for i2, b2 in enumerate(bc["bcr2"]):
                for i3, b3 in enumerate(bc["bc3_tails"]):
                    combo = b1 + b2 + b3
                    wl.write(combo + "\n")
                    tr.write(f"{combo}\t{cell_id(i1, i2, i3)}\t{i1 + 1}\t{i2 + 1}\t{i3 + 1}\n")

    print(f"{len(bc['bcr1'])} x {len(bc['bcr2'])} x {len(bc['bc3_tails'])} = {n} cell barcodes "
          f"-> {args.whitelist_out}, {args.translate_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
