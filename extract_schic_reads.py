#!/usr/bin/env python3
# Author: Lissa Cruz-Saavedra
# Date: 24-09-2026
"""
extract_schic_reads.py  (step 2)

Finds the barcode linker anywhere in R1 and R2, assigns the cell, removes
the linker (including its mirror copy) and writes the genomic pieces on
each side as single-end segments for mapping.

  read = [arm A] + GATC-BC1-BC2-BC3 + mirror BC2-BC1 + [arm B]

For each read pair:
  - a mate with a junction gives two segments, armA and armB
  - a mate without a junction is kept whole ("whole") if its partner
    found the cell barcode (rescue)
  - pairs with no junction in either mate are dropped
  - segments shorter than --min-arm-len are dropped

Output header format (read by summarize_segments.py and segments_to_pairs.py):
    @<read_id>_<R1|R2>_<armA|armB|whole> CB:Z:cell_XX_YY_ZZ

Usage:
    python3 extract_schic_reads.py --r1 R1.fastq.gz --r2 R2.fastq.gz \
        --bcr1 bcr1.txt --bcr2 bcr2.txt --bcr3 bcr3.txt \
        --out-fastq sample_segments.fastq.gz --stats-out sample_extract_stats.tsv
"""
import argparse
import sys
from collections import Counter

from schic_utils import (JunctionFinder, add_barcode_args, cell_id, iter_fastq,
                         load_barcodes, open_text)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    add_barcode_args(ap)
    ap.add_argument("--out-fastq", required=True, help="single-end segments (.fastq.gz)")
    ap.add_argument("--stats-out", required=True, help="metric<TAB>count table")
    ap.add_argument("--min-arm-len", type=int, default=20, help="shortest segment kept [20]")
    ap.add_argument("--n-reads", type=int, default=None, help="only the first N pairs (testing)")
    args = ap.parse_args()

    finder = JunctionFinder(load_barcodes(args.bcr1, args.bcr2, args.bcr3))
    stats = Counter()
    out = open_text(args.out_fastq, "wt")

    def emit(read_id, mate, arm, seq, qual, cell):
        if len(seq) < args.min_arm_len:
            stats[f"{arm}_too_short"] += 1
            return
        out.write(f"@{read_id}_{mate}_{arm} CB:Z:{cell}\n{seq}\n+\n{qual}\n")
        stats[f"{arm}_emitted"] += 1

    def scan(mate, seq):
        hit = finder.find(seq)
        if hit is None:
            stats[f"{mate}_no_junction"] += 1
            return None
        start, b_start, i1, i2, i3, mirror = hit
        stats[f"{mate}_junction_found"] += 1
        stats[f"{mate}_mirror_{mirror}"] += 1
        return cell_id(i1, i2, i3), start, b_start

    r2_iter = iter_fastq(args.r2, args.n_reads)
    for n, (id1, s1, q1) in enumerate(iter_fastq(args.r1, args.n_reads), 1):
        try:
            id2, s2, q2 = next(r2_iter)
        except StopIteration:
            sys.exit(f"ERROR: R2 has fewer reads than R1 (stopped at pair {n})")
        if id1 != id2:
            sys.exit(f"ERROR: R1/R2 out of sync at pair {n}: {id1} vs {id2}")
        stats["pairs_total"] += 1

        hit1, hit2 = scan("R1", s1), scan("R2", s2)
        if hit1 and hit2:
            stats["pairs_junction_both_mates"] += 1
            if hit1[0] != hit2[0]:
                stats["pairs_cell_disagree_between_mates"] += 1
        cell = (hit1 or hit2 or (None,))[0]
        if cell is None:
            stats["pairs_no_barcode_either_mate"] += 1
            continue

        for mate, hit, seq, qual in (("R1", hit1, s1, q1), ("R2", hit2, s2, q2)):
            if hit:
                _, start, b_start = hit
                emit(id1, mate, "armA", seq[:start], qual[:start], cell)
                emit(id1, mate, "armB", seq[b_start:], qual[b_start:], cell)
            else:
                emit(id1, mate, "whole", seq, qual, cell)

        if n % 5_000_000 == 0:
            print(f"  ... {n:,} pairs", file=sys.stderr, flush=True)

    if next(r2_iter, None) is not None and args.n_reads is None:
        sys.exit("ERROR: R2 has more reads than R1")
    out.close()

    with open(args.stats_out, "w") as fh:
        fh.write("metric\tcount\n")
        for k in sorted(stats):
            fh.write(f"{k}\t{stats[k]}\n")

    total = max(stats["pairs_total"], 1)
    print("Done.", file=sys.stderr)
    for k in sorted(stats):
        print(f"  {k:<36} {stats[k]:>12,}  ({100 * stats[k] / total:.4f}% of pairs)",
              file=sys.stderr)


if __name__ == "__main__":
    main()
