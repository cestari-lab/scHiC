#!/usr/bin/env python3
# Author: Lissa Cruz-Saavedra
# Date: 24-09-2026
"""
summarize_segments.py  (step 4)

Reads the segments FASTQ from extract_schic_reads.py and reports, per
original read pair, whether the junction was found in R1 only, R2 only or
both mates (and whether both mates agree on the cell). Optionally writes a
per-cell table of read pairs and segments.

Usage:
    python3 summarize_segments.py --segments sample_segments.fastq.gz \
        --per-cell-out sample_segments_per_cell.tsv
"""
import argparse
from collections import Counter, defaultdict

from schic_utils import open_text


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--segments", required=True)
    ap.add_argument("--per-cell-out", default=None)
    args = ap.parse_args()

    # read_id -> mate -> list of (arm, cell)
    pairs = defaultdict(lambda: defaultdict(list))
    with open_text(args.segments) as fh:
        for i, line in enumerate(fh):
            if i % 4 != 0:
                continue
            name, _, tag = line.rstrip("\n")[1:].partition(" CB:Z:")
            read_id, mate, arm = name.rsplit("_", 2)
            pairs[read_id][mate].append((arm, tag))

    c = Counter()
    per_cell_pairs, per_cell_segments = Counter(), Counter()
    for mates in pairs.values():
        junction = {m: {cell for arm, cell in segs if arm != "whole"}
                    for m, segs in mates.items()}
        r1, r2 = bool(junction.get("R1")), bool(junction.get("R2"))
        if r1 and r2:
            c["both"] += 1
            c["agree" if junction["R1"] == junction["R2"] else "disagree"] += 1
        elif r1:
            c["r1_only"] += 1
        elif r2:
            c["r2_only"] += 1
        else:
            c["segments_without_junction"] += 1  # e.g. junction arms too short
        cell = next(cell for segs in mates.values() for _, cell in segs)
        per_cell_pairs[cell] += 1
        per_cell_segments[cell] += sum(len(s) for s in mates.values())
        for segs in mates.values():
            for arm, _ in segs:
                c[f"seg_{arm}"] += 1

    total = len(pairs)
    print(f"Read pairs with >= 1 segment: {total:,}")
    print(f"  junction in R1 only:        {c['r1_only']:,}")
    print(f"  junction in R2 only:        {c['r2_only']:,}")
    print(f"  junction in both mates:     {c['both']:,}")
    if c["both"]:
        print(f"    same cell in both mates:  {c['agree']:,}")
        print(f"    different cells:          {c['disagree']:,}")
    if c["segments_without_junction"]:
        print(f"  only 'whole' segments kept: {c['segments_without_junction']:,}")
    print(f"Segments: armA {c['seg_armA']:,} | armB {c['seg_armB']:,} | whole {c['seg_whole']:,}")
    print(f"Cells with >= 1 read pair: {len(per_cell_pairs):,}")
    for cell, n in per_cell_pairs.most_common(10):
        print(f"  {cell}  {n:,} pairs")

    if args.per_cell_out:
        with open(args.per_cell_out, "w") as fh:
            fh.write("cell_id\tread_pairs\tsegments\n")
            for cell, n in per_cell_pairs.most_common():
                fh.write(f"{cell}\t{n}\t{per_cell_segments[cell]}\n")


if __name__ == "__main__":
    main()
