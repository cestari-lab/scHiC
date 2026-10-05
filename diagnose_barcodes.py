#!/usr/bin/env python3
# Author: Lissa Cruz-Saavedra
# Date: 24-09-2026
"""
diagnose_barcodes.py  (step 3)

Barcode QC on a sample of raw reads. For each mate it reports:
  - reads with Illumina adapter readthrough (adapter-dimer signature)
  - reads with each barcode round anywhere, next to the rate expected by
    pure chance (a 6 bp BC3 tail matches ~50% of 150 bp random reads)
  - complete junctions (GATC-BC1-BC2-BC3, <= 1 mismatch per piece),
    where they sit in the read, whether the mirror copy follows, and the
    most frequent cells
  - junction completeness: of reads with an exact BC1, how many have BC2
    right after it, and of those how many have BC3

Real barcodes: rates well above chance and many complete junctions.
Rates close to chance: the reads contain little or no barcoded library.

Usage:
    python3 diagnose_barcodes.py --r1 R1.fastq.gz [--r2 R2.fastq.gz] \
        --bcr1 bcr1.txt --bcr2 bcr2.txt --bcr3 bcr3.txt --n-reads 2000000
"""
import argparse
import json
from collections import Counter

from schic_utils import (TRUSEQ_READTHROUGH, JunctionFinder, add_barcode_args,
                         cell_id, chance_rate, iter_fastq, load_barcodes)

POSITION_BINS = [(0, 0), (1, 20), (21, 50), (51, 100), (101, 10 ** 9)]


def pct(a, b):
    return 100 * a / max(b, 1)


def diagnose(path, n_reads, bc, finder):
    rounds = {"BC1": bc["bcr1"], "BC2": bc["bcr2"], "BC3": bc["bc3_tails"]}
    bc2_set, bc3_set = set(bc["bcr2"]), set(bc["bc3_tails"])
    c = Counter()
    positions, mirrors, cells = Counter(), Counter(), Counter()
    total_len = 0

    for _, seq, _ in iter_fastq(path, n_reads):
        c["reads"] += 1
        total_len += len(seq)
        c["adapter"] += TRUSEQ_READTHROUGH in seq
        c["gatc"] += "GATC" in seq
        for name, barcodes in rounds.items():
            if any(b in seq for b in barcodes):
                c[name] += 1
        for b1 in bc["bcr1"]:  # junction completeness, exact matches
            p = seq.find(b1)
            if p != -1:
                c["bc1_exact"] += 1
                if seq[p + 10:p + 24] in bc2_set:
                    c["bc1_bc2"] += 1
                    if seq[p + 24:p + 30] in bc3_set:
                        c["bc1_bc2_bc3"] += 1
                break
        hit = finder.find(seq)
        if hit:
            start, _, i1, i2, i3, mirror = hit
            c["junction"] += 1
            positions[start] += 1
            mirrors[mirror] += 1
            cells[cell_id(i1, i2, i3)] += 1

    n = c["reads"]
    mean_len = total_len / max(n, 1)
    result = {
        "file": path, "reads": n, "mean_length": round(mean_len, 1),
        "adapter_readthrough_pct": pct(c["adapter"], n),
        "contains_GATC_pct": pct(c["gatc"], n),
        "rounds": {name: {"length": len(bs[0]), "observed_pct": pct(c[name], n),
                          "chance_pct": 100 * chance_rate(len(bs[0]), len(bs), mean_len)}
                   for name, bs in rounds.items()},
        "junction_pct": pct(c["junction"], n),
        "junctions": c["junction"],
        "completeness": {"reads_with_exact_BC1": c["bc1_exact"],
                         "then_BC2": c["bc1_bc2"], "then_BC3": c["bc1_bc2_bc3"]},
        "junction_position_bins": {
            f"{lo}-{hi}" if hi < 10 ** 9 else f">{lo - 1}":
                sum(v for p, v in positions.items() if lo <= p <= hi)
            for lo, hi in POSITION_BINS},
        "top_positions": positions.most_common(5),
        "mirror": dict(mirrors),
        "top_cells": cells.most_common(10),
        "distinct_cells": len(cells),
    }
    return result


def report(r):
    n = r["reads"]
    print(f"\n=== {r['file']} ===")
    print(f"Reads sampled: {n:,} (mean length {r['mean_length']} bp)")
    print(f"Adapter readthrough:  {r['adapter_readthrough_pct']:8.3f}%")
    print(f"Contain GATC:         {r['contains_GATC_pct']:8.3f}%")
    print(f"\n{'round':<6}{'len':>4}{'observed':>12}{'chance':>12}{'obs/chance':>12}")
    for name, d in r["rounds"].items():
        ratio = d["observed_pct"] / d["chance_pct"] if d["chance_pct"] else float("nan")
        print(f"{name:<6}{d['length']:>4}{d['observed_pct']:>11.4f}%"
              f"{d['chance_pct']:>11.4f}%{ratio:>12.1f}")
    comp = r["completeness"]
    print(f"\nExact BC1 found: {comp['reads_with_exact_BC1']:,} | "
          f"+BC2 next: {comp['then_BC2']:,} | +BC3 next: {comp['then_BC3']:,}")
    print(f"Complete junctions (<=1 mismatch/piece): {r['junctions']:,} "
          f"({r['junction_pct']:.4f}%) in {r['distinct_cells']} cells")
    if r["junctions"]:
        print("Junction start position:",
              ", ".join(f"{k}: {v}" for k, v in r["junction_position_bins"].items()))
        print("Mirror after junction:  ",
              ", ".join(f"{k}: {v}" for k, v in sorted(r["mirror"].items())))
        print("Top cells:              ",
              ", ".join(f"{k} ({v})" for k, v in r["top_cells"]))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", default=None)
    add_barcode_args(ap)
    ap.add_argument("--n-reads", type=int, default=1_000_000, help="reads sampled per mate")
    ap.add_argument("--json-out", default=None, help="optional machine-readable output")
    args = ap.parse_args()

    bc = load_barcodes(args.bcr1, args.bcr2, args.bcr3)
    finder = JunctionFinder(bc)
    results = [diagnose(p, args.n_reads, bc, finder) for p in (args.r1, args.r2) if p]
    for r in results:
        report(r)
    if args.json_out:
        with open(args.json_out, "w") as fh:
            json.dump(results, fh, indent=2)


if __name__ == "__main__":
    main()
