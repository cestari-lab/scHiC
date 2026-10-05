#!/usr/bin/env python3
# Author: Lissa Cruz-Saavedra
# Date: 24-09-2026
"""
Builds a small synthetic dataset with known answers (50 read pairs per case):

  case1  R1: junction at position 0 + forward mirror + arm B    R2: random
  case2  R1: arm A + junction + reverse-complement mirror + arm B
  case3  R1: random    R2: arm A + junction + arm B (rescue of R1 as 'whole')
  case4  random / random                           (dropped)
  case5  R1: arm A + junction with 1 mismatch in BC2 + arm B (still matched)
  case6  adapter dimer / random                    (dropped)
  case7  junction in both mates, same cell, forward mirror

Usage: python3 tests/make_test_data.py --out tests/data
"""
import argparse
import gzip
import json
import os
import random

random.seed(11)
COMP = str.maketrans("ACGT", "TGCA")
rc = lambda s: s.translate(COMP)[::-1]
rnd = lambda n: "".join(random.choice("ACGT") for _ in range(n))
READTHROUGH = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
N = 50


def unique(make, n, min_dist=3):
    """n distinct barcodes, each >= min_dist substitutions from the others,
    so that 1-mismatch correction is never ambiguous."""
    out = []
    while len(out) < n:
        b = make()
        if all(sum(x != y for x, y in zip(b, o)) >= min_dist for o in out):
            out.append(b)
    return sorted(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="tests/data")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    # barcodes without internal GATC so the tests are unambiguous
    ok = lambda s: "GATC" not in s[1:]
    bcr1 = unique(lambda: "GATC" + rnd(6), 20)
    bcr2 = unique(lambda: next(s for s in iter(lambda: "TTGG" + rnd(6) + "AGAG", None) if ok(s)), 20)
    bcr3 = unique(lambda: rnd(6) + "CTCT", 20)
    for name, bs in (("bcr1", bcr1), ("bcr2", bcr2), ("bcr3", bcr3)):
        with open(f"{a.out}/{name}.txt", "w") as fh:
            fh.write("\n".join(bs) + "\n")

    def module(i1, i2, i3, bc2=None):
        return bcr1[i1] + (bc2 or bcr2[i2]) + rc(bcr3[i3])[4:]

    fwd = lambda i1, i2: bcr2[i2] + bcr1[i1]
    rev = lambda i1, i2: rc(bcr2[i2]) + rc(bcr1[i1])
    expected = {}

    def cases(k):
        i1, i2, i3 = random.randrange(20), random.randrange(20), random.randrange(20)
        if k == 1:
            arm_b = rnd(60)
            expected[f"case1_{i}"] = arm_b
            return module(i1, i2, i3) + fwd(i1, i2) + arm_b, rnd(120)
        if k == 2:
            return rnd(40) + module(i1, i2, i3) + rev(i1, i2) + rnd(50), rnd(120)
        if k == 3:
            return rnd(150), rnd(30) + module(i1, i2, i3) + rnd(60)
        if k == 4:
            return rnd(150), rnd(150)
        if k == 5:
            b2 = list(bcr2[i2])
            b2[6] = "A" if b2[6] != "A" else "C"
            return rnd(25) + module(i1, i2, i3, "".join(b2)) + rnd(60), rnd(120)
        if k == 6:
            return rnd(10) + READTHROUGH + rnd(40), rnd(84)
        return (rnd(30) + module(i1, i2, i3) + fwd(i1, i2) + rnd(40),
                rnd(25) + module(i1, i2, i3) + fwd(i1, i2) + rnd(45))

    with gzip.open(f"{a.out}/test_R1.fastq.gz", "wt") as r1, \
         gzip.open(f"{a.out}/test_R2.fastq.gz", "wt") as r2:
        for k in range(1, 8):
            for i in range(N):
                s1, s2 = cases(k)
                r1.write(f"@case{k}_{i} 1:N:0:ACGT\n{s1}\n+\n{'I' * len(s1)}\n")
                r2.write(f"@case{k}_{i} 2:N:0:ACGT\n{s2}\n+\n{'I' * len(s2)}\n")

    with open(f"{a.out}/expected_armB.json", "w") as fh:
        json.dump(expected, fh)
    print(f"Wrote {7 * N} read pairs to {a.out}")


if __name__ == "__main__":
    main()
