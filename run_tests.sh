#!/bin/bash
# Runs every step on synthetic data and checks the answers.
# Usage (from the repository root): bash tests/run_tests.sh
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
S="${REPO}/scripts"; D="${REPO}/tests/data"; O="${REPO}/tests/out"
rm -rf "$O"; mkdir -p "$O"
python3 "${REPO}/tests/make_test_data.py" --out "$D"
BC="--bcr1 $D/bcr1.txt --bcr2 $D/bcr2.txt --bcr3 $D/bcr3.txt"

echo "== 1 build_whitelist"
python3 $S/build_whitelist.py $BC --whitelist-out $O/wl.txt --translate-out $O/tr.tsv
test "$(wc -l < $O/wl.txt)" -eq 8000

echo "== 2 extract_schic_reads"
python3 $S/extract_schic_reads.py --r1 $D/test_R1.fastq.gz --r2 $D/test_R2.fastq.gz $BC \
    --out-fastq $O/segments.fastq.gz --stats-out $O/extract_stats.tsv 2> $O/extract.log

echo "== 3 diagnose_barcodes"
python3 $S/diagnose_barcodes.py --r1 $D/test_R1.fastq.gz --r2 $D/test_R2.fastq.gz $BC \
    --n-reads 1000 --json-out $O/diagnose.json > $O/diagnose.txt

echo "== 4 summarize_segments"
python3 $S/summarize_segments.py --segments $O/segments.fastq.gz \
    --per-cell-out $O/per_cell.tsv > $O/summary.txt

echo "== 5 segments_to_pairs (fake alignments; bwa itself is not tested)"
python3 - "$O" << 'PY'
import gzip, random, sys
o = sys.argv[1]; random.seed(3)
with gzip.open(f"{o}/segments.fastq.gz", "rt") as fh, open(f"{o}/fake.sam", "w") as out:
    out.write("@SQ\tSN:chr1\tLN:1000000\n@SQ\tSN:chr2\tLN:500000\n")
    lines = fh.read().splitlines()
    for i in range(0, len(lines), 4):
        name, tag = lines[i][1:].split(" ")
        seq = lines[i + 1]
        chrom, flag = random.choice(["chr1", "chr2"]), random.choice([0, 16])
        out.write(f"{name}\t{flag}\t{chrom}\t{random.randint(1, 400000)}\t60\t{len(seq)}M\t*\t0\t0\t{seq}\t*\t{tag}\n")
PY
python3 $S/segments_to_pairs.py --sam $O/fake.sam --out-pairs $O/test.pairs.gz \
    --per-cell-dir $O/per_cell_pairs --stats-out $O/contacts.tsv 2> $O/pairs.log

echo "== checks"
python3 - "$O" "$D" << 'PY'
import csv, gzip, json, sys
from collections import Counter
from itertools import combinations
o, d = sys.argv[1], sys.argv[2]
st = {r["metric"]: int(r["count"]) for r in csv.DictReader(open(f"{o}/extract_stats.tsv"), delimiter="\t")}
exp = {"pairs_total": 350, "R1_junction_found": 200, "R2_junction_found": 100,
       "pairs_no_barcode_either_mate": 100, "pairs_junction_both_mates": 50,
       "R1_mirror_forward": 100, "R1_mirror_revcomp": 50, "R2_mirror_forward": 50,
       "whole_emitted": 200}
bad = {k: (st.get(k, 0), v) for k, v in exp.items() if st.get(k, 0) != v}
assert not bad, f"extract stats (got, expected): {bad}"
assert "pairs_cell_disagree_between_mates" not in st

seqs = {}
with gzip.open(f"{o}/segments.fastq.gz", "rt") as fh:
    lines = fh.read().splitlines()
for i in range(0, len(lines), 4):
    seqs[lines[i][1:].split()[0]] = lines[i + 1]
for rid, arm_b in json.load(open(f"{d}/expected_armB.json")).items():
    assert seqs[f"{rid}_R1_armB"] == arm_b, f"mirror not removed correctly in {rid}"

summ = open(f"{o}/summary.txt").read()
for s in ("junction in R1 only:        150", "junction in R2 only:        50",
          "junction in both mates:     50", "same cell in both mates:  50"):
    assert s in summ, f"summary missing '{s}'"

groups = Counter(k.rsplit("_", 2)[0] for k in seqs)
expected_contacts = sum(len(list(combinations(range(n), 2))) for n in groups.values())
with gzip.open(f"{o}/test.pairs.gz", "rt") as fh:
    body = [l for l in fh if not l.startswith("#")]
assert len(body) == expected_contacts, (len(body), expected_contacts)
diag = json.load(open(f"{o}/diagnose.json"))
assert diag[0]["junctions"] == 200 and diag[1]["junctions"] == 100
print(f"ALL TESTS PASSED ({len(body)} contacts)")
PY
