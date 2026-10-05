# scHiC

Processing for single-cell Hi-C with a three-round combinatorial barcode
linker: find the barcode junction inside each read, assign the cell, cut the
linker out, map the genomic pieces and build per-cell contacts.

## Read structure

The linker sits **between** the two genomic sides of a ligation, so it can
start anywhere in a read, not only at position 0. It is palindromic
(BC1-BC2-BC3-BC2-BC1):

```
[arm A] GATC+BC1tail  BC2  BC3tail  | BC2  BC1(+GATC) [arm B]
        '--- bcr1 ---''bcr2''rc(bcr3)[4:]' '-- mirror copy --'
```

| Round | File | Length | Rule |
|-------|------|--------|------|
| BC1 | `bcr1.txt` | 10 bp | starts with `GATC` (Sau3AI overhang) |
| BC2 | `bcr2.txt` | 14 bp | ends with `AGAG` |
| BC3 | `bcr3.txt` | 10 bp | ends with `CTCT`; anneals to `AGAG`, is filled in with biotin-dNTPs, so only `revcomp(bcr3)[4:]` (6 bp) is in the read |

The 30 bp forward module (GATC + BC1 tail + BC2 + BC3 tail) identifies the
cell. Each piece is matched with at most 1 mismatch. Cells are named
`cell_<BC1>_<BC2>_<BC3>` using 1-based line numbers in the barcode files.

## Requirements

- Python 3.8+ (standard library only)
- `bwa` and `samtools` for mapping (step 5)
- SLURM for the provided job scripts (optional)

## Setup

```bash
git clone https://github.com/<user>/scHiC.git && cd scHiC
cp config.example.sh config.sh     # edit paths, genome and module loads
mkdir -p logs
```

`config.sh` is git-ignored, so your paths are never committed.

## Workflow

| Step | Run | Output |
|------|-----|--------|
| 1. Whitelist | `python3 scripts/build_whitelist.py` | `chromap_whitelist.txt`, `barcode_translate.tsv` |
| 2. Extract | `sbatch slurm/02_extract.sh` | `<sample>_segments.fastq.gz`, `<sample>_extract_stats.tsv` |
| 3. Diagnose | `sbatch slurm/03_diagnose.sh` | `<sample>_diagnose.txt`, `.json` |
| 4. Summarize | `python3 scripts/summarize_segments.py` | text report, `<sample>_segments_per_cell.tsv` |
| 5. Map + contacts | `sbatch slurm/05_map_segments.sh` | BAMs, `<sample>.pairs.gz`, `per_cell_pairs/`, `<sample>_contacts_per_cell.tsv` |

Steps 1 and 4 take seconds and can run on a login node:

```bash
source config.sh
python3 scripts/build_whitelist.py --bcr1 $BCR1 --bcr2 $BCR2 --bcr3 $BCR3 \
    --whitelist-out $OUT_DIR/chromap_whitelist.txt \
    --translate-out $OUT_DIR/barcode_translate.tsv

python3 scripts/summarize_segments.py --segments $OUT_DIR/${SAMPLE}_segments.fastq.gz \
    --per-cell-out $OUT_DIR/${SAMPLE}_segments_per_cell.tsv
```

Every script has `--help`. To try step 2 on a subset first, add
`--n-reads 1000000` to `extract_schic_reads.py`.

### 1. build_whitelist.py
All BC1 x BC2 x BC3 combinations and a table from each 30 bp barcode to its
cell ID. It warns if two barcodes of a round differ at fewer than 3
positions, because reads with a mismatch between them cannot be assigned.

### 2. extract_schic_reads.py
For each read pair:
- a mate with a junction gives two segments, `armA` (before the linker) and
  `armB` (after the linker and its mirror copy);
- a mate without a junction is kept `whole` if its partner found the cell;
- pairs with no junction in either mate are dropped;
- segments shorter than `--min-arm-len` (default 20 bp) are dropped.

The mirror copy is searched in both orientations (forward BC2+BC1 and
reverse complement); the stats file reports which one was found
(`R1_mirror_forward`, `R1_mirror_revcomp`, `R1_mirror_none`, ...).

Segment headers: `@<read_id>_<R1|R2>_<armA|armB|whole> CB:Z:<cell_id>`

### 3. diagnose_barcodes.py
QC on a sample of raw reads: adapter readthrough, each round's hit rate
against the rate expected by chance, junction completeness (BC1 → BC2 →
BC3), junction positions, mirror orientation and top cells.

### 4. summarize_segments.py
Per read pair: junction found in R1 only, R2 only or both, whether both
mates agree on the cell, and read pairs per cell.

### 5. Mapping and contacts
`05_map_segments.sh` maps segments single-end with `bwa mem -C` (keeps
the `CB:Z:` tag), keeps alignments with MAPQ ≥ `MIN_MAPQ` and drops
unmapped, secondary and supplementary ones. It writes a coordinate-sorted
BAM for IGV and QC, then runs `segments_to_pairs.py` on the unsorted BAM.

`segments_to_pairs.py` pairs every uniquely mapped segment of a read pair
with the others:
- `junction`: armA × armB of the same mate (the two sides of one ligation)
- `mate`: any other combination

Output is a 4DN `.pairs` file (upper triangle, unsorted) with columns
`readID chrom1 pos1 chrom2 pos2 strand1 strand2 cell_id contact_type`,
where positions are 5' ends. `per_cell_pairs/` has one file per cell.
These can go to `pairtools sort`/`dedup` and `cooler cload pairs`.

## Interpreting QC

- **Chance hits.** The expected-by-chance rate is `1 - (1 - n/4^k)^(L-k+1)`
  for `n` barcodes of length `k` in reads of length `L`. For 20 barcodes and
  150 bp reads: about 0.27% for BC1 (10 bp), 0.001% for BC2 (14 bp) and 50%
  for the 6 bp BC3 tail. A BC3 hit alone means nothing; only complete
  junctions assign cells.
- **Adapter dimers.** A high adapter-readthrough rate means many reads have
  no insert.
- **Real signal** looks like complete junctions well above chance,
  concentrated in a limited set of cells, and often confirmed in both mates.

## Testing

```bash
bash tests/run_tests.sh
```

Builds 350 synthetic read pairs with known answers (junction at position 0
and mid-read, both mirror orientations, rescue of a mate without junction,
a 1-mismatch barcode, adapter dimers, junctions in both mates) and checks
every step. Mapping is simulated with fake alignments because bwa is not
needed for the test.

## Layout

```
scHiC/
├── README.md  config.example.sh  environment.yml  .gitignore
├── barcodes/README.md
├── scripts/
│   ├── schic_utils.py            shared barcode logic
│   ├── build_whitelist.py        step 1
│   ├── extract_schic_reads.py    step 2
│   ├── diagnose_barcodes.py      step 3
│   ├── summarize_segments.py     step 4
│   └── segments_to_pairs.py      step 5 (contacts)
├── slurm/
│   ├── 02_extract.sh
│   ├── 03_diagnose.sh
│   └── 05_map_segments.sh
└── tests/
    ├── make_test_data.py
    └── run_tests.sh
```
