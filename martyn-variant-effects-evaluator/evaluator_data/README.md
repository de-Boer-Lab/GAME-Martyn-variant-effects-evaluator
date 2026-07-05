# `evaluator_data/` — inputs and sequence generation

This directory holds everything the Martyn Variant Effects Evaluator reads at **runtime**, as well as the raw data + scripts to recreate the input data. At run time it is bind-mounted into the container
at `/evaluator_data`. See the [top-level README](../README.md) for how the container consumes these
files; this document covers **how they are made**.

---

## Contents

```
evaluator_data/
├── parse_Engreitz_data.py         Script to create Evaluator input data (included as a reference)
├── gene_reference_sequences.csv   Reference promoter sequences, one row per gene
├── File Specification_...png       Upstream Variant-EFFECTS column spec (reference)
├── THP1/
│   ├── *_VariantEffectsFile.tsv    Raw per-assay measurement tables (PPIF only)
│   ├── all_THP1.tsv                Merged measurements (Evaluator input)
│   └── variant_sequences.csv       Per-variant sequences (Evaluator input)
└── Jurkat/
    ├── *_VariantEffectsFile.tsv    Raw per-assay tables (PPIF + IL2RA)
    ├── all_Jurkat.tsv              Merged measurements (Evaluator input)
    └── variant_sequences.csv       Per-variant sequences (Evaluator input)
```

---


### `gene_reference_sequences.csv` — comma-separated
| Column | Type | Notes |
|--------|------|-------|
| `gene` | str | `PPIF` or `IL2RA`; used as the sequence key in the request |
| `promoter_length` | int | Width of the prediction range (`PPIF 200`, `IL2RA 184`) |
| `strand` | str | `positive` / `negative` (documentary; not read by the container) |
| `sequence` | str | 2 Mb reference sequence (±1 Mb around the promoter start) |

### `{cell}/variant_sequences.csv` — **tab-separated** 
This is the merged measurement table (`all_{cell}.tsv`) with a `variant_sequence` column appended,
so it carries every measurement column plus a leading `Unnamed: 0` index and **both** coordinate
systems (GRCh38 `variant`/`chr`/`pos` from the source spec, and the added hg19 `VariantID_h19`).
`data_loader.create_json_from_tsv()` reads it with `sep='\t'` but consumes only three columns:
| Column | Notes |
|--------|-------|
| `VariantID_h19` | `chrom:pos:REF>ALT` (hg19), e.g. `chr10:81107165:G>GAACGGAGC` (insertions supported); request key and merge key |
| `gene_symbol` | `PPIF` / `IL2RA`; selects the reference sequence and prediction range |
| `variant_sequence` | ~2 Mb sequence with the variant spliced in (reverse-complemented if negative strand) |

### `{cell}/all_{cell}.tsv` — tab-separated
Read at the metrics step for the **measured** effect sizes. Required columns:
`VariantID_h19`, `gene_symbol`, `log2_fold_change` (the measured log2 fold-change; corresponds to
the upstream spec's `log2fc`).

---

## The Variant-EFFECTS source files

Each `*_VariantEffectsFile.tsv` is one screen (element × cell type). The upstream column spec
(`File Specification_VariantEffectsFiles.png`) defines, per the paper: a SPDI `variant`, `chr`,
`pos` (GRCh38), `ref`, `alt`, `effect_allele`, `other_allele`, `gene` (Ensembl ID), `gene_symbol`,
optional `transcript`, `effect_size` (fractional change, range `[−1, ∞)`), `log2fc`,
`p_nominal_nlog10`, `fdr_nlog10`, `fdr_method`, and `power` `[0,1]`.

**This evaluator keys off `VariantID_h19` (hg19), not the GRCh38 SPDI `variant`.** Using the hg19 ID
keeps variant coordinates consistent with the hg19 promoter coordinates from the paper and `hg19.fa` used to build the sequences.

Screens present:
- **THP1/** — PPIF only: promoter tiling, enhancer tiling, 5′ splice pool, promoter 8-mer insertion,
  promoter Enformer edits.
- **Jurkat/** — PPIF (promoter Enformer edits, 8-mer insertion ± stimulation) **and** IL2RA
  (promoter tiling, stimulation).

---

## `parse_Engreitz_data.py` — data generation steps

Run **once**, outside the container (needs `pysam`, `biopython`, `pandas`, and a local `hg19.fa`).
Absolute paths are hard-coded — edit them for your environment first.

**Step 1 — merge (currently commented out).** `merge_variants_from_tsv()` concatenates the per-assay
`*_VariantEffectsFile.tsv` files in a cell-type folder into `all_{cell}.tsv`. It was run once and
left commented; uncomment to regenerate the merged tables.

**Step 2 — reference sequences.** For each promoter, fetch ±1 Mb around the start coordinate from
`hg19.fa` (a 2 Mb window). Write `gene_reference_sequences.csv` with `PPIF` on the forward strand and
`IL2RA` reverse-complemented (see strand handling below).

Promoter anchors (hg19; coordinates from paper):
- **PPIF** — `chr10:81,107,026`, positive strand, 220 bp target region (81,107,026–81,107,246).
- **IL2RA** — `chr10:6,104,468`, negative strand, 184 bp target region (6,104,468–6,104,652).

**Step 3 — variant sequences.** For each row of `all_{cell}.tsv`:
1. Parse `VariantID_h19` into `chrom, pos, ref, alt`.
2. Select the gene's 2 Mb backbone; locate the variant at index `pos − (promoter_start − 1e6)`.
3. Insert `alt` (replacing `ref`); warn if the backbone base ≠ `ref`.
4. If the gene is on the negative strand, reverse-complement the whole sequence.
5. Append as `variant_sequence` and write `{cell}/variant_sequences.csv`.

---

## Strand handling (intended design)

| Gene | Strand | Backbone for variant insertion | After insertion | Reference stored |
|------|--------|-------------------------------|-----------------|------------------|
| PPIF | + | forward 2 Mb | (unchanged) | forward |
| IL2RA | − | **forward** 2 Mb | reverse-complement whole sequence | reverse-complement |

For IL2RA the variant is inserted using forward-strand coordinates (so the REF-allele check runs
against the forward reference), and only then is the whole sequence reverse-complemented — so the
stored reference and every variant sequence share the same (negative-strand) orientation. Prediction
ranges stay `start → start+length` **without flipping**: because the reference and its variants are
reverse-complemented identically, the same range measures the same location in both, and the
variant-vs-reference log2FC is computed consistently. 

---


## Regenerating the inputs

```bash
# 1. edit the hard-coded paths and hg19.fa location in parse_Engreitz_data.py
# 2. (optional) uncomment merge_variants_from_tsv() calls to rebuild all_{cell}.tsv
conda activate <env-with-pysam-biopython-pandas>
python3 parse_Engreitz_data.py
# -> writes gene_reference_sequences.csv and {THP1,Jurkat}/variant_sequences.csv
```
