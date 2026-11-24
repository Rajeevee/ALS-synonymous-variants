# ALS-synonymous-variants

Code and variant-level data for:

**Sivadasan R. "Silent no more? Synonymous variants in ALS genes resolved by prediction, patient brain RNA and minigene splicing."** (submitted)

## Overview

We catalogued 344 synonymous variants across 21 ALS-associated genes from ClinVar (GRCh38), scored each on four prediction axes (splicing, RNA structure, microRNA binding, canonical splice-site strength), and integrated predictions with patient brain RNA-seq from the NYGC ALS Consortium (468 people, GEO GSE137810) and an Alzheimer's disease cohort (GSE125583, 289 people). Three variants were tested in minigene splicing assays.

## Repository structure

```
curation_scoring/   ClinVar curation, SpliceAI, MaxEntScan, RNAfold, miRanda
rna_structure/      mRNA window extraction and structure visualisation
rna_genotyping/     RNA-based genotyping and allele-specific expression (R)
figures/            All main and supplementary figure scripts (Python)
data/               Variant-level scores and ClinVar counts (no individual data)
envs/               Conda environment specifications
```

## Pipeline

| Step | Script | Input | Output |
|---|---|---|---|
| 1. ClinVar landscape | `curation_scoring/pie_counts.py` | ClinVar GRCh38 VCF | SNV counts by consequence × class (Fig. 1a–c) |
| 2. Curation on MANE Select | `curation_scoring/syn_recheck.py curate` | ClinVar VCF, GRCh38.p14 FASTA, GENCODE v44 GTF | 344 verified synonymous VUS/conflicting variants |
| 3. SpliceAI | `curation_scoring/run_recheck.sh` | step 2 VCF | SpliceAI Δ scores |
| 4. Scoring | `curation_scoring/syn_recheck.py score` | steps 2–3, miRBase 22 | canonical-site MaxEntScan, RNAfold ΔΔG, miRanda Δsites |
| 5. Structure windows | `rna_structure/extract_windows.py` | step 2 output | WT/variant 201-nt mRNA windows + MFE structures |
| 6. RNA genotyping (ALS) | `rna_genotyping/SYN_02_rna_genotype_recheck.R` | NYGC ALS BAMs (GSE137810) | per-sample genotypes, carriers, allele-specific expression |
| 7. Sample identity | `rna_genotyping/SYN_03_sample_identity.R` | same BAMs, common ClinVar SNPs | sample-to-person map (concordance ≥0.97) |
| 8. RNA genotyping (AD) | `rna_genotyping/SYN_04_AD_genotype.R` | GSE125583 BAMs | carriers in Alzheimer's/control brains |
| 9. Figures | `figures/*.py` | outputs above | Figs 1–9, Supp. Figs 2–4 |

`run_recheck.sh` runs steps 2–4 end to end. Set the paths to reference files at the top of each script before running.

## Composite score (0–9)

| Axis | Threshold | Points |
|---|---|---|
| SpliceAI | ≥0.8 / ≥0.5 / ≥0.2 | 3 / 2 / 1 |
| Canonical splice-site ΔMaxEntScan | ≥3 / ≥1 units loss | 2 / 1 |
| ΔΔG (mRNA context) | ≥3 / ≥1 kcal mol⁻¹ | 2 / 1 |
| ΔmiRNA sites | ≥10 / ≥5 | 2 / 1 |

Tiers: **high** ≥5 (9 variants), **medium** 3–4 (32), **low** 1–2 (189), **minimal** 0 (114).

## Data

- `data/SuppData1_variant_scores_344.csv` — all 344 variants with ClinVar accession, MANE cDNA change, SpliceAI Δ score, canonical-site ΔMaxEntScan, ΔΔG (mRNA and pre-mRNA contexts), ΔmiRNA sites, composite score and tier.
- `data/SuppData2_ClinVar_SNV_counts.tsv` — ClinVar SNVs in 21 ALS genes by consequence and clinical significance class.

Individual-level genotypes and clinical data from the NYGC ALS Consortium are **not** included; they are available from the consortium under its data-access agreement.

## Software environments

See `envs/README.md` for conda specifications.

Key versions: Python 3.10, SpliceAI 1.3.1, TensorFlow-CPU 2.11, ViennaRNA 2.7.2, miRanda 3.3a, R 4.x, Rsamtools, GenomicRanges.

## Citation

If you use this code or data, please cite the paper above.

## Licence

Code: MIT. Data files (`data/`): CC BY 4.0.
