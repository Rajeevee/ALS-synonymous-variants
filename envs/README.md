# Software environments

All steps run in conda environments. Create them as follows.

## Environments

| env | used for | key packages |
|---|---|---|
| `mes` | curation, MaxEntScan, RNAfold, miRanda | python 3.10, pyfaidx, maxentpy, viennarna 2.7.2, miranda 3.3a |
| `spliceai` | SpliceAI scoring | python 3.9, spliceai 1.3.1, tensorflow-cpu 2.11, numpy<2 |
| `mirdb` | miRBase 22 mature sequences | R, miRBaseVersions.db, RSQLite |
| `apa` | RNA-based genotyping | R 4.x, Rsamtools, GenomicRanges, data.table, readxl |
| `figures` | plotting | python 3.10, matplotlib, numpy, scipy, viennarna, openpyxl |

## Setup commands

```bash
# Main scoring environment
conda create -n mes -c conda-forge -c bioconda python=3.10 pyfaidx viennarna miranda
conda activate mes
pip install maxentpy

# SpliceAI (requires older numpy)
conda create -n spliceai -c conda-forge python=3.9
conda activate spliceai
pip install spliceai tensorflow-cpu==2.11 "numpy<2" "setuptools<81"

# miRBase (R)
conda create -n mirdb -c conda-forge -c bioconda r-base bioconductor-mirbaseversions.db r-rsqlite
conda activate mirdb

# RNA genotyping (R)
conda create -n apa -c conda-forge -c bioconda r-base bioconductor-rsamtools bioconductor-genomicranges r-data.table r-readxl
conda activate apa

# Figures
conda create -n figures -c conda-forge python=3.10 matplotlib numpy scipy viennarna openpyxl
conda activate figures
```

## Reference files required (not included — download separately)

| File | Source |
|---|---|
| GRCh38.p14 genome FASTA | Ensembl / NCBI |
| GENCODE v44 GTF | gencodegenes.org |
| ClinVar GRCh38 VCF | ftp.ncbi.nlm.nih.gov/pub/clinvar/ |
| miRBase 22 mature.fa | mirbase.org |
| SpliceAI pre-scored VCF (optional) | Illumina / Zenodo |
