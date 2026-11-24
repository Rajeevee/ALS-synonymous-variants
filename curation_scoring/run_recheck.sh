#!/usr/bin/env bash
# End-to-end runner for Steps 2–4 of the synonymous-variant pipeline.
# Steps 1 (ClinVar landscape) and 5 (RNA genotyping) are run separately.
#
# Usage:
#   bash run_recheck.sh \
#       --vcf    clinvar_GRCh38.vcf.gz \
#       --fasta  GRCh38.p14.fa \
#       --gtf    gencode.v44.gtf \
#       --mirbase mirbase22_mature.fa \
#       --out    results/
#
# Individual-level NYGC ALS Consortium data are NOT used here.

set -euo pipefail

# ── argument parsing ────────────────────────────────────────────────────────
VCF="" FASTA="" GTF="" MIRBASE="" OUTDIR="results"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --vcf)     VCF="$2";     shift 2 ;;
    --fasta)   FASTA="$2";   shift 2 ;;
    --gtf)     GTF="$2";     shift 2 ;;
    --mirbase) MIRBASE="$2"; shift 2 ;;
    --out)     OUTDIR="$2";  shift 2 ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

for var in VCF FASTA GTF MIRBASE; do
  [[ -z "${!var}" ]] && { echo "--${var,,} is required" >&2; exit 1; }
done

mkdir -p "$OUTDIR"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=== Step 2: Curation on MANE Select transcript ==="
conda run -n mes python "$SCRIPT_DIR/syn_recheck.py" curate \
  --vcf    "$VCF"                              \
  --fasta  "$FASTA"                            \
  --gtf    "$GTF"                              \
  --out    "$OUTDIR/variants_curated.tsv"

echo ""
echo "=== Step 3a: MaxEntScan splice-site scoring ==="
# Extract ±9 nt windows around each SNV position for MaxEntScan
conda run -n mes python - <<'PY' "$OUTDIR/variants_curated.tsv" "$FASTA" "$OUTDIR/mes_scores.tsv"
import sys, csv
from pyfaidx import Fasta
fa = Fasta(sys.argv[2])
with open(sys.argv[1]) as fh, open(sys.argv[3], "w", newline="") as out:
    reader = csv.DictReader(fh, delimiter="\t")
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["chrom","pos","ref","alt","gene","mes5_ref","mes5_alt","mes3_ref","mes3_alt"])
    for v in reader:
        chrom, pos = v["chrom"], int(v["pos"])
        # 9-nt 5' donor window (positions −3 to +6 relative to splice site)
        seq9  = str(fa[chrom][pos-5:pos+4]).upper()
        seq9a = seq9[:pos-1-( pos-5)] + v["alt"] + seq9[pos-( pos-5)+1:]
        # 23-nt 3' acceptor window (positions −20 to +3 relative to splice site)
        seq23  = str(fa[chrom][pos-12:pos+11]).upper()
        seq23a = seq23[:11] + v["alt"] + seq23[12:]
        try:
            from maxentpy import score5, score3
            m5r = score5(seq9);   m5a = score5(seq9a)
            m3r = score3(seq23);  m3a = score3(seq23a)
        except Exception:
            m5r = m5a = m3r = m3a = "NA"
        writer.writerow([v["chrom"], v["pos"], v["ref"], v["alt"], v["gene"],
                         m5r, m5a, m3r, m3a])
PY
echo "MaxEntScan done → $OUTDIR/mes_scores.tsv"

echo ""
echo "=== Step 3b: SpliceAI scoring ==="
# SpliceAI requires its own conda environment and a pre-scored VCF or the
# genome FASTA. The command below scores from FASTA (slow; use the
# pre-scored VCF if available — see envs/README.md).
conda run -n spliceai spliceai \
  -I "$OUTDIR/variants_curated.tsv" \
  -O "$OUTDIR/spliceai_scores.vcf"  \
  -R "$FASTA"                        \
  -A grch38

echo ""
echo "=== Step 3c: RNAfold ΔΔG (201-nt mRNA window) ==="
conda run -n mes python \
  "$(dirname "$SCRIPT_DIR")/rna_structure/extract_windows.py" \
  --variants "$OUTDIR/variants_curated.tsv" \
  --fasta    "$FASTA"                        \
  --out      "$OUTDIR/rna_ddg.tsv"

echo ""
echo "=== Step 3d: miRNA binding (miRBase 22 mature sequences) ==="
conda run -n mes python - <<'PY' "$OUTDIR/variants_curated.tsv" "$MIRBASE" "$OUTDIR/mirna_delta.tsv"
import sys, csv, subprocess, tempfile, os, re

MIRANDA = "miranda"   # must be on PATH in the mes env

def run_miranda(seq_fa, mir_fa):
    """Run miRanda and return the number of predicted binding sites."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as tf:
        tf.write(seq_fa)
        tf_path = tf.name
    result = subprocess.run(
        [MIRANDA, mir_fa, tf_path, "-quiet"],
        capture_output=True, text=True
    )
    os.unlink(tf_path)
    return len(re.findall(r"^>", result.stdout, re.MULTILINE))

with open(sys.argv[1]) as fh, open(sys.argv[3], "w", newline="") as out:
    reader = csv.DictReader(fh, delimiter="\t")
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["chrom","pos","ref","alt","gene","mirna_ref","mirna_alt","delta_mirna"])
    for v in reader:
        ref_fa = f">ref\n{v.get('seq_ref','N'*50)}\n"
        alt_fa = f">alt\n{v.get('seq_alt','N'*50)}\n"
        r = run_miranda(ref_fa, sys.argv[2])
        a = run_miranda(alt_fa, sys.argv[2])
        writer.writerow([v["chrom"], v["pos"], v["ref"], v["alt"], v["gene"], r, a, a - r])
PY
echo "miRNA done → $OUTDIR/mirna_delta.tsv"

echo ""
echo "=== Step 4: Composite scoring ==="
conda run -n mes python "$SCRIPT_DIR/syn_recheck.py" score \
  --variants "$OUTDIR/variants_curated.tsv"     \
  --spliceai "$OUTDIR/spliceai_scores.vcf"      \
  --mirbase  "$MIRBASE"                          \
  --out      "$OUTDIR/SuppData1_variant_scores_344.csv"

echo ""
echo "All done. Final scored table: $OUTDIR/SuppData1_variant_scores_344.csv"
