"""
RNA structure — Step 3c
Extracts 201-nt mRNA windows centred on each variant, folds them with
RNAfold (ViennaRNA), and reports ΔΔG (alt MFE − ref MFE).

Usage:
    conda activate mes
    python extract_windows.py \
        --variants variants_curated.tsv \
        --fasta    GRCh38.p14.fa \
        --out      rna_ddg.tsv

Required packages: pyfaidx, viennarna (conda env: mes)
"""

import argparse
import csv
import sys

WIN = 100   # ±100 nt around the variant → 201-nt window


def fold_mfe(seq: str) -> float:
    """Return RNAfold minimum free energy for seq."""
    import RNA
    md = RNA.md()
    fc = RNA.fold_compound(seq, md)
    _structure, mfe = fc.mfe()
    return mfe


def rc(seq: str) -> str:
    """Reverse-complement a DNA/RNA sequence (DNA bases)."""
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variants", required=True)
    ap.add_argument("--fasta",    required=True)
    ap.add_argument("--out",      required=True)
    args = ap.parse_args()

    from pyfaidx import Fasta
    fa = Fasta(args.fasta)

    with open(args.variants) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    out_rows = []
    for v in rows:
        chrom = v["chrom"]
        pos   = int(v["pos"])   # 1-based
        ref   = v["ref"]
        alt   = v["alt"]
        gene  = v["gene"]

        start = max(1, pos - WIN)
        end   = pos + WIN
        seq   = str(fa[chrom][start - 1 : end]).upper()

        # Offset of the variant within the window
        offset = pos - start
        if seq[offset] != ref:
            print(f"[WARN] {chrom}:{pos} ref mismatch: expected {ref}, got {seq[offset]}",
                  file=sys.stderr)

        seq_alt = seq[:offset] + alt + seq[offset + 1:]

        # Convert to RNA (T→U)
        rna_ref = seq.replace("T", "U")
        rna_alt = seq_alt.replace("T", "U")

        try:
            mfe_ref = fold_mfe(rna_ref)
            mfe_alt = fold_mfe(rna_alt)
            ddg     = round(mfe_alt - mfe_ref, 2)
        except Exception as exc:
            print(f"[WARN] RNAfold failed for {chrom}:{pos}: {exc}", file=sys.stderr)
            mfe_ref = mfe_alt = ddg = "NA"

        out_rows.append({
            "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
            "gene": gene,
            "mfe_ref": mfe_ref, "mfe_alt": mfe_alt, "ddg": ddg,
        })

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    print(f"Written: {args.out}")


if __name__ == "__main__":
    main()
