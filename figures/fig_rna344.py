"""
Figure 5 — RNA-based genotyping in the NYGC ALS cohort (344 variants)
Allelic imbalance plot: RNA ALT allele fraction vs. expected fraction from
RNA-seq pileup for carriers of each variant.

NOTE: Individual-level genotype and clinical data from the NYGC ALS
Consortium are NOT included in this repository and are NOT read by this
script from the repository. The input files (rna_genotypes_long.tsv,
manifest_identity_ok.tsv) must be generated locally using the
rna_genotyping/ scripts and are excluded via .gitignore.

Usage:
    conda activate figures
    python fig_rna344.py \
        --pileup  rna_genotypes_long.tsv \
        --manifest manifest_identity_ok.tsv \
        --scores  SuppData1_variant_scores_344.csv \
        --out     Fig5_rna_genotyping.pdf
"""

import argparse
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np

TIER_COLOURS = {
    "high":    "#b2182b",
    "medium":  "#ef8a62",
    "low":     "#fddbc7",
    "minimal": "#d9d9d9",
}


def load_csv(path, delimiter="\t"):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter=delimiter))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pileup",   required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--scores",   required=True)
    ap.add_argument("--out",      required=True)
    args = ap.parse_args()

    pileup   = load_csv(args.pileup)
    manifest = load_csv(args.manifest)
    scores   = {(r["gene"], r["pos"], r["ref"], r["alt"]): r
                for r in load_csv(args.scores, delimiter=",")}

    ok_samples = {r["sample_id"] for r in manifest}

    # Keep only identity-verified carriers
    carriers = [r for r in pileup
                if r["sample_id"] in ok_samples
                and r.get("carried", "FALSE") == "TRUE"]

    if not carriers:
        print("[WARN] No carriers found after filtering — check inputs.")
        return

    afs   = [float(r["af"]) for r in carriers]
    genes = [r["gene"] for r in carriers]
    tiers = [scores.get((r["gene"], r["pos"], r["ref"], r["alt"]),
                        {}).get("tier", "minimal") for r in carriers]
    colours = [TIER_COLOURS.get(t, "#aaaaaa") for t in tiers]

    fig, ax = plt.subplots(figsize=(6, 5))

    # Jittered strip plot: allele fraction per carrier observation
    jitter = np.random.default_rng(42).uniform(-0.3, 0.3, len(afs))
    ax.scatter(jitter, afs, c=colours, s=15, alpha=0.7, linewidths=0)

    ax.axhline(0.5, ls="--", lw=0.8, color="#444444", label="Expected (het)")
    ax.set_xlim(-0.6, 0.6)
    ax.set_xticks([])
    ax.set_ylim(0, 1)
    ax.set_ylabel("RNA ALT allele fraction", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("Allelic imbalance — ALS cohort carriers", fontsize=9)

    from matplotlib.patches import Patch
    legend_elements = [
        mlines.Line2D([], [], ls="--", color="#444444", label="Expected (het, 0.5)"),
    ] + [Patch(facecolor=TIER_COLOURS[t], label=t) for t in
         ["high", "medium", "low", "minimal"]]
    ax.legend(handles=legend_elements, fontsize=7, frameon=False)

    plt.tight_layout()
    plt.savefig(args.out, dpi=300, bbox_inches="tight")
    print(f"Written: {args.out}")


if __name__ == "__main__":
    main()
