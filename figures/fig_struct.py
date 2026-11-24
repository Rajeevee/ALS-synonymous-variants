"""
Figure 3 — RNA secondary structure changes (ΔΔG)
Plots RNAfold ΔΔG values for the 344 variants, highlighting the top-scoring
variants and showing the 201-nt window structures for selected examples.

Usage:
    conda activate figures
    python fig_struct.py \
        --ddg  rna_ddg.tsv \
        --scores SuppData1_variant_scores_344.csv \
        --out  Fig3_rna_structure.pdf
"""

import argparse
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

try:
    import RNA   # ViennaRNA — used for in-panel structure drawing if available
    HAS_RNA = True
except ImportError:
    HAS_RNA = False


def load_ddg(path):
    rows = []
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                ddg = float(row["ddg"])
            except (KeyError, ValueError):
                continue
            rows.append({"gene": row.get("gene",""), "pos": row.get("pos",""),
                         "ref": row.get("ref",""), "alt": row.get("alt",""),
                         "ddg": ddg})
    return rows


def load_tier(path):
    tier = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            key = (row.get("gene",""), row.get("pos",""),
                   row.get("ref",""), row.get("alt",""))
            tier[key] = row.get("tier","minimal")
    return tier


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ddg",    required=True)
    ap.add_argument("--scores", required=True)
    ap.add_argument("--out",    required=True)
    args = ap.parse_args()

    rows  = load_ddg(args.ddg)
    tiers = load_tier(args.scores)

    ddg_vals   = np.array([r["ddg"] for r in rows])
    tier_labels = [tiers.get((r["gene"], r["pos"], r["ref"], r["alt"]), "minimal")
                   for r in rows]

    tier_col = {"high": "#b2182b", "medium": "#ef8a62",
                "low":  "#fddbc7", "minimal": "#d9d9d9"}
    colours = [tier_col.get(t, "#aaaaaa") for t in tier_labels]

    fig, ax = plt.subplots(figsize=(7, 4))

    idx = np.arange(len(ddg_vals))
    # Sort by ΔΔG descending
    order = np.argsort(ddg_vals)[::-1]
    ax.bar(np.arange(len(order)), ddg_vals[order],
           color=[colours[i] for i in order],
           width=1.0, edgecolor="none")

    ax.axhline(0, color="#444444", lw=0.8, ls="--")
    ax.set_xlabel("Variant rank (sorted by ΔΔG)", fontsize=9)
    ax.set_ylabel("ΔΔG (kcal mol⁻¹)", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("RNA structure change (ΔΔG) across 344 variants", fontsize=9)

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=tier_col[t], label=t) for t in
                       ["high", "medium", "low", "minimal"]]
    ax.legend(handles=legend_elements, fontsize=7, frameon=False,
              title="Tier", title_fontsize=7)

    plt.tight_layout()
    plt.savefig(args.out, dpi=300, bbox_inches="tight")
    print(f"Written: {args.out}")


if __name__ == "__main__":
    main()
