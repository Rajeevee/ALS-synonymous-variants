"""
Figure 1 — ClinVar variant landscape in ALS genes
Pie/bar charts of SNV counts by consequence and clinical significance
(panel a: synonymous vs other consequence; panel b: synonymous by clinsig tier).

Usage:
    conda activate figures
    python fig_landscape.py \
        --counts SuppData2_ClinVar_SNV_counts.tsv \
        --out    Fig1_landscape.pdf
"""

import argparse
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

# ── colour palette (accessible, perceptually ordered) ───────────────────────
CLINSIG_COLOURS = {
    "Benign":                   "#4393c3",
    "Likely benign":            "#92c5de",
    "Uncertain significance":   "#f4a582",
    "Conflicting classifications": "#d6604d",
    "Likely pathogenic":        "#b2182b",
    "Pathogenic":               "#67001f",
}

CONS_COLOURS = {
    "synonymous":  "#b2df8a",
    "missense":    "#1f78b4",
    "nonsense":    "#e31a1c",
    "frameshift":  "#ff7f00",
    "splice_site": "#6a3d9a",
    "intronic":    "#a6cee3",
    "other":       "#cab2d6",
}

CLINSIG_ORDER = [
    "Benign", "Likely benign", "Uncertain significance",
    "Conflicting classifications", "Likely pathogenic", "Pathogenic",
]

CONSEQUENCE_ORDER = [
    "synonymous", "missense", "nonsense", "frameshift",
    "splice_site", "intronic", "other",
]


def load_counts(path):
    """Load SuppData2 TSV into {consequence: {clinsig: count}}."""
    data = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            cons = row["consequence"]
            data[cons] = {cs: int(row.get(cs, 0)) for cs in CLINSIG_ORDER}
    return data


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--out",    required=True)
    args = ap.parse_args()

    data = load_counts(args.counts)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    fig.subplots_adjust(wspace=0.4)

    # Panel a — total SNVs by consequence (horizontal bar)
    ax = axes[0]
    totals = {c: sum(data.get(c, {}).values()) for c in CONSEQUENCE_ORDER}
    y = np.arange(len(CONSEQUENCE_ORDER))
    bars = ax.barh(y, [totals[c] for c in CONSEQUENCE_ORDER],
                   color=[CONS_COLOURS[c] for c in CONSEQUENCE_ORDER],
                   height=0.6, edgecolor="white", linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels([c.replace("_", " ").capitalize() for c in CONSEQUENCE_ORDER],
                       fontsize=8)
    ax.set_xlabel("Number of ClinVar SNVs", fontsize=9)
    ax.set_title("a", loc="left", fontweight="bold", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))

    # Panel b — synonymous SNVs by clinical significance (stacked bar)
    ax = axes[1]
    syn = data.get("synonymous", {cs: 0 for cs in CLINSIG_ORDER})
    bottoms = 0
    for cs in CLINSIG_ORDER:
        count = syn.get(cs, 0)
        ax.bar(0, count, bottom=bottoms, width=0.4,
               color=CLINSIG_COLOURS[cs], label=cs,
               edgecolor="white", linewidth=0.5)
        bottoms += count
    ax.set_xlim(-0.5, 0.5)
    ax.set_xticks([])
    ax.set_ylabel("Synonymous SNVs (21 ALS genes)", fontsize=9)
    ax.set_title("b", loc="left", fontweight="bold", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(fontsize=7, frameon=False, loc="upper right",
              bbox_to_anchor=(2.2, 1.0))

    plt.savefig(args.out, dpi=300, bbox_inches="tight")
    print(f"Written: {args.out}")


if __name__ == "__main__":
    main()
