"""
Figure 2 — Composite score distribution and tier breakdown
Scatter / violin of composite scores across the 344 synonymous VUS/conflicting
variants; colour-coded by tier (high/medium/low/minimal).

Usage:
    conda activate figures
    python fig_comp2.py \
        --scores SuppData1_variant_scores_344.csv \
        --out    Fig2_composite_scores.pdf
"""

import argparse
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

TIER_COLOURS = {
    "high":    "#b2182b",
    "medium":  "#ef8a62",
    "low":     "#fddbc7",
    "minimal": "#d9d9d9",
}

TIER_ORDER = ["high", "medium", "low", "minimal"]


def load_scores(path):
    rows = []
    with open(path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            try:
                score = float(row["composite_score"])
            except (KeyError, ValueError):
                continue
            rows.append({"score": score, "tier": row.get("tier", "minimal"),
                         "gene": row.get("gene", "")})
    return rows


def jitter(n, width=0.3, seed=42):
    rng = np.random.default_rng(seed)
    return rng.uniform(-width, width, n)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", required=True)
    ap.add_argument("--out",    required=True)
    args = ap.parse_args()

    rows = load_scores(args.scores)
    scores = np.array([r["score"] for r in rows])
    tiers  = [r["tier"] for r in rows]
    colours = [TIER_COLOURS.get(t, "#aaaaaa") for t in tiers]

    fig, axes = plt.subplots(1, 2, figsize=(9, 4), gridspec_kw={"width_ratios": [2, 1]})
    fig.subplots_adjust(wspace=0.3)

    # Panel a — strip plot of all 344 scores
    ax = axes[0]
    x = jitter(len(scores))
    ax.scatter(x, scores, c=colours, s=12, alpha=0.75, linewidths=0,
               zorder=2)
    # Tier threshold lines
    for thresh, label in [(5, "high ≥5"), (3, "medium ≥3"), (1, "low ≥1")]:
        ax.axhline(thresh, ls="--", lw=0.8, color="#666666", zorder=1)
        ax.text(0.32, thresh + 0.08, label, fontsize=7, color="#444444", va="bottom")
    ax.set_xlim(-0.5, 0.5)
    ax.set_xticks([])
    ax.set_ylabel("Composite score (0–9)", fontsize=9)
    ax.set_yticks(range(0, 10))
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("a", loc="left", fontweight="bold", fontsize=10)

    # Panel b — bar chart of tier counts
    ax = axes[1]
    tier_counts = {t: tiers.count(t) for t in TIER_ORDER}
    y = np.arange(len(TIER_ORDER))
    ax.barh(y, [tier_counts[t] for t in TIER_ORDER],
            color=[TIER_COLOURS[t] for t in TIER_ORDER],
            height=0.6, edgecolor="#888888", linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(TIER_ORDER, fontsize=9)
    ax.set_xlabel("Number of variants", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("b", loc="left", fontweight="bold", fontsize=10)
    for i, t in enumerate(TIER_ORDER):
        ax.text(tier_counts[t] + 1, i, str(tier_counts[t]), va="center", fontsize=8)

    plt.savefig(args.out, dpi=300, bbox_inches="tight")
    print(f"Written: {args.out}")


if __name__ == "__main__":
    main()
