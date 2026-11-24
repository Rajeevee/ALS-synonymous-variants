"""
Figure 6 / Supplementary Figure 4 — Alzheimer's disease brain cohort
Carrier frequency comparison (ALS cohort vs. AD brain, GSE125583).

Usage:
    conda activate figures
    python fig_ad.py \
        --als_carriers   carriers_als_summary.tsv \
        --ad_carriers    ad_carrier_summary.tsv \
        --scores         SuppData1_variant_scores_344.csv \
        --out            Fig6_AD_comparison.pdf

NOTE: carriers_als_summary.tsv contains only variant-level aggregate
carrier frequencies (not individual genotypes) and is NOT excluded from
the repository. ad_carrier_summary.tsv is generated locally by
SYN_04_AD_genotype.R.
"""

import argparse
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
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
    ap.add_argument("--als_carriers", required=True)
    ap.add_argument("--ad_carriers",  required=True)
    ap.add_argument("--scores",       required=True)
    ap.add_argument("--out",          required=True)
    args = ap.parse_args()

    als = {(r["gene"], r["pos"], r["ref"], r["alt"]): float(r.get("carrier_freq", 0))
           for r in load_csv(args.als_carriers)}
    ad  = {(r["gene"], r["pos"], r["ref"], r["alt"]): float(r.get("n_carrier", 0)) /
                                                       max(float(r.get("n_tested", 1)), 1)
           for r in load_csv(args.ad_carriers) if r.get("diagnosis", "") == "AD"}

    scores = {(r["gene"], r["pos"], r["ref"], r["alt"]): r
              for r in load_csv(args.scores, delimiter=",")}

    # Variants covered in both cohorts
    common_keys = set(als) & set(ad)

    als_freq = np.array([als[k] for k in common_keys])
    ad_freq  = np.array([ad[k]  for k in common_keys])
    tiers    = [scores.get(k, {}).get("tier", "minimal") for k in common_keys]
    colours  = [TIER_COLOURS.get(t, "#aaaaaa") for t in tiers]

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    fig.subplots_adjust(wspace=0.4)

    # Panel a — scatter: ALS carrier freq vs AD carrier freq
    ax = axes[0]
    ax.scatter(als_freq, ad_freq, c=colours, s=18, alpha=0.75, linewidths=0.3,
               edgecolors="#888888")
    lim_max = max(als_freq.max(), ad_freq.max()) * 1.05
    ax.plot([0, lim_max], [0, lim_max], ls="--", lw=0.8, color="#444444",
            label="Equal frequency")
    ax.set_xlabel("Carrier frequency — ALS cohort", fontsize=9)
    ax.set_ylabel("Carrier frequency — Alzheimer's cohort", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("a", loc="left", fontweight="bold", fontsize=10)

    from matplotlib.patches import Patch
    import matplotlib.lines as mlines
    legend_els = [
        mlines.Line2D([], [], ls="--", color="#444444", label="Equal frequency"),
    ] + [Patch(facecolor=TIER_COLOURS[t], label=t) for t in
         ["high", "medium", "low", "minimal"]]
    ax.legend(handles=legend_els, fontsize=7, frameon=False)

    # Panel b — per-diagnosis breakdown for carried variants (≥2 carriers)
    ax = axes[1]
    ad_long = load_csv(args.ad_carriers)
    carried_keys = {(r["gene"], r["pos"], r["ref"], r["alt"])
                    for r in ad_long if int(r.get("n_carrier", 0)) >= 2}

    if carried_keys:
        # For each carried variant, build an x-tick label and stack bars by diagnosis
        diagnoses = sorted({r["diagnosis"] for r in ad_long})
        diag_col  = dict(zip(diagnoses, plt.cm.Set2.colors[:len(diagnoses)]))
        x = np.arange(len(carried_keys))
        carried_list = sorted(carried_keys)

        bottom = np.zeros(len(carried_list))
        for diag in diagnoses:
            counts = []
            for k in carried_list:
                hit = next((r for r in ad_long
                            if (r["gene"], r["pos"], r["ref"], r["alt"]) == k
                            and r["diagnosis"] == diag), None)
                counts.append(int(hit["n_carrier"]) if hit else 0)
            ax.bar(x, counts, bottom=bottom, width=0.6,
                   label=diag, color=diag_col[diag], edgecolor="white", linewidth=0.3)
            bottom += np.array(counts)

        ax.set_xticks(x)
        ax.set_xticklabels([f"{k[0]}\n{k[2]}>{k[3]}" for k in carried_list],
                           fontsize=6, rotation=45, ha="right")
        ax.set_ylabel("Carriers", fontsize=9)
        ax.legend(fontsize=7, frameon=False)
    else:
        ax.text(0.5, 0.5, "No variants carried by ≥2 people",
                ha="center", va="center", transform=ax.transAxes, fontsize=8)

    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("b", loc="left", fontweight="bold", fontsize=10)

    plt.savefig(args.out, dpi=300, bbox_inches="tight")
    print(f"Written: {args.out}")


if __name__ == "__main__":
    main()
