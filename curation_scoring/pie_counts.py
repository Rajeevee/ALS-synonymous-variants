"""
Step 1 — ClinVar landscape
Counts ClinVar SNVs in ALS genes by consequence and clinical significance class.
Produces SuppData2_ClinVar_SNV_counts.tsv and the counts behind Fig. 1a–c.

Usage:
    python pie_counts.py --vcf clinvar_GRCh38.vcf.gz --genes gene_list.txt \
                         --gtf gencode.v44.gtf --out SuppData2_ClinVar_SNV_counts.tsv

Required packages (conda env: mes):
    pyfaidx, pysam
"""

import argparse
import csv
import gzip
import re
from collections import defaultdict

ALS_GENES = [
    "SOD1","ALS2","SETX","SPG11","FUS","VAPB","ANG","TARDBP","FIG4",
    "OPTN","ATXN2","UBQLN2","SIGMAR1","HNRNPA1","MATR3","TUBA4A",
    "NEK1","C21orf2","DCTN1","SQSTM1","VCP",
]

CLINSIG_ORDER = [
    "Benign", "Likely benign", "Uncertain significance",
    "Conflicting classifications", "Likely pathogenic", "Pathogenic",
]

CONSEQUENCE_ORDER = [
    "synonymous", "missense", "nonsense", "frameshift",
    "splice_site", "intronic", "other",
]


def clinsig_class(raw):
    r = raw.lower()
    if "pathogenic" in r and "likely" in r:
        return "Likely pathogenic"
    if "pathogenic" in r:
        return "Pathogenic"
    if "likely benign" in r:
        return "Likely benign"
    if "benign" in r:
        return "Benign"
    if "conflict" in r:
        return "Conflicting classifications"
    if "uncertain" in r or "vus" in r:
        return "Uncertain significance"
    return "other"


def consequence_class(mc):
    mc = mc.lower()
    if "synonymous" in mc:
        return "synonymous"
    if "missense" in mc:
        return "missense"
    if "nonsense" in mc or "stop_gained" in mc:
        return "nonsense"
    if "frameshift" in mc:
        return "frameshift"
    if "splice" in mc:
        return "splice_site"
    if "intron" in mc:
        return "intronic"
    return "other"


def parse_info(info_str):
    d = {}
    for field in info_str.split(";"):
        if "=" in field:
            k, v = field.split("=", 1)
            d[k] = v
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    counts = defaultdict(lambda: defaultdict(int))  # counts[consequence][clinsig]

    opener = gzip.open if args.vcf.endswith(".gz") else open
    with opener(args.vcf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            info = parse_info(parts[7])
            gene = info.get("GENEINFO", "").split(":")[0]
            if gene not in ALS_GENES:
                continue
            # SNVs only
            ref, alt = parts[3], parts[4]
            if len(ref) != 1 or len(alt) != 1:
                continue
            mc = info.get("MC", "")
            cs_raw = info.get("CLNSIG", "uncertain_significance")
            cons = consequence_class(mc)
            cs = clinsig_class(cs_raw)
            counts[cons][cs] += 1

    with open(args.out, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["consequence"] + CLINSIG_ORDER)
        for cons in CONSEQUENCE_ORDER:
            row = [counts[cons].get(cs, 0) for cs in CLINSIG_ORDER]
            writer.writerow([cons] + row)

    print(f"Written: {args.out}")


if __name__ == "__main__":
    main()
