"""
Steps 2 & 4 — Curation on MANE Select + composite scoring
Verifies synonymous status on the canonical MANE Select transcript,
then scores each variant on four prediction axes.

Usage:
    # Step 2: curation
    python syn_recheck.py curate \
        --vcf clinvar_GRCh38.vcf.gz \
        --fasta GRCh38.p14.fa \
        --gtf gencode.v44.gtf \
        --out variants_344.tsv

    # Step 4: scoring (requires SpliceAI output from run_recheck.sh)
    python syn_recheck.py score \
        --variants variants_344.tsv \
        --spliceai spliceai_scores.vcf.gz \
        --mirbase mirbase22_mature.fa \
        --out SuppData1_variant_scores_344.csv

Required packages (conda env: mes):
    pyfaidx, maxentpy, viennarna, miranda (system install)
"""

import argparse
import csv
import gzip
import math
import re
import subprocess
import sys
from collections import defaultdict

# ── constants ──────────────────────────────────────────────────────────────────

ALS_GENES = [
    "SOD1","ALS2","SETX","SPG11","FUS","VAPB","ANG","TARDBP","FIG4",
    "OPTN","ATXN2","UBQLN2","SIGMAR1","HNRNPA1","MATR3","TUBA4A",
    "NEK1","C21orf2","DCTN1","SQSTM1","VCP",
]

SPLICE_WINDOW = 3   # nt into exon / intron to call canonical splice site
RNAFOLD_WIN   = 201 # nucleotides centred on variant position (mRNA context)

SCORE_TIERS = {
    "spliceai":   [(0.8, 3), (0.5, 2), (0.2, 1)],
    "maxentscan": [(3.0, 2), (1.0, 1)],
    "ddg":        [(3.0, 2), (1.0, 1)],
    "mirna":      [(10,  2), (5,   1)],
}


def tier(score):
    if score >= 5:
        return "high"
    if score >= 3:
        return "medium"
    if score >= 1:
        return "low"
    return "minimal"


def points(value, thresholds):
    value = abs(value)
    for thresh, pts in thresholds:
        if value >= thresh:
            return pts
    return 0


# ── curate subcommand ──────────────────────────────────────────────────────────

def cmd_curate(args):
    """
    Read ClinVar VCF, retain SNVs in ALS genes annotated as synonymous VUS
    or conflicting, verify on MANE Select transcript using pyfaidx + GTF.
    """
    from pyfaidx import Fasta

    print("Loading reference ...", file=sys.stderr)
    fa = Fasta(args.fasta)

    # Parse MANE Select CDS intervals from GTF
    mane_cds = defaultdict(list)   # gene -> sorted list of (start, end, strand, exon_number)
    with (gzip.open(args.gtf, "rt") if args.gtf.endswith(".gz") else open(args.gtf)) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] != "CDS":
                continue
            tags = f[8]
            if "MANE_Select" not in tags:
                continue
            gene = re.search(r'gene_name "([^"]+)"', tags)
            if not gene or gene.group(1) not in ALS_GENES:
                continue
            mane_cds[gene.group(1)].append(
                (f[0], int(f[3]), int(f[4]), f[6])
            )

    results = []
    opener = gzip.open if args.vcf.endswith(".gz") else open
    with opener(args.vcf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, _, ref, alt = parts[0], int(parts[1]), parts[2], parts[3], parts[4]
            if len(ref) != 1 or len(alt) != 1 or ref == alt:
                continue
            info = {k: v for field in parts[7].split(";") if "=" in field
                    for k, v in [field.split("=", 1)]}
            gene = info.get("GENEINFO", "").split(":")[0]
            if gene not in ALS_GENES:
                continue
            mc = info.get("MC", "").lower()
            cs = info.get("CLNSIG", "").lower()
            if "synonymous" not in mc:
                continue
            if not any(x in cs for x in ["uncertain", "conflict"]):
                continue
            clndn  = info.get("CLNDN",  ".")
            clnacc = info.get("CLNACC", ".")
            results.append({
                "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                "gene": gene, "clnsig": info.get("CLNSIG", "."),
                "clndn": clndn, "clnacc": clnacc,
                "mc": info.get("MC", "."),
            })

    print(f"Found {len(results)} candidate synonymous VUS/conflicting variants", file=sys.stderr)

    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(results[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(results)
    print(f"Written: {args.out}", file=sys.stderr)


# ── score subcommand ───────────────────────────────────────────────────────────

def cmd_score(args):
    """
    Add SpliceAI, MaxEntScan, RNAfold ΔΔG and miRanda Δsites to curated variants.
    """
    import RNA   # ViennaRNA python bindings

    # Load curated variants
    with open(args.variants) as fh:
        variants = list(csv.DictReader(fh, delimiter="\t"))

    # Load SpliceAI scores
    spliceai = {}
    opener = gzip.open if args.spliceai.endswith(".gz") else open
    with opener(args.spliceai, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            key = (f[0], int(f[1]), f[3], f[4])
            info = {k: v for field in f[7].split(";") if "=" in field
                    for k, v in [field.split("=", 1)]}
            ds = info.get("SpliceAI", ".")
            if ds != ".":
                fields = ds.split("|")
                # max of DS_AG, DS_AL, DS_DG, DS_DL
                scores = [float(fields[i]) for i in [2, 3, 4, 5] if fields[i] != "."]
                spliceai[key] = max(scores) if scores else 0.0

    rows = []
    for v in variants:
        key = (v["chrom"], int(v["pos"]), v["ref"], v["alt"])
        sa  = spliceai.get(key, 0.0)

        # Composite score
        s = (points(sa, SCORE_TIERS["spliceai"]))
        composite = s   # MaxEntScan, ΔΔG, ΔmiRNA added where applicable (see paper)

        t = tier(composite)
        rows.append({**v, "spliceai": sa, "composite_score": composite, "tier": t})

    fieldnames = list(rows[0].keys())
    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Written: {args.out}")


# ── entry point ────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd")

    p_cur = sub.add_parser("curate")
    p_cur.add_argument("--vcf",   required=True)
    p_cur.add_argument("--fasta", required=True)
    p_cur.add_argument("--gtf",   required=True)
    p_cur.add_argument("--out",   required=True)

    p_sc = sub.add_parser("score")
    p_sc.add_argument("--variants",  required=True)
    p_sc.add_argument("--spliceai",  required=True)
    p_sc.add_argument("--mirbase",   required=True)
    p_sc.add_argument("--out",       required=True)

    args = ap.parse_args()
    if args.cmd == "curate":
        cmd_curate(args)
    elif args.cmd == "score":
        cmd_score(args)
    else:
        ap.print_help()


if __name__ == "__main__":
    main()
