# SYN_03_sample_identity.R
# Step 5b — Sample identity check (RNA-seq ↔ DNA WGS concordance)
#
# Compares allele fractions at common SNP positions between RNA-seq pileup
# and WGS genotypes to confirm that each RNA-seq BAM matches its declared
# WGS sample. Flags discordant pairs for exclusion.
#
# NOTE: Individual-level genotypes derived from the NYGC ALS Consortium
# are NOT included in this repository. This script reads and writes
# per-sample data locally; output files are excluded via .gitignore.
#
# Usage:
#   conda activate apa
#   Rscript SYN_03_sample_identity.R \
#       --rna_pileup  rna_genotypes_long.tsv \
#       --wgs_vcf     wgs_snps_common.vcf.gz \
#       --manifest    manifest_bams.tsv \
#       --out_ok      manifest_identity_ok.tsv \
#       --out_fail    manifest_identity_fail.tsv

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(Rsamtools)
  library(GenomicRanges)
})

opt_list <- list(
  make_option("--rna_pileup", type = "character", help = "Long-format RNA pileup TSV"),
  make_option("--wgs_vcf",    type = "character", help = "WGS SNP VCF (gzipped)"),
  make_option("--manifest",   type = "character", help = "BAM manifest TSV"),
  make_option("--out_ok",     type = "character", default = "manifest_identity_ok.tsv"),
  make_option("--out_fail",   type = "character", default = "manifest_identity_fail.tsv"),
  make_option("--r2_thresh",  type = "numeric",   default = 0.9,
              help = "Minimum Pearson r^2 for concordance [default: %default]"),
  make_option("--min_snps",   type = "integer",   default = 50L,
              help = "Minimum overlapping SNPs for a valid check [default: %default]")
)
opt <- parse_args(OptionParser(option_list = opt_list))

manifest  <- fread(opt$manifest,   sep = "\t")
rna_pileup <- fread(opt$rna_pileup, sep = "\t")

# ── parse WGS genotypes from VCF ────────────────────────────────────────────
# Read WGS VCF header to get sample columns
vcf_header_lines <- system(paste("zcat", opt$wgs_vcf, "| grep '^#CHROM'"), intern = TRUE)
vcf_cols <- strsplit(vcf_header_lines, "\t")[[1]]
sample_cols <- vcf_cols[-(1:9)]

vcf_dt <- fread(cmd = paste("zcat", opt$wgs_vcf, "| grep -v '^#'"),
                col.names = vcf_cols, sep = "\t")

# Keep biallelic SNPs only
vcf_dt <- vcf_dt[nchar(REF) == 1 & nchar(ALT) == 1]

# Parse GT fields — assume GT is first FORMAT field
gt_mat <- as.data.table(lapply(sample_cols, function(s) {
  gt <- sapply(strsplit(vcf_dt[[s]], ":"), `[`, 1)
  # Convert 0/0→0, 0/1→0.5, 1/1→1
  dplyr::case_when(
    gt %in% c("0/0","0|0") ~ 0,
    gt %in% c("0/1","1/0","0|1","1|0") ~ 0.5,
    gt %in% c("1/1","1|1") ~ 1,
    TRUE ~ NA_real_
  )
}))
setnames(gt_mat, sample_cols)
vcf_dt <- cbind(vcf_dt[, .(CHROM, POS, REF, ALT)], gt_mat)

# ── compute concordance per BAM sample ──────────────────────────────────────
results <- vector("list", nrow(manifest))

for (i in seq_len(nrow(manifest))) {
  sid <- manifest$sample_id[i]

  # RNA allele fraction at WGS SNP positions
  rna_sid <- rna_pileup[sample_id == sid & !is.na(af)]
  rna_sid[, key := paste(chrom, pos, ref, alt, sep = ":")]

  vcf_sid <- vcf_dt[, .(CHROM, POS, REF, ALT, wgs_dosage = get(sid))]
  vcf_sid[, key := paste(CHROM, POS, REF, ALT, sep = ":")]
  vcf_sid <- vcf_sid[!is.na(wgs_dosage)]

  common <- merge(rna_sid[, .(key, rna_af = af)],
                  vcf_sid[, .(key, wgs_dosage)], by = "key")

  if (nrow(common) < opt$min_snps) {
    r2 <- NA_real_
    n_snps <- nrow(common)
    pass   <- FALSE
  } else {
    r2     <- cor(common$rna_af, common$wgs_dosage, use = "complete.obs")^2
    n_snps <- nrow(common)
    pass   <- r2 >= opt$r2_thresh
  }

  results[[i]] <- data.table(sample_id = sid, r2 = round(r2, 4),
                             n_snps = n_snps, pass = pass)
}

res_dt <- rbindlist(results)
manifest_ok   <- merge(manifest, res_dt[pass == TRUE],  by = "sample_id")
manifest_fail <- merge(manifest, res_dt[pass == FALSE | is.na(pass)], by = "sample_id")

fwrite(manifest_ok,   opt$out_ok,   sep = "\t")
fwrite(manifest_fail, opt$out_fail, sep = "\t")

message(sprintf("Identity check: %d passed, %d failed/insufficient SNPs",
                nrow(manifest_ok), nrow(manifest_fail)))
if (nrow(manifest_fail) > 0) {
  message("Failed samples:")
  print(manifest_fail[, .(sample_id, r2, n_snps)])
}
