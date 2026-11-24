# SYN_02_rna_genotype_recheck.R
# Step 5a — RNA-based genotyping of synonymous SNVs from NYGC ALS RNA-seq BAMs
#
# Counts reads bearing REF and ALT alleles at each variant position using
# Rsamtools pileup(), computes allelic fractions, and writes a long-format
# summary table.
#
# NOTE: Individual-level genotypes and clinical data derived from the NYGC ALS
# Consortium are NOT included in this repository and are NOT written by this
# script. This script processes BAM files locally; the output rna_genotypes_long.tsv
# is excluded from the repository via .gitignore.
#
# Usage:
#   conda activate apa
#   Rscript SYN_02_rna_genotype_recheck.R \
#       --variants variants_scored.csv \
#       --bamdir   /path/to/bams \
#       --manifest manifest_bams.tsv \
#       --out      rna_genotypes_long.tsv
#
# Required packages (conda env: apa):
#   Rsamtools, GenomicRanges, data.table, optparse

suppressPackageStartupMessages({
  library(optparse)
  library(Rsamtools)
  library(GenomicRanges)
  library(data.table)
})

opt_list <- list(
  make_option("--variants", type = "character", help = "Scored variant CSV (SuppData1)"),
  make_option("--bamdir",   type = "character", help = "Directory containing BAM files"),
  make_option("--manifest", type = "character", help = "TSV: sample_id, bam_filename, diagnosis"),
  make_option("--out",      type = "character", default = "rna_genotypes_long.tsv",
              help = "Output long-format TSV [default: %default]"),
  make_option("--min_dp",   type = "integer",   default = 10L,
              help = "Minimum total depth to report a genotype [default: %default]"),
  make_option("--min_af",   type = "numeric",   default = 0.1,
              help = "Minimum ALT allele fraction to call carried [default: %default]")
)
opt <- parse_args(OptionParser(option_list = opt_list))

# ── load inputs ──────────────────────────────────────────────────────────────
variants <- fread(opt$variants)
manifest <- fread(opt$manifest, sep = "\t")

required_cols <- c("chrom", "pos", "ref", "alt", "gene")
if (!all(required_cols %in% names(variants))) {
  stop("variants file missing columns: ",
       paste(setdiff(required_cols, names(variants)), collapse = ", "))
}

# Build GRanges for pileup
gr_vars <- GRanges(
  seqnames = variants$chrom,
  ranges   = IRanges(start = variants$pos, end = variants$pos),
  ref      = variants$ref,
  alt      = variants$alt,
  gene     = variants$gene
)

# Pileup parameters — count each nucleotide per position
pup_params <- PileupParam(
  max_depth              = 1e5,
  min_base_quality       = 20,
  min_mapq               = 20,
  distinguish_strands    = FALSE,
  distinguish_nucleotides = TRUE
)

# ── process each BAM ─────────────────────────────────────────────────────────
results <- list()

for (i in seq_len(nrow(manifest))) {
  sample_id <- manifest$sample_id[i]
  bam_path  <- file.path(opt$bamdir, manifest$bam_filename[i])
  diagnosis <- manifest$diagnosis[i]

  if (!file.exists(bam_path)) {
    warning("BAM not found, skipping: ", bam_path)
    next
  }

  bf  <- BamFile(bam_path)
  pup <- pileup(bf, scanBamParam = ScanBamParam(which = gr_vars),
                pileupParam = pup_params)

  # Aggregate per position × nucleotide
  pup_dt <- as.data.table(pup)
  pup_dt[, pos := pos]   # already 1-based from Rsamtools

  # Join to variant table
  var_dt <- as.data.table(variants)[, .(chrom, pos, ref, alt, gene)]
  pup_dt <- merge(pup_dt, var_dt,
                  by.x = c("seqnames", "pos"),
                  by.y = c("chrom",    "pos"),
                  all.x = TRUE)

  pup_wide <- dcast(pup_dt, seqnames + pos + ref + alt + gene ~ nucleotide,
                    value.var = "count", fill = 0L)

  # Compute depth and allele fractions
  nuc_cols <- intersect(c("A", "C", "G", "T"), names(pup_wide))
  pup_wide[, depth := rowSums(.SD), .SDcols = nuc_cols]
  pup_wide[, ref_count := mapply(function(r, dp) {
    if (r %in% nuc_cols) .SD[[r]] else 0L }, ref, depth)]
  pup_wide[, alt_count := mapply(function(a, dp) {
    if (a %in% nuc_cols) .SD[[a]] else 0L }, alt, depth)]
  pup_wide[, af := ifelse(depth >= opt$min_dp, alt_count / depth, NA_real_)]
  pup_wide[, carried := (!is.na(af)) & (af >= opt$min_af)]
  pup_wide[, sample_id := sample_id]
  pup_wide[, diagnosis := diagnosis]

  results[[i]] <- pup_wide[, .(seqnames, pos, ref, alt, gene,
                                depth, ref_count, alt_count, af,
                                carried, sample_id, diagnosis)]
}

out_dt <- rbindlist(results)
setnames(out_dt, "seqnames", "chrom")
fwrite(out_dt, opt$out, sep = "\t", na = "NA")
message(sprintf("Written %d rows to %s", nrow(out_dt), opt$out))
