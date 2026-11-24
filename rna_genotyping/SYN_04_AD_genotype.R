# SYN_04_AD_genotype.R
# Step 5c — RNA-based genotyping in the Alzheimer's disease brain RNA-seq cohort
#            (GSE125583; 289 donors; Mathys et al. 2019 Nature)
#
# Reuses the pileup approach of SYN_02 on the GSE125583 pseudo-bulk BAMs.
# Outputs carrier counts and allele fractions per variant per diagnosis group.
# No individual-level clinical data from the NYGC ALS Consortium are used here.
#
# Usage:
#   conda activate apa
#   Rscript SYN_04_AD_genotype.R \
#       --variants variants_scored.csv \
#       --bamdir   /path/to/gse125583/bams \
#       --manifest manifest_gse125583.tsv \
#       --out      ad_carrier_summary.tsv

suppressPackageStartupMessages({
  library(optparse)
  library(Rsamtools)
  library(GenomicRanges)
  library(data.table)
})

opt_list <- list(
  make_option("--variants",   type = "character", help = "Scored variant CSV"),
  make_option("--bamdir",     type = "character", help = "Directory of GSE125583 BAMs"),
  make_option("--manifest",   type = "character", help = "TSV: sample_id, bam_filename, diagnosis"),
  make_option("--out",        type = "character", default = "ad_carrier_summary.tsv"),
  make_option("--min_dp",     type = "integer",   default = 10L),
  make_option("--min_af",     type = "numeric",   default = 0.1,
              help = "ALT AF threshold for calling a carrier [default: %default]")
)
opt <- parse_args(OptionParser(option_list = opt_list))

variants <- fread(opt$variants)
manifest <- fread(opt$manifest, sep = "\t")

gr_vars <- GRanges(
  seqnames = variants$chrom,
  ranges   = IRanges(start = variants$pos, end = variants$pos),
  ref      = variants$ref,
  alt      = variants$alt,
  gene     = variants$gene
)

pup_params <- PileupParam(
  max_depth              = 1e5,
  min_base_quality       = 20,
  min_mapq               = 20,
  distinguish_strands    = FALSE,
  distinguish_nucleotides = TRUE
)

all_pileups <- list()

for (i in seq_len(nrow(manifest))) {
  sample_id <- manifest$sample_id[i]
  bam_path  <- file.path(opt$bamdir, manifest$bam_filename[i])
  diagnosis <- manifest$diagnosis[i]

  if (!file.exists(bam_path)) {
    warning("BAM not found, skipping: ", bam_path)
    next
  }

  bf  <- BamFile(bam_path)
  pup <- as.data.table(pileup(bf,
           scanBamParam = ScanBamParam(which = gr_vars),
           pileupParam  = pup_params))

  var_dt <- as.data.table(variants)[, .(chrom, pos, ref, alt, gene)]
  pup    <- merge(pup, var_dt,
                  by.x = c("seqnames", "pos"),
                  by.y = c("chrom",    "pos"), all.x = TRUE)

  pup_w <- dcast(pup, seqnames + pos + ref + alt + gene ~ nucleotide,
                 value.var = "count", fill = 0L)
  nuc_cols <- intersect(c("A","C","G","T"), names(pup_w))
  pup_w[, depth     := rowSums(.SD), .SDcols = nuc_cols]
  pup_w[, alt_count := mapply(function(a) if (a %in% nuc_cols) .SD[[a]] else 0L, alt)]
  pup_w[, af        := ifelse(depth >= opt$min_dp, alt_count / depth, NA_real_)]
  pup_w[, carried   := (!is.na(af)) & (af >= opt$min_af)]
  pup_w[, sample_id := sample_id]
  pup_w[, diagnosis := diagnosis]

  all_pileups[[i]] <- pup_w[, .(seqnames, pos, ref, alt, gene,
                                 depth, alt_count, af, carried,
                                 sample_id, diagnosis)]
}

long_dt <- rbindlist(all_pileups)
setnames(long_dt, "seqnames", "chrom")

# Summarise: carrier counts and mean AF per variant × diagnosis
summary_dt <- long_dt[, .(
    n_tested  = .N,
    n_carrier = sum(carried, na.rm = TRUE),
    mean_af   = round(mean(af[carried], na.rm = TRUE), 4)
  ), by = .(chrom, pos, ref, alt, gene, diagnosis)]

fwrite(summary_dt, opt$out, sep = "\t", na = "NA")
message(sprintf("Written %d rows to %s", nrow(summary_dt), opt$out))
