# get_mirbase.R
# Retrieve miRBase 22 mature sequences for hsa-miRNAs using miRBaseVersions.db
# and write them to FASTA for use with miRanda.
#
# Usage:
#   conda activate mirdb
#   Rscript get_mirbase.R --out mirbase22_mature.fa
#
# Output: FASTA file of all hsa mature miRNA sequences from miRBase release 22.

suppressPackageStartupMessages({
  library(optparse)
  library(miRBaseVersions.db)
  library(RSQLite)
})

opt_list <- list(
  make_option("--out", type = "character", default = "mirbase22_mature.fa",
              help = "Output FASTA file [default: %default]")
)
opt <- parse_args(OptionParser(option_list = opt_list))

# Connect to the miRBase 22 database bundled with miRBaseVersions.db
db <- miRBaseVersions.db(release = 22)

# Pull mature sequences for Homo sapiens
sql <- "
  SELECT m.name, m.sequence
  FROM   mature         AS m
  JOIN   mirna          AS mi ON mi.id = m.mirna_id
  WHERE  mi.organism = 'hsa'
  AND    m.sequence IS NOT NULL
"
matures <- dbGetQuery(db, sql)
dbDisconnect(db)

if (nrow(matures) == 0) {
  stop("No hsa mature sequences found — check that miRBaseVersions.db is installed correctly.")
}

# Write FASTA
fa_lines <- unlist(lapply(seq_len(nrow(matures)), function(i) {
  c(paste0(">", matures$name[i]), matures$sequence[i])
}))
writeLines(fa_lines, opt$out)

message(sprintf("Written %d hsa mature miRNA sequences to %s", nrow(matures), opt$out))
