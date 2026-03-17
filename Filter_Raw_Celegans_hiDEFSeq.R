#!/usr/bin/env Rscript
#
# Usage:
#   Rscript Filtering_Raw_Celegans_hiDEFSeq.R <input.tdv> [reference.fa]
#
# Outputs (in current directory):
#   duplex_variants.tsv   -> duplex-supported high-confidence de novos
#   ss_variants.tsv       -> single-strand high-confidence events
#   diagnostics.txt       -> summary stats and filter counts
#
# Dependencies: data.table, dplyr, stringr, Biostrings, Rsamtools (Biostrings/Rsamtools only if reference.fa provided)
#
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1) stop("Provide input TDV file path as first argument.\nUsage: Rscript <script.R> <input.tdv> [reference.fa]")

tdv_path <- args[1]
fasta_path <- ifelse(length(args) >= 2, args[2], NA)

# ---------------------------
# Configurable thresholds
# ---------------------------
min_mapq <- 60                    # keep only MAPQ == 60
allow_multi_map_matches <- 1      # num_other_ref_matches should equal this (1 means unique)
min_dist_to_end <- 50             # read-level: distance from end (rms_dist_to_end)
max_mismatch_per_read <- 10       # read-level: exclude very mismatchy reads
max_rms_indel_count <- 5          # read-level: exclude reads with many indels (tune)
min_site_depth <- 50              # site-level: minimum depth to consider (total reads supporting site)
min_alt_reads <- 3                # site-level: minimum total alt reads
min_alt_in_each_strand_for_duplex <- 1  # duplex requires >= this on each duplex strand
min_distinct_barcodes <- 1        # require support from at least N distinct barcodes / molecules (set 2 for stricter)
cluster_bp <- 5                   # remove variants with another variant within +/- cluster_bp
homopolymer_len <- 6              # mask variants in >= this homopolymer length
# VAF upper bound to remove common variants (we look for very rare events)
max_vaf_allowed <- 0.05           # remove high-frequency variants (likely recurrent/systematic)
# ---------------------------

# helper to parse "x|y" fields into two numeric columns
split_pair <- function(x) {
  # return matrix n x 2
  sapply(x, function(s) {
    if (is.na(s) || s == "") return(c(NA_integer_, NA_integer_))
    parts <- str_split(s, "\\|")[[1]]
    # sometimes counts have more than 2? take first two
    parts <- parts[1:2]
    parts[is.na(parts)] <- "0"
    as.integer(parts)
  }) %>% t()
}

cat("Reading TDV (using data.table::fread). This may take a while for large files...\n")
dt <- fread(tdv_path, sep = "\t", header = TRUE, data.table = FALSE)

# remove leading # from first column name if present
colnames(dt)[1] <- sub("^#", "", colnames(dt)[1])

if(nrow(dt) == 0) stop("No rows read from TDV. Check file and header.")

# Inspect columns we know; if names differ adapt here
# Expected columns from your head(): chrom,pos,alignment_pos,ref,alt,read_id,barcode,ref_mapq,num_other_ref_matches,coverage,ref_count,alt_count,deletion_count,...
# We'll parse fields that are pipe-separated (many fields like ref_count, alt_count etc).
cat("Parsing strand-paired fields...\n")
# Fields with x|y representation - adjust this vector if more/less columns exist
pair_cols <- c("ref_count","alt_count","deletion_count",
               "rms_ref_matches_10x","rms_alt_matches_10x",
               "rms_ref_indel_10x","rms_alt_indel_10x",
               "rms_indel_count","rms_mismatch_per_read",
               "rms_subread_len","rms_dist_to_end")

# Only keep ones that exist
pair_cols <- intersect(pair_cols, colnames(dt))

# create columns with suffix _f (_r for strand2)
for (col in pair_cols) {
  mat <- split_pair(dt[[col]])
  # names
  fcol <- paste0(col, "_s1")
  rcol <- paste0(col, "_s2")
  # If the column originally numeric without pipe, treat as total and split as NA/NA
  if(ncol(mat) == 0) {
    dt[[fcol]] <- NA_integer_
    dt[[rcol]] <- NA_integer_
  } else {
    dt[[fcol]] <- as.integer(mat[,1])
    dt[[rcol]] <- as.integer(mat[,2])
  }
}

# ensure numeric conversion for some single-valued columns
numcols <- c("ref_mapq","num_other_ref_matches","coverage","rms_mismatch_per_read")
numcols <- intersect(numcols, colnames(dt))
for (cname in numcols) dt[[cname]] <- as.numeric(dt[[cname]])

# Some fields might be combined (like ref_count column showing "11|10" etc). We created ref_count_s1 and _s2.
# For alt counts we expect alt_count_s1 and alt_count_s2
if(!("alt_count_s1" %in% colnames(dt))) stop("alt_count parsing failed - ensure 'alt_count' column exists in the TDV and uses x|y format.")

# Basic read-level filters
cat("Applying read-level filters...\n")
initial_n_reads <- nrow(dt)

read_filter <- rep(TRUE, nrow(dt))
# MAPQ
if("ref_mapq" %in% colnames(dt)) {
  read_filter <- read_filter & (dt$ref_mapq >= min_mapq)
}
# unique mapping
if("num_other_ref_matches" %in% colnames(dt)) {
  read_filter <- read_filter & (dt$num_other_ref_matches == allow_multi_map_matches)
}
# distance to end (we have per-strand; require both strands? We'll require the maximum of the two)
if("rms_dist_to_end_s1" %in% colnames(dt) & "rms_dist_to_end_s2" %in% colnames(dt)) {
  dt$rms_dist_to_end_max <- pmax(dt$rms_dist_to_end_s1, dt$rms_dist_to_end_s2, na.rm = TRUE)
  read_filter <- read_filter & (dt$rms_dist_to_end_max >= min_dist_to_end)
} else if("rms_dist_to_end" %in% colnames(dt)) {
  read_filter <- read_filter & (dt$rms_dist_to_end >= min_dist_to_end)
}
# mismatch rate (use per-row summary column if present otherwise per strand max)
if("rms_mismatch_per_read" %in% colnames(dt)) {
  read_filter <- read_filter & (dt$rms_mismatch_per_read <= max_mismatch_per_read)
} else if("rms_mismatch_per_read_s1" %in% colnames(dt)) {
  dt$rms_mm_max <- pmax(dt$rms_mismatch_per_read_s1, dt$rms_mismatch_per_read_s2, na.rm = TRUE)
  read_filter <- read_filter & (dt$rms_mm_max <= max_mismatch_per_read)
}
# indel count per-read
if("rms_indel_count_s1" %in% colnames(dt) & "rms_indel_count_s2" %in% colnames(dt)) {
  dt$rms_indel_max <- pmax(dt$rms_indel_count_s1, dt$rms_indel_count_s2, na.rm = TRUE)
  read_filter <- read_filter & (dt$rms_indel_max <= max_rms_indel_count)
}

dt_filtered_reads <- dt[read_filter, ]
filtered_n_reads <- nrow(dt_filtered_reads)

cat(sprintf("Reads: initial=%d , after_read_filters=%d (removed %d)\n",
            initial_n_reads, filtered_n_reads, initial_n_reads - filtered_n_reads))

# ---------------------------
# Aggregate to site-level
# ---------------------------
cat("Aggregating to site-level (chrom,pos,ref,alt)...\n")
# create site key
dt_filtered_reads$site_key <- paste(dt_filtered_reads$chrom, dt_filtered_reads$pos, dt_filtered_reads$ref, dt_filtered_reads$alt, sep = ":")

# compute per-row alt totals and per-strand alt counts (they already parsed)
dt_filtered_reads$alt_total_row <- rowSums(cbind(dt_filtered_reads$alt_count_s1, dt_filtered_reads$alt_count_s2), na.rm = TRUE)
dt_filtered_reads$alt_s1 <- dt_filtered_reads$alt_count_s1
dt_filtered_reads$alt_s2 <- dt_filtered_reads$alt_count_s2

# We'll treat barcode/read_id as molecule identifiers if available
if(!("barcode" %in% colnames(dt_filtered_reads))) {
  dt_filtered_reads$barcode <- dt_filtered_reads$read_id
}

# site aggregation
site_dt <- dt_filtered_reads %>%
  group_by(site_key, chrom, pos, ref, alt) %>%
  summarise(
    site_depth = n(),                              # number of reads covering the site (after read filters)
    alt_reads = sum(alt_total_row, na.rm = TRUE),  # total alt-supporting read fragments
    alt_s1 = sum(alt_s1, na.rm = TRUE),
    alt_s2 = sum(alt_s2, na.rm = TRUE),
    distinct_barcodes = n_distinct(barcode),
    mean_mapq = mean(ref_mapq, na.rm = TRUE),
    .groups = "drop"
  ) %>% as.data.frame()

# compute VAF using site_depth as denominator: here site_depth is count of rows (reads),
# but if you want per-base coverage you may want to use a field 'coverage' if that denotes coverage at site.
if("coverage" %in% colnames(dt_filtered_reads)) {
  # coverage per-row is the same for that read? More correct: use sum of coverage? ambiguous.
  # We'll compute VAF using alt_reads / (alt_reads + ref_reads) if ref_count parsed available
  if("ref_count_s1" %in% colnames(dt_filtered_reads)) {
    # compute per-site total ref support
    ref_support_site <- dt_filtered_reads %>%
      group_by(site_key) %>%
      summarise(ref_reads = sum(ref_count_s1 + ref_count_s2, na.rm = TRUE), .groups = "drop")
    site_dt <- merge(site_dt, ref_support_site, by = "site_key", all.x = TRUE)
    site_dt$site_depth_bases <- site_dt$alt_reads + site_dt$ref_reads
    site_dt$vaf <- site_dt$alt_reads / pmax(site_dt$site_depth_bases, 1)
  } else {
    # fallback: use alt_reads / site_depth (reads)
    site_dt$vaf <- site_dt$alt_reads / pmax(site_dt$site_depth, 1)
    site_dt$site_depth_bases <- NA
  }
} else {
  site_dt$vaf <- site_dt$alt_reads / pmax(site_dt$site_depth, 1)
  site_dt$site_depth_bases <- NA
}

cat(sprintf("Number of aggregated sites: %d\n", nrow(site_dt)))

# ---------------------------
# Site-level filters for de novo (initial)
# ---------------------------
cat("Applying site-level candidate filters (depth, alt reads, vaf upper bound)...\n")
pre_site_n <- nrow(site_dt)
site_filter <- with(site_dt,
                    (site_depth >= min_site_depth) &
                      (alt_reads >= min_alt_reads) &
                      (vaf <= max_vaf_allowed) &
                      (distinct_barcodes >= min_distinct_barcodes)
)
site_candidates <- site_dt[site_filter, ]
cat(sprintf("Sites passing basic site filters: %d (from %d)\n", nrow(site_candidates), pre_site_n))

# ---------------------------
# Classify duplex vs single-strand
# ---------------------------
cat("Classifying duplex-supported vs single-strand events...\n")
site_candidates$duplex_supported <- (site_candidates$alt_s1 >= min_alt_in_each_strand_for_duplex) &
  (site_candidates$alt_s2 >= min_alt_in_each_strand_for_duplex)

site_candidates$single_strand <- ( (site_candidates$alt_s1 > 0 & site_candidates$alt_s2 == 0) |
                                     (site_candidates$alt_s2 > 0 & site_candidates$alt_s1 == 0) )

# sanity: some sites might be both (if counts >= threshold on both) - then duplex_supported TRUE
# keep both classification columns to allow separate analyses

# ---------------------------
# Remove clustered variants (within +/- cluster_bp)
# ---------------------------
cat("Removing clustered variants within +/- ", cluster_bp, " bp\n", sep = "")
site_candidates <- site_candidates %>% arrange(chrom, pos) %>% as.data.frame()

# function to mark clusters
site_candidates$cluster_flag <- FALSE
if(nrow(site_candidates) > 1) {
  for(chr in unique(site_candidates$chrom)) {
    sidx <- which(site_candidates$chrom == chr)
    poses <- site_candidates$pos[sidx]
    if(length(poses) <= 1) next
    # compute distance to previous and next
    dprev <- c(Inf, diff(poses))
    dnext <- c(diff(poses), Inf)
    cluster_here <- (dprev <= cluster_bp) | (dnext <= cluster_bp)
    site_candidates$cluster_flag[sidx] <- cluster_here
  }
}
cat(sprintf("Clustered sites flagged: %d\n", sum(site_candidates$cluster_flag)))

# ---------------------------
# Optional: mask homopolymers using reference FASTA if provided
# ---------------------------
mask_hp <- rep(FALSE, nrow(site_candidates))
if(!is.na(fasta_path) && file.exists(fasta_path)) {
  cat("Reference FASTA provided; checking homopolymers (requires Biostrings & Rsamtools)...\n")
  suppressPackageStartupMessages({
    library(Rsamtools)
    library(Biostrings)
  })
  fa <- FaFile(fasta_path)
  open(fa)
  # get window +/- (homopolymer_len) around site and test for homopolymer run covering the variant position
  halfwin <- homopolymer_len
  for(i in seq_len(nrow(site_candidates))) {
    chr <- as.character(site_candidates$chrom[i])
    pos <- as.integer(site_candidates$pos[i])
    start <- max(1, pos - halfwin)
    end <- pos + halfwin
    # try/catch because some chromosomes may not exist in FASTA
    ok <- TRUE
    seq_region <- tryCatch({
      getSeq(fa, chr, start = start, end = end)
    }, error = function(e) {
      ok <<- FALSE
      DNAString("")
    })
    if(!ok) {
      mask_hp[i] <- FALSE; next
    }
    s <- as.character(seq_region)
    # compute longest homopolymer run that covers central position
    # find run lengths
    chars <- strsplit(s, "")[[1]]
    runs <- rle(chars)
    # rebuild run positions
    ends <- cumsum(runs$lengths)
    starts <- ends - runs$lengths + 1
    # map centre position in window to index in s
    centre_idx <- pos - start + 1
    # find which run contains centre_idx
    which_run <- which( (starts <= centre_idx) & (ends >= centre_idx) )
    if(length(which_run) == 1) {
      runlen <- runs$lengths[which_run]
      if(runlen >= homopolymer_len) mask_hp[i] <- TRUE
    }
  }
  close(fa)
  cat(sprintf("Sites in homopolymers >=%d: %d\n", homopolymer_len, sum(mask_hp)))
} else {
  cat("No FASTA provided or file not found. Skipping homopolymer masking.\n")
}

site_candidates$homopolymer_flag <- mask_hp

# ---------------------------
# Build final filtered sets
# ---------------------------
# Duplex high-confidence: duplex_supported TRUE, not clustered, not homopolymer
duplex_hc <- site_candidates %>%
  filter(duplex_supported == TRUE,
         cluster_flag == FALSE,
         homopolymer_flag == FALSE)

# Single-strand high-confidence: single_strand TRUE, not clustered, not homopolymer
ss_hc <- site_candidates %>%
  filter(single_strand == TRUE,
         cluster_flag == FALSE,
         homopolymer_flag == FALSE)

# Additional optional strictness: require at least 2 distinct barcodes for duplex as extra confidence
if(min_distinct_barcodes >= 2) {
  duplex_hc <- duplex_hc %>% filter(distinct_barcodes >= min_distinct_barcodes)
  ss_hc <- ss_hc %>% filter(distinct_barcodes >= min_distinct_barcodes)
}

cat(sprintf("Final duplex HC sites: %d\n", nrow(duplex_hc)))
cat(sprintf("Final single-strand HC sites: %d\n", nrow(ss_hc)))

# ---------------------------
# Output results
# ---------------------------
cat("Writing output files...\n")
write.table(duplex_hc, file = "duplex_variants.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ss_hc, file = "ss_variants.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# diagnostics
diag <- c(
  paste0("Input TDV: ", tdv_path),
  paste0("Reads initial: ", initial_n_reads),
  paste0("Reads after read-filters: ", filtered_n_reads),
  paste0("Sites aggregated: ", nrow(site_dt)),
  paste0("Sites after basic site filters: ", nrow(site_candidates)),
  paste0("Clustered sites flagged: ", sum(site_candidates$cluster_flag)),
  paste0("Homopolymer masked sites (if FASTA provided): ", sum(site_candidates$homopolymer_flag)),
  paste0("Final duplex HC sites: ", nrow(duplex_hc)),
  paste0("Final single-strand HC sites: ", nrow(ss_hc)),
  paste0("Thresholds used: min_mapq=", min_mapq,
         ", min_dist_to_end=", min_dist_to_end,
         ", max_mismatch_per_read=", max_mismatch_per_read,
         ", min_site_depth=", min_site_depth,
         ", min_alt_reads=", min_alt_reads,
         ", max_vaf_allowed=", max_vaf_allowed,
         ", cluster_bp=", cluster_bp,
         ", homopolymer_len=", homopolymer_len)
)
writeLines(diag, con = "diagnostics.txt")
cat("Done. Outputs: duplex_variants.tsv , ss_variants.tsv , diagnostics.txt\n")
