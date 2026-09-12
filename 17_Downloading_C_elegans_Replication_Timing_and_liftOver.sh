#!/usr/bin/env bash
set -euo pipefail

#SBATCH --job-name=17_Downloading_C_elegans_Replication_Timing_and_liftOver
#SBATCH --output=/home/ocdm0351/DPhil/logs/%x_%A.log
#SBATCH --error=/home/ocdm0351/DPhil/logs/%x_%A.err

module load R
module load deepTools/3.5.2-foss-2022a
module load Anaconda3
source activate UCSC_liftOver

# ============================================================
# C. elegans replication timing track
#
# Source:
#   GSE140804 Repli-seq WT Early and Late tracks
#
# Processing:
#   1. Download Early and Late processed Repli-seq tracks
#   2. Build continuous Early - Late replication-timing signal
#   3. Lift WBcel215 / UCSC ce10 -> WBcel235 / UCSC ce11
#   4. Leave ONLY the final WBcel235 bedGraph in R_Data
#
# Final output:
#   C_elegans_WBcel235_replication_timing_Early_minus_Late.bedgraph
#
# Final chromosome names:
#   chrI chrII chrIII chrIV chrV chrX
#
# Requirements:
#   bash, curl, wget, gunzip, Rscript, data.table, liftOver
# ============================================================

OUTDIR="${1:-/home/ocdm0351/DPhil/R_Data}"
BASE_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140804/suppl"

EARLY="GSE140804_Repliseq_WT_Early_final.bed.gz"
LATE="GSE140804_Repliseq_WT_Late_final.bed.gz"

CHAIN_DIR="/home/ocdm0351/DPhil/liftingOver/chain_files/C_elegans"
CHAIN="$CHAIN_DIR/ce10ToCe11.over.chain"

LOG_DIR="/home/ocdm0351/DPhil/logs"

mkdir -p "$OUTDIR" "$CHAIN_DIR" "$LOG_DIR"

# ------------------------------------------------------------
# Check dependencies
# ------------------------------------------------------------

command -v curl >/dev/null 2>&1 || {
    echo "ERROR: curl is required." >&2
    exit 1
}

command -v wget >/dev/null 2>&1 || {
    echo "ERROR: wget is required." >&2
    exit 1
}

command -v gunzip >/dev/null 2>&1 || {
    echo "ERROR: gunzip is required." >&2
    exit 1
}

command -v Rscript >/dev/null 2>&1 || {
    echo "ERROR: Rscript is required." >&2
    exit 1
}

command -v liftOver >/dev/null 2>&1 || {
    echo "ERROR: UCSC liftOver is required and must be on PATH." >&2
    exit 1
}

# ------------------------------------------------------------
# Temporary working directory
# Nothing here will remain in R_Data.
# ------------------------------------------------------------

TMPDIR=$(mktemp -d "${OUTDIR}/.C_elegans_replication_timing_XXXXXX")
trap 'rm -rf "$TMPDIR"' EXIT

EARLY_PATH="$TMPDIR/$EARLY"
LATE_PATH="$TMPDIR/$LATE"
PRELIFT="$TMPDIR/C_elegans_replication_timing_Early_minus_Late.ce10.bedgraph"
UCSC_INPUT="$TMPDIR/C_elegans_replication_timing_Early_minus_Late.ce10.ucsc.bedgraph"
LIFTED_TMP="$TMPDIR/C_elegans_WBcel235_replication_timing_Early_minus_Late.bedgraph"
UNMAPPED="$LOG_DIR/unmapped_C_elegans_replication_timing_Early_minus_Late.bedgraph"
FINAL="$OUTDIR/C_elegans_WBcel235_replication_timing_Early_minus_Late.bedgraph"
QC="$LOG_DIR/C_elegans_replication_timing_QC.tsv"

# ------------------------------------------------------------
# Download WT Early and Late Repli-seq tracks
# ------------------------------------------------------------

echo "Downloading WT Early Repli-seq..."
curl -L --fail --retry 3 \
    -o "$EARLY_PATH" \
    "$BASE_URL/$EARLY"

echo "Downloading WT Late Repli-seq..."
curl -L --fail --retry 3 \
    -o "$LATE_PATH" \
    "$BASE_URL/$LATE"

# ------------------------------------------------------------
# Build continuous Early - Late track in temporary space
# ------------------------------------------------------------

echo
echo "Building continuous replication-timing track..."

Rscript - "$TMPDIR" "$EARLY" "$LATE" "$PRELIFT" "$QC" <<'RSCRIPT'

args <- commandArgs(trailingOnly = TRUE)

tmpdir <- args[[1]]
early_file <- file.path(tmpdir, args[[2]])
late_file <- file.path(tmpdir, args[[3]])
prelift_file <- args[[4]]
qc_file <- args[[5]]

if (!requireNamespace("data.table", quietly = TRUE)) {
    stop(
        "R package 'data.table' is required. Install with install.packages('data.table')."
    )
}

library(data.table)

read_track <- function(path) {

    x <- fread(
        cmd = paste("gzip -dc", shQuote(path)),
        header = FALSE,
        data.table = TRUE
    )

    # Remove UCSC track header lines.
    x <- x[!grepl("^track", V1)]

    if (ncol(x) < 4L) {
        stop("Expected at least 4 columns in: ", path)
    }

    setnames(x, 1:4, c("chrom", "start", "end", "signal"))

    x <- x[, .(
        chrom = as.character(chrom),
        start = as.integer(start),
        end = as.integer(end),
        signal = as.numeric(signal)
    )]

    x <- x[
        !is.na(chrom) &
        !is.na(start) &
        !is.na(end) &
        !is.na(signal)
    ]

    # GEO file uses 1-based inclusive coordinates.
    # Convert to 0-based, half-open bedGraph coordinates.
    x[, start := start - 1L]

    x
}

early <- read_track(early_file)
late <- read_track(late_file)

cat("Early bins:", nrow(early), "\n")
cat("Late bins:", nrow(late), "\n")

# Check that the genomic bins match exactly.
bins_e <- unique(early[, .(chrom, start, end)])
bins_l <- unique(late[, .(chrom, start, end)])

setkey(bins_e, chrom, start, end)
setkey(bins_l, chrom, start, end)

if (!identical(bins_e, bins_l)) {
    stop("Early and Late tracks do not have identical genomic bins.")
}

# Merge signals.
rt <- merge(
    early[, .(chrom, start, end, early_signal = signal)],
    late[, .(chrom, start, end, late_signal = signal)],
    by = c("chrom", "start", "end"),
    all = FALSE
)

# Continuous relative replication timing.
rt[, replication_timing := early_signal - late_signal]

# Validate chromosome names.
chrom_order <- c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrX")
rt[, chrom_order_num := match(chrom, chrom_order)]

if (any(is.na(rt$chrom_order_num))) {
    bad <- unique(rt[is.na(chrom_order_num), chrom])
    stop(
        "Unexpected chromosome names: ",
        paste(bad, collapse = ", ")
    )
}

setorder(rt, chrom_order_num, start, end)

out <- rt[, .(chrom, start, end, replication_timing)]

# Write temporary pre-liftover track.
fwrite(
    out,
    prelift_file,
    sep = "\t",
    col.names = FALSE,
    quote = FALSE
)

# QC is written to logs, not R_Data, so it cannot be picked up by the pipeline.
qc <- data.table(
    metric = c(
        "early_bins",
        "late_bins",
        "output_bins",
        "chromosomes",
        "bin_size_min",
        "bin_size_max",
        "score_min",
        "score_max"
    ),
    value = c(
        nrow(early),
        nrow(late),
        nrow(out),
        length(unique(out$chrom)),
        min(out$end - out$start),
        max(out$end - out$start),
        min(out$replication_timing),
        max(out$replication_timing)
    )
)

fwrite(qc, qc_file, sep = "\t", quote = FALSE)

cat("Temporary pre-liftover track:", prelift_file, "\n")
cat("QC:", qc_file, "\n")

RSCRIPT

# ------------------------------------------------------------
# Download ce10 -> ce11 liftOver chain
# ------------------------------------------------------------

echo
echo "Downloading ce10 -> ce11 liftOver chain..."

wget -q -O - \
    "https://hgdownload.soe.ucsc.edu/goldenPath/ce10/liftOver/ce10ToCe11.over.chain.gz" \
    | gunzip -c > "$CHAIN"

# ------------------------------------------------------------
# Input is already using UCSC chromosome names:
# chrI chrII chrIII chrIV chrV chrX
# ------------------------------------------------------------

cp "$PRELIFT" "$UCSC_INPUT"

# ------------------------------------------------------------
# LiftOver WBcel215 / ce10 -> WBcel235 / ce11
# ------------------------------------------------------------

echo
echo "Running UCSC liftOver..."

echo "Input:  $UCSC_INPUT"
echo "Output: $FINAL"

liftOver \
    "$UCSC_INPUT" \
    "$CHAIN" \
    "$LIFTED_TMP" \
    "$UNMAPPED"

# ------------------------------------------------------------
# Sort and write ONLY the final bedGraph to R_Data
# ------------------------------------------------------------

sort \
    -k1,1 \
    -k2,2n \
    -k3,3n \
    "$LIFTED_TMP" \
    > "$FINAL"

# ------------------------------------------------------------
# Final validation
# ------------------------------------------------------------

if [[ ! -s "$FINAL" ]]; then
    echo "ERROR: final WBcel235 bedGraph is empty." >&2
    echo "Check: $UNMAPPED" >&2
    exit 1
fi

echo
echo "============================================================"
echo "Done"
echo "============================================================"
echo "Final WBcel235 bedGraph:"
echo "  $FINAL"
echo
echo "QC:"
echo "  $QC"
echo
echo "Unmapped intervals:"
echo "  $UNMAPPED"
echo
echo "The temporary source/pre-liftover files will be removed."
echo "No other replication-timing bedGraph will remain in R_Data."
echo "============================================================"
