#!/usr/bin/env bash
set -euo pipefail

#SBATCH --job-name=17_Downloading_C_elegans_Recombination_Rates
#SBATCH --output=/home/ocdm0351/DPhil/logs/%x_%A.log
#SBATCH --error=/home/ocdm0351/DPhil/logs/%x_%A.err

# ============================================================
# C. elegans WT sex-specific recombination rate -> WBcel235
#
# Source:
#   Bush et al. 2026, PLOS Genetics 10.1371/journal.pgen.1012237
#   S1 WT crossover CSV
#
# Source assembly:
#   Libuda/DLW N2 Bristol, GCA_028201515.1
#
# Target assembly:
#   Local WBcel235 / WS235 FASTA
#
# Workflow:
#   1. Download S1 WT crossover calls
#   2. Download Libuda/DLW N2 Bristol assembly
#   3. Align Libuda assembly -> local WBcel235 with minimap2
#   4. Lift UNIQUE source intervals with paftools.js
#   5. Map lifted intervals back to every original S1 event sharing that interval
#   6. Split oocyte / spermatocyte events
#   7. Calculate 200-kb windows with 100-kb step
#   8. Output cM/Mb bedGraphs using the paper's successfully mapped
#      WT meiosis denominators: 297 oocytes and 300 spermatocytes
#
# Requirements:
#   - minimap2 module
#   - Rscript in pandoc_env
#   - R package data.table
#   - local k8 + paftools.js
# ============================================================

module load minimap2
module load Anaconda3
source activate pandoc_env

OUTDIR="${1:-/home/ocdm0351/DPhil/R_Data}"
LOG_DIR="/home/ocdm0351/DPhil/logs"
WORKDIR="$OUTDIR/.C_elegans_recombination_work"

mkdir -p "$OUTDIR" "$LOG_DIR" "$WORKDIR"

# ------------------------------------------------------------
# Fixed paths
# ------------------------------------------------------------

WBCEL235_FASTA="/home/ocdm0351/.conda/envs/snpEff/share/snpeff-5.2-2/data/C_elegans/WS235/sequences.fa"
K8="$HOME/DPhil/software/minimap2_helpers/k8"
PAFTOOLS="$HOME/DPhil/software/minimap2_helpers/paftools.js"

S1_URL="https://journals.plos.org/plosgenetics/article/file?type=supplementary&id=10.1371/journal.pgen.1012237.s001"
S1_FILE="$WORKDIR/S1_WT_crossover_data.csv"

LIBUDA_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/028/201/515/GCA_028201515.1_ASM2820151v1/GCA_028201515.1_ASM2820151v1_genomic.fna.gz"
LIBUDA_GZ="$WORKDIR/GCA_028201515.1_DLW_N2_Bristol_genomic.fna.gz"
LIBUDA_FASTA="$WORKDIR/GCA_028201515.1_DLW_N2_Bristol_genomic.fna"

PAF="$WORKDIR/Libuda_to_WBcel235.paf"
PAF_SORTED="$WORKDIR/Libuda_to_WBcel235.sorted.paf"

# Unique intervals are what we actually lift. This avoids ambiguity when the
# same crossover coordinates occur in more than one S1 sample.
LIFT_INPUT="$WORKDIR/S1_unique_source_intervals_for_liftover.bed"
LIFT_OUTPUT="$WORKDIR/S1_unique_source_intervals_WBcel235.bed"

# All-event metadata maps each original S1 event to its source interval key.
META_FILE="$WORKDIR/S1_WT_event_metadata.tsv"

OOCYTE_MIDPOINTS="$WORKDIR/oocyte_crossover_midpoints_WBcel235.bed"
SPERM_MIDPOINTS="$WORKDIR/spermatocyte_crossover_midpoints_WBcel235.bed"

OOCYTE_OUT="$OUTDIR/C_elegans_WBcel235_oocyte_Recombination.bedgraph"
SPERM_OUT="$OUTDIR/C_elegans_WBcel235_spermatocyte_Recombination.bedgraph"

QC_OUT="$LOG_DIR/C_elegans_recombination_QC.tsv"
UNMAPPED="$LOG_DIR/unmapped_C_elegans_WT_crossover_intervals.bed"

# Paper-reported successfully called WT meioses.
# The paper sequenced 300 oocyte-derived and 310 spermatocyte-derived samples,
# with crossover calls successfully obtained for 297 and 300, respectively.
N_OOCYTE_MEIOSES=297
N_SPERMATOCYTE_MEIOSES=300

# WBcel235 chromosome lengths from the local reference headers.
declare -A CHR_LEN=(
    [I]=15072434
    [II]=15279421
    [III]=13783801
    [IV]=17493829
    [V]=20924180
    [X]=17718942
)

# ------------------------------------------------------------
# Check dependencies
# ------------------------------------------------------------

command -v minimap2 >/dev/null 2>&1 || {
    echo "ERROR: minimap2 not found on PATH." >&2
    exit 1
}

command -v curl >/dev/null 2>&1 || {
    echo "ERROR: curl not found on PATH." >&2
    exit 1
}

command -v gzip >/dev/null 2>&1 || {
    echo "ERROR: gzip not found on PATH." >&2
    exit 1
}

command -v sort >/dev/null 2>&1 || {
    echo "ERROR: sort not found on PATH." >&2
    exit 1
}

command -v Rscript >/dev/null 2>&1 || {
    echo "ERROR: Rscript not found. Check pandoc_env." >&2
    exit 1
}

[[ -x "$K8" ]] || {
    echo "ERROR: k8 is not executable: $K8" >&2
    exit 1
}

[[ -r "$PAFTOOLS" ]] || {
    echo "ERROR: paftools.js is not readable: $PAFTOOLS" >&2
    exit 1
}

[[ -f "$WBCEL235_FASTA" ]] || {
    echo "ERROR: WBcel235 FASTA not found: $WBCEL235_FASTA" >&2
    exit 1
}

Rscript -e 'quit(status = if (requireNamespace("data.table", quietly = TRUE)) 0 else 1)' || {
    echo "ERROR: R package data.table is not available in pandoc_env." >&2
    exit 1
}

echo "Dependencies found."
echo "  minimap2: $(command -v minimap2)"
echo "  Rscript:  $(command -v Rscript)"
echo "  k8:       $K8"
echo "  paftools: $PAFTOOLS"
echo ""
echo "Using local WBcel235 FASTA:"
echo "  $WBCEL235_FASTA"

# ------------------------------------------------------------
# Remove stale outputs from previous runs
# ------------------------------------------------------------

rm -f "$OOCYTE_OUT" "$SPERM_OUT"
rm -f "$OOCYTE_MIDPOINTS" "$SPERM_MIDPOINTS"
rm -f "$LIFT_INPUT" "$LIFT_OUTPUT" "$META_FILE"
rm -f "$UNMAPPED" "$UNMAPPED.ambiguous.tsv" "$UNMAPPED.partial.tsv"
rm -f "$WORKDIR/paftools_liftover.err"

# ------------------------------------------------------------
# Download S1 WT crossover CSV
# ------------------------------------------------------------

echo
echo "Downloading WT crossover S1 file..."

curl -L --fail --retry 3 \
    -o "$S1_FILE" \
    "$S1_URL"

# ------------------------------------------------------------
# Download Libuda/DLW N2 Bristol assembly
# ------------------------------------------------------------

echo
echo "Downloading Libuda/DLW N2 Bristol assembly..."

curl -L --fail --retry 3 \
    -o "$LIBUDA_GZ" \
    "$LIBUDA_URL"

gzip -dc "$LIBUDA_GZ" > "$LIBUDA_FASTA"

# ------------------------------------------------------------
# Show assembly headers
# ------------------------------------------------------------

echo
echo "Source assembly headers:"
grep '^>' "$LIBUDA_FASTA" | head -6

echo
echo "Target assembly headers:"
grep '^>' "$WBCEL235_FASTA" | head -7

# ------------------------------------------------------------
# Build assembly-to-assembly alignment
#
# -c        emit cg\:Z: CIGAR tag required by paftools liftover
# -x asm5   assembly-to-assembly preset
# --cs      retain cs\:Z: alignment string for diagnostics
# ------------------------------------------------------------

echo
echo "Building assembly-to-assembly alignment with minimap2..."
echo "Target: $WBCEL235_FASTA"
echo "Query:  $LIBUDA_FASTA"

minimap2 -cx asm5 --cs \
    "$WBCEL235_FASTA" \
    "$LIBUDA_FASTA" \
    > "$PAF"

if [[ ! -s "$PAF" ]]; then
    echo "ERROR: minimap2 produced an empty PAF." >&2
    exit 1
fi

if ! grep -m1 -q 'cg\:Z:' "$PAF"; then
    echo "ERROR: PAF does not contain cg\:Z: tags required by paftools liftover." >&2
    exit 1
fi

# ------------------------------------------------------------
# Prepare S1 data.
#
# IMPORTANT: only UNIQUE source intervals are lifted. Every original S1 event
# is retained in META_FILE and will be assigned the lifted coordinates of its
# source interval after liftover. This makes duplicated coordinates safe.
# ------------------------------------------------------------

echo
echo "Preparing WT crossover intervals for liftover..."

Rscript - "$S1_FILE" "$LIFT_INPUT" "$META_FILE" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
s1_path <- args[[1]]
lift_input <- args[[2]]
meta_path <- args[[3]]

library(data.table)

s1 <- fread(s1_path, data.table = TRUE)

required <- c("sample", "chromosome", "start", "stop", "germ_cell", "genotype")
missing <- setdiff(required, names(s1))
if (length(missing) > 0L) {
    stop("S1 file is missing columns: ", paste(missing, collapse = ", "))
}

s1 <- s1[
    genotype == "WT" &
    germ_cell %in% c("oocyte", "spermatocyte")
]

if (nrow(s1) == 0L) {
    stop("No WT oocyte/spermatocyte crossover rows found in S1.")
}

# Libuda assembly FASTA names corresponding to the S1 chromosome labels.
source_chrom_map <- c(
    I = "CP116366.1",
    II = "CP116367.1",
    III = "CP116368.1",
    IV = "CP116369.1",
    V = "CP116370.1",
    X = "CP116371.1"
)

s1[, source_chrom := unname(source_chrom_map[as.character(chromosome)])]
s1[, bed_start := as.integer(start) - 1L]
s1[, bed_end := as.integer(stop)]

bad <- s1[
    is.na(source_chrom) |
    is.na(bed_start) |
    is.na(bed_end) |
    bed_start < 0L |
    bed_end <= bed_start
]

if (nrow(bad) > 0L) {
    stop("Found invalid/unrecognised crossover intervals: ", nrow(bad))
}

# This is exactly the name paftools.js writes into its lifted BED output.
s1[, source_key := paste(source_chrom, bed_start, bed_end, sep = "_")]

# Stable event IDs for QC only.
s1[, event_id := sprintf("event_%07d", seq_len(.N))]

meta <- s1[, .(
    event_id,
    sample_id = as.character(sample),
    cell_type = as.character(germ_cell),
    chromosome = as.character(chromosome),
    source_chrom,
    original_start = as.integer(start),
    original_stop = as.integer(stop),
    bed_start,
    bed_end,
    source_key
)]

# One row per original S1 event.
fwrite(meta, meta_path, sep = "\t", quote = FALSE)

# One row per UNIQUE source interval for liftover.
unique_bed <- unique(s1[, .(source_chrom, bed_start, bed_end, source_key)])
unique_bed <- unique_bed[order(source_chrom, bed_start, bed_end)]

# BED4: chrom, start, end, source_key. paftools will replace the name with
# the source coordinate name, which is exactly source_key.
fwrite(
    unique_bed[, .(source_chrom, bed_start, bed_end, source_key)],
    lift_input,
    sep = "\t",
    col.names = FALSE,
    quote = FALSE
)

cat("Input WT crossover events:", nrow(meta), "\n")
cat("Unique source intervals to lift:", nrow(unique_bed), "\n")
cat("Unique oocyte samples represented in S1:", uniqueN(meta[cell_type == "oocyte", sample_id]), "\n")
cat("Unique spermatocyte samples represented in S1:", uniqueN(meta[cell_type == "spermatocyte", sample_id]), "\n")
RSCRIPT

# ------------------------------------------------------------
# Sort PAF for paftools
# ------------------------------------------------------------

echo
echo "Sorting assembly PAF for paftools.js..."
sort -k6,6 -k8,8n "$PAF" > "$PAF_SORTED"

# ------------------------------------------------------------
# Lift unique source intervals
# ------------------------------------------------------------

echo
echo "Lifting unique crossover intervals to WBcel235 with paftools.js..."

if ! "$K8" "$PAFTOOLS" liftover \
    "$PAF_SORTED" \
    "$LIFT_INPUT" \
    > "$LIFT_OUTPUT" \
    2> "$WORKDIR/paftools_liftover.err"; then

    echo "ERROR: paftools liftover failed." >&2
    cat "$WORKDIR/paftools_liftover.err" >&2 || true
    exit 1
fi

if [[ ! -s "$LIFT_OUTPUT" ]]; then
    echo "ERROR: liftover produced no mapped intervals." >&2
    cat "$WORKDIR/paftools_liftover.err" >&2 || true
    exit 1
fi

# ------------------------------------------------------------
# Parse lifted BED and restore all S1 events
# ------------------------------------------------------------

echo
echo "Inspecting lifted BED format..."
head -5 "$LIFT_OUTPUT"

echo
echo "Mapping lifted coordinates back to all S1 crossover events..."

Rscript - "$LIFT_OUTPUT" "$META_FILE" "$OOCYTE_MIDPOINTS" "$SPERM_MIDPOINTS" "$UNMAPPED" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
lifted_path <- args[[1]]
meta_path <- args[[2]]
oocyte_out <- args[[3]]
sperm_out <- args[[4]]
unmapped_path <- args[[5]]

library(data.table)

# paftools output is standard BED6:
# target_chrom, target_start, target_end, source_key, score, strand
lifted <- fread(
    lifted_path,
    header = FALSE,
    fill = TRUE,
    data.table = TRUE
)

if (ncol(lifted) < 6L) {
    stop("Expected 6 columns in paftools lifted BED; found ", ncol(lifted))
}

lifted <- lifted[, 1:6]
setnames(
    lifted,
    c("target_chrom", "target_start", "target_end", "lifted_name", "score", "strand")
)

lifted[, target_start := as.integer(target_start)]
lifted[, target_end := as.integer(target_end)]
lifted[, lifted_name := as.character(lifted_name)]

lifted <- lifted[
    !is.na(target_chrom) &
    !is.na(target_start) &
    !is.na(target_end) &
    !is.na(lifted_name) &
    target_end > target_start
]

meta <- fread(meta_path, data.table = TRUE)

required_meta <- c(
    "event_id", "sample_id", "cell_type", "source_key",
    "source_chrom", "bed_start", "bed_end"
)
missing_meta <- setdiff(required_meta, names(meta))
if (length(missing_meta) > 0L) {
    stop("Metadata file missing columns: ", paste(missing_meta, collapse = ", "))
}

# paftools may append _t5 or _t3 to names for partial mappings.
lifted[, source_key := sub("_(t5|t3)$", "", lifted_name)]
lifted[, partial_mapping := grepl("_(t5|t3)$", lifted_name)]

# We lifted unique source intervals, so every clean source key should map to
# exactly one target interval. Anything mapping zero or multiple times is not
# safe to use.
lift_counts <- lifted[, .(n_lifted = .N), by = source_key]
meta_counts <- unique(meta[, .(source_key)])

count_check <- merge(
    meta_counts,
    lift_counts,
    by = "source_key",
    all.x = TRUE
)
count_check[is.na(n_lifted), n_lifted := 0L]
count_check[, ambiguous := n_lifted != 1L]

ambiguous_keys <- count_check[count_check$ambiguous == TRUE, source_key]

cat("Original S1 events:", nrow(meta), "\n")
cat("Unique source intervals:", nrow(meta_counts), "\n")
cat("Lifted unique intervals:", nrow(lifted), "\n")
cat("Ambiguous/unmapped source keys:", length(ambiguous_keys), "\n")

if (length(ambiguous_keys) > 0L) {
    fwrite(
        count_check[count_check$ambiguous == TRUE],
        paste0(unmapped_path, ".ambiguous.tsv"),
        sep = "\t",
        quote = FALSE
    )
}

# Only fully mapped, non-partial, one-to-one source intervals.
good_lifted <- lifted[
    !lifted$source_key %in% ambiguous_keys &
    lifted$partial_mapping == FALSE
]

# There can be multiple original S1 events at the same source interval.
# That is intentional: a single lifted coordinate is inherited by each event.
joined <- merge(
    meta,
    good_lifted[, .(source_key, target_chrom, target_start, target_end)],
    by = "source_key",
    all.x = TRUE,
    sort = FALSE
)

matched <- !is.na(joined$target_chrom) &
           !is.na(joined$target_start) &
           !is.na(joined$target_end)

cat("Matched S1 events:", sum(matched), "\n")
cat("Unmatched S1 events:", sum(!matched), "\n")

unmatched <- joined[matched == FALSE]
if (nrow(unmatched) > 0L) {
    fwrite(
        unmatched[, .(event_id, sample_id, cell_type, source_key)],
        unmapped_path,
        sep = "\t",
        quote = FALSE
    )
}

joined <- joined[matched == TRUE]

if (nrow(joined) == 0L) {
    stop("No S1 events survived the source-to-WBcel235 mapping.")
}

# Midpoint of the lifted WBcel235 interval.
joined[, midpoint := floor((target_start + target_end) / 2)]

valid_chr <- c("I", "II", "III", "IV", "V", "X")
joined <- joined[as.character(target_chrom) %in% valid_chr]

# Store one line per original S1 event so duplicates are retained correctly.
oocyte <- joined[
    cell_type == "oocyte",
    .(
        chrom = as.character(target_chrom),
        midpoint = as.integer(midpoint),
        event_id = as.character(event_id)
    )
]

sperm <- joined[
    cell_type == "spermatocyte",
    .(
        chrom = as.character(target_chrom),
        midpoint = as.integer(midpoint),
        event_id = as.character(event_id)
    )
]

fwrite(oocyte, oocyte_out, sep = "\t", col.names = FALSE, quote = FALSE)
fwrite(sperm, sperm_out, sep = "\t", col.names = FALSE, quote = FALSE)

cat("Final mapped oocyte events:", nrow(oocyte), "\n")
cat("Final mapped spermatocyte events:", nrow(sperm), "\n")
RSCRIPT

# ------------------------------------------------------------
# Calculate NON-OVERLAPPING recombination bins
#
# Purpose:
#   The downstream pipeline intersects signal tracks with gene bodies
#   and takes the mean signal per gene. Therefore the recombination
#   track must not contain overlapping windows: overlapping 200-kb
#   sliding windows would cause the same genomic region to be counted
#   multiple times and create artificial intermediate values.
#
# We therefore use non-overlapping 200-kb bins directly from the
# underlying crossover midpoint events. This keeps the same spatial
# scale as the paper's 200-kb windows, while producing a conventional
# non-overlapping bedGraph suitable for gene-level mean summaries.
# ------------------------------------------------------------

echo
echo "Calculating NON-OVERLAPPING 200-kb recombination bins..."

Rscript - "$OOCYTE_MIDPOINTS" "$SPERM_MIDPOINTS" "$OOCYTE_OUT" "$SPERM_OUT" "$QC_OUT" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
oocyte_path <- args[[1]]
sperm_path <- args[[2]]
oocyte_out <- args[[3]]
sperm_out <- args[[4]]
qc_out <- args[[5]]

library(data.table)

# WBcel235 chromosome lengths.
chrom_lengths <- c(
    I = 15072434L,
    II = 15279421L,
    III = 13783801L,
    IV = 17493829L,
    V = 20924180L,
    X = 17718942L
)

# Use the same 200-kb spatial scale as the source analysis, but make
# the bins non-overlapping so they can safely be averaged over genes.
BIN_SIZE <- 200000L

# One crossover per chromosome corresponds to a 50 cM genetic map length
# in C. elegans. Convert crossover counts per successfully called meiosis
# into cM/Mb using the actual bin width in Mb.
N_OOCYTE_MEIOSES <- 297L
N_SPERMATOCYTE_MEIOSES <- 300L

read_events <- function(path) {
    if (!file.exists(path) || file.info(path)$size == 0) {
        return(data.table(chrom = character(), midpoint = integer(), event_id = character()))
    }

    x <- fread(path, header = FALSE, fill = TRUE, data.table = TRUE)
    if (ncol(x) < 3L) {
        stop("Unexpected crossover midpoint format in ", path)
    }

    x <- x[, 1:3]
    setnames(x, c("chrom", "midpoint", "event_id"))
    x[, chrom := as.character(chrom)]
    x[, midpoint := as.integer(midpoint)]
    x[, event_id := as.character(event_id)]

    x[
        chrom %in% names(chrom_lengths) &
        !is.na(midpoint) &
        midpoint >= 0L
    ]
}

make_bins <- function(chr, chr_len) {
    starts <- seq.int(0L, chr_len - 1L, by = BIN_SIZE)
    ends <- pmin(starts + BIN_SIZE, chr_len)

    data.table(
        chrom = chr,
        start = starts,
        end = ends
    )[end > start]
}

all_bins <- rbindlist(
    lapply(names(chrom_lengths), function(chr) make_bins(chr, chrom_lengths[[chr]]))
)
all_bins[, bin_id := .I]

calc_rates <- function(events, n_meioses) {
    b <- copy(all_bins)
    b[, count := 0L]

    if (nrow(events) > 0L) {
        # Each crossover midpoint belongs to exactly one non-overlapping bin.
        e <- copy(events)
        e[, bin_id := floor(midpoint / BIN_SIZE) + 1L]

        counts <- e[, .(count = .N), by = .(chrom, bin_id)]

        if (nrow(counts) > 0L) {
            b[counts, on = .(chrom, bin_id), count := i.count]
        }
    }

    # cM/Mb = (COs / successfully called meioses) * 50 cM / bin length (Mb).
    b[, width_mb := (end - start) / 1e6]
    b[, rate := (count / n_meioses) * 50 / width_mb]
    b[, chrom_order := match(chrom, names(chrom_lengths))]
    setorder(b, chrom_order, start)
    b[, chrom_order := NULL]
    b[, width_mb := NULL]
    b
}

oocyte <- read_events(oocyte_path)
sperm <- read_events(sperm_path)

oocyte_rate <- calc_rates(oocyte, N_OOCYTE_MEIOSES)
sperm_rate <- calc_rates(sperm, N_SPERMATOCYTE_MEIOSES)

oocyte_out_dt <- oocyte_rate[, .(
    chrom = paste0("chr", chrom),
    start,
    end,
    rate
)]

sperm_out_dt <- sperm_rate[, .(
    chrom = paste0("chr", chrom),
    start,
    end,
    rate
)]

fwrite(oocyte_out_dt, oocyte_out, sep = "\t", col.names = FALSE, quote = FALSE)
fwrite(sperm_out_dt, sperm_out, sep = "\t", col.names = FALSE, quote = FALSE)

qc <- data.table(
    metric = c(
        "oocyte_successfully_called_meioses",
        "spermatocyte_successfully_called_meioses",
        "oocyte_mapped_crossover_events",
        "spermatocyte_mapped_crossover_events",
        "oocyte_bins",
        "spermatocyte_bins",
        "bin_size_bp",
        "oocyte_max_cM_per_Mb",
        "spermatocyte_max_cM_per_Mb"
    ),
    value = c(
        N_OOCYTE_MEIOSES,
        N_SPERMATOCYTE_MEIOSES,
        nrow(oocyte),
        nrow(sperm),
        nrow(oocyte_out_dt),
        nrow(sperm_out_dt),
        BIN_SIZE,
        max(oocyte_out_dt$rate, na.rm = TRUE),
        max(sperm_out_dt$rate, na.rm = TRUE)
    )
)

fwrite(qc, qc_out, sep = "\t", quote = FALSE)

cat("Oocyte successfully called meioses:", N_OOCYTE_MEIOSES, "\n")
cat("Spermatocyte successfully called meioses:", N_SPERMATOCYTE_MEIOSES, "\n")
cat("Oocyte mapped crossover events:", nrow(oocyte), "\n")
cat("Spermatocyte mapped crossover events:", nrow(sperm), "\n")
cat("Oocyte non-overlapping bins:", nrow(oocyte_out_dt), "\n")
cat("Spermatocyte non-overlapping bins:", nrow(sperm_out_dt), "\n")
cat("Bin size (bp):", BIN_SIZE, "\n")
cat("Wrote:", oocyte_out, "\n")
cat("Wrote:", sperm_out, "\n")
RSCRIPT

# ------------------------------------------------------------
# Final checks
# ------------------------------------------------------------

echo
echo "Final output check:"

[[ -s "$OOCYTE_OUT" ]] || {
    echo "ERROR: oocyte bedGraph is empty." >&2
    exit 1
}

[[ -s "$SPERM_OUT" ]] || {
    echo "ERROR: spermatocyte bedGraph is empty." >&2
    exit 1
}

head -3 "$OOCYTE_OUT"
head -3 "$SPERM_OUT"

echo
echo "============================================================"
echo "Done"
echo "============================================================"
echo "Final WBcel235 tracks:"
echo "  $OOCYTE_OUT"
echo "  $SPERM_OUT"
echo "QC:"
echo "  $QC_OUT"
echo "Unmapped / ambiguous intervals:"
echo "  $UNMAPPED"
echo "============================================================"

