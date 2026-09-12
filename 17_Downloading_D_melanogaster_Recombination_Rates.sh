#!/usr/bin/env bash
set -euo pipefail

#SBATCH --job-name=17_Downloading_C_elegans_Replication_Timing_and_liftOver
#SBATCH --output=/home/ocdm0351/DPhil/logs/%x_%A.log
#SBATCH --error=/home/ocdm0351/DPhil/logs/%x_%A.err

module load deepTools/3.5.2-foss-2022a
module load Anaconda3
source activate pandoc_env

URL='https://comeron.lab.uiowa.edu/sites/comeron.lab.uiowa.edu/files/2023-01/Dmelanogaster-Crossover%20Rates-All-Chromosomes-100kb.R5%20%26%20R6.xlsx'
OUTDIR="${1:-/home/ocdm0351/DPhil/R_Data}"
XLSX="$OUTDIR/Dmelanogaster-Crossover-Rates-R6.xlsx"
BEDGRAPH="$OUTDIR/D_melanogaster_female_Recombination_Rates.bedgraph"

mkdir -p "$OUTDIR"

command -v curl >/dev/null 2>&1 || { echo 'ERROR: curl is required.' >&2; exit 1; }
command -v Rscript >/dev/null 2>&1 || { echo 'ERROR: Rscript is required (R installation).' >&2; exit 1; }

curl -L --fail --retry 3 -o "$XLSX" "$URL"

Rscript - "$XLSX" "$BEDGRAPH" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
xlsx <- args[1]
out <- args[2]

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop(
    "The R package 'readxl' is required. Install it once with:\n",
    "  install.packages('readxl')",
    call. = FALSE
  )
}

library(readxl)

sheets <- excel_sheets(xlsx)
sheet <- sheets[tolower(trimws(sheets)) == "r6 recomb melanogaster"]
if (length(sheet) != 1) {
  stop("Could not find the R6 worksheet. Found: ", paste(sheets, collapse = ", "))
}

x <- read_excel(xlsx, sheet = sheet, col_names = FALSE, .name_repair = "minimal")

# Locate the header row.
header_row <- which(apply(x, 1, function(z) {
  vals <- trimws(ifelse(is.na(z), "", as.character(z)))
  any(grepl("Midpoint 100kb window", vals, fixed = TRUE)) &&
    any(grepl("c (cM/Mb/female meiosis)", vals, fixed = TRUE), na.rm = TRUE)
}))[1]

if (is.na(header_row)) stop("Could not locate the R6 data header row.")

header <- trimws(ifelse(is.na(x[header_row, ]), "", as.character(x[header_row, ])))

chrom_idx <- which(header == "Chromosome /Arm")
blocks <- lapply(chrom_idx, function(i) {
  mids <- which(grepl("^Midpoint 100kb window", header))
  rates <- which(grepl("^c \\(cM/Mb/female meiosis\\)$", header))
  mid <- mids[mids > i][1]
  rate <- rates[rates > i][1]
  if (is.na(mid) || is.na(rate)) return(NULL)
  c(chrom = i, midpoint = mid, rate = rate)
})
blocks <- Filter(Negate(is.null), blocks)

if (!length(blocks)) stop("Could not identify chromosome/rate blocks.")

out_rows <- list()
k <- 0L
if (header_row < nrow(x)) {
  for (b in blocks) {
    ci <- as.integer(b[["chrom"]])
    mi <- as.integer(b[["midpoint"]])
    ri <- as.integer(b[["rate"]])

    chrom <- trimws(as.character(x[(header_row + 1):nrow(x), ci][[1]]))
    midpoint <- suppressWarnings(as.numeric(x[(header_row + 1):nrow(x), mi][[1]]))
    rate <- suppressWarnings(as.numeric(x[(header_row + 1):nrow(x), ri][[1]]))

    keep <- !is.na(chrom) & chrom != "" & !is.na(midpoint) & !is.na(rate)
    if (any(keep)) {
      d <- data.frame(chrom = chrom[keep], midpoint = midpoint[keep], rate = rate[keep], stringsAsFactors = FALSE)
      d$chrom <- sub("^Chr", "", d$chrom, ignore.case = TRUE)
      k <- k + 1L
      out_rows[[k]] <- d
    }
  }
}

d <- do.call(rbind, out_rows)
d <- d[!duplicated(d[c("chrom", "midpoint", "rate")]), ]

chrom_order <- c("X", "2L", "2R", "3L", "3R", "4")
d$chrom_rank <- match(d$chrom, chrom_order)
d$chrom_rank[is.na(d$chrom_rank)] <- 999L
d <- d[order(d$chrom_rank, d$chrom, d$midpoint), ]

# Convert midpoint measurements to non-overlapping intervals.
# Boundaries are halfway between adjacent R6 midpoints. This preserves the
# coordinate-lifted R6 map without creating overlapping 100-kb windows.
lines <- character(0)
for (ch in unique(d$chrom)) {
  z <- d[d$chrom == ch, c("midpoint", "rate")]
  z <- z[order(z$midpoint), ]
  if (nrow(z) == 1L) {
    start <- max(0, round(z$midpoint[1] - 50000))
    end <- round(z$midpoint[1] + 50000)
    lines <- c(lines, paste(ch, start, end, format(z$rate[1], trim = TRUE, scientific = FALSE), sep = "\t"))
    next
  }

  bounds <- numeric(nrow(z) + 1L)
  bounds[1] <- max(0, round(z$midpoint[1] - (z$midpoint[2] - z$midpoint[1]) / 2))
  bounds[nrow(z) + 1L] <- round(z$midpoint[nrow(z)] + (z$midpoint[nrow(z)] - z$midpoint[nrow(z) - 1L]) / 2)
  if (nrow(z) > 1L) {
    bounds[2:nrow(z)] <- round((z$midpoint[-nrow(z)] + z$midpoint[-1]) / 2)
  }

  for (i in seq_len(nrow(z))) {
    if (bounds[i + 1L] > bounds[i]) {
      lines <- c(lines, paste(ch, bounds[i], bounds[i + 1L], format(z$rate[i], trim = TRUE, scientific = FALSE), sep = "\t"))
    }
  }
}

writeLines(lines, out)
cat("Wrote", out, "\n")
RSCRIPT

printf 'Done. Output: %s\n' "$BEDGRAPH"
