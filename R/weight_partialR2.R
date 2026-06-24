#!/usr/bin/env Rscript
###############################################################################
# weight_partialR2.R  --  TRACE-CpG mode C
#
# Re-weight / re-rank an EWAS result by how well the peripheral tissue tracks
# brain at each CpG, using the partial R2 in Data S2 (all CpGs):
#
#     weighted_t = stat * sqrt(Partial_R2_<tissue>)       # sign preserved
#
# The weight sqrt(partial R2) is always positive, so weighted_t keeps the sign
# (direction) of the input statistic. Rows are ranked by |weighted_t|
# (strongest first) but the value returned is signed.
#
# Usage:
#   Rscript R/weight_partialR2.R \
#       --ewas   my_ewas.csv \
#       --s2     Data_S2_TRACE_partialR2_ALL_CpGs.csv.gz \
#       --tissue Blood \
#       --stat   t \
#       --out    my_ewas_weighted.csv
#
# Inputs
#   --ewas   CSV with 'CpG_ID' and the statistic column named by --stat
#   --s2     TRACE Data S2 (wide; .csv or .csv.gz)
#   --tissue one of Blood | Buccal | Saliva
#   --stat   a SIGNED statistic that carries direction (e.g. t, z, logFC)
#   --out    output CSV path
#
# Output: the input EWAS columns unchanged, plus Partial_R2_<tissue> and the
#         signed 'weighted_t', sorted by descending |weighted_t|. No new p-value
#         is computed (this is a re-ranking, not a test); keep using your EWAS
#         p-values for significance. CpGs without a partial R2 in S2 get NA
#         weighted_t and sink to the bottom (each CpG handled independently).
#
# Requires: data.table.
###############################################################################

suppressPackageStartupMessages(library(data.table))

parse_args <- function(args) {
  out <- list(); i <- 1
  while (i <= length(args)) {
    key <- args[i]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    if (i == length(args) || startsWith(args[i + 1], "--"))
      stop("Missing value for ", key)
    out[[sub("^--", "", key)]] <- args[i + 1]; i <- i + 2
  }
  out
}
a <- parse_args(commandArgs(trailingOnly = TRUE))
req <- c("ewas", "s2", "tissue", "stat", "out")
miss <- setdiff(req, names(a))
if (length(miss))
  stop("Missing required argument(s): ", paste0("--", miss, collapse = " "))
if (!a$tissue %in% c("Blood", "Buccal", "Saliva"))
  stop("--tissue must be Blood, Buccal, or Saliva (got: ", a$tissue, ")")

pr2_col <- paste0("Partial_R2_", a$tissue)

ewas <- fread(a$ewas, showProgress = FALSE)
if (!"CpG_ID" %in% names(ewas)) stop("--ewas must contain a 'CpG_ID' column.")
if (!a$stat %in% names(ewas))
  stop("--stat column '", a$stat, "' not found in --ewas.")

s2 <- fread(a$s2, showProgress = FALSE, select = c("CpG_ID", pr2_col))

out <- merge(ewas, s2[, .(CpG_ID, Partial_R2 = get(pr2_col))],
             by = "CpG_ID", all.x = TRUE, sort = FALSE)
setnames(out, "Partial_R2", pr2_col)

stat_vec <- as.numeric(out[[a$stat]])
out[, weighted_t := stat_vec * sqrt(get(pr2_col))]          # signed weight
out[, .abs_wt := abs(weighted_t)]
setorder(out, -.abs_wt, na.last = TRUE)
out[, .abs_wt := NULL]

fwrite(out, a$out)
n_w <- sum(!is.na(out$weighted_t))
cat(sprintf("[weight_partialR2] %s: weighted %d of %d CpGs (weight = |%s| * sqrt(%s)) -> %s\n",
            a$tissue, n_w, nrow(out), a$stat, pr2_col, a$out))
