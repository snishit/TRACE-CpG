#!/usr/bin/env Rscript
###############################################################################
# filter_bess.R  --  TRACE-CpG mode A
#
# Restrict an EWAS result to the TRACE BESS for a given peripheral tissue,
# using the is_BESS flag in Data S2.
#
# Usage:
#   Rscript R/filter_bess.R \
#       --ewas   my_ewas.csv \
#       --s2     Data_S2_TRACE_partialR2_ALL_CpGs.csv.gz \
#       --tissue Blood \
#       --out    my_ewas_BESS.csv
#
# Inputs
#   --ewas   CSV with a 'CpG_ID' column (any other columns are kept as-is)
#   --s2     TRACE Data S2 (wide; .csv or .csv.gz)
#   --tissue one of Blood | Buccal | Saliva
#   --out    output CSV path
#
# Output: the rows of --ewas whose CpG is a BESS in the chosen tissue,
#         with the matching Partial_R2_<tissue> appended.
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
req <- c("ewas", "s2", "tissue", "out")
miss <- setdiff(req, names(a))
if (length(miss))
  stop("Missing required argument(s): ", paste0("--", miss, collapse = " "))
if (!a$tissue %in% c("Blood", "Buccal", "Saliva"))
  stop("--tissue must be Blood, Buccal, or Saliva (got: ", a$tissue, ")")

bess_col <- paste0("is_BESS_",   a$tissue)
pr2_col  <- paste0("Partial_R2_", a$tissue)

ewas <- fread(a$ewas, showProgress = FALSE)
if (!"CpG_ID" %in% names(ewas)) stop("--ewas must contain a 'CpG_ID' column.")

s2 <- fread(a$s2, showProgress = FALSE,
            select = c("CpG_ID", bess_col, pr2_col))
if (!all(c(bess_col, pr2_col) %in% names(s2)))
  stop("Data S2 is missing ", bess_col, " / ", pr2_col)

bess_ids <- s2[get(bess_col) == TRUE, CpG_ID]
out <- merge(ewas, s2[, .(CpG_ID, Partial_R2 = get(pr2_col))],
             by = "CpG_ID", all.x = FALSE, sort = FALSE)
out <- out[CpG_ID %in% bess_ids]
setnames(out, "Partial_R2", pr2_col)

fwrite(out, a$out)
cat(sprintf("[filter_bess] %s: %d of %d EWAS CpGs are BESS -> %s\n",
            a$tissue, nrow(out), nrow(ewas), a$out))
