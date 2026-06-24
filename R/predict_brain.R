#!/usr/bin/env Rscript
###############################################################################
# predict_brain.R  --  TRACE-CpG mode B
#
# Predict brain DNA methylation (beta values) from a peripheral tissue, using
# the TRACE BESS coefficients in Data S1 and the fixed variance-rescaling
# coefficient (Scale).
#
# Model (M-value scale), fit per CpG with brain ~ peripheral + age + sex + REs:
#     pred = Coef_Intercept
#          + Coef_Peripheral_M * peripheral_M
#          + Coef_Age          * age
#          + Coef_Sex          * sex            (sex: Male = 1, Female = 0)
#
# Variance rescaling (fixed from the discovery cohort):
#     corrected = rowMean(pred) + (pred - rowMean(pred)) * Scale
#     (Scale = sqrt(var_target / var_pred_discovery); CpGs absent from the
#      coefficient table are left uncorrected, i.e. Scale = 1.)
#
#   The rescaling spreads each CpG's predicted values around their mean ACROSS
#   THE INPUT SAMPLES. The Scale factor is fixed (from the discovery cohort),
#   but the centre is the mean of the samples you provide. Consequences:
#     * This is a cohort-level correction: provide multiple samples (tens or
#       more recommended) for a stable result.
#     * Corrected absolute values depend slightly on the set of samples
#       supplied, so the same individual may differ a little between, e.g., a
#       10-sample run and a 100-sample run. The spread (relative ordering and
#       variance) is unaffected; only the absolute level shifts, and it
#       stabilises as the sample size grows.
#     * With a SINGLE sample the correction is inactive (the value equals its
#       own mean, so corrected == pred). Single-sample input is not recommended.
#
# Output: predicted BRAIN beta values (corrected), one row per BESS CpG that is
#         present in the input beta matrix.
#
# Usage:
#   Rscript R/predict_brain.R \
#       --tissue Blood \
#       --beta   peripheral_beta.csv \
#       --pheno  phenotype.csv \
#       --s1     Data_S1_TRACE_model_coefficients.csv.gz \
#       --out    predicted_brain_beta.csv
#
# Inputs
#   --beta   CSV: first column = CpG_ID, remaining columns = samples (beta in [0,1])
#   --pheno  CSV: first column = sample ID (matching --beta column names),
#                 plus columns 'Age' (numeric) and 'Sex'
#                 ('M'/'F' or 'Male'/'Female', case-insensitive)
#   --s1     TRACE Data S1 (confirmed BESS coefficients; .csv or .csv.gz)
#   --tissue one of Blood | Buccal | Saliva  (selects the rows of S1)
#   --out    output CSV path
#
# Notes
#   - Age is required; any missing Age stops the run.
#   - The variance rescaling is a cohort-level correction centred on the mean of
#     the input samples; use multiple samples (single-sample input leaves the
#     correction inactive). See "Variance rescaling" above.
#   - Each CpG is predicted independently. BESS CpGs absent from the input beta
#     matrix are omitted from the output (they do not affect other CpGs).
#   - Requires: data.table, lumi (Bioconductor).
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(lumi)
})

## ---- minimal --key value argument parser ----------------------------------
parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    if (i == length(args) || startsWith(args[i + 1], "--"))
      stop("Missing value for ", key)
    out[[sub("^--", "", key)]] <- args[i + 1]
    i <- i + 2
  }
  out
}
a <- parse_args(commandArgs(trailingOnly = TRUE))
req <- c("tissue", "beta", "pheno", "s1", "out")
miss <- setdiff(req, names(a))
if (length(miss))
  stop("Missing required argument(s): ", paste0("--", miss, collapse = " "))
if (!a$tissue %in% c("Blood", "Buccal", "Saliva"))
  stop("--tissue must be Blood, Buccal, or Saliva (got: ", a$tissue, ")")

## ---- 1. coefficients (Data S1), restricted to the chosen tissue -----------
cat(sprintf("[S1] reading %s ; tissue = %s\n", a$s1, a$tissue))
s1 <- fread(a$s1, showProgress = FALSE)
need <- c("CpG_ID", "Tissue", "Coef_Intercept", "Coef_Peripheral_M",
          "Coef_Age", "Coef_Sex", "Scale")
miss <- setdiff(need, names(s1))
if (length(miss)) stop("Data S1 is missing column(s): ", paste(miss, collapse = ", "))
coef <- s1[Tissue == a$tissue, ..need]
if (nrow(coef) == 0) stop("No rows in S1 for tissue ", a$tissue)
coef[!is.finite(Scale), Scale := 1]               # uncorrected if Scale absent
cat(sprintf("[S1] %d BESS CpGs for %s\n", nrow(coef), a$tissue))

## ---- 2. peripheral beta -> M ----------------------------------------------
cat(sprintf("[BETA] reading %s\n", a$beta))
beta <- fread(a$beta, data.table = FALSE)
rownames(beta) <- beta[[1]]; beta <- beta[, -1, drop = FALSE]
bmat <- as.matrix(beta)
if (any(bmat < 0, na.rm = TRUE) || any(bmat > 1, na.rm = TRUE))
  stop("Input --beta has values outside [0,1]; supply beta values, not M-values.")
bmat[bmat >= 1] <- 1 - 1e-6
bmat[bmat <= 0] <- 1e-6
Mmat <- lumi::beta2m(bmat)
rm(beta, bmat); gc(verbose = FALSE)

## ---- 3. phenotype: align samples, build age/sex vectors -------------------
cat(sprintf("[PHENO] reading %s\n", a$pheno))
ph <- fread(a$pheno, data.table = FALSE)
rownames(ph) <- ph[[1]]
if (!all(c("Age", "Sex") %in% names(ph)))
  stop("--pheno must contain columns 'Age' and 'Sex'.")

common <- intersect(colnames(Mmat), rownames(ph))
if (length(common) == 0)
  stop("No samples shared between --beta columns and --pheno rows.")
Mmat <- Mmat[, common, drop = FALSE]
ph   <- ph[common, , drop = FALSE]
cat(sprintf("[PHENO] %d samples used\n", length(common)))

age_vec <- suppressWarnings(as.numeric(ph$Age))
if (any(is.na(age_vec)))
  stop("Age is missing or non-numeric for ", sum(is.na(age_vec)),
       " sample(s); Age is required for every sample.")

sex_raw <- toupper(trimws(as.character(ph$Sex)))
sex_vec <- ifelse(sex_raw %in% c("M", "MALE", "1"), 1,
           ifelse(sex_raw %in% c("F", "FEMALE", "0"), 0, NA_real_))
if (any(is.na(sex_vec)))
  stop("Unrecognised Sex value(s): ",
       paste(unique(ph$Sex[is.na(sex_vec)]), collapse = ", "),
       " (use M/F or Male/Female).")

## ---- 4. predict on M scale (BESS CpGs present in the input) ---------------
present <- coef$CpG_ID[coef$CpG_ID %in% rownames(Mmat)]
cat(sprintf("[PRED] %d / %d BESS CpGs available in input; %d omitted (absent from input)\n",
            length(present), nrow(coef), nrow(coef) - length(present)))

cf <- coef[match(present, CpG_ID)]
Mp <- Mmat[present, , drop = FALSE]

pred <- cf$Coef_Intercept +
        cf$Coef_Peripheral_M * Mp +
        outer(cf$Coef_Age, age_vec) +
        outer(cf$Coef_Sex, sex_vec)
rownames(pred) <- present; colnames(pred) <- common

## ---- 5. variance rescaling (corrected) ------------------------------------
row_mean <- rowMeans(pred, na.rm = TRUE)
corrected <- row_mean + (pred - row_mean) * cf$Scale

## ---- 6. back to beta; output only the BESS CpGs present in the input -------
beta_out <- lumi::m2beta(corrected)        # rows = present, cols = common

res <- data.table(CpG_ID = present)
res <- cbind(res, as.data.table(beta_out))
fwrite(res, a$out)
cat(sprintf("[OUT] wrote %s  (%d CpGs x %d samples)\n",
            a$out, nrow(res), length(common)))
