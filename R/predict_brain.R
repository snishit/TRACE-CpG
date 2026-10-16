#!/usr/bin/env Rscript
###############################################################################
# predict_brain.R  --  TRACE-CpG
#
# Estimate brain DNA methylation (beta values) from blood or buccal epithelium
# at brain epigenome shared sites (BrESS), using the model coefficients and the
# fixed variance-rescaling factor (Scale) in Data S2.
#
# Model (M-value scale), fit per CpG on the IMAGE-CpG2 discovery data with
#     brain_M ~ peripheral_M + Age + Sex + (1|PatientID) + (1|Region)
#     pred = Coef_Intercept
#          + Coef_Peripheral_M * peripheral_M
#          + Coef_Age          * age
#          + Coef_Sex          * sex            (sex: Male = 1, Female = 0)
#
# Variance rescaling (factor fixed from the discovery cohort):
#     corrected = rowMean(pred) + (pred - rowMean(pred)) * Scale
#
#   The rescaling spreads each CpG's predicted values around their mean ACROSS
#   THE INPUT SAMPLES. Scale is fixed, but the centre is the mean of the samples
#   you provide. Consequences:
#     * This is a cohort-level correction: provide multiple samples (tens or
#       more recommended) for a stable result.
#     * Corrected absolute values depend slightly on the set of samples
#       supplied; the relative ordering and variance are unaffected, and the
#       level stabilises as the sample size grows.
#     * With a SINGLE sample the correction is inactive (corrected == pred).
#
# Output: estimated BRAIN beta values (corrected), one row per BrESS that is
#         present in the input beta matrix.
#
# Usage:
#   Rscript R/predict_brain.R \
#       --tissue Blood \
#       --beta   peripheral_beta.csv \
#       --pheno  phenotype.csv \
#       --s2     Data_S2_TRACE_BrESS_models.csv.gz \
#       --out    predicted_brain_beta.csv
#
# Inputs
#   --beta   CSV: first column = CpG_ID, remaining columns = samples (beta in [0,1])
#   --pheno  CSV: first column = sample ID (matching --beta column names),
#                 plus columns 'Age' (numeric, years) and 'Sex'
#                 ('M'/'F' or 'Male'/'Female', case-insensitive)
#   --s2     TRACE Data S2 (BrESS model parameters; .csv or .csv.gz)
#   --tissue Blood | Buccal  (selects the rows of S2)
#   --out    output CSV path
#
# Notes
#   - Age is required; any missing Age stops the run.
#   - Each CpG is predicted independently. BrESS absent from the input beta
#     matrix are omitted from the output (they do not affect other CpGs).
#   - Requires: data.table.
###############################################################################

suppressPackageStartupMessages(library(data.table))

## beta <-> M (same definitions as lumi::beta2m / lumi::m2beta)
beta2m <- function(b) log2(b / (1 - b))
m2beta <- function(m) 2^m / (2^m + 1)

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
req <- c("tissue", "beta", "pheno", "s2", "out")
miss <- setdiff(req, names(a))
if (length(miss))
  stop("Missing required argument(s): ", paste0("--", miss, collapse = " "))
if (!a$tissue %in% c("Blood", "Buccal"))
  stop("--tissue must be Blood or Buccal (got: ", a$tissue, ")")

## ---- 1. model parameters (Data S2), restricted to the chosen tissue -------
cat(sprintf("[S2] reading %s ; tissue = %s\n", a$s2, a$tissue))
s2 <- fread(a$s2, showProgress = FALSE)
need <- c("CpG_ID", "Tissue", "Coef_Intercept", "Coef_Peripheral_M",
          "Coef_Age", "Coef_Sex", "Scale")
miss <- setdiff(need, names(s2))
if (length(miss)) stop("Data S2 is missing column(s): ", paste(miss, collapse = ", "))
coef <- s2[Tissue == a$tissue, ..need]
if (nrow(coef) == 0) stop("No rows in S2 for tissue ", a$tissue)
cat(sprintf("[S2] %d BrESS for %s\n", nrow(coef), a$tissue))

## ---- 2. peripheral beta -> M ----------------------------------------------
cat(sprintf("[BETA] reading %s\n", a$beta))
beta <- fread(a$beta, data.table = FALSE)
rownames(beta) <- beta[[1]]; beta <- beta[, -1, drop = FALSE]
bmat <- as.matrix(beta)
if (any(bmat < 0, na.rm = TRUE) || any(bmat > 1, na.rm = TRUE))
  stop("Input --beta has values outside [0,1]; supply beta values, not M-values.")
bmat[bmat >= 1] <- 1 - 1e-6
bmat[bmat <= 0] <- 1e-6
Mmat <- beta2m(bmat)
rm(beta, bmat); invisible(gc(verbose = FALSE))

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

## ---- 4. predict on M scale (BrESS present in the input) -------------------
present <- coef$CpG_ID[coef$CpG_ID %in% rownames(Mmat)]
cat(sprintf("[PRED] %d / %d BrESS available in input; %d omitted (absent from input)\n",
            length(present), nrow(coef), nrow(coef) - length(present)))
if (length(present) == 0) stop("None of the BrESS for ", a$tissue, " are in the input.")

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

## ---- 6. back to beta; output only the BrESS present in the input ---------
beta_out <- m2beta(corrected)

res <- data.table(CpG_ID = present)
res <- cbind(res, as.data.table(beta_out))
fwrite(res, a$out)
cat(sprintf("[OUT] wrote %s  (%d CpGs x %d samples)\n",
            a$out, nrow(res), length(common)))
