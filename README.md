# TRACE-CpG

**TRACE-CpG** is a resource for relating peripheral-tissue DNA methylation (DNAm)
to brain DNAm. It provides, for blood, buccal epithelium, and saliva:

1. **BESS** (Brain–peripheral Epigenetic Surrogate Sites): CpGs whose peripheral
   DNAm reliably tracks brain DNAm.
2. **Prediction coefficients** to estimate brain DNAm from a peripheral sample.
3. **Partial R² values** for *every* analyzed CpG, quantifying how much of the
   brain signal each peripheral CpG captures.

Three ways to use it (the three "modes"):

| Mode | Question | Needs |
|------|----------|-------|
| **A** Filter   | "Which of my CpGs are reliable brain surrogates?"        | `is_BESS_<tissue>` (Data S2) |
| **B** Predict  | "What is the predicted brain DNAm for my samples?"       | coefficients + `Scale` (Data S1) |
| **C** Weight   | "Re-rank my EWAS by brain relevance."                    | `Partial_R2_<tissue>` (Data S2) |

---

## How BESS were defined

A two-stage procedure:

1. **Candidate selection** by leave-one-dataset-out cross-validation (LODO-CV).
2. **Confirmation** by a permutation test; CpGs with **FDR < 0.05** become BESS.

Confirmed BESS per tissue: **blood 155,658**, **buccal 151,214**, **saliva 67,955**.

---

## Data files

### `Data_S1_TRACE_model_coefficients.csv.gz` — confirmed BESS, prediction model
One row per BESS CpG × tissue (long format). Model:
`brain_M ~ peripheral_M + Age + Sex + (1|PatientID) + (1|Region)` (ML fit).

| Column | Meaning |
|--------|---------|
| `CpG_ID` | Illumina CpG identifier |
| `Tissue` | Blood / Buccal / Saliva |
| `Partial_R2` | marginal R² gained by the peripheral term (R2_Marginal − R2m_Reduced) |
| `R2_Conditional`, `R2_Marginal`, `R2m_Reduced` | mixed-model R² components |
| `Coef_Intercept` | model intercept (β₀) |
| `Coef_Peripheral_M` | peripheral-M coefficient (β_peripheral) |
| `Coef_Age` | age coefficient (β_age) |
| `Coef_Sex` | sex coefficient (β_sex); **Sex coded Male = 1, Female = 0** |
| `Var_PatientID`, `Var_Region`, `Var_Residual` | variance components |
| `Permutation_p_value`, `FDR` | permutation significance (BESS = FDR < 0.05) |
| `Scale` | fixed variance-rescaling factor (√(var_target / var_pred_discovery)) |

### `Data_S2_TRACE_partialR2_ALL_CpGs.csv.gz` — all CpGs, partial R² + BESS flag
One row per CpG (wide format), 687,873 CpGs.

| Column | Meaning |
|--------|---------|
| `CpG_ID` | Illumina CpG identifier |
| `Partial_R2_Blood`, `Partial_R2_Buccal`, `Partial_R2_Saliva` | partial R² (NA if the CpG was not modeled in that tissue) |
| `is_BESS_Blood`, `is_BESS_Buccal`, `is_BESS_Saliva` | TRUE if the CpG is a confirmed BESS (permutation FDR < 0.05) |

The partial R² in Data S1 (BESS) and Data S2 (all CpGs) come from the **same
model and definition** and agree for shared CpGs.

---

## Prediction model (mode B)

On the M-value scale, per CpG:

```
pred = Coef_Intercept
     + Coef_Peripheral_M * peripheral_M
     + Coef_Age          * age
     + Coef_Sex          * sex        # Sex: Male = 1, Female = 0
```

Then a fixed variance rescaling (parameters frozen from the discovery cohort):

```
corrected = rowMean(pred) + (pred - rowMean(pred)) * Scale
```

where `rowMean(pred)` is the mean predicted value of that CpG across the samples
being predicted, and `Scale` is the per-CpG column in Data S1 (CpGs absent from
the table are left uncorrected, i.e. Scale = 1). Predicted values are returned as
**beta** values.

The rescaling is a **cohort-level** correction: the `Scale` factor is fixed (from
the discovery cohort), but the centre is the mean of *your* input samples. As a
result, corrected absolute values depend slightly on the set of samples you
supply — the same individual can differ a little between, say, a 10-sample and a
100-sample run. Only the absolute level shifts; the relative ordering and the
variance are unaffected, and the level stabilises as the sample size grows.
Provide multiple samples (tens or more recommended). With a **single sample** the
correction is inactive (the value equals its own mean, so `corrected = pred`).

---

## Usage

All scripts take `--key value` arguments and run with `Rscript`.
Run them from the repository root (so that `R/...` paths resolve).

### Mode B — predict brain beta from a peripheral tissue
```
Rscript R/predict_brain.R \
    --tissue Blood \
    --beta   peripheral_beta.csv \
    --pheno  phenotype.csv \
    --s1     Data_S1_TRACE_model_coefficients.csv.gz \
    --out    predicted_brain_beta.csv
```
- `--beta`: CSV, first column `CpG_ID`, remaining columns are samples (beta in [0,1]).
- `--pheno`: CSV, first column = sample ID (matching the beta columns), plus
  `Age` (numeric) and `Sex` (`M`/`F` or `Male`/`Female`).
- **Age is required**; any missing Age stops the run.
- Each CpG is predicted independently; BESS CpGs absent from the input beta
  matrix are omitted from the output (they do not affect other CpGs).
- The variance rescaling is a cohort-level correction centred on the mean of
  the input samples — use multiple samples (see "Prediction model" above).

### Mode A — keep only BESS CpGs in an EWAS
```
Rscript R/filter_bess.R \
    --ewas   my_ewas.csv \
    --s2     Data_S2_TRACE_partialR2_ALL_CpGs.csv.gz \
    --tissue Blood \
    --out    my_ewas_BESS.csv
```
`--ewas` must contain a `CpG_ID` column; all other columns are preserved.

### Mode C — re-rank an EWAS by brain relevance
```
Rscript R/weight_partialR2.R \
    --ewas   my_ewas.csv \
    --s2     Data_S2_TRACE_partialR2_ALL_CpGs.csv.gz \
    --tissue Blood \
    --stat   t \
    --out    my_ewas_weighted.csv
```
Adds a signed `weighted_t = stat * sqrt(Partial_R2_<tissue>)` and sorts by
`|weighted_t|` (strongest first). Because `sqrt(partial R2)` is always positive,
`weighted_t` keeps the direction (sign) of your statistic.
`--stat` must be a **signed** quantity that carries direction (e.g. `t`, `z`,
`logFC`). No new p-value is produced — this is a re-ranking, not a test, so keep
using your EWAS p-values for significance.

Each command above is self-contained — copy it and replace the input/output file
names with your own.

---

## Quick test (with the included toy inputs)

The repository ships three small **synthetic** input files so you can check that
all three scripts run end to end. The beta values, phenotypes, and statistics are
artificial (no real individual-level data), but the `CpG_ID`s are 40 CpGs that
exist in Data S1 / Data S2 (25 BESS and 15 non-BESS), so the scripts produce
non-empty output when combined with the data files.

These tests use the same `Data_S1_*.csv.gz` / `Data_S2_*.csv.gz` as normal use.
Place the toy files and the data files in the repository root and run from there.

Toy files:
- `toy_beta.csv` — peripheral beta values (40 CpGs × 10 samples)
- `toy_pheno.csv` — phenotype (`SampleID`, `Age`, `Sex`), 10 samples
- `toy_ewas.csv` — EWAS results (`CpG_ID`, `t`, `P.Value`), 40 CpGs

### Mode B — predict
```
Rscript R/predict_brain.R \
    --tissue Blood \
    --beta   toy_beta.csv \
    --pheno  toy_pheno.csv \
    --s1     Data_S1_TRACE_model_coefficients.csv.gz \
    --out    toy_pred_out.csv
```
Expected: `[OUT] wrote toy_pred_out.csv  (25 CpGs x 10 samples)`. The output has
**25 rows** (the 25 toy BESS present in Data S1); CpGs absent from the input
are not written. All values are beta in [0,1]. (Tip: blanking an `Age` in
`toy_pheno.csv` and re-running should stop with `Age is missing ...`.)

### Mode A — filter
```
Rscript R/filter_bess.R \
    --ewas   toy_ewas.csv \
    --s2     Data_S2_TRACE_partialR2_ALL_CpGs.csv.gz \
    --tissue Blood \
    --out    toy_filtered_out.csv
```
Expected: `[filter_bess] Blood: 25 of 40 EWAS CpGs are BESS ...`. The output has
**25 rows** (BESS only) with a `Partial_R2_Blood` column added.

### Mode C — weight
```
Rscript R/weight_partialR2.R \
    --ewas   toy_ewas.csv \
    --s2     Data_S2_TRACE_partialR2_ALL_CpGs.csv.gz \
    --tissue Blood \
    --stat   t \
    --out    toy_weighted_out.csv
```
Expected: `[weight_partialR2] Blood: weighted 40 of 40 CpGs ...`. The output has
**40 rows** with a signed `weighted_t` column, sorted by `|weighted_t|`
(descending); any CpG without a partial R² gets `NA` and sinks to the bottom.

If all three run as above, the same commands work on your own data — just
replace the `toy_*` files with your own (the `.gz` data files are read as-is).

**Troubleshooting**
- `could not find function "beta2m"` → install `lumi` (see Requirements).
- `Missing required argument(s)` → check the flag spelling (`--tissue`, etc.).
- `--tissue must be Blood, Buccal, or Saliva` → use one of those three values.
- If your path contains spaces, quote the whole folder when you `cd`.

---

## Requirements

- R (≥ 4.0) with **data.table**.
- Mode B additionally requires **lumi** (Bioconductor) for beta↔M conversion:
  ```r
  install.packages("data.table")
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("lumi")
  ```

Input methylation must be **beta** values on the Illumina EPIC platform; the
scripts convert to M-values internally (epsilon clipping at 1e-6).

---

## License

- **Code** (`R/`): MIT License — see [`LICENSE`](LICENSE).
- **Data** (`Data_S1_*.csv.gz`, `Data_S2_*.csv.gz`): CC-BY 4.0 — see
  [`DATA_LICENSE`](DATA_LICENSE).

## Citation

See [`CITATION.cff`](CITATION.cff). (Article and DOI to be added on publication.)
