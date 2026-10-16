# TRACE-CpG

**TRACE-CpG** (Trans-tissue Regression Approximation of Central Epigenome from
Correlated peripheral & Gray-matter) relates peripheral-tissue DNA methylation
(DNAm) to brain DNAm. For blood and buccal epithelium it provides:

1. **BrESS** (Brain Epigenome Shared Sites): CpGs at which peripheral DNAm
   reproducibly captures inter-individual variation in brain DNAm.
2. **Prediction models** to estimate brain DNAm from a peripheral sample at BrESS.

A web version for users who prefer not to run R is available at
https://epigen.shinozakilab.net/.

---

## How BrESS were defined

Models were built on IMAGE-CpG2, paired brain and peripheral EPIC-array data from
three datasets (DB6, DB7, DB8). For each CpG, brain DNAm was predicted from a model
with age and sex, with and without peripheral DNAm, in leave-one-dataset-out
cross-validation (train on two datasets, test on the third). **R²gain** is the
increase in out-of-sample R² for brain DNAm when peripheral DNAm is added.

A CpG is a **BrESS** when R²gain > 0.2 in at least two of the three cycles.
This was evaluated for 687,873 CpGs shared between EPIC v1 and EPIC v2.

| | BrESS |
|---|---|
| Blood | 121,167 |
| Buccal epithelium | 109,385 |
| Both (common BrESS) | 76,208 |

---

## Data files

### `Data_S1_LODO_CV_all_CpGs.csv.gz` — cross-validation results, all CpGs
One row per CpG (wide format), 687,873 CpGs.

| Column | Meaning |
|--------|---------|
| `CpG_ID` | Illumina CpG identifier |
| `Blood_R2gain_DB8`, `Blood_R2gain_DB7`, `Blood_R2gain_DB6` | R²gain in the cycle tested on the named dataset (trained on the other two) |
| `Blood_n_cycles_R2gain_gt0.2` | number of cycles with R²gain > 0.2 (0–3) |
| `Blood_BrESS` | TRUE if R²gain > 0.2 in at least two cycles |
| `Buccal_…` | same columns for buccal epithelium |
| `Common_BrESS` | TRUE if the CpG is a BrESS in both blood and buccal epithelium |

### `Data_S2_TRACE_BrESS_models.csv.gz` — prediction models, BrESS only
One row per BrESS × tissue (long format); blood 121,167 rows, buccal 109,385 rows.
Model, fit on the full IMAGE-CpG2 data (DB6–DB8) by maximum likelihood:
`brain_M ~ peripheral_M + Age + Sex + (1|PatientID) + (1|Region)`.

| Column | Meaning |
|--------|---------|
| `CpG_ID` | Illumina CpG identifier |
| `Tissue` | Blood / Buccal |
| `R2_Marginal`, `R2_Conditional` | marginal and conditional R² of the full model |
| `R2m_Reduced` | marginal R² of the reduced model (age and sex only) |
| `Coef_Intercept` | intercept (β₀) |
| `Coef_Peripheral_M` | peripheral-M coefficient |
| `Coef_Age` | age coefficient (years) |
| `Coef_Sex` | sex coefficient; **Sex coded Male = 1, Female = 0** |
| `Var_PatientID`, `Var_Region`, `Var_Residual` | random-effect variance components |
| `Scale` | fixed variance-rescaling factor, √(σ²_brain / σ²_pBrain) from the discovery data |

---

## Prediction model

On the M-value scale, per CpG:

```
pred = Coef_Intercept
     + Coef_Peripheral_M * peripheral_M
     + Coef_Age          * age
     + Coef_Sex          * sex        # Sex: Male = 1, Female = 0
```

Because a linear prediction has smaller inter-individual variance than measured
brain DNAm, a fixed variance rescaling is then applied:

```
corrected = rowMean(pred) + (pred - rowMean(pred)) * Scale
```

where `rowMean(pred)` is the mean predicted value of that CpG across the samples
being predicted. Predicted values are returned as **beta** values.

The rescaling is a **cohort-level** correction: `Scale` is fixed from the discovery
data, but the centre is the mean of *your* input samples. Corrected absolute values
therefore depend slightly on the set of samples supplied; the relative ordering and
the variance are unaffected, and the level stabilises as the sample size grows.
Provide multiple samples (tens or more recommended). With a **single sample** the
correction is inactive (`corrected = pred`).

Predicted values are statistical estimates, not measurements of brain DNAm, and
are provided only for BrESS. Because age and sex are predictors, testing the
predicted values for association with age or sex is circular and should be avoided.

---

## Usage

Run from the repository root (so that `R/...` paths resolve).

```
Rscript R/predict_brain.R \
    --tissue Blood \
    --beta   peripheral_beta.csv \
    --pheno  phenotype.csv \
    --s2     Data_S2_TRACE_BrESS_models.csv.gz \
    --out    predicted_brain_beta.csv
```
- `--tissue`: `Blood` or `Buccal`.
- `--beta`: CSV, first column `CpG_ID`, remaining columns are samples (beta in [0,1]).
  Input should be EPIC v1 or EPIC v2 data, background-corrected and normalized
  from IDAT files, with poor-quality samples and probes removed. For EPIC v2,
  collapse replicate (suffix) probes to their v1-equivalent CpG identifiers.
- `--pheno`: CSV, first column = sample ID (matching the beta columns), plus
  `Age` (numeric, years) and `Sex` (`M`/`F` or `Male`/`Female`).
- **Age is required**; any missing Age stops the run.
- Each CpG is predicted independently; BrESS absent from the input are omitted
  from the output.

### Checking which of your CpGs are BrESS

No script is needed; match your CpG IDs against Data S1. For example, in R:

```r
library(data.table)
s1   <- fread("Data_S1_LODO_CV_all_CpGs.csv.gz")
ewas <- fread("my_ewas.csv")                      # a table with a CpG_ID column
ewas_bress <- ewas[CpG_ID %in% s1[Blood_BrESS == TRUE, CpG_ID]]
```

---

## Quick test (with the included toy inputs)

Two small **synthetic** input files are included. The beta values and phenotypes
are artificial (no real individual-level data), but the `CpG_ID`s are 40 real CpGs,
21 of which are blood BrESS and 18 buccal BrESS.

- `toy_beta.csv` — peripheral beta values (40 CpGs × 10 samples)
- `toy_pheno.csv` — phenotype (`SampleID`, `Age`, `Sex`), 10 samples

```
Rscript R/predict_brain.R \
    --tissue Blood \
    --beta   toy_beta.csv \
    --pheno  toy_pheno.csv \
    --s2     Data_S2_TRACE_BrESS_models.csv.gz \
    --out    toy_pred_out.csv
```
Expected: `[OUT] wrote toy_pred_out.csv  (21 CpGs x 10 samples)`, all values beta
in [0,1]. With `--tissue Buccal` the output has 18 rows.

**Troubleshooting**
- `Missing required argument(s)` → check the flag spelling (`--tissue`, etc.).
- `--tissue must be Blood or Buccal` → use one of those two values.
- If your path contains spaces, quote the whole folder when you `cd`.

---

## Requirements

R (≥ 4.0) with **data.table**:
```r
install.packages("data.table")
```

---

## License

- **Code** (`R/`): MIT License — see [`LICENSE`](LICENSE).
- **Data** (`Data_S1_*.csv.gz`, `Data_S2_*.csv.gz`): CC-BY 4.0 — see
  [`DATA_LICENSE`](DATA_LICENSE).

## Citation

See [`CITATION.cff`](CITATION.cff). (Article and DOI to be added on publication.)
