# CDNCH-DID

R analysis code for continuous heat exposure and cardiovascular disease using an extended continuous-treatment difference-in-differences framework.

## Scripts

| Script | Purpose |
| --- | --- |
| `01_run_did.R` | Primary models M0–M5, subgroup analyses, nonlinear exposure–response relationships, attributable burden, and robustness analyses. |
| `09_run_population_sensitivity.R` | Population-weighted exposure sensitivity analysis and comparison with the original M5 model. |
| `check_population_model_setup.R` | Check eligible samples and reproduce the original mortality CDNCH coefficient. |
| `editorial/validate_attribution.R` | Independently verify nonlinear predictions and attributable burden calculations. |

## Requirements

The scripts were developed with R 4.5.2. Install the required packages:

```r
install.packages(c("data.table", "sandwich"))
```

`splines` is included with R.

## Input data and execution

Input data are not included. Place the analysis panel at `data/input_panel.csv`. Required variables and sample restrictions are defined in `01_run_did.R`.

Run the primary analysis from the repository directory:

```sh
Rscript 01_run_did.R
```

The primary script creates the `outputs/` directory and the reference results used by the supplementary scripts. The primary script also accepts an input path through the `DID_INPUT_FILE` environment variable.

After completing the primary analysis:

```sh
Rscript check_population_model_setup.R
```

The population-weighted analysis additionally requires the externally prepared `outputs/population_exposure/country_year_exposures.csv`. The upstream population-raster preprocessing workflow is not included in this R-code repository.

```sh
Rscript 09_run_population_sensitivity.R
```

The attribution validation script retains the original local project path `D:/DID`; update its path references for a different checkout before running it. It uses the original main-script initialization and refits models for independent verification.

## Model overview

Main count-outcome models use a quasi-Poisson log link, a population offset, exposure increments of one IQR, and country-clustered robust inference. The core M5 model adjusts for seven covariates and superregion-specific natural-spline time trends. Consult the scripts for the full model definitions and output conventions.

## Validation

All four R scripts passed syntax parsing with R 4.5.2 before upload. This check does not rerun the statistical analyses.
