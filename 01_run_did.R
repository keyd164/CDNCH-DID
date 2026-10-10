#!/usr/bin/env Rscript

# Continuous heat exposure and CVD: revised continuous-treatment DID analysis
# Version: V1.0 | Author: dy | Implementation updated: 2026-09-19

required_packages <- c("data.table", "splines", "sandwich")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required R packages: ", paste(missing_packages, collapse = ", "),
    ". Install them before running the analysis.", call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(splines)
  library(sandwich)
})

options(encoding = "UTF-8", stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg) == 1L) sub("^--file=", "", script_arg) else "01_run_did.R"
project_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
if (!dir.exists(project_dir)) project_dir <- getwd()

input_default <- file.path(project_dir, "data", "input_panel.csv")
input_file <- Sys.getenv("DID_INPUT_FILE", unset = input_default)
output_dir <- file.path(project_dir, "outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) stop("Input file not found: ", input_file, call. = FALSE)

heat_vars <- c(
  "absolute_cumulative_heat", "relative_cumulative_heat",
  "mixed_cumulative_heat", "daytime_cumulative_heat",
  "nighttime_cumulative_heat", "compound_day_night_cumulative_heat",
  "annual_mean_humidex", "annual_mean_tmean"
)

heat_labels <- c(
  absolute_cumulative_heat = "绝对阈值累积热",
  relative_cumulative_heat = "相对阈值累积热",
  mixed_cumulative_heat = "混合阈值累积热",
  daytime_cumulative_heat = "昼间累积热",
  nighttime_cumulative_heat = "夜间累积热",
  compound_day_night_cumulative_heat = "昼夜复合累积热",
  annual_mean_humidex = "年均 humidex",
  annual_mean_tmean = "年均平均温度"
)

heat_units <- c(
  absolute_cumulative_heat = "°C·days",
  relative_cumulative_heat = "°C·days",
  mixed_cumulative_heat = "°C·days",
  daytime_cumulative_heat = "°C·days",
  nighttime_cumulative_heat = "°C·days",
  compound_day_night_cumulative_heat = "°C·days",
  annual_mean_humidex = "humidex",
  annual_mean_tmean = "°C"
)

outcome_labels <- c(`1` = "死亡", `2` = "DALYs", `5` = "患病", `6` = "发病")
outcome_burden_units <- c(
  `1` = "Deaths", `2` = "DALYs", `5` = "Prevalent cases", `6` = "Incident cases"
)
expected_measures <- c(1L, 2L, 5L, 6L)
model_ids <- c("M0", "M1", "M2", "M3a", "M3b", "M4", "M5")
core_model_id <- "M5"
nonlinear_model_ids <- core_model_id

covariate_vars <- c(
  "mean_value", "SEV_agestd_air_hap", "SEV_agestd_air_pm",
  "SEV_agestd_diet_fruit", "SEV_agestd_diet_veg",
  "SEV_agestd_metab_bmi", "mean_LDL"
)
covariate_labels <- c(
  mean_value = "SDI",
  SEV_agestd_air_hap = "Household air-pollution SEV",
  SEV_agestd_air_pm = "Ambient particulate-matter SEV",
  SEV_agestd_diet_fruit = "Low-fruit-intake SEV",
  SEV_agestd_diet_veg = "Low-vegetable-intake SEV",
  SEV_agestd_metab_bmi = "High-BMI SEV",
  mean_LDL = "Mean LDL"
)

selected_columns <- c(
  "measure_id", "measure_name", "metric_id", "metric_name",
  "location_id", "location", "sex_id", "sex_name", "age_id", "age_name",
  "year", "val", "lower", "upper", "population", "superregion", "region_wb", covariate_vars,
  "fid", "iso_a3", heat_vars
)

message("Reading selected columns from: ", input_file)
panel <- fread(input_file, select = selected_columns, encoding = "UTF-8", showProgress = interactive())

required_columns <- c(
  "measure_id", "metric_id", "age_id", "sex_id", "location_id", "year",
  "superregion", "region_wb", "val", "lower", "upper", "population", "fid",
  covariate_vars, heat_vars
)
missing_columns <- setdiff(required_columns, names(panel))
if (length(missing_columns) > 0L) {
  stop("Input is missing columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}

assert_unique_panel <- function(dat, stage_label) {
  duplicates <- dat[, .N, by = .(location_id, year, measure_id)][N > 1L, .N]
  if (duplicates > 0L) stop("Duplicate location-year-outcome rows in ", stage_label, ".", call. = FALSE)
}

add_quality_flags <- function(dat) {
  dat <- copy(dat)
  dat[, invalid_outcome := !is.finite(val) | val < 0]
  dat[, invalid_outcome_bounds := !is.finite(lower) | !is.finite(upper) |
    lower < 0 | lower > val | upper < val
  ]
  dat[, invalid_any_covariate := !Reduce(`&`, lapply(.SD, is.finite)), .SDcols = covariate_vars]
  dat[, invalid_population := !is.finite(population) | population <= 0]
  dat[, missing_superregion := is.na(superregion) | trimws(superregion) == ""]
  dat[, missing_region_wb := is.na(region_wb) | trimws(region_wb) == ""]
  dat[, invalid_any_heat := !Reduce(`&`, lapply(.SD, is.finite)), .SDcols = heat_vars]
  dat[, complete_case := !(
    invalid_outcome | invalid_outcome_bounds | invalid_any_covariate | invalid_population |
      missing_superregion | missing_region_wb | invalid_any_heat
  )]
  dat
}

# Prespecified main sample: count, all ages, both sexes, four outcomes.
eligible <- panel[
  metric_id == 1L & age_id == 22L & sex_id == 3L &
    measure_id %in% expected_measures
]
assert_unique_panel(eligible, "main analysis")
eligible <- add_quality_flags(eligible)
analysis_data <- copy(eligible[complete_case == TRUE])
quality_columns <- c(
  "invalid_outcome", "invalid_outcome_bounds", "invalid_any_covariate", "invalid_population",
  "missing_superregion", "missing_region_wb", "invalid_any_heat", "complete_case"
)
analysis_data[, (quality_columns) := NULL]

if (!setequal(unique(analysis_data$measure_id), expected_measures)) {
  stop("Not all four outcomes remain after filtering.", call. = FALSE)
}

balance <- analysis_data[, .(
  n_years = uniqueN(year), first_year = min(year), last_year = max(year)
), by = .(measure_id, location_id)]
if (balance[, any(n_years != 31L | first_year != 1990L | last_year != 2020L)]) {
  stop("The complete-case panel is not balanced over 1990-2020.", call. = FALSE)
}

count_location_years <- function(dat) dat[, uniqueN(paste(location_id, year, sep = ":"))]
sample_flow <- rbindlist(list(
  data.table(
    stage = "Input", step = 1L, rule = "输入文件全部行",
    rows = nrow(panel), location_years = count_location_years(panel),
    locations = panel[, uniqueN(location_id)]
  ),
  data.table(
    stage = "A", step = 2L, rule = "数量、全部年龄、两性合计、4类结局",
    rows = nrow(eligible), location_years = count_location_years(eligible),
    locations = eligible[, uniqueN(location_id)]
  ),
  data.table(
    stage = "A", step = 3L,
    rule = "结局及其输入边界有效，且全部协变量、人口、GBD/世行区域和全部热暴露完整",
    rows = nrow(analysis_data), location_years = count_location_years(analysis_data),
    locations = analysis_data[, uniqueN(location_id)]
  )
), use.names = TRUE)

excluded_locations <- unique(eligible[complete_case == FALSE, .(
  location_id, location, missing_superregion, missing_region_wb, invalid_outcome,
  invalid_outcome_bounds,
  invalid_any_covariate, invalid_population, invalid_any_heat
)])
setorder(excluded_locations, location_id)

# Common exposure scaling reference across every model and outcome.
scaling_panel <- unique(
  analysis_data[, c("location_id", "year", heat_vars), with = FALSE],
  by = c("location_id", "year")
)
exposure_scaling <- rbindlist(lapply(heat_vars, function(v) {
  x <- scaling_panel[[v]]
  data.table(
    heat_variable = v,
    heat_label = unname(heat_labels[[v]]),
    unit = unname(heat_units[[v]]),
    center_median = median(x), iqr = IQR(x),
    q25 = unname(quantile(x, 0.25)), q75 = unname(quantile(x, 0.75)),
    minimum = min(x), maximum = max(x), n_location_years = length(x),
    scaling_reference = "Common complete-case sample"
  )
}))
if (exposure_scaling[, any(!is.finite(iqr) | iqr <= 0)]) {
  stop("At least one heat exposure has a zero or invalid IQR.", call. = FALSE)
}

# Supplement Table 1 source: study-period medians by GBD superregion.
descriptive_heat_covariate_vars <- c(
  "daytime_cumulative_heat", "nighttime_cumulative_heat",
  "compound_day_night_cumulative_heat", covariate_vars
)
descriptive_base_panel <- unique(
  analysis_data[
    measure_id == 1L,
    c("location_id", "year", "superregion", descriptive_heat_covariate_vars),
    with = FALSE
  ],
  by = c("location_id", "year")
)
descriptive_outcome_vars <- c(
  "death_burden", "dalys_burden", "prevalence_burden", "incidence_burden"
)
descriptive_outcome_panel <- dcast(
  analysis_data[, .(location_id, year, superregion, measure_id, val)],
  location_id + year + superregion ~ measure_id,
  value.var = "val"
)
if (!all(as.character(expected_measures) %in% names(descriptive_outcome_panel))) {
  stop("The descriptive panel is missing one or more prespecified CVD outcomes.", call. = FALSE)
}
setnames(
  descriptive_outcome_panel,
  as.character(expected_measures),
  descriptive_outcome_vars
)
descriptive_panel <- merge(
  descriptive_base_panel,
  descriptive_outcome_panel,
  by = c("location_id", "year", "superregion"),
  all = FALSE,
  sort = FALSE
)
descriptive_vars <- c(
  "daytime_cumulative_heat", "nighttime_cumulative_heat",
  "compound_day_night_cumulative_heat", "mean_value",
  descriptive_outcome_vars,
  setdiff(covariate_vars, "mean_value")
)
descriptive_superregion_medians <- descriptive_panel[
  , c(
    list(n_countries = uniqueN(location_id), n_country_years = .N),
    lapply(.SD, median)
  ),
  by = superregion,
  .SDcols = descriptive_vars
]
setorder(descriptive_superregion_medians, superregion)

model_definitions <- data.table(
  model_id = model_ids,
  outcome_rows = "metric_id=1; age_id=22; sex_id=3",
  family = "quasi-Poisson log",
  effect_scale = "log rate ratio per IQR",
  population_offset = TRUE,
  covariates = c(
    "None", rep(paste(covariate_vars, collapse = " + "), length(model_ids) - 1L)
  ),
  year_fixed_effects = model_ids %in% c("M2", "M4"),
  region_wb_fixed_effects = model_ids == "M3a",
  superregion_fixed_effects = model_ids %in% c("M3b", "M4"),
  regional_spline_df3 = model_ids == "M5",
  regional_spline_interaction = ifelse(
    model_ids == "M5", "ns(year, df=3):factor(superregion)", NA_character_
  ),
  model_role = ifelse(model_ids == core_model_id, "Core continuous DID", "Comparison")
)

cluster_vcov <- function(model, cluster_id) {
  vcovCL(model, cluster = cluster_id, type = "HC1", cadjust = TRUE)
}

cluster_inference <- function(model, cluster_id, coefficient = "heat_iqr") {
  clustered_vcov <- cluster_vcov(model, cluster_id)
  beta <- unname(coef(model)[coefficient])
  se <- unname(sqrt(diag(clustered_vcov))[coefficient])
  n_clusters <- uniqueN(cluster_id)
  inferential_df <- n_clusters - 1L
  critical_value <- qt(0.975, df = inferential_df)
  list(
    beta = beta, se = se,
    ci_low = beta - critical_value * se,
    ci_high = beta + critical_value * se,
    p_value = 2 * pt(abs(beta / se), df = inferential_df, lower.tail = FALSE),
    n_clusters = n_clusters, inferential_df = inferential_df
  )
}

scale_heat <- function(dat, heat_var) {
  scale_row <- exposure_scaling[heat_variable == heat_var]
  dat <- copy(dat)
  dat[, heat_iqr := (get(heat_var) - scale_row$center_median) / scale_row$iqr]
  dat
}

model_rhs <- function(exposure_term, model_id) {
  covariate_rhs <- paste(covariate_vars, collapse = " + ")
  switch(
    model_id,
    M0 = exposure_term,
    M1 = paste(exposure_term, "+", covariate_rhs),
    M2 = paste(exposure_term, "+", covariate_rhs, "+ factor(year)"),
    M3a = paste(exposure_term, "+", covariate_rhs, "+ factor(region_wb)"),
    M3b = paste(exposure_term, "+", covariate_rhs, "+ factor(superregion)"),
    M4 = paste(
      exposure_term, "+", covariate_rhs,
      "+ factor(year) + factor(superregion)"
    ),
    M5 = paste(
      exposure_term, "+", covariate_rhs,
      "+ ns(year, df = 3):factor(superregion)"
    ),
    stop("Unknown model id: ", model_id, call. = FALSE)
  )
}

fit_linear_model <- function(dat, heat_var, model_id) {
  dat <- scale_heat(dat, heat_var)
  model <- glm(
    as.formula(paste0(
      "val ~ ", model_rhs("heat_iqr", model_id),
      " + offset(log(population))"
    )),
    family = quasipoisson(link = "log"), data = dat,
    control = glm.control(maxit = 100L, epsilon = 1e-8)
  )
  inference <- cluster_inference(model, dat$location_id)
  scale_row <- exposure_scaling[heat_variable == heat_var]
  result <- data.table(
    measure_id = unique(dat$measure_id),
    outcome = unname(outcome_labels[[as.character(unique(dat$measure_id))]]),
    heat_variable = heat_var, heat_label = unname(heat_labels[[heat_var]]),
    increment = scale_row$iqr, increment_unit = "IQR", model_id = model_id,
    family = "quasi-Poisson log", effect_scale = "log rate ratio per IQR",
    beta = inference$beta, std_error_clustered = inference$se,
    ci_low = inference$ci_low, ci_high = inference$ci_high,
    rate_ratio = exp(inference$beta), rr_ci_low = exp(inference$ci_low),
    rr_ci_high = exp(inference$ci_high), p_value = inference$p_value,
    n_observations = nobs(model), n_countries = inference$n_clusters,
    inference_df = inference$inferential_df
  )
  diagnostic <- data.table(
    analysis_type = "linear", measure_id = unique(dat$measure_id),
    heat_variable = heat_var, model_id = model_id, family = "quasi-Poisson log",
    converged = isTRUE(model$converged), iterations = model$iter,
    rank = model$rank, df_residual = df.residual(model),
    dispersion = summary(model)$dispersion, deviance = deviance(model),
    n_observations = nobs(model), n_countries = inference$n_clusters
  )
  list(result = result, diagnostic = diagnostic)
}

message("Fitting 8 exposures x 4 outcomes x 7 prespecified model specifications...")
fits <- list()
fit_index <- 0L
for (measure in expected_measures) {
  dat <- analysis_data[measure_id == measure]
  for (heat_var in heat_vars) {
    for (model_id in model_ids) {
      fit_index <- fit_index + 1L
      fits[[fit_index]] <- fit_linear_model(dat, heat_var, model_id)
    }
  }
}

results <- rbindlist(lapply(fits, `[[`, "result"), use.names = TRUE)
diagnostics <- rbindlist(lapply(fits, `[[`, "diagnostic"), use.names = TRUE)
results[, p_adjust_bh := p.adjust(p_value, method = "BH"), by = model_id]
results[, `:=`(
  heat_order = match(heat_variable, heat_vars),
  outcome_order = match(measure_id, expected_measures),
  model_order = match(model_id, model_ids)
)]
setorder(results, heat_order, outcome_order, model_order)
results[, c("heat_order", "outcome_order", "model_order") := NULL]

# Additional reviewer-requested country and year two-way fixed-effects model.
# This is reported separately from M0--M5 because it changes the identifying
# variation from regional trends to within-country changes net of common year
# shocks. The outcome, covariates, offset, heat scaling, and clustered
# inference otherwise match the main analysis.
fit_country_year_fe_model <- function(dat, heat_var) {
  dat <- scale_heat(dat, heat_var)
  covariate_rhs <- paste(covariate_vars, collapse = " + ")
  model <- glm(
    as.formula(paste0(
      "val ~ heat_iqr + ", covariate_rhs,
      " + factor(location_id) + factor(year) + offset(log(population))"
    )),
    family = quasipoisson(link = "log"), data = dat,
    control = glm.control(maxit = 100L, epsilon = 1e-8)
  )
  inference <- cluster_inference(model, dat$location_id)
  scale_row <- exposure_scaling[heat_variable == heat_var]
  result <- data.table(
    specification_id = "country_year_twfe",
    specification_label = "Country and year two-way fixed effects",
    measure_id = unique(dat$measure_id),
    outcome = unname(outcome_labels[[as.character(unique(dat$measure_id))]]),
    heat_variable = heat_var,
    heat_label = unname(heat_labels[[heat_var]]),
    increment = scale_row$iqr,
    increment_unit = "Primary-sample IQR",
    family = "quasi-Poisson log",
    effect_scale = "log rate ratio per primary-sample IQR",
    covariates = paste(covariate_vars, collapse = " + "),
    country_fixed_effects = TRUE,
    year_fixed_effects = TRUE,
    population_offset = TRUE,
    beta = inference$beta,
    std_error_clustered = inference$se,
    ci_low = inference$ci_low,
    ci_high = inference$ci_high,
    rate_ratio = exp(inference$beta),
    rr_ci_low = exp(inference$ci_low),
    rr_ci_high = exp(inference$ci_high),
    p_value = inference$p_value,
    n_observations = nobs(model),
    n_countries = inference$n_clusters,
    analysis_year_first = min(dat$year),
    analysis_year_last = max(dat$year),
    inference_df = inference$inferential_df
  )
  diagnostic <- data.table(
    analysis_type = "country_year_twfe",
    measure_id = unique(dat$measure_id),
    heat_variable = heat_var,
    model_id = "country_year_twfe",
    family = "quasi-Poisson log",
    converged = isTRUE(model$converged),
    iterations = model$iter,
    rank = model$rank,
    df_residual = df.residual(model),
    dispersion = summary(model)$dispersion,
    deviance = deviance(model),
    n_observations = nobs(model),
    n_countries = inference$n_clusters
  )
  list(result = result, diagnostic = diagnostic)
}

message("Fitting reviewer-requested country and year two-way fixed-effects models...")
twfe_fits <- list()
twfe_index <- 0L
for (measure in expected_measures) {
  dat <- analysis_data[measure_id == measure]
  for (heat_var in heat_vars) {
    twfe_index <- twfe_index + 1L
    twfe_fits[[twfe_index]] <- fit_country_year_fe_model(dat, heat_var)
  }
}
country_year_twfe_results <- rbindlist(
  lapply(twfe_fits, `[[`, "result"), use.names = TRUE
)
country_year_twfe_diagnostics <- rbindlist(
  lapply(twfe_fits, `[[`, "diagnostic"), use.names = TRUE
)
country_year_twfe_results[, p_adjust_bh := p.adjust(p_value, method = "BH")]
country_year_twfe_results[, `:=`(
  heat_order = match(heat_variable, heat_vars),
  outcome_order = match(measure_id, expected_measures)
)]
setorder(country_year_twfe_results, heat_order, outcome_order)
country_year_twfe_results[, c("heat_order", "outcome_order") := NULL]

# Indicator-specific mutual-adjustment analysis requested after review.
# CDNCH and RCH are each scaled by their own primary-sample IQR. Single-
# exposure estimates come from core M5; this block fits their joint core-M5
# model so that before/after mutual adjustment remains directly comparable.
mutual_heat_vars <- c(
  CDNCH = "compound_day_night_cumulative_heat",
  RCH = "relative_cumulative_heat"
)
mutual_indicator_correlation <- cor(
  scaling_panel[[mutual_heat_vars[["CDNCH"]]]],
  scaling_panel[[mutual_heat_vars[["RCH"]]]]
)
mutual_pairwise_vif <- 1 / (1 - mutual_indicator_correlation^2)

# Exposure-specific VIFs in the full joint-M5 design matrix. The population
# offset is not part of the design matrix and therefore does not enter VIF.
mutual_vif_data <- copy(analysis_data[measure_id == expected_measures[[1L]]])
mutual_vif_data[, cdnch_iqr := (
  get(mutual_heat_vars[["CDNCH"]]) -
    exposure_scaling[
      heat_variable == mutual_heat_vars[["CDNCH"]]
    ]$center_median
) / exposure_scaling[
  heat_variable == mutual_heat_vars[["CDNCH"]]
]$iqr]
mutual_vif_data[, rch_iqr := (
  get(mutual_heat_vars[["RCH"]]) -
    exposure_scaling[
      heat_variable == mutual_heat_vars[["RCH"]]
    ]$center_median
) / exposure_scaling[
  heat_variable == mutual_heat_vars[["RCH"]]
]$iqr]
mutual_vif_design <- model.matrix(
  as.formula(paste0(
    "~ cdnch_iqr + rch_iqr + ",
    paste(covariate_vars, collapse = " + "),
    " + ns(year, df = 3):factor(superregion)"
  )),
  data = mutual_vif_data
)
calculate_design_vif <- function(design_matrix, coefficient) {
  response <- design_matrix[, coefficient]
  other_columns <- setdiff(
    colnames(design_matrix), c("(Intercept)", coefficient)
  )
  auxiliary_fit <- lm.fit(
    cbind(`(Intercept)` = 1, design_matrix[, other_columns, drop = FALSE]),
    response
  )
  r_squared <- 1 - sum(auxiliary_fit$residuals^2) /
    sum((response - mean(response))^2)
  list(r_squared = r_squared, vif = 1 / (1 - r_squared))
}
cdnch_vif <- calculate_design_vif(mutual_vif_design, "cdnch_iqr")
rch_vif <- calculate_design_vif(mutual_vif_design, "rch_iqr")
cdnch_rch_collinearity_diagnostics <- data.table(
  n_country_years = nrow(mutual_vif_data),
  pearson_correlation = mutual_indicator_correlation,
  pairwise_vif = mutual_pairwise_vif,
  cdnch_full_m5_r_squared = cdnch_vif$r_squared,
  cdnch_full_m5_vif = cdnch_vif$vif,
  rch_full_m5_r_squared = rch_vif$r_squared,
  rch_full_m5_vif = rch_vif$vif,
  vif_design = paste0(
    "CDNCH + RCH + 7 covariates + ",
    "ns(year, df=3):factor(superregion)"
  )
)

fit_cd_rch_joint_model <- function(dat) {
  dat <- copy(dat)
  cd_scale <- exposure_scaling[
    heat_variable == mutual_heat_vars[["CDNCH"]]
  ]
  rch_scale <- exposure_scaling[
    heat_variable == mutual_heat_vars[["RCH"]]
  ]
  dat[, cdnch_iqr := (
    get(mutual_heat_vars[["CDNCH"]]) - cd_scale$center_median
  ) / cd_scale$iqr]
  dat[, rch_iqr := (
    get(mutual_heat_vars[["RCH"]]) - rch_scale$center_median
  ) / rch_scale$iqr]
  model <- glm(
    as.formula(paste0(
      "val ~ ", model_rhs("cdnch_iqr + rch_iqr", core_model_id),
      " + offset(log(population))"
    )),
    family = quasipoisson(link = "log"), data = dat,
    control = glm.control(maxit = 100L, epsilon = 1e-8)
  )
  clustered_vcov <- cluster_vcov(model, dat$location_id)
  n_clusters <- uniqueN(dat$location_id)
  inferential_df <- n_clusters - 1L
  critical_value <- qt(0.975, df = inferential_df)

  coefficient_rows <- data.table(
    indicator = c("CDNCH", "RCH"),
    heat_variable = unname(mutual_heat_vars[c("CDNCH", "RCH")]),
    coefficient = c("cdnch_iqr", "rch_iqr"),
    increment = c(cd_scale$iqr, rch_scale$iqr)
  )
  coefficient_rows[, `:=`(
    measure_id = unique(dat$measure_id),
    outcome = unname(outcome_labels[[as.character(unique(dat$measure_id))]]),
    model_group = "CDNCH+RCH",
    model_id = "M5_CDNCH_RCH_joint",
    family = "quasi-Poisson log",
    effect_scale = "log rate ratio per indicator-specific primary-sample IQR",
    n_observations = nobs(model),
    n_countries = n_clusters,
    analysis_year_first = min(dat$year),
    analysis_year_last = max(dat$year),
    inference_df = inferential_df
  )]
  coefficient_rows[, beta := unname(coef(model)[coefficient])]
  coefficient_rows[, std_error_clustered := unname(
    sqrt(diag(clustered_vcov))[coefficient]
  )]
  coefficient_rows[, `:=`(
    ci_low = beta - critical_value * std_error_clustered,
    ci_high = beta + critical_value * std_error_clustered,
    p_value = 2 * pt(
      abs(beta / std_error_clustered), df = inferential_df,
      lower.tail = FALSE
    )
  )]
  coefficient_rows[, `:=`(
    rate_ratio = exp(beta),
    rr_ci_low = exp(ci_low),
    rr_ci_high = exp(ci_high)
  )]
  setcolorder(coefficient_rows, c(
    "measure_id", "outcome", "indicator", "heat_variable", "increment",
    "model_group", "model_id", "family", "effect_scale", "beta",
    "std_error_clustered", "ci_low", "ci_high", "rate_ratio", "rr_ci_low",
    "rr_ci_high", "p_value", "n_observations", "n_countries",
    "analysis_year_first", "analysis_year_last", "inference_df", "coefficient"
  ))
  diagnostic <- data.table(
    analysis_type = "m5_cd_rch_joint",
    measure_id = unique(dat$measure_id),
    heat_variable = "compound_day_night_cumulative_heat + relative_cumulative_heat",
    model_id = "M5_CDNCH_RCH_joint",
    family = "quasi-Poisson log",
    converged = isTRUE(model$converged),
    iterations = model$iter,
    rank = model$rank,
    df_residual = df.residual(model),
    dispersion = summary(model)$dispersion,
    deviance = deviance(model),
    n_observations = nobs(model),
    n_countries = n_clusters
  )
  list(result = coefficient_rows, diagnostic = diagnostic)
}

message("Fitting core-M5 joint CDNCH + RCH mutual-adjustment models...")
mutual_adjustment_fits <- lapply(expected_measures, function(measure) {
  fit_cd_rch_joint_model(analysis_data[measure_id == measure])
})
cdnch_rch_joint_results <- rbindlist(
  lapply(mutual_adjustment_fits, `[[`, "result"), use.names = TRUE
)
cdnch_rch_joint_diagnostics <- rbindlist(
  lapply(mutual_adjustment_fits, `[[`, "diagnostic"), use.names = TRUE
)
cdnch_rch_joint_results[, p_adjust_bh := p.adjust(p_value, method = "BH")]
setorder(cdnch_rch_joint_results, measure_id, indicator)

# Supplementary core-M5 subgroup analyses: rate outcomes by sex and age group.
# These models use the same country-year sample and exposure scaling as the
# primary analysis, but use the M5 adjustment set with OLS (not a count model),
# so there is no population offset. Age-standardized rates are excluded because
# they do not represent a distinct age subgroup.
subgroup_metric_id <- 3L
subgroup_sex_ids <- c(1L, 2L)
subgroup_heat_vars <- c("compound_day_night_cumulative_heat")
subgroup_sample_keys <- unique(
  scaling_panel[, .(location_id, year)],
  by = c("location_id", "year")
)
subgroup_derived_age_id <- 2070L
subgroup_observed_age_ids <- c(158L, 26L)
subgroup_age_ids <- c(158L, subgroup_derived_age_id, 26L)
subgroup_interaction_age_ids <- c(158L, subgroup_derived_age_id, 26L)
subgroup_age_display_labels <- c(
  `158` = "<20 years", `2070` = "20–70 years (derived)", `26` = ">70 years"
)
subgroup_sex_display_labels <- c(`1` = "Male", `2` = "Female")
subgroup_outcome_display_labels <- c(
  `1` = "Death", `2` = "DALYs", `5` = "Prevalence", `6` = "Incidence"
)
subgroup_heat_display_labels <- c(
  absolute_cumulative_heat = "Absolute-threshold cumulative heat",
  relative_cumulative_heat = "Relative-threshold cumulative heat",
  mixed_cumulative_heat = "Mixed-threshold cumulative heat",
  daytime_cumulative_heat = "Daytime cumulative heat",
  nighttime_cumulative_heat = "Nighttime cumulative heat",
  compound_day_night_cumulative_heat = "Compound day-night cumulative heat",
  annual_mean_humidex = "Annual mean humidex",
  annual_mean_tmean = "Annual mean temperature"
)
if (!all(subgroup_sex_ids %in% unique(panel[metric_id == subgroup_metric_id, sex_id]))) {
  stop("Rate subgroup analysis requires both male and female sex strata.", call. = FALSE)
}
if (!all(subgroup_observed_age_ids %in% unique(panel[metric_id == subgroup_metric_id, age_id]))) {
  stop("An observed subgroup age required by the design is unavailable.", call. = FALSE)
}
if (!all(as.character(subgroup_age_ids) %in% names(subgroup_age_display_labels))) {
  stop("A subgroup age identifier lacks a display label.", call. = FALSE)
}

# Derive the requested middle-age rate from mutually exclusive burden counts
# and age-specific effective population denominators. Direct subtraction of
# rates is invalid, so each component population is recovered exactly from
# count / rate * 100,000 before taking the residual.
derived_component_age_ids <- c(22L, 158L, 26L)
derived_keys <- c("measure_id", "sex_id", "location_id", "year", "age_id")
derived_source <- panel[
  metric_id %in% c(1L, subgroup_metric_id) &
    age_id %in% derived_component_age_ids &
    sex_id %in% subgroup_sex_ids &
    measure_id %in% expected_measures
]
derived_source <- derived_source[
  subgroup_sample_keys,
  on = .(location_id, year),
  nomatch = 0L
]
derived_duplicates <- derived_source[, .N, by = c(derived_keys, "metric_id")][N > 1L, .N]
if (derived_duplicates > 0L) {
  stop("Duplicate rows prevent derivation of the 20–70-year subgroup.", call. = FALSE)
}
derived_counts <- derived_source[metric_id == 1L, c(derived_keys, "val"), with = FALSE]
setnames(derived_counts, "val", "count_value")
derived_rates <- derived_source[metric_id == subgroup_metric_id, c(derived_keys, "val"), with = FALSE]
setnames(derived_rates, "val", "rate_value")
derived_components <- merge(
  derived_counts, derived_rates, by = derived_keys, all = TRUE, sort = FALSE
)
expected_component_rows <- length(expected_measures) * length(subgroup_sex_ids) *
  nrow(subgroup_sample_keys) * length(derived_component_age_ids)
if (nrow(derived_components) != expected_component_rows ||
    derived_components[, any(!is.finite(count_value) | count_value < 0 |
      !is.finite(rate_value) | rate_value <= 0)]) {
  stop("Count/rate components for the 20–70-year subgroup are incomplete or invalid.", call. = FALSE)
}
derived_components[, implied_population := count_value / rate_value * 100000]
if (derived_components[, any(!is.finite(implied_population) | implied_population <= 0)]) {
  stop("Invalid age-specific implied population in the 20–70-year derivation.", call. = FALSE)
}
population_consistency <- derived_components[, .(
  relative_spread = (max(implied_population) - min(implied_population)) / mean(implied_population)
), by = .(sex_id, location_id, year, age_id)]
if (population_consistency[, any(relative_spread > 1e-8)]) {
  stop("Count/rate-implied age populations disagree across outcomes.", call. = FALSE)
}

derived_wide <- dcast(
  derived_components,
  measure_id + sex_id + location_id + year ~ age_id,
  value.var = c("count_value", "implied_population")
)
derived_wide[, `:=`(
  count_all = count_value_22,
  count_under20 = count_value_158,
  count_over70 = count_value_26,
  population_all = implied_population_22,
  population_under20 = implied_population_158,
  population_over70 = implied_population_26,
  derived_count = count_value_22 - count_value_158 - count_value_26,
  derived_population = implied_population_22 - implied_population_158 - implied_population_26
)]
derived_wide[, derived_rate := derived_count / derived_population * 100000]
if (derived_wide[, any(!is.finite(derived_count) | derived_count < 0 |
    !is.finite(derived_population) | derived_population <= 0 |
    !is.finite(derived_rate) | derived_rate < 0)]) {
  stop("The derived 20–70-year count, population, or rate is invalid.", call. = FALSE)
}

derived_template <- copy(panel[
  metric_id == subgroup_metric_id & age_id == 22L &
    sex_id %in% subgroup_sex_ids & measure_id %in% expected_measures
])
derived_template <- derived_template[
  subgroup_sample_keys,
  on = .(location_id, year),
  nomatch = 0L
]
derived_template <- merge(
  derived_template,
  derived_wide[, .(
    measure_id, sex_id, location_id, year, count_all, count_under20, count_over70,
    population_all, population_under20, population_over70,
    derived_count, derived_population, derived_rate
  )],
  by = c("measure_id", "sex_id", "location_id", "year"),
  all.x = TRUE,
  sort = FALSE
)
derived_template[, `:=`(
  age_id = subgroup_derived_age_id,
  age_name = unname(subgroup_age_display_labels[[as.character(subgroup_derived_age_id)]]),
  val = derived_rate,
  age_derivation = "All ages - <20 years - >70 years; rate from residual count / residual implied population"
)]
subgroup_age_derivation_audit <- derived_template[, .(
  measure_id, sex_id, location_id, location, year,
  count_all, count_under20, count_over70, derived_count,
  population_all, population_under20, population_over70, derived_population,
  derived_rate, age_derivation
)]

add_rate_subgroup_quality_flags <- function(dat) {
  dat <- copy(dat)
  dat[, invalid_outcome := !is.finite(val) | val < 0]
  dat[, invalid_any_covariate := !Reduce(`&`, lapply(.SD, is.finite)), .SDcols = covariate_vars]
  dat[, missing_superregion := is.na(superregion) | trimws(superregion) == ""]
  dat[, missing_region_wb := is.na(region_wb) | trimws(region_wb) == ""]
  dat[, invalid_any_heat := !Reduce(`&`, lapply(.SD, is.finite)), .SDcols = subgroup_heat_vars]
  dat[, invalid_derived_age := fifelse(
    age_id == subgroup_derived_age_id,
    !is.finite(derived_count) | derived_count < 0 |
      !is.finite(derived_population) | derived_population <= 0,
    FALSE
  )]
  dat[, complete_case := !(invalid_outcome | invalid_any_covariate |
    missing_superregion | missing_region_wb | invalid_any_heat | invalid_derived_age)]
  dat
}

subgroup_observed <- panel[
  metric_id == subgroup_metric_id &
    age_id %in% subgroup_observed_age_ids &
    sex_id %in% subgroup_sex_ids &
    measure_id %in% expected_measures
]
subgroup_observed <- subgroup_observed[
  subgroup_sample_keys,
  on = .(location_id, year),
  nomatch = 0L
]
subgroup_eligible <- rbindlist(list(subgroup_observed, derived_template), use.names = TRUE, fill = TRUE)
subgroup_duplicates <- subgroup_eligible[, .N,
  by = .(location_id, year, measure_id, age_id, sex_id)
][N > 1L, .N]
if (subgroup_duplicates > 0L) {
  stop("Duplicate location-year-outcome-age-sex rows in rate subgroup analysis.", call. = FALSE)
}
subgroup_eligible <- add_rate_subgroup_quality_flags(subgroup_eligible)
subgroup_data <- copy(subgroup_eligible[complete_case == TRUE])
subgroup_quality_columns <- c(
  "invalid_outcome", "invalid_any_covariate", "missing_superregion",
  "missing_region_wb", "invalid_any_heat", "invalid_derived_age", "complete_case"
)
subgroup_data[, (subgroup_quality_columns) := NULL]

subgroup_balance <- subgroup_data[, .(
  n_years = uniqueN(year), first_year = min(year), last_year = max(year)
), by = .(measure_id, age_id, sex_id, location_id)]
if (subgroup_balance[, any(n_years != 31L | first_year != 1990L | last_year != 2020L)]) {
  stop("The rate subgroup panel is not balanced over 1990-2020.", call. = FALSE)
}
subgroup_cell_counts <- subgroup_data[, .(
  n_observations = .N, n_countries = uniqueN(location_id)
), by = .(measure_id, age_id, sex_id)]
expected_subgroup_observations <- nrow(subgroup_sample_keys)
expected_subgroup_countries <- subgroup_sample_keys[, uniqueN(location_id)]
if (nrow(subgroup_cell_counts) != length(expected_measures) * length(subgroup_age_ids) *
    length(subgroup_sex_ids) ||
    subgroup_cell_counts[, any(
      n_observations != expected_subgroup_observations |
        n_countries != expected_subgroup_countries
    )]) {
  stop("A subgroup cell does not match the primary-analysis country-year sample.", call. = FALSE)
}

subgroup_scaling_panel <- unique(
  subgroup_data[, c("location_id", "year", subgroup_heat_vars), with = FALSE],
  by = c("location_id", "year")
)
if (nrow(subgroup_scaling_panel) != nrow(subgroup_sample_keys) ||
    subgroup_scaling_panel[, uniqueN(location_id)] != expected_subgroup_countries) {
  stop("The subgroup scaling panel does not match the primary-analysis sample.", call. = FALSE)
}
subgroup_exposure_scaling <- copy(
  exposure_scaling[heat_variable %in% subgroup_heat_vars]
)
subgroup_exposure_scaling[, scaling_reference :=
  "Primary-analysis complete-case country-year panel (shared with subgroup analysis)"
]
if (subgroup_exposure_scaling[, any(!is.finite(iqr) | iqr <= 0)]) {
  stop("At least one subgroup heat exposure has a zero or invalid IQR.", call. = FALSE)
}

scale_subgroup_heat <- function(dat, heat_var) {
  scale_row <- subgroup_exposure_scaling[heat_variable == heat_var]
  dat <- copy(dat)
  dat[, heat_iqr := (get(heat_var) - scale_row$center_median) / scale_row$iqr]
  dat
}

fit_subgroup_m5_model <- function(dat, heat_var) {
  dat <- scale_subgroup_heat(dat, heat_var)
  model <- lm(
    as.formula(paste0("val ~ ", model_rhs("heat_iqr", core_model_id))),
    data = dat
  )
  inference <- cluster_inference(model, dat$location_id)
  scale_row <- subgroup_exposure_scaling[heat_variable == heat_var]
  result <- data.table(
    metric_id = subgroup_metric_id,
    metric = "Rate",
    measure_id = unique(dat$measure_id),
    outcome = unname(outcome_labels[[as.character(unique(dat$measure_id))]]),
    age_id = unique(dat$age_id),
    age_group = unique(dat$age_name),
    age_derivation = if (unique(dat$age_id) == subgroup_derived_age_id) {
      unique(na.omit(dat$age_derivation))
    } else {
      "Observed rate row"
    },
    sex_id = unique(dat$sex_id),
    sex = unique(dat$sex_name),
    heat_variable = heat_var,
    heat_label = unname(heat_labels[[heat_var]]),
    increment = scale_row$iqr,
    increment_unit = "Common primary-analysis IQR",
    model_id = core_model_id,
    family = "OLS",
    effect_scale = "Rate difference per common primary-analysis IQR",
    adjustment_terms = paste(covariate_vars, collapse = " + "),
    regional_trend = "ns(year, df=3):factor(superregion)",
    population_offset = FALSE,
    beta = inference$beta,
    std_error_clustered = inference$se,
    ci_low = inference$ci_low,
    ci_high = inference$ci_high,
    p_value = inference$p_value,
    n_observations = nobs(model),
    n_countries = inference$n_clusters,
    analysis_year_first = min(dat$year),
    analysis_year_last = max(dat$year),
    inference_df = inference$inferential_df
  )
  diagnostic <- data.table(
    analysis_type = "subgroup_rate_m5",
    metric_id = subgroup_metric_id,
    measure_id = unique(dat$measure_id),
    age_id = unique(dat$age_id),
    sex_id = unique(dat$sex_id),
    heat_variable = heat_var,
    model_id = core_model_id,
    family = "OLS",
    converged = TRUE,
    iterations = NA_integer_,
    rank = model$rank,
    df_residual = df.residual(model),
    dispersion = NA_real_,
    deviance = sum(residuals(model)^2),
    n_observations = nobs(model),
    n_countries = inference$n_clusters
  )
  list(result = result, diagnostic = diagnostic)
}

message("Fitting core M5 rate-outcome subgroup models by sex and age group...")
subgroup_fits <- list()
subgroup_index <- 0L
for (measure in expected_measures) {
  for (age in subgroup_age_ids) {
    for (sex in subgroup_sex_ids) {
      dat <- subgroup_data[measure_id == measure & age_id == age & sex_id == sex]
      for (heat_var in subgroup_heat_vars) {
        subgroup_index <- subgroup_index + 1L
        subgroup_fits[[subgroup_index]] <- fit_subgroup_m5_model(dat, heat_var)
      }
    }
  }
}
subgroup_m5_results <- rbindlist(lapply(subgroup_fits, `[[`, "result"), use.names = TRUE)
subgroup_m5_diagnostics <- rbindlist(lapply(subgroup_fits, `[[`, "diagnostic"), use.names = TRUE)
subgroup_m5_results[, p_adjust_bh := p.adjust(p_value, method = "BH"),
  by = .(measure_id, age_id, sex_id)
]
subgroup_m5_results[, `:=`(
  age_order = match(age_id, subgroup_age_ids),
  heat_order = match(heat_variable, subgroup_heat_vars)
)]
setorder(subgroup_m5_results, measure_id, age_order, sex_id, heat_order)
subgroup_m5_results[, c("age_order", "heat_order") := NULL]
subgroup_m5_summary <- subgroup_m5_results[, .(
  tests = .N,
  raw_p_lt_0_05 = sum(p_value < 0.05),
  bh_q_lt_0_05 = sum(p_adjust_bh < 0.05),
  n_observations = unique(n_observations),
  n_countries = unique(n_countries),
  first_year = unique(analysis_year_first),
  last_year = unique(analysis_year_last)
), by = .(measure_id, outcome, age_id, age_group, sex_id, sex)]
subgroup_m5_summary[, age_order := match(age_id, subgroup_age_ids)]
setorder(subgroup_m5_summary, measure_id, age_order, sex_id)
subgroup_m5_summary[, age_order := NULL]

# Formal effect-modification tests use the three mutually exclusive age groups
# (<20, derived 20–70, >70) and both sexes in one pooled model per outcome.
cluster_wald_test <- function(model, cluster_id, coefficient_names, test_label) {
  if (length(coefficient_names) == 0L) {
    stop("No coefficients selected for interaction test: ", test_label, call. = FALSE)
  }
  vc <- cluster_vcov(model, cluster_id)
  b <- coef(model)[coefficient_names]
  v <- vc[coefficient_names, coefficient_names, drop = FALSE]
  if (any(!is.finite(b)) || any(!is.finite(v))) {
    stop("Non-finite interaction estimate or covariance: ", test_label, call. = FALSE)
  }
  q <- length(b)
  wald_chisq <- as.numeric(t(b) %*% solve(v, b))
  n_clusters <- uniqueN(cluster_id)
  df_den <- n_clusters - 1L
  wald_f <- wald_chisq / q
  list(
    wald_chisq = wald_chisq,
    wald_f = wald_f,
    df_num = q,
    df_den = df_den,
    p_value = pf(wald_f, df1 = q, df2 = df_den, lower.tail = FALSE),
    n_clusters = n_clusters,
    coefficient_names = paste(coefficient_names, collapse = " | ")
  )
}

# Recover the heat-IQR slope for any age-sex cell from either a conventional
# age-by-sex interaction model or the fully stratified stacked model below.
subgroup_heat_slope_vector <- function(model, age_id, sex_id) {
  coefficient_names <- names(coef(model))
  slope_vector <- setNames(rep(0, length(coefficient_names)), coefficient_names)
  heat_terms <- coefficient_names[
    vapply(strsplit(coefficient_names, ":", fixed = TRUE), function(parts) {
      "heat_iqr" %in% parts
    }, logical(1))
  ]
  expected_age_term <- paste0("age_factor", age_id)
  expected_sex_term <- paste0("sex_factor", sex_id)
  expected_stratum_term <- paste0("stratum_factor", age_id, "_", sex_id)
  for (term in heat_terms) {
    parts <- strsplit(term, ":", fixed = TRUE)[[1L]]
    age_terms <- parts[grepl("^age_factor", parts)]
    sex_terms <- parts[grepl("^sex_factor", parts)]
    stratum_terms <- parts[grepl("^stratum_factor", parts)]
    age_applies <- length(age_terms) == 0L || expected_age_term %in% age_terms
    sex_applies <- length(sex_terms) == 0L || expected_sex_term %in% sex_terms
    stratum_applies <- length(stratum_terms) == 0L ||
      expected_stratum_term %in% stratum_terms
    if (age_applies && sex_applies && stratum_applies) slope_vector[[term]] <- 1
  }
  if (sum(abs(slope_vector)) == 0) {
    stop("Unable to recover an age-sex heat slope from the interaction model.", call. = FALSE)
  }
  slope_vector
}

cluster_linear_contrast <- function(model, clustered_vcov, contrast_vector, df_den) {
  b <- coef(model)
  estimate <- sum(contrast_vector * b)
  variance <- as.numeric(t(contrast_vector) %*% clustered_vcov %*% contrast_vector)
  if (!is.finite(estimate) || !is.finite(variance) || variance <= 0) {
    stop("Invalid estimate or variance in subgroup linear contrast.", call. = FALSE)
  }
  se <- sqrt(variance)
  t_value <- estimate / se
  critical_value <- qt(0.975, df = df_den)
  list(
    estimate = estimate,
    std_error_clustered = se,
    ci_low = estimate - critical_value * se,
    ci_high = estimate + critical_value * se,
    t_value = t_value,
    df_den = df_den,
    p_value = 2 * pt(abs(t_value), df = df_den, lower.tail = FALSE)
  )
}

cluster_joint_contrast <- function(model, clustered_vcov, contrast_matrix, df_den) {
  b <- coef(model)
  contrast_estimates <- as.numeric(contrast_matrix %*% b)
  contrast_covariance <- contrast_matrix %*% clustered_vcov %*% t(contrast_matrix)
  if (any(!is.finite(contrast_estimates)) || any(!is.finite(contrast_covariance))) {
    stop("Invalid estimate or covariance in subgroup joint contrast.", call. = FALSE)
  }
  q <- nrow(contrast_matrix)
  wald_chisq <- as.numeric(
    t(contrast_estimates) %*% solve(contrast_covariance, contrast_estimates)
  )
  wald_f <- wald_chisq / q
  list(
    contrast_estimates = contrast_estimates,
    wald_chisq = wald_chisq,
    wald_f = wald_f,
    df_num = q,
    df_den = df_den,
    p_value = pf(wald_f, df1 = q, df2 = df_den, lower.tail = FALSE)
  )
}

fit_subgroup_interaction_model <- function(dat) {
  dat <- scale_subgroup_heat(dat, subgroup_heat_vars[[1L]])
  dat[, age_factor := factor(age_id, levels = subgroup_interaction_age_ids)]
  dat[, sex_factor := factor(sex_id, levels = subgroup_sex_ids)]
  dat[, stratum_factor := factor(
    paste(age_id, sex_id, sep = "_"),
    levels = as.vector(outer(
      subgroup_interaction_age_ids, subgroup_sex_ids, paste, sep = "_"
    ))
  )]

  # Interact every adjustment term with the age-sex stratum so each heat slope
  # is identical to its separately fitted subgroup M5 estimate. This prevents
  # the contrast from depending on an implausible common covariate or regional
  # time-trend coefficient across populations with very different rate scales.
  model <- lm(
    as.formula(paste0(
      "val ~ 0 + stratum_factor + stratum_factor:heat_iqr + ",
      paste0("stratum_factor:", covariate_vars, collapse = " + "),
      " + stratum_factor:ns(year, df = 3):factor(superregion)"
    )),
    data = dat
  )
  clustered_vcov <- cluster_vcov(model, dat$location_id)
  df_den <- uniqueN(dat$location_id) - 1L
  cell_slope_vectors <- list()
  for (age in subgroup_interaction_age_ids) {
    for (sex in subgroup_sex_ids) {
      cell_key <- paste(age, sex, sep = "_")
      cell_slope_vectors[[cell_key]] <- subgroup_heat_slope_vector(model, age, sex)
    }
  }
  get_cell_slope <- function(age, sex) {
    cell_slope_vectors[[paste(age, sex, sep = "_")]]
  }

  age_contrast_matrix <- rbind(
    male_middle_minus_under20 = get_cell_slope(subgroup_derived_age_id, 1L) -
      get_cell_slope(158L, 1L),
    male_older_minus_under20 = get_cell_slope(26L, 1L) - get_cell_slope(158L, 1L),
    female_middle_minus_under20 = get_cell_slope(subgroup_derived_age_id, 2L) -
      get_cell_slope(158L, 2L),
    female_older_minus_under20 = get_cell_slope(26L, 2L) - get_cell_slope(158L, 2L)
  )
  sex_contrast_matrix <- rbind(
    under20_male_minus_female = get_cell_slope(158L, 1L) - get_cell_slope(158L, 2L),
    middle_male_minus_female = get_cell_slope(subgroup_derived_age_id, 1L) -
      get_cell_slope(subgroup_derived_age_id, 2L),
    older_male_minus_female = get_cell_slope(26L, 1L) - get_cell_slope(26L, 2L)
  )
  age_sex_contrast_matrix <- rbind(
    middle_vs_under20_sex_difference = sex_contrast_matrix[2L, ] -
      sex_contrast_matrix[1L, ],
    older_vs_under20_sex_difference = sex_contrast_matrix[3L, ] -
      sex_contrast_matrix[1L, ]
  )
  tests <- list(
    Age = cluster_joint_contrast(model, clustered_vcov, age_contrast_matrix, df_den),
    Sex = cluster_joint_contrast(model, clustered_vcov, sex_contrast_matrix, df_den),
    `Age × sex` = cluster_joint_contrast(
      model, clustered_vcov, age_sex_contrast_matrix, df_den
    )
  )
  tested_constraints <- c(
    Age = paste(rownames(age_contrast_matrix), collapse = " | "),
    Sex = paste(rownames(sex_contrast_matrix), collapse = " | "),
    `Age × sex` = paste(rownames(age_sex_contrast_matrix), collapse = " | ")
  )
  test_table <- rbindlist(lapply(names(tests), function(test_name) {
    test <- tests[[test_name]]
    data.table(
      measure_id = unique(dat$measure_id),
      outcome = unname(outcome_labels[[as.character(unique(dat$measure_id))]]),
      heat_variable = subgroup_heat_vars[[1L]],
      model_id = core_model_id,
      interaction_test = test_name,
      wald_chisq = test$wald_chisq,
      wald_f = test$wald_f,
      df_num = test$df_num,
      df_den = test$df_den,
      p_value = test$p_value,
      tested_coefficients = unname(tested_constraints[[test_name]]),
      n_observations = nobs(model),
      n_countries = uniqueN(dat$location_id),
      age_groups = paste(
        unname(subgroup_age_display_labels[as.character(subgroup_interaction_age_ids)]),
        collapse = " | "
      ),
      sexes = "Male | Female",
      effect_scale = "Formal clustered Wald F test of CDNCH slope heterogeneity"
    )
  }))
  age_label <- function(age) {
    unname(subgroup_age_display_labels[[as.character(age)]])
  }
  sex_label <- function(sex) {
    unname(subgroup_sex_display_labels[[as.character(sex)]])
  }

  # Prespecified pairwise age-slope contrasts: older adults versus each of the
  # two other mutually exclusive age groups, separately for men and women.
  age_contrast_specs <- rbindlist(lapply(subgroup_sex_ids, function(sex) {
    data.table(
      contrast_family = "Older vs other age group",
      lhs_age_id = 26L,
      lhs_sex_id = sex,
      rhs_age_id = c(158L, subgroup_derived_age_id),
      rhs_sex_id = sex,
      comparison_order = c(1L, 2L)
    )
  }))
  # Within-age sex contrasts for every mutually exclusive age group.
  sex_contrast_specs <- data.table(
    contrast_family = "Male vs female within age group",
    lhs_age_id = subgroup_interaction_age_ids,
    lhs_sex_id = 1L,
    rhs_age_id = subgroup_interaction_age_ids,
    rhs_sex_id = 2L,
    comparison_order = seq_along(subgroup_interaction_age_ids)
  )
  contrast_specs <- rbindlist(
    list(age_contrast_specs, sex_contrast_specs), use.names = TRUE
  )
  contrast_table <- rbindlist(lapply(seq_len(nrow(contrast_specs)), function(i) {
    spec <- contrast_specs[i]
    lhs_vector <- get_cell_slope(spec$lhs_age_id, spec$lhs_sex_id)
    rhs_vector <- get_cell_slope(spec$rhs_age_id, spec$rhs_sex_id)
    contrast_vector <- lhs_vector - rhs_vector
    inference <- cluster_linear_contrast(
      model, clustered_vcov, contrast_vector, df_den
    )
    lhs_slope <- sum(lhs_vector * coef(model))
    rhs_slope <- sum(rhs_vector * coef(model))
    comparison <- if (spec$contrast_family == "Older vs other age group") {
      paste0(
        age_label(spec$lhs_age_id), " vs ", age_label(spec$rhs_age_id),
        " within ", sex_label(spec$lhs_sex_id)
      )
    } else {
      paste0(
        sex_label(spec$lhs_sex_id), " vs ", sex_label(spec$rhs_sex_id),
        " within ", age_label(spec$lhs_age_id)
      )
    }
    data.table(
      measure_id = unique(dat$measure_id),
      outcome = unname(outcome_labels[[as.character(unique(dat$measure_id))]]),
      heat_variable = subgroup_heat_vars[[1L]],
      model_id = core_model_id,
      contrast_family = spec$contrast_family,
      comparison = comparison,
      lhs_age_id = spec$lhs_age_id,
      lhs_age_group = age_label(spec$lhs_age_id),
      lhs_sex_id = spec$lhs_sex_id,
      lhs_sex = sex_label(spec$lhs_sex_id),
      rhs_age_id = spec$rhs_age_id,
      rhs_age_group = age_label(spec$rhs_age_id),
      rhs_sex_id = spec$rhs_sex_id,
      rhs_sex = sex_label(spec$rhs_sex_id),
      lhs_heat_slope = lhs_slope,
      rhs_heat_slope = rhs_slope,
      slope_difference = inference$estimate,
      std_error_clustered = inference$std_error_clustered,
      ci_low = inference$ci_low,
      ci_high = inference$ci_high,
      t_value = inference$t_value,
      df_den = inference$df_den,
      p_value = inference$p_value,
      comparison_order = spec$comparison_order,
      n_observations = nobs(model),
      n_countries = uniqueN(dat$location_id),
      effect_scale = "Difference in rate change per heat IQR (lhs minus rhs)"
    )
  }))

  # A two-degree-of-freedom test asks whether the older-adult heat slope equals
  # both the <20 and the 20–70 slopes within each sex.
  older_joint_tests <- rbindlist(lapply(subgroup_sex_ids, function(sex) {
    older_vector <- get_cell_slope(26L, sex)
    under20_vector <- get_cell_slope(158L, sex)
    middle_vector <- get_cell_slope(subgroup_derived_age_id, sex)
    contrast_matrix <- rbind(
      older_minus_under20 = older_vector - under20_vector,
      older_minus_middle = older_vector - middle_vector
    )
    test <- cluster_joint_contrast(model, clustered_vcov, contrast_matrix, df_den)
    data.table(
      measure_id = unique(dat$measure_id),
      outcome = unname(outcome_labels[[as.character(unique(dat$measure_id))]]),
      heat_variable = subgroup_heat_vars[[1L]],
      model_id = core_model_id,
      sex_id = sex,
      sex = sex_label(sex),
      joint_test = paste0(">70 years vs <20 and 20–70 years within ", sex_label(sex)),
      older_minus_under20 = test$contrast_estimates[[1L]],
      older_minus_middle = test$contrast_estimates[[2L]],
      wald_chisq = test$wald_chisq,
      wald_f = test$wald_f,
      df_num = test$df_num,
      df_den = test$df_den,
      p_value = test$p_value,
      n_observations = nobs(model),
      n_countries = uniqueN(dat$location_id),
      effect_scale = "Joint clustered Wald F test of older-adult slope equality"
    )
  }))
  diagnostic <- data.table(
    analysis_type = "subgroup_age_sex_interaction_m5",
    measure_id = unique(dat$measure_id),
    heat_variable = subgroup_heat_vars[[1L]],
    model_id = core_model_id,
    family = "OLS",
    converged = TRUE,
    rank = model$rank,
    df_residual = df.residual(model),
    deviance = sum(residuals(model)^2),
    n_observations = nobs(model),
    n_countries = uniqueN(dat$location_id)
  )
  list(
    tests = test_table,
    contrasts = contrast_table,
    older_joint_tests = older_joint_tests,
    diagnostic = diagnostic
  )
}

message("Fitting formal age, sex, and age-by-sex interaction tests...")
subgroup_interaction_fits <- lapply(expected_measures, function(measure) {
  fit_subgroup_interaction_model(subgroup_data[
    measure_id == measure & age_id %in% subgroup_interaction_age_ids
  ])
})
subgroup_interaction_tests <- rbindlist(
  lapply(subgroup_interaction_fits, `[[`, "tests"), use.names = TRUE
)
subgroup_interaction_diagnostics <- rbindlist(
  lapply(subgroup_interaction_fits, `[[`, "diagnostic"), use.names = TRUE
)
subgroup_difference_contrasts <- rbindlist(
  lapply(subgroup_interaction_fits, `[[`, "contrasts"), use.names = TRUE
)
subgroup_older_joint_tests <- rbindlist(
  lapply(subgroup_interaction_fits, `[[`, "older_joint_tests"), use.names = TRUE
)
subgroup_interaction_tests[, p_adjust_bh_all := p.adjust(p_value, method = "BH")]
subgroup_interaction_tests[, p_adjust_bh_within_test := p.adjust(p_value, method = "BH"),
  by = interaction_test
]
subgroup_interaction_tests[, interaction_order := match(
  interaction_test, c("Age", "Sex", "Age × sex")
)]
setorder(subgroup_interaction_tests, measure_id, interaction_order)
subgroup_interaction_tests[, interaction_order := NULL]
subgroup_difference_contrasts[, p_adjust_bh_all := p.adjust(p_value, method = "BH")]
subgroup_difference_contrasts[, p_adjust_bh_within_family := p.adjust(
  p_value, method = "BH"
), by = contrast_family]
subgroup_difference_contrasts[, family_order := match(
  contrast_family,
  c("Older vs other age group", "Male vs female within age group")
)]
setorder(subgroup_difference_contrasts, measure_id, family_order, comparison_order)
subgroup_difference_contrasts[, c("family_order", "comparison_order") := NULL]
subgroup_older_joint_tests[, p_adjust_bh_all := p.adjust(p_value, method = "BH")]
subgroup_older_joint_tests[, p_adjust_bh_within_sex := p.adjust(
  p_value, method = "BH"
), by = sex_id]
setorder(subgroup_older_joint_tests, measure_id, sex_id)

# Prespecified core-M5 sensitivity analyses. All estimates retain the primary
# sample's heat centering and IQR so coefficients remain directly comparable.
sensitivity_definitions <- data.table(
  sensitivity_id = c("trim_1_99", "lag1", "lag01_mean", "exclude_2020"),
  sensitivity_label = c(
    "Trimmed 1st–99th percentiles", "Lag 1 year", "Mean of lag 0–1",
    "Excluding 2020 (COVID-19 shock)"
  ),
  sensitivity_order = 1:4,
  analysis_role = c(rep("Prespecified", 3L), "Reviewer-requested")
)

prepare_sensitivity_data <- function(dat, heat_var, sensitivity_id) {
  dat <- copy(dat)
  scale_row <- exposure_scaling[heat_variable == heat_var]
  trim_lower <- NA_real_
  trim_upper <- NA_real_

  if (sensitivity_id == "trim_1_99") {
    trim_lower <- unname(quantile(scaling_panel[[heat_var]], 0.01))
    trim_upper <- unname(quantile(scaling_panel[[heat_var]], 0.99))
    dat <- dat[get(heat_var) >= trim_lower & get(heat_var) <= trim_upper]
    dat[, sensitivity_heat := get(heat_var)]
    exposure_definition <- "Observed heat after excluding values outside the primary-sample 1st–99th percentiles"
  } else if (sensitivity_id == "lag1") {
    setorderv(dat, c("location_id", "year"))
    dat[, sensitivity_heat := shift(get(heat_var), type = "lag"), by = location_id]
    exposure_definition <- "Heat in the preceding calendar year assigned to the current outcome year"
  } else if (sensitivity_id == "lag01_mean") {
    setorderv(dat, c("location_id", "year"))
    dat[, lag_heat := shift(get(heat_var), type = "lag"), by = location_id]
    dat[, sensitivity_heat := (get(heat_var) + lag_heat) / 2]
    dat[, lag_heat := NULL]
    exposure_definition <- "Mean of current-year and preceding-year heat assigned to the current outcome year"
  } else if (sensitivity_id == "exclude_2020") {
    dat <- dat[year != 2020L]
    dat[, sensitivity_heat := get(heat_var)]
    exposure_definition <- "Observed current-year heat after excluding outcome year 2020 to remove the COVID-19 shock"
  } else {
    stop("Unknown sensitivity analysis: ", sensitivity_id, call. = FALSE)
  }

  dat <- dat[is.finite(sensitivity_heat)]
  dat[, heat_iqr := (sensitivity_heat - scale_row$center_median) / scale_row$iqr]
  attr(dat, "trim_lower") <- trim_lower
  attr(dat, "trim_upper") <- trim_upper
  attr(dat, "exposure_definition") <- exposure_definition
  dat
}

fit_sensitivity_model <- function(dat, heat_var, sensitivity_id) {
  dat <- prepare_sensitivity_data(dat, heat_var, sensitivity_id)
  model <- glm(
    as.formula(paste0(
      "val ~ ", model_rhs("heat_iqr", core_model_id),
      " + offset(log(population))"
    )),
    family = quasipoisson(link = "log"), data = dat,
    control = glm.control(maxit = 100L, epsilon = 1e-8)
  )
  inference <- cluster_inference(model, dat$location_id)
  scale_row <- exposure_scaling[heat_variable == heat_var]
  scenario_id <- sensitivity_id
  sensitivity_row <- sensitivity_definitions[sensitivity_id == scenario_id]
  result <- data.table(
    sensitivity_id = sensitivity_id,
    sensitivity_label = sensitivity_row$sensitivity_label,
    sensitivity_order = sensitivity_row$sensitivity_order,
    analysis_role = sensitivity_row$analysis_role,
    measure_id = unique(dat$measure_id),
    outcome = unname(outcome_labels[[as.character(unique(dat$measure_id))]]),
    heat_variable = heat_var,
    heat_label = unname(heat_labels[[heat_var]]),
    increment = scale_row$iqr,
    increment_unit = "Primary-sample IQR",
    model_id = core_model_id,
    family = "quasi-Poisson log",
    effect_scale = "log rate ratio per primary-sample IQR",
    exposure_definition = attr(dat, "exposure_definition"),
    trim_lower = attr(dat, "trim_lower"),
    trim_upper = attr(dat, "trim_upper"),
    beta = inference$beta,
    std_error_clustered = inference$se,
    ci_low = inference$ci_low,
    ci_high = inference$ci_high,
    rate_ratio = exp(inference$beta),
    rr_ci_low = exp(inference$ci_low),
    rr_ci_high = exp(inference$ci_high),
    p_value = inference$p_value,
    n_observations = nobs(model),
    n_countries = inference$n_clusters,
    analysis_year_first = min(dat$year),
    analysis_year_last = max(dat$year),
    inference_df = inference$inferential_df
  )
  diagnostic <- data.table(
    analysis_type = "sensitivity_m5",
    sensitivity_id = sensitivity_id,
    measure_id = unique(dat$measure_id),
    heat_variable = heat_var,
    model_id = core_model_id,
    family = "quasi-Poisson log",
    converged = isTRUE(model$converged),
    iterations = model$iter,
    rank = model$rank,
    df_residual = df.residual(model),
    dispersion = summary(model)$dispersion,
    deviance = deviance(model),
    n_observations = nobs(model),
    n_countries = inference$n_clusters
  )
  list(result = result, diagnostic = diagnostic)
}

message("Fitting four core-M5 sensitivity analyses, including exclusion of 2020...")
sensitivity_fits <- list()
sensitivity_index <- 0L
for (sensitivity_id in sensitivity_definitions$sensitivity_id) {
  for (measure in expected_measures) {
    dat <- analysis_data[measure_id == measure]
    for (heat_var in heat_vars) {
      sensitivity_index <- sensitivity_index + 1L
      sensitivity_fits[[sensitivity_index]] <- fit_sensitivity_model(
        dat, heat_var, sensitivity_id
      )
    }
  }
}
sensitivity_results <- rbindlist(
  lapply(sensitivity_fits, `[[`, "result"), use.names = TRUE
)
sensitivity_diagnostics <- rbindlist(
  lapply(sensitivity_fits, `[[`, "diagnostic"), use.names = TRUE
)
sensitivity_results[, p_adjust_bh := p.adjust(p_value, method = "BH"), by = sensitivity_id]
setorder(sensitivity_results, sensitivity_order, heat_variable, measure_id)
sensitivity_summary <- sensitivity_results[, .(
  tests = .N,
  raw_p_lt_0_05 = sum(p_value < 0.05),
  bh_q_lt_0_05 = sum(p_adjust_bh < 0.05),
  n_observations_min = min(n_observations),
  n_observations_max = max(n_observations),
  n_countries_min = min(n_countries),
  n_countries_max = max(n_countries),
  first_year = min(analysis_year_first),
  last_year = max(analysis_year_last)
), by = .(sensitivity_id, sensitivity_label, sensitivity_order, analysis_role)]
setorder(sensitivity_summary, sensitivity_order)
diagnostics <- rbindlist(
  list(
    diagnostics, country_year_twfe_diagnostics,
    cdnch_rch_joint_diagnostics, sensitivity_diagnostics
  ),
  use.names = TRUE, fill = TRUE
)

# Nonlinear compound day-night heat curves: 4 df natural spline.
standard_nonlinear_rr_definition <- paste0(
  "RR(h)=exp{f(h)-f(h_ref)}; h_ref is the global median of CDNCH"
)
standard_nonlinear_rr_source <- paste0(
  "Same outcome-specific core-M5 4-df natural-spline exposure-response ",
  "function used in Figure 3"
)

# Propagate uncertainty in the fitted nonlinear exposure-response function to
# cumulative attributable burden. Draws use the country-clustered covariance
# matrix of the four spline coefficients while keeping the observed GBD point
# burden fixed. GBD lower/upper inputs remain a separate boundary sensitivity.
attributable_uncertainty_draws <- 10000L
attributable_uncertainty_seed_base <- 20260919L

simulate_attributable_burden_ci <- function(
    basis_difference, beta, covariance, observed_burden,
    n_draws, seed, chunk_size = 500L) {
  covariance <- (covariance + t(covariance)) / 2
  covariance_eigen <- eigen(covariance, symmetric = TRUE)
  covariance_tolerance <- max(abs(covariance_eigen$values)) * 1e-10
  if (min(covariance_eigen$values) < -covariance_tolerance) {
    stop("Spline covariance matrix is not positive semidefinite.", call. = FALSE)
  }
  covariance_root <- covariance_eigen$vectors %*%
    diag(sqrt(pmax(covariance_eigen$values, 0)), nrow = length(beta))

  set.seed(seed)
  normal_draws <- matrix(
    rnorm(length(beta) * n_draws), nrow = length(beta), ncol = n_draws
  )
  coefficient_draws <- sweep(
    covariance_root %*% normal_draws, 1L, beta, FUN = "+"
  )
  burden_draws <- numeric(n_draws)
  for (draw_start in seq.int(1L, n_draws, by = chunk_size)) {
    draw_index <- draw_start:min(draw_start + chunk_size - 1L, n_draws)
    log_rr_draw <- basis_difference %*%
      coefficient_draws[, draw_index, drop = FALSE]
    # Put the matrix first so pmax preserves its dimensions for row-weighting.
    attributable_fraction_draw <- pmax(1 - exp(-log_rr_draw), 0)
    burden_draws[draw_index] <- colSums(sweep(
      attributable_fraction_draw, 1L, observed_burden, FUN = "*"
    ))
  }
  if (any(!is.finite(burden_draws)) || any(burden_draws < 0)) {
    stop("Invalid attributable-burden simulation draw.", call. = FALSE)
  }
  as.numeric(quantile(
    burden_draws, probs = c(0.025, 0.975), names = FALSE, type = 8
  ))
}

spline_basis_contrast <- function(basis_template, evaluation_iqr, reference_iqr) {
  evaluated_basis <- as.matrix(predict(basis_template, evaluation_iqr))
  reference_basis <- as.numeric(predict(basis_template, reference_iqr))
  sweep(evaluated_basis, 2L, reference_basis, FUN = "-")
}

fit_compound_spline <- function(dat, model_id, linear_p_value) {
  heat_var <- "compound_day_night_cumulative_heat"
  dat <- scale_heat(dat, heat_var)
  observed_basis <- ns(dat$heat_iqr, df = 4)
  basis_names <- paste0("heat_ns", seq_len(ncol(observed_basis)))
  dat[, (basis_names) := as.data.table(observed_basis)]
  model <- glm(
    as.formula(paste0(
      "val ~ ", model_rhs(paste(basis_names, collapse = " + "), model_id),
      " + offset(log(population))"
    )),
    family = quasipoisson(link = "log"), data = dat,
    control = glm.control(maxit = 100L, epsilon = 1e-8)
  )
  vc <- cluster_vcov(model, dat$location_id)
  b <- coef(model)[basis_names]
  v_sub <- vc[basis_names, basis_names, drop = FALSE]
  wald_stat <- as.numeric(t(b) %*% solve(v_sub, b))
  overall_p <- pchisq(wald_stat, df = length(basis_names), lower.tail = FALSE)

  scale_row <- exposure_scaling[heat_variable == heat_var]
  lower_heat <- unname(quantile(dat[[heat_var]], 0.01))
  upper_heat <- unname(quantile(dat[[heat_var]], 0.99))
  # CDNCH is an anomalous-heat metric, so use the global median rather than
  # zero exposure as the RR reference.
  reference_heat <- scale_row$center_median
  heat_original <- sort(unique(c(
    seq(lower_heat, upper_heat, length.out = 199L), reference_heat
  )))
  grid_scaled <- (heat_original - scale_row$center_median) / scale_row$iqr
  reference_scaled <- (reference_heat - scale_row$center_median) / scale_row$iqr
  grid_basis_difference <- spline_basis_contrast(
    observed_basis, grid_scaled, reference_scaled
  )
  critical_value <- qt(0.975, df = uniqueN(dat$location_id) - 1L)

  log_rr <- as.numeric(grid_basis_difference %*% b)
  log_rr_variance <- rowSums(
    (grid_basis_difference %*% v_sub) * grid_basis_difference
  )
  log_rr_se <- sqrt(pmax(0, log_rr_variance))
  nonlinear_curve_rr <- exp(log_rr)
  reference_index <- which.min(abs(heat_original - reference_heat))
  if (abs(nonlinear_curve_rr[reference_index] - 1) > 1e-10) {
    stop("Standard nonlinear exposure-response RR is not 1 at the median reference exposure.", call. = FALSE)
  }

  # Core-M5 attributable burden for all four prespecified count outcomes. The
  # RR is evaluated from the exact same standard nonlinear exposure-response
  # function used for Figure 3. The global median CDNCH (RR=1) is the
  # counterfactual and
  # RR <= 1 is assigned zero attributable burden. Each
  # outcome retains its native burden unit.
  measure_id_value <- unique(dat$measure_id)
  attribution_basis_difference <- spline_basis_contrast(
    observed_basis, dat$heat_iqr, reference_scaled
  )
  nonlinear_attribution_rr <- exp(as.numeric(attribution_basis_difference %*% b))
  attributable_fraction <- pmax(
    0, (nonlinear_attribution_rr - 1) / nonlinear_attribution_rr
  )
  attributable_burden_ci <- simulate_attributable_burden_ci(
    basis_difference = attribution_basis_difference,
    beta = b,
    covariance = v_sub,
    observed_burden = dat$val,
    n_draws = attributable_uncertainty_draws,
    seed = attributable_uncertainty_seed_base + measure_id_value
  )
  attributable_model_uncertainty <- data.table(
    measure_id = measure_id_value,
    outcome = unname(outcome_labels[[as.character(measure_id_value)]]),
    burden_unit = unname(outcome_burden_units[[as.character(measure_id_value)]]),
    model_id = model_id,
    heat_variable = heat_var,
    cumulative_attributable_burden = sum(dat$val * attributable_fraction),
    cumulative_attributable_burden_model_ci_low = attributable_burden_ci[[1L]],
    cumulative_attributable_burden_model_ci_high = attributable_burden_ci[[2L]],
    cumulative_attributable_fraction =
      sum(dat$val * attributable_fraction) / sum(dat$val),
    cumulative_attributable_fraction_model_ci_low =
      attributable_burden_ci[[1L]] / sum(dat$val),
    cumulative_attributable_fraction_model_ci_high =
      attributable_burden_ci[[2L]] / sum(dat$val),
    simulation_draws = attributable_uncertainty_draws,
    simulation_seed = attributable_uncertainty_seed_base + measure_id_value,
    uncertainty_method = paste0(
      "Parametric multivariate-normal simulation of the four spline ",
      "coefficients using their country-clustered covariance matrix; ",
      "GBD point burden held fixed"
    ),
    uncertainty_scope = paste0(
      "Model-based 95% CI only; excludes GBD burden-input uncertainty and ",
      "counterfactual/specification uncertainty"
    )
  )
  attributable <- dat[, .(
    measure_id = measure_id_value,
    outcome = unname(outcome_labels[[as.character(measure_id_value)]]),
    burden_unit = unname(outcome_burden_units[[as.character(measure_id_value)]]),
    location_id, location, year, superregion, region_wb,
    exposure = get(heat_var),
    observed_burden_lower_input = lower,
    observed_burden = val,
    observed_burden_upper_input = upper
  )]
  attributable[, `:=`(
    model_id = model_id,
    heat_variable = heat_var,
    reference_heat = reference_heat,
    reference_heat_unit = unname(heat_units[[heat_var]]),
    relative_risk = nonlinear_attribution_rr,
    rr_definition = standard_nonlinear_rr_definition,
    rr_source = standard_nonlinear_rr_source,
    attributable_fraction = attributable_fraction,
    attributable_burden_lower_input = observed_burden_lower_input * attributable_fraction,
    attributable_burden = observed_burden * attributable_fraction,
    attributable_burden_upper_input = observed_burden_upper_input * attributable_fraction,
    counterfactual_definition = paste0(
      "Global median CDNCH (",
      formatC(reference_heat, format = "f", digits = 2),
      " °C·days) where RR(h)=1"
    ),
    safe_range_definition = "RR <= 1 assigned attributable fraction 0",
    input_boundary_definition = paste0(
      "GBD row-level lower/upper burden inputs substituted while holding ",
      "the point-estimate fitted spline-RR function fixed"
    )
  )]

  curve <- data.table(
    measure_id = unique(dat$measure_id),
    outcome = unname(outcome_labels[[as.character(unique(dat$measure_id))]]),
    model_id = model_id, heat_variable = heat_var,
    heat_unit = unname(heat_units[[heat_var]]), heat_value = heat_original,
    heat_iqr = grid_scaled,
    reference_heat = reference_heat,
    reference_heat_iqr = reference_scaled,
    log_rr = log_rr,
    relative_risk = nonlinear_curve_rr,
    rr_ci_low = exp(log_rr - critical_value * log_rr_se),
    rr_ci_high = exp(log_rr + critical_value * log_rr_se),
    rr_definition = standard_nonlinear_rr_definition,
    rr_source = standard_nonlinear_rr_source,
    effect_scale = "Standard nonlinear exposure-response rate ratio",
    linear_p_value = linear_p_value,
    linear_significant = linear_p_value < 0.05,
    spline_overall_p = overall_p,
    prediction_range = "Heat at observed 1st through 99th percentiles"
  )
  test <- data.table(
    measure_id = unique(dat$measure_id),
    outcome = unname(outcome_labels[[as.character(unique(dat$measure_id))]]),
    model_id = model_id, linear_p_value = linear_p_value,
    linear_significant = linear_p_value < 0.05, spline_df = 4L,
    spline_wald_chisq = wald_stat, spline_overall_p = overall_p,
    reference_heat = reference_heat,
    reference_heat_unit = unname(heat_units[[heat_var]]),
    rr_definition = standard_nonlinear_rr_definition,
    rr_source = standard_nonlinear_rr_source,
    effect_scale = "Standard nonlinear exposure-response rate ratio",
    n_observations = nobs(model), n_countries = uniqueN(dat$location_id),
    converged = isTRUE(model$converged)
  )
  diagnostic <- data.table(
    analysis_type = "nonlinear_ns4", measure_id = unique(dat$measure_id),
    heat_variable = heat_var, model_id = model_id, family = "quasi-Poisson log",
    converged = isTRUE(model$converged), iterations = model$iter,
    rank = model$rank, df_residual = df.residual(model),
    dispersion = summary(model)$dispersion, deviance = deviance(model),
    n_observations = nobs(model), n_countries = uniqueN(dat$location_id)
  )
  list(
    curve = curve, test = test, diagnostic = diagnostic,
    attributable = attributable,
    attributable_model_uncertainty = attributable_model_uncertainty
  )
}

message("Fitting core M5 4 df nonlinear compound day-night heat models...")
nonlinear_fits <- list()
nonlinear_index <- 0L
for (measure in expected_measures) {
  dat <- analysis_data[measure_id == measure]
  for (model_id in nonlinear_model_ids) {
    nonlinear_index <- nonlinear_index + 1L
    current_model <- model_id
    linear_p <- results[
      measure_id == measure &
        heat_variable == "compound_day_night_cumulative_heat" &
        model_id == current_model,
      p_value
    ]
    nonlinear_fits[[nonlinear_index]] <- fit_compound_spline(dat, model_id, linear_p)
  }
}
nonlinear_curves <- rbindlist(lapply(nonlinear_fits, `[[`, "curve"), use.names = TRUE)
nonlinear_tests <- rbindlist(lapply(nonlinear_fits, `[[`, "test"), use.names = TRUE)
nonlinear_diagnostics <- rbindlist(lapply(nonlinear_fits, `[[`, "diagnostic"), use.names = TRUE)
diagnostics <- rbindlist(list(diagnostics, nonlinear_diagnostics), use.names = TRUE, fill = TRUE)
attributable_burden <- rbindlist(
  Filter(function(x) nrow(x) > 0L, lapply(nonlinear_fits, `[[`, "attributable")),
  use.names = TRUE
)
attributable_model_uncertainty <- rbindlist(
  lapply(nonlinear_fits, `[[`, "attributable_model_uncertainty"),
  use.names = TRUE
)
if (nrow(attributable_burden) != 24552L ||
    attributable_burden[, any(!is.finite(attributable_burden) | attributable_burden < 0)] ||
    attributable_burden[, any(
      !is.finite(attributable_burden_lower_input) |
        !is.finite(attributable_burden_upper_input) |
        attributable_burden_lower_input < 0 |
        attributable_burden_upper_input < attributable_burden_lower_input
    )] ||
    any(attributable_burden[, .N, by = measure_id]$N != 6138L)) {
  stop("Core-M5 four-outcome attributable-burden output is incomplete or invalid.", call. = FALSE)
}
attributable_burden_summary <- attributable_burden[, .(
  countries = uniqueN(location_id),
  observed_burden_lower_input = sum(observed_burden_lower_input),
  observed_burden = sum(observed_burden),
  observed_burden_upper_input = sum(observed_burden_upper_input),
  attributable_burden_lower_input = sum(attributable_burden_lower_input),
  attributable_burden = sum(attributable_burden),
  attributable_burden_upper_input = sum(attributable_burden_upper_input),
  attributable_fraction_lower_input = sum(attributable_burden_lower_input) /
    sum(observed_burden_lower_input),
  attributable_fraction = sum(attributable_burden) / sum(observed_burden),
  attributable_fraction_upper_input = sum(attributable_burden_upper_input) /
    sum(observed_burden_upper_input)
), by = .(
  measure_id, outcome, burden_unit, year, model_id, heat_variable,
  reference_heat, reference_heat_unit, rr_definition, rr_source
)]
setorder(attributable_burden_summary, measure_id, year)

# Table 2 source: one overall row for each prespecified burden outcome.
attributable_burden_overall_summary <- attributable_burden_summary[, .(
  model_id = unique(model_id),
  heat_variable = unique(heat_variable),
  reference_heat = unique(reference_heat),
  reference_heat_unit = unique(reference_heat_unit),
  rr_definition = unique(rr_definition),
  rr_source = unique(rr_source),
  first_year = min(year),
  last_year = max(year),
  countries = max(countries),
  cumulative_observed_burden_lower_input = sum(observed_burden_lower_input),
  cumulative_observed_burden = sum(observed_burden),
  cumulative_observed_burden_upper_input = sum(observed_burden_upper_input),
  cumulative_attributable_burden_lower_input = sum(attributable_burden_lower_input),
  cumulative_attributable_burden = sum(attributable_burden),
  cumulative_attributable_burden_upper_input = sum(attributable_burden_upper_input),
  cumulative_attributable_fraction_lower_input =
    sum(attributable_burden_lower_input) / sum(observed_burden_lower_input),
  cumulative_attributable_fraction = sum(attributable_burden) / sum(observed_burden),
  cumulative_attributable_fraction_upper_input =
    sum(attributable_burden_upper_input) / sum(observed_burden_upper_input),
  annual_attributable_fraction_min = min(attributable_fraction),
  annual_attributable_fraction_max = max(attributable_fraction)
), by = .(measure_id, outcome, burden_unit)]
attributable_burden_overall_summary <- merge(
  attributable_burden_overall_summary,
  attributable_model_uncertainty[, .(
    measure_id,
    cumulative_attributable_burden_model_ci_low,
    cumulative_attributable_burden_model_ci_high,
    cumulative_attributable_fraction_model_ci_low,
    cumulative_attributable_fraction_model_ci_high,
    model_uncertainty_draws = simulation_draws,
    model_uncertainty_seed = simulation_seed,
    model_uncertainty_method = uncertainty_method,
    model_uncertainty_scope = uncertainty_scope
  )],
  by = "measure_id", all.x = TRUE, sort = FALSE
)
if (attributable_burden_overall_summary[, any(
  cumulative_attributable_burden_model_ci_low > cumulative_attributable_burden |
    cumulative_attributable_burden_model_ci_high < cumulative_attributable_burden |
    !is.finite(cumulative_attributable_burden_model_ci_low) |
    !is.finite(cumulative_attributable_burden_model_ci_high)
)]) {
  stop("Invalid model-based attributable-burden confidence interval.", call. = FALSE)
}
attributable_burden_overall_summary[, outcome_order := match(measure_id, expected_measures)]
setorder(attributable_burden_overall_summary, outcome_order)
attributable_burden_overall_summary[, outcome_order := NULL]

# Backward-compatible mortality-only files retain their original column names.
attributable_deaths <- attributable_burden[measure_id == 1L, .(
  location_id, location, year, superregion, region_wb, exposure,
  observed_deaths_lower_input = observed_burden_lower_input,
  observed_deaths = observed_burden, model_id, heat_variable, reference_heat,
  observed_deaths_upper_input = observed_burden_upper_input,
  reference_heat_unit, relative_risk, rr_definition, rr_source,
  attributable_fraction,
  attributable_deaths_lower_input = attributable_burden_lower_input,
  attributable_deaths = attributable_burden, counterfactual_definition,
  attributable_deaths_upper_input = attributable_burden_upper_input,
  safe_range_definition, input_boundary_definition
)]
attributable_deaths_annual_summary <- attributable_burden_summary[measure_id == 1L, .(
  year, model_id, heat_variable, reference_heat, reference_heat_unit, countries,
  rr_definition, rr_source,
  observed_deaths_lower_input = observed_burden_lower_input,
  observed_deaths = observed_burden, attributable_deaths = attributable_burden,
  observed_deaths_upper_input = observed_burden_upper_input,
  attributable_deaths_lower_input = attributable_burden_lower_input,
  attributable_deaths_upper_input = attributable_burden_upper_input,
  attributable_fraction_lower_input, attributable_fraction,
  attributable_fraction_upper_input
)]

# Source data for 1990–2020 mean all-age, both-sex outcome-rate maps.
map_outcome_labels <- c(
  `1` = "Death", `2` = "DALYs", `5` = "Prevalence", `6` = "Incidence"
)
outcome_map_data <- panel[
  metric_id == 3L & age_id == 22L & sex_id == 3L &
    year >= 1990L & year <= 2020L &
    measure_id %in% expected_measures & is.finite(val) & val >= 0,
  .(
    value = mean(val),
    years_observed = uniqueN(year),
    first_year = min(year),
    last_year = max(year)
  ),
  by = .(measure_id, location_id, location, fid, iso_a3)
]
outcome_map_data[, `:=`(
  outcome = unname(map_outcome_labels[as.character(measure_id)]),
  metric_id = 3L,
  metric_name = "Rate",
  age_id = 22L,
  age_name = "All ages",
  sex_id = 3L,
  sex_name = "Both sexes",
  period = "1990–2020 mean",
  unit = "All-age, both-sex rate per 100,000 population"
)]
setcolorder(outcome_map_data, c(
  "measure_id", "outcome", "metric_id", "metric_name", "age_id", "age_name",
  "sex_id", "sex_name", "location_id", "location", "fid", "iso_a3",
  "period", "first_year", "last_year", "years_observed", "value", "unit"
))
if (outcome_map_data[, anyDuplicated(paste(measure_id, location_id, sep = ":"))] > 0L) {
  stop("Duplicate country-outcome rows in the 1990–2020 mean outcome map source.", call. = FALSE)
}
if (outcome_map_data[, any(years_observed != 31L)]) {
  stop("Outcome map means do not all contain 31 annual observations.", call. = FALSE)
}

# Mean 1990–2020 heat indicators for the main-model country sample.
heat_map_vars <- c(
  "daytime_cumulative_heat", "nighttime_cumulative_heat",
  "compound_day_night_cumulative_heat"
)
map_heat_labels <- c(
  daytime_cumulative_heat = "Daytime cumulative heat",
  nighttime_cumulative_heat = "Nighttime cumulative heat",
  compound_day_night_cumulative_heat = "Compound day–night cumulative heat"
)
heat_panel <- unique(
  analysis_data[, c(
    "location_id", "location", "fid", "iso_a3", "year",
    heat_map_vars
  ), with = FALSE],
  by = c("location_id", "year")
)
indicator_map_data <- rbindlist(lapply(heat_map_vars, function(v) {
  valid <- heat_panel[year >= 1990L & year <= 2020L & is.finite(get(v))]
  means <- valid[, .(
    value = mean(get(v)),
    years_observed = uniqueN(year),
    first_year = min(year),
    last_year = max(year)
  ), by = .(location_id, location, fid, iso_a3)]
  means[, .(
    location_id, location, fid, iso_a3,
    period = "1990–2020 mean",
    first_year, last_year, years_observed,
    indicator = v,
    indicator_label = unname(map_heat_labels[[v]]),
    value,
    unit = unname(heat_units[[v]])
  )]
}), use.names = TRUE)
if (indicator_map_data[, any(years_observed != 31L)]) {
  stop("Heat-indicator map means do not all contain 31 annual observations.", call. = FALSE)
}

format_number <- function(x) {
  ifelse(
    abs(x) >= 100, formatC(x, format = "f", digits = 1),
    ifelse(abs(x) >= 10, formatC(x, format = "f", digits = 2),
      ifelse(abs(x) < 0.001 & x != 0, formatC(x, format = "e", digits = 2),
             formatC(x, format = "f", digits = 3)))
  )
}
format_p <- function(x) ifelse(x < 0.001, "<0.001", formatC(x, format = "f", digits = 3))
format_percent <- function(x, digits = 2L) paste0(formatC(100 * x, format = "f", digits = digits), "%")

markdown_subgroup_m5_results <- function(x) {
  rows <- x[, .(
    subgroup = paste0(
      unname(subgroup_age_display_labels[as.character(age_id)]), " × ",
      unname(subgroup_sex_display_labels[as.character(sex_id)])
    ),
    outcome = unname(subgroup_outcome_display_labels[as.character(measure_id)]),
    heat = unname(subgroup_heat_display_labels[heat_variable]),
    beta = format_number(beta),
    ci = paste0("[", format_number(ci_low), ", ", format_number(ci_high), "]"),
    p = format_p(p_value),
    q = format_p(p_adjust_bh),
    n = format(n_observations, big.mark = ","),
    countries = as.character(n_countries)
  )]
  c(
    "| Age group × sex | Outcome | Heat exposure | Rate difference per IQR | 95% CI | P value | BH q | Observations | Countries |",
    "|---|---|---|---:|---:|---:|---:|---:|---:|",
    apply(rows, 1L, function(row) paste0("| ", paste(row, collapse = " | "), " |")),
    "",
    paste0(
      "All models use OLS for rates, the core M5 covariates and superregion-specific 3 df natural-spline time trends. ",
      "Population is not used as an offset. The 20–70-year rate is derived as (all-age count - <20 count - >70 count) / ",
      "(all-age implied population - <20 implied population - >70 implied population) × 100,000. ",
      "Only compound day-night cumulative heat is analysed; therefore the BH q value equals the P value within each subgroup."
    )
  )
}

markdown_subgroup_m5_summary <- function(x) {
  rows <- x[, .(
    subgroup = paste0(
      unname(subgroup_age_display_labels[as.character(age_id)]), " × ",
      unname(subgroup_sex_display_labels[as.character(sex_id)])
    ),
    outcome = unname(subgroup_outcome_display_labels[as.character(measure_id)]),
    tests = as.character(tests),
    raw_p = as.character(raw_p_lt_0_05),
    bh_q = as.character(bh_q_lt_0_05),
    n = format(n_observations, big.mark = ","),
    countries = as.character(n_countries),
    years = paste0(first_year, "–", last_year)
  )]
  c(
    "| Age group × sex | Outcome | Tests | Raw P<0.05 | BH q<0.05 | Observations/model | Countries | Years |",
    "|---|---|---:|---:|---:|---:|---:|---|",
    apply(rows, 1L, function(row) paste0("| ", paste(row, collapse = " | "), " |")),
    "",
    "Rate-outcome subgroup analysis of compound day-night cumulative heat: male and female strata within the mutually exclusive <20-year, derived 20–70-year, and >70-year groups."
  )
}

markdown_subgroup_interaction_tests <- function(x) {
  rows <- x[, .(
    outcome = unname(subgroup_outcome_display_labels[as.character(measure_id)]),
    test = interaction_test,
    f_test = paste0(
      formatC(wald_f, format = "f", digits = 3), " (",
      df_num, ", ", df_den, ")"
    ),
    p = format_p(p_value),
    q_all = format_p(p_adjust_bh_all),
    q_family = format_p(p_adjust_bh_within_test),
    n = format(n_observations, big.mark = ","),
    countries = as.character(n_countries)
  )]
  c(
    "| Outcome | Interaction test | Clustered F (df1, df2) | P value | BH q (all 12) | BH q (within test) | Observations | Countries |",
    "|---|---|---:|---:|---:|---:|---:|---:|",
    apply(rows, 1L, function(row) paste0("| ", paste(row, collapse = " | "), " |")),
    "",
    paste0(
      "Formal effect-modification tests use one fully stratified stacked OLS core-M5 model per outcome. ",
      "Heat, covariate, and superregion-specific time-trend coefficients are allowed to vary by age-sex stratum, ",
      "so every cell heat slope reproduces its separately fitted subgroup M5 estimate. Age tests jointly compare ",
      "age slopes within sex; sex tests jointly compare male and female slopes within age; age × sex tests whether ",
      "the sex difference varies by age. Age groups are mutually exclusive: ",
      "<20 years, derived 20–70 years, and >70 years. ",
      "Inference uses country-clustered HC1 covariance and F reference degrees of freedom based on ",
      unique(x$n_countries), " countries."
    )
  )
}

markdown_subgroup_difference_contrasts <- function(x) {
  rows <- x[, .(
    outcome = unname(subgroup_outcome_display_labels[as.character(measure_id)]),
    contrast = comparison,
    lhs_slope = format_number(lhs_heat_slope),
    rhs_slope = format_number(rhs_heat_slope),
    difference_ci = paste0(
      format_number(slope_difference), " [",
      format_number(ci_low), ", ", format_number(ci_high), "]"
    ),
    p = format_p(p_value),
    q_all = format_p(p_adjust_bh_all),
    q_family = format_p(p_adjust_bh_within_family),
    n = format(n_observations, big.mark = ","),
    countries = as.character(n_countries)
  )]
  c(
    "| Outcome | Formal contrast | LHS heat slope | RHS heat slope | Slope difference [95% CI] | P value | BH q (all 28) | BH q (within family) | Observations | Countries |",
    "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    apply(rows, 1L, function(row) paste0("| ", paste(row, collapse = " | "), " |")),
    "",
    paste0(
      "Contrasts come from the same fully stratified stacked OLS core-M5 model used for formal effect-modification tests. ",
      "Each age-sex stratum has its own heat, covariate, and superregion-specific time-trend coefficients. ",
      "Older-age contrasts are >70 years minus <20 years or derived 20–70 years within sex. ",
      "Sex contrasts are male minus female within age group. A positive difference means the LHS heat slope is larger. ",
      "Inference uses country-clustered HC1 covariance with ", unique(x$df_den),
      " reference degrees of freedom."
    )
  )
}

markdown_subgroup_older_joint_tests <- function(x) {
  rows <- x[, .(
    outcome = unname(subgroup_outcome_display_labels[as.character(measure_id)]),
    sex = unname(subgroup_sex_display_labels[as.character(sex_id)]),
    old_minus_under20 = format_number(older_minus_under20),
    old_minus_middle = format_number(older_minus_middle),
    f_test = paste0(
      formatC(wald_f, format = "f", digits = 3), " (",
      df_num, ", ", df_den, ")"
    ),
    p = format_p(p_value),
    q_all = format_p(p_adjust_bh_all),
    q_sex = format_p(p_adjust_bh_within_sex),
    n = format(n_observations, big.mark = ","),
    countries = as.character(n_countries)
  )]
  c(
    "| Outcome | Sex | >70 minus <20 | >70 minus 20–70 | Clustered F (df1, df2) | P value | BH q (all 8) | BH q (within sex) | Observations | Countries |",
    "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    apply(rows, 1L, function(row) paste0("| ", paste(row, collapse = " | "), " |")),
    "",
    paste0(
      "Each two-degree-of-freedom Wald test evaluates whether the >70-year heat slope equals both other mutually exclusive age-group slopes within one sex. ",
      "The two displayed differences are >70 years minus <20 years and >70 years minus derived 20–70 years. ",
      "Inference uses country-clustered HC1 covariance with 200 reference degrees of freedom."
    )
  )
}

markdown_table <- function(x, note) {
  rows <- x[, .(
    comparison = paste0(heat_label, " × ", outcome),
    coefficient = format_number(beta),
    ci = paste0("[", format_number(ci_low), ", ", format_number(ci_high), "]"),
    p = format_p(p_value), model = model_id
  )]
  c(
    "| 热指标 × 结局 | 系数 β（每 IQR） | 95% CI | P 值 | 模型编号 |",
    "|---|---:|---:|---:|---|",
    paste0("| ", rows$comparison, " | ", rows$coefficient, " | ", rows$ci,
           " | ", rows$p, " | ", rows$model, " |"),
    "", note, ""
  )
}

write_markdown_utf8 <- function(lines, path) {
  fwrite(data.table(line = lines), path, sep = "\t", col.names = FALSE, quote = FALSE, bom = FALSE)
}

table_note <- paste0(
  "注：M0、M1、M2、M3a、M3b、M4和M5均使用相同的计数结局样本和ln(population) offset；",
  "β是每增加1个IQR的对数率比；每模型n=6,138，包含198个国家。",
  "除M0外，比较模型均调整7个预先规定协变量；M5为核心连续DID。",
  "95% CI和P值使用国家聚类稳健标准误（HC1）；BH校正在每个模型的32项检验内进行。"
)
core_table_note <- paste0(
  "注：M5（核心连续DID）使用计数结局样本和ln(population) offset，",
  "β是每增加1个IQR的对数率比；n=6,138，包含198个国家。",
  "模型调整7个预先规定协变量和超级区域特异3 df年份自然样条。",
  "95% CI和P值使用国家聚类稳健标准误（HC1）。"
)

markdown_descriptive_table <- function(x) {
  columns <- c("superregion", "n_countries", "n_country_years", descriptive_vars)
  headers <- c(
    "GBD superregion", "Countries", "Country-years",
    "Daytime heat (°C·days)", "Nighttime heat (°C·days)",
    "Compound heat (°C·days)", unname(covariate_labels[covariate_vars])
  )
  formatted <- copy(x)[, ..columns]
  formatted[, n_countries := as.character(n_countries)]
  formatted[, n_country_years := as.character(n_country_years)]
  numeric_columns <- setdiff(columns, c("superregion", "n_countries", "n_country_years"))
  formatted[, (numeric_columns) := lapply(.SD, format_number), .SDcols = numeric_columns]
  c(
    paste0("| ", paste(headers, collapse = " | "), " |"),
    paste0("|", paste(rep("---:", length(headers)), collapse = "|"), "|"),
    apply(formatted, 1L, function(row) paste0("| ", paste(row, collapse = " | "), " |")),
    "",
    "Note: Values are study-period medians calculated from the main complete-case country-year panel."
  )
}

markdown_sensitivity_summary <- function(x) {
  rows <- x[, .(
    sensitivity = sensitivity_label,
    tests = as.character(tests),
    raw_p = as.character(raw_p_lt_0_05),
    bh_q = as.character(bh_q_lt_0_05),
    observations = ifelse(
      n_observations_min == n_observations_max,
      format(n_observations_min, big.mark = ","),
      paste0(format(n_observations_min, big.mark = ","), "–", format(n_observations_max, big.mark = ","))
    ),
    countries = ifelse(
      n_countries_min == n_countries_max,
      as.character(n_countries_min),
      paste0(n_countries_min, "–", n_countries_max)
    ),
    outcome_years = paste0(first_year, "–", last_year)
  )]
  c(
    "| Sensitivity analysis | Tests | Raw P<0.05 | BH q<0.05 | Observations/model | Countries | Outcome years |",
    "|---|---:|---:|---:|---:|---:|---|",
    apply(rows, 1L, function(row) paste0("| ", paste(row, collapse = " | "), " |")),
    "",
    "All models use core M5 adjustment and the primary-sample heat IQR."
  )
}

markdown_attributable_summary <- function(x) {
  get_row <- function(measure) x[measure_id == measure]
  format_fraction <- function(measure) {
    row <- get_row(measure)
    format_percent(row$cumulative_attributable_fraction)
  }
  format_burden <- function(measure) {
    row <- get_row(measure)
    paste0(
      format(round(row$cumulative_attributable_burden), big.mark = ",", scientific = FALSE),
      " ", row$burden_unit
    )
  }
  format_range <- function(measure) {
    row <- get_row(measure)
    paste0(
      format_percent(row$annual_attributable_fraction_min), "–",
      format_percent(row$annual_attributable_fraction_max)
    )
  }
  rows <- data.table(
    statistic = c("累计归因分数", "累计归因负担（估计）", "年度归因分数范围（1990–2020）"),
    death = c(format_fraction(1L), format_burden(1L), format_range(1L)),
    incidence = c(format_fraction(6L), format_burden(6L), format_range(6L)),
    dalys = c(format_fraction(2L), format_burden(2L), format_range(2L)),
    prevalence = c(format_fraction(5L), format_burden(5L), format_range(5L))
  )
  c(
    "| 归因指标 | 死亡 | 发病 | DALYs | 患病 |",
    "|---|---:|---|---|---|",
    apply(rows, 1L, function(row) paste0("| ", paste(row, collapse = " | "), " |")),
    "",
    paste0(
      "注：四种结局均基于各自核心 M5 的 CDNCH 4 df 非线性 RR；以全球中位数的 CDNCH 为反事实，",
      "而非零暴露（CDNCH 是异常热指标）；RR≤1 时归因分数置零。累计归因分数按累计归因负担/累计观察负担计算；归因负担保留各结局原始单位。"
    )
  )
}

# Prespecified presentation tables. These objects contain exactly the columns
# and ordering defined in Section 4 of the V1.0 analysis plan.
presentation_heat_labels <- c(
  absolute_cumulative_heat = "Absolute-threshold cumulative heat",
  relative_cumulative_heat = "Relative-threshold cumulative heat",
  mixed_cumulative_heat = "Mixed-threshold cumulative heat",
  daytime_cumulative_heat = "Daytime cumulative heat",
  nighttime_cumulative_heat = "Nighttime cumulative heat",
  compound_day_night_cumulative_heat = "Compound day–night cumulative heat",
  annual_mean_humidex = "Annual mean humidex",
  annual_mean_tmean = "Annual mean temperature"
)
presentation_outcome_labels <- c(
  `1` = "Death", `2` = "DALYs", `5` = "Prevalence", `6` = "Incidence"
)
presentation_age_labels <- c(
  `158` = "<20 years", `2070` = "20–70 years", `26` = ">70 years"
)
presentation_sex_labels <- c(`1` = "Male", `2` = "Female")

table_0_did_results <- results[, .(
  `Heat indicator ×outcome` = paste0(
    unname(presentation_heat_labels[heat_variable]), " × ",
    unname(presentation_outcome_labels[as.character(measure_id)])
  ),
  `β` = format_number(beta),
  `CI` = paste0("[", format_number(ci_low), ", ", format_number(ci_high), "]"),
  `P value` = format_p(p_value),
  `Model` = model_id
)]

# Reviewer-requested fixed-effect and COVID-19 sensitivity reporting tables.
covid_exclude_2020_results <- copy(
  sensitivity_results[sensitivity_id == "exclude_2020"]
)

country_year_twfe_table <- country_year_twfe_results[, .(
  `Heat indicator` = unname(presentation_heat_labels[heat_variable]),
  `Outcome` = unname(presentation_outcome_labels[as.character(measure_id)]),
  `β per IQR` = format_number(beta),
  `95% CI` = paste0("[", format_number(ci_low), ", ", format_number(ci_high), "]"),
  `Rate ratio` = format_number(rate_ratio),
  `RR 95% CI` = paste0("[", format_number(rr_ci_low), ", ", format_number(rr_ci_high), "]"),
  `P value` = format_p(p_value),
  `BH q` = format_p(p_adjust_bh),
  `Observations` = format(n_observations, big.mark = ","),
  `Countries` = as.character(n_countries)
)]

m5_excluding_2020_table <- covid_exclude_2020_results[, .(
  `Heat indicator` = unname(presentation_heat_labels[heat_variable]),
  `Outcome` = unname(presentation_outcome_labels[as.character(measure_id)]),
  `β per IQR` = format_number(beta),
  `95% CI` = paste0("[", format_number(ci_low), ", ", format_number(ci_high), "]"),
  `Rate ratio` = format_number(rate_ratio),
  `RR 95% CI` = paste0("[", format_number(rr_ci_low), ", ", format_number(rr_ci_high), "]"),
  `P value` = format_p(p_value),
  `BH q` = format_p(p_adjust_bh),
  `Observations` = format(n_observations, big.mark = ","),
  `Countries` = as.character(n_countries)
)]

m5_comparison <- results[model_id == core_model_id, .(
  measure_id, heat_variable,
  m5_beta = beta, m5_ci_low = ci_low, m5_ci_high = ci_high,
  m5_p_value = p_value, m5_p_adjust_bh = p_adjust_bh,
  m5_n_observations = n_observations, m5_n_countries = n_countries
)]
twfe_comparison <- country_year_twfe_results[, .(
  measure_id, heat_variable,
  twfe_beta = beta, twfe_ci_low = ci_low, twfe_ci_high = ci_high,
  twfe_p_value = p_value, twfe_p_adjust_bh = p_adjust_bh,
  twfe_n_observations = n_observations, twfe_n_countries = n_countries
)]
covid_comparison <- covid_exclude_2020_results[, .(
  measure_id, heat_variable,
  no2020_beta = beta, no2020_ci_low = ci_low, no2020_ci_high = ci_high,
  no2020_p_value = p_value, no2020_p_adjust_bh = p_adjust_bh,
  no2020_n_observations = n_observations, no2020_n_countries = n_countries
)]
reviewer_robustness_comparison <- Reduce(
  function(x, y) merge(x, y, by = c("measure_id", "heat_variable"), all = FALSE),
  list(m5_comparison, twfe_comparison, covid_comparison)
)
reviewer_robustness_comparison[, `:=`(
  heat_order = match(heat_variable, heat_vars),
  outcome_order = match(measure_id, expected_measures),
  m5_twfe_same_sign = sign(m5_beta) == sign(twfe_beta),
  m5_no2020_same_sign = sign(m5_beta) == sign(no2020_beta)
)]
setorder(reviewer_robustness_comparison, heat_order, outcome_order)
reviewer_robustness_comparison[, c("heat_order", "outcome_order") := NULL]

reviewer_robustness_table <- reviewer_robustness_comparison[, .(
  `Heat indicator` = unname(presentation_heat_labels[heat_variable]),
  `Outcome` = unname(presentation_outcome_labels[as.character(measure_id)]),
  `Original M5 β (95% CI)` = paste0(
    format_number(m5_beta), " [", format_number(m5_ci_low), ", ",
    format_number(m5_ci_high), "]"
  ),
  `Original M5 P` = format_p(m5_p_value),
  `Country + year FE β (95% CI)` = paste0(
    format_number(twfe_beta), " [", format_number(twfe_ci_low), ", ",
    format_number(twfe_ci_high), "]"
  ),
  `Country + year FE P` = format_p(twfe_p_value),
  `M5 excluding 2020 β (95% CI)` = paste0(
    format_number(no2020_beta), " [", format_number(no2020_ci_low), ", ",
    format_number(no2020_ci_high), "]"
  ),
  `M5 excluding 2020 P` = format_p(no2020_p_value)
)]

cdnch_m5_twfe_no2020_table <- reviewer_robustness_comparison[
  heat_variable == "compound_day_night_cumulative_heat",
  .(
    `Outcome` = unname(
      presentation_outcome_labels[as.character(measure_id)]
    ),
    `Primary M5 β (95% CI)` = paste0(
      formatC(m5_beta, format = "f", digits = 3), " [",
      formatC(m5_ci_low, format = "f", digits = 3), ", ",
      formatC(m5_ci_high, format = "f", digits = 3), "]"
    ),
    `P value` = format_p(m5_p_value),
    `Country + year TWFE β (95% CI)` = paste0(
      formatC(twfe_beta, format = "f", digits = 3), " [",
      formatC(twfe_ci_low, format = "f", digits = 3), ", ",
      formatC(twfe_ci_high, format = "f", digits = 3), "]"
    ),
    `P value ` = format_p(twfe_p_value),
    `M5 excluding 2020 β (95% CI)` = paste0(
      formatC(no2020_beta, format = "f", digits = 3), " [",
      formatC(no2020_ci_low, format = "f", digits = 3), ", ",
      formatC(no2020_ci_high, format = "f", digits = 3), "]"
    ),
    `P value  ` = format_p(no2020_p_value)
  )
]

# Three indicator-model groups: CDNCH-only, RCH-only, and CDNCH+RCH.
cdnch_rch_single_results <- copy(results[
  model_id == core_model_id & heat_variable %in% unname(mutual_heat_vars)
])
cdnch_rch_single_results[, indicator := fifelse(
  heat_variable == mutual_heat_vars[["CDNCH"]], "CDNCH", "RCH"
)]
cdnch_rch_single_results[, `:=`(
  model_group = paste0(indicator, "-only"),
  analysis_year_first = 1990L,
  analysis_year_last = 2020L,
  coefficient = "heat_iqr"
)]
cdnch_rch_single_results <- cdnch_rch_single_results[, .(
  measure_id, outcome, indicator, heat_variable, increment,
  model_group, model_id, family, effect_scale, beta,
  std_error_clustered, ci_low, ci_high, rate_ratio, rr_ci_low, rr_ci_high,
  p_value, p_adjust_bh, n_observations, n_countries,
  analysis_year_first, analysis_year_last, inference_df, coefficient
)]
cdnch_rch_model_results_long <- rbindlist(
  list(cdnch_rch_single_results, cdnch_rch_joint_results),
  use.names = TRUE, fill = TRUE
)
cdnch_rch_model_results_long[, model_group_order := match(
  model_group, c("CDNCH-only", "RCH-only", "CDNCH+RCH")
)]
cdnch_rch_model_results_long[, indicator_order := match(
  indicator, c("CDNCH", "RCH")
)]
setorder(
  cdnch_rch_model_results_long,
  measure_id, model_group_order, indicator_order
)
cdnch_rch_model_results_long[, c(
  "model_group_order", "indicator_order"
) := NULL]

cdnch_single_compare <- cdnch_rch_single_results[
  indicator == "CDNCH",
  .(
    measure_id, outcome,
    cdnch_single_beta = beta,
    cdnch_single_ci_low = ci_low,
    cdnch_single_ci_high = ci_high
  )
]
rch_single_compare <- cdnch_rch_single_results[
  indicator == "RCH",
  .(
    measure_id,
    rch_single_beta = beta,
    rch_single_ci_low = ci_low,
    rch_single_ci_high = ci_high
  )
]
cdnch_joint_compare <- cdnch_rch_joint_results[
  indicator == "CDNCH",
  .(
    measure_id,
    cdnch_joint_beta = beta,
    cdnch_joint_ci_low = ci_low,
    cdnch_joint_ci_high = ci_high
  )
]
rch_joint_compare <- cdnch_rch_joint_results[
  indicator == "RCH",
  .(
    measure_id,
    rch_joint_beta = beta,
    rch_joint_ci_low = ci_low,
    rch_joint_ci_high = ci_high
  )
]
cdnch_rch_mutual_adjustment_comparison <- Reduce(
  function(x, y) merge(x, y, by = "measure_id", all = FALSE),
  list(
    cdnch_single_compare, rch_single_compare,
    cdnch_joint_compare, rch_joint_compare
  )
)
cdnch_rch_mutual_adjustment_comparison[, `:=`(
  cdnch_beta_change = cdnch_joint_beta - cdnch_single_beta,
  cdnch_attenuation_percent = 100 * (
    1 - abs(cdnch_joint_beta) / abs(cdnch_single_beta)
  ),
  cdnch_change_abs_beta_percent = 100 * (
    abs(cdnch_joint_beta) - abs(cdnch_single_beta)
  ) / abs(cdnch_single_beta),
  cdnch_sign_changed = sign(cdnch_joint_beta) != sign(cdnch_single_beta),
  rch_beta_change = rch_joint_beta - rch_single_beta,
  rch_attenuation_percent = 100 * (
    1 - abs(rch_joint_beta) / abs(rch_single_beta)
  ),
  rch_change_abs_beta_percent = 100 * (
    abs(rch_joint_beta) - abs(rch_single_beta)
  ) / abs(rch_single_beta),
  rch_sign_changed = sign(rch_joint_beta) != sign(rch_single_beta)
)]
setorder(cdnch_rch_mutual_adjustment_comparison, measure_id)

format_attenuation <- function(x) paste0(
  formatC(x, format = "f", digits = 1), "%"
)
cdnch_rch_mutual_adjustment_table <- cdnch_rch_mutual_adjustment_comparison[, .(
  `Outcome` = unname(presentation_outcome_labels[as.character(measure_id)]),
  `CDNCH-only β (95% CI)` = paste0(
    format_number(cdnch_single_beta), " [",
    format_number(cdnch_single_ci_low), ", ",
    format_number(cdnch_single_ci_high), "]"
  ),
  `after add RCH β (95% CI)` = paste0(
    format_number(cdnch_joint_beta), " [",
    format_number(cdnch_joint_ci_low), ", ",
    format_number(cdnch_joint_ci_high), "]"
  ),
  `CDNCH Change in \\|β\\|` = format_attenuation(
    cdnch_change_abs_beta_percent
  ),
  `RCH-only β (95% CI)` = paste0(
    format_number(rch_single_beta), " [",
    format_number(rch_single_ci_low), ", ",
    format_number(rch_single_ci_high), "]"
  ),
  `after add CDNCH β (95% CI)` = paste0(
    format_number(rch_joint_beta), " [",
    format_number(rch_joint_ci_low), ", ",
    format_number(rch_joint_ci_high), "]"
  ),
  `RCH Change in \\|β\\|` = format_attenuation(rch_change_abs_beta_percent)
)]

supplement_table_1_superregion <- descriptive_superregion_medians[, .(
  `GBDregion` = superregion,
  `DCH` = format_number(daytime_cumulative_heat),
  `NCH` = format_number(nighttime_cumulative_heat),
  `CDNCH` = format_number(compound_day_night_cumulative_heat),
  `SDI` = format_number(mean_value),
  `Deaths` = format_number(death_burden),
  `Incidence` = format_number(incidence_burden),
  `DALYs` = format_number(dalys_burden),
  `Prevalence` = format_number(prevalence_burden),
  `Household air pollution` = format_number(SEV_agestd_air_hap),
  `Ambient particulate matter pollution` = format_number(SEV_agestd_air_pm),
  `Fruit intake` = format_number(SEV_agestd_diet_fruit),
  `Vegetable intake` = format_number(SEV_agestd_diet_veg),
  `BMI` = format_number(SEV_agestd_metab_bmi),
  `LDL` = format_number(mean_LDL)
)]

table_1_subgroup_significant <- subgroup_m5_results[p_value < 0.05, .(
  `Outcome` = unname(presentation_outcome_labels[as.character(measure_id)]),
  `age` = unname(presentation_age_labels[as.character(age_id)]),
  `sex` = unname(presentation_sex_labels[as.character(sex_id)]),
  `rate difference` = format_number(beta),
  `P value` = format_p(p_value)
)]

supplement_table_2_older_differences <- subgroup_difference_contrasts[
  contrast_family == "Older vs other age group",
  .(
    `Outcome` = unname(presentation_outcome_labels[as.character(measure_id)]),
    `age` = ">70 years",
    `sex` = unname(presentation_sex_labels[as.character(lhs_sex_id)]),
    `reference group` = unname(presentation_age_labels[as.character(rhs_age_id)]),
    `rate difference` = format_number(slope_difference),
    `P value` = format_p(p_value)
  )
]

supplement_table_3_sex_differences <- subgroup_difference_contrasts[
  contrast_family == "Male vs female within age group",
  .(
    `Outcome` = unname(presentation_outcome_labels[as.character(measure_id)]),
    `age` = unname(presentation_age_labels[as.character(lhs_age_id)]),
    `male−female rate difference` = format_number(slope_difference),
    `P value` = format_p(p_value)
  )
]

get_attributable_row <- function(measure) {
  attributable_burden_overall_summary[measure_id == measure]
}
presentation_attributable_fraction <- function(measure) {
  row <- get_attributable_row(measure)
  paste0(
    format_percent(row$cumulative_attributable_fraction),
    " (95% CI: ",
    format_percent(row$cumulative_attributable_fraction_model_ci_low),
    "–",
    format_percent(row$cumulative_attributable_fraction_model_ci_high),
    ")"
  )
}
presentation_attributable_burden <- function(measure) {
  row <- get_attributable_row(measure)
  paste0(
    "≈", formatC(
      row$cumulative_attributable_burden / 1e6,
      format = "f", digits = 1, big.mark = ","
    ), " million (95% CI: ", formatC(
      row$cumulative_attributable_burden_model_ci_low / 1e6,
      format = "f", digits = 1, big.mark = ","
    ), "–", formatC(
      row$cumulative_attributable_burden_model_ci_high / 1e6,
      format = "f", digits = 1, big.mark = ","
    ), " million)"
  )
}
presentation_attributable_input_bounds <- function(measure) {
  row <- get_attributable_row(measure)
  paste0(
    "≈", formatC(
      row$cumulative_attributable_burden_lower_input / 1e6,
      format = "f", digits = 1
    ), "–", formatC(
      row$cumulative_attributable_burden_upper_input / 1e6,
      format = "f", digits = 1
    ), " million"
  )
}
presentation_attributable_range <- function(measure) {
  row <- get_attributable_row(measure)
  paste0(
    format_percent(row$annual_attributable_fraction_min), "–",
    format_percent(row$annual_attributable_fraction_max)
  )
}
table_2_attributable_summary <- data.table(
  `Attributable measure` = c(
    "Cumulative attributable fraction (95% CI)",
    "Excess attributable number (95% CI)",
    "Attributable fraction range (1990–2020)"
  ),
  `Deaths` = c(
    presentation_attributable_fraction(1L),
    presentation_attributable_burden(1L),
    presentation_attributable_range(1L)
  ),
  `Incidence` = c(
    presentation_attributable_fraction(6L),
    presentation_attributable_burden(6L),
    presentation_attributable_range(6L)
  ),
  `DALYS` = c(
    presentation_attributable_fraction(2L),
    presentation_attributable_burden(2L),
    presentation_attributable_range(2L)
  ),
  `Prevalence` = c(
    presentation_attributable_fraction(5L),
    presentation_attributable_burden(5L),
    presentation_attributable_range(5L)
  )
)
setnames(table_2_attributable_summary, "Attributable measure", "")

attributable_input_bounds_summary <- attributable_burden_overall_summary[, .(
  measure_id, outcome, burden_unit, model_id, heat_variable,
  reference_heat, reference_heat_unit, first_year, last_year, countries,
  attributable_burden_lower_input = cumulative_attributable_burden_lower_input,
  attributable_burden_point = cumulative_attributable_burden,
  attributable_burden_upper_input = cumulative_attributable_burden_upper_input,
  attributable_fraction_lower_input = cumulative_attributable_fraction_lower_input,
  attributable_fraction_point = cumulative_attributable_fraction,
  attributable_fraction_upper_input = cumulative_attributable_fraction_upper_input,
  boundary_method = paste0(
    "Row-level GBD lower/upper burden inputs substituted with the fitted ",
    "point-estimate M5 spline-RR held fixed; deterministic input-boundary ",
    "sensitivity, not a total 95% uncertainty interval"
  )
)]

attributable_model_uncertainty_summary <- attributable_burden_overall_summary[, .(
  measure_id, outcome, burden_unit, model_id, heat_variable,
  reference_heat, reference_heat_unit, first_year, last_year, countries,
  attributable_burden_point = cumulative_attributable_burden,
  attributable_burden_model_ci_low =
    cumulative_attributable_burden_model_ci_low,
  attributable_burden_model_ci_high =
    cumulative_attributable_burden_model_ci_high,
  attributable_fraction_point = cumulative_attributable_fraction,
  attributable_fraction_model_ci_low =
    cumulative_attributable_fraction_model_ci_low,
  attributable_fraction_model_ci_high =
    cumulative_attributable_fraction_model_ci_high,
  simulation_draws = model_uncertainty_draws,
  simulation_seed = model_uncertainty_seed,
  uncertainty_method = model_uncertainty_method,
  uncertainty_scope = model_uncertainty_scope
)]

attributable_model_uncertainty_table <- attributable_model_uncertainty_summary[
  order(match(measure_id, c(1L, 6L, 2L, 5L))), .(
  `Outcome` = unname(presentation_outcome_labels[as.character(measure_id)]),
  `Excess attributable number (95% CI)` = paste0(
    "≈", formatC(
      attributable_burden_point / 1e6,
      format = "f", digits = 1, big.mark = ","
    ), " million (", formatC(
      attributable_burden_model_ci_low / 1e6,
      format = "f", digits = 1, big.mark = ","
    ), "–", formatC(
      attributable_burden_model_ci_high / 1e6,
      format = "f", digits = 1, big.mark = ","
    ), " million)"
  ),
  `Model-based AF (95% CI)` = paste0(
    format_percent(attributable_fraction_point), " (",
    format_percent(attributable_fraction_model_ci_low), "–",
    format_percent(attributable_fraction_model_ci_high), ")"
  )
)]

attributable_model_uncertainty_report <- c(
  "# Model-based uncertainty for cumulative attributable burden",
  "",
  paste0(
    "The 95% confidence intervals were obtained from ",
    format(attributable_uncertainty_draws, big.mark = ","),
    " parametric multivariate-normal draws of the four core-M5 spline ",
    "coefficients using their country-clustered covariance matrix. For each ",
    "draw, RR, harmful-only AF, and cumulative attributable burden were ",
    "recalculated over all country-years."
  ),
  "",
  paste0(
    "These intervals quantify fitted spline-coefficient uncertainty with the ",
    "GBD point burden held fixed. They exclude GBD burden-input uncertainty, ",
    "counterfactual uncertainty, truncation-rule uncertainty, and broader ",
    "model-specification uncertainty. The GBD lower/upper input-boundary ",
    "analysis is reported separately and must not be interpreted as a 95% CI."
  )
)

attributable_input_bounds_table <- attributable_input_bounds_summary[, .(
  `Outcome` = unname(presentation_outcome_labels[as.character(measure_id)]),
  `Unit` = burden_unit,
  `Lower-input burden` = paste0(
    "≈", formatC(
      attributable_burden_lower_input / 1e6,
      format = "f", digits = 1
    ), " million"
  ),
  `Point estimate` = paste0(
    "≈", formatC(
      attributable_burden_point / 1e6,
      format = "f", digits = 1
    ), " million"
  ),
  `Upper-input burden` = paste0(
    "≈", formatC(
      attributable_burden_upper_input / 1e6,
      format = "f", digits = 1
    ), " million"
  ),
  `Lower / point / upper AF` = paste0(
    format_percent(attributable_fraction_lower_input, digits = 1L), " / ",
    format_percent(attributable_fraction_point, digits = 1L), " / ",
    format_percent(attributable_fraction_upper_input, digits = 1L)
  )
)]

presentation_sensitivity_labels <- c(
  trim_1_99 = "Excluding exposures outside the 1st–99th percentiles",
  lag1 = "lag 1",
  lag01_mean = "lag 0–1 mean",
  exclude_2020 = "Excluding 2020 (COVID-19 shock)"
)
supplement_table_4_sensitivity <- sensitivity_summary[, .(
  `Sensitivity analysis` = unname(presentation_sensitivity_labels[sensitivity_id]),
  `Tests` = as.character(tests),
  `P<0.05` = as.character(raw_p_lt_0_05),
  `BH q<0.05` = as.character(bh_q_lt_0_05),
  `Observations per model` = ifelse(
    n_observations_min == n_observations_max,
    format(n_observations_min, big.mark = ","),
    paste0(
      format(n_observations_min, big.mark = ","), "–",
      format(n_observations_max, big.mark = ",")
    )
  ),
  `Countries` = ifelse(
    n_countries_min == n_countries_max,
    as.character(n_countries_min),
    paste0(n_countries_min, "–", n_countries_max)
  ),
  `Outcome years` = paste0(first_year, "–", last_year)
)]

assert_presentation_table <- function(x, expected_names, expected_rows, label) {
  if (!identical(names(x), expected_names)) {
    stop(label, " column schema does not match the V1.0 plan.", call. = FALSE)
  }
  if (nrow(x) != expected_rows) {
    stop(label, " row count does not match the V1.0 plan.", call. = FALSE)
  }
  if (any(vapply(x, anyNA, logical(1)))) {
    stop(label, " contains missing values.", call. = FALSE)
  }
}

presentation_tables <- list(
  table_0_did_results = table_0_did_results,
  supplement_table_1_superregion_medians = supplement_table_1_superregion,
  table_1_subgroup_significant = table_1_subgroup_significant,
  supplement_table_2_older_slope_differences = supplement_table_2_older_differences,
  supplement_table_3_sex_slope_differences = supplement_table_3_sex_differences,
  table_2_attributable_summary = table_2_attributable_summary,
  supplement_table_4_sensitivity_summary = supplement_table_4_sensitivity
)
presentation_schemas <- list(
  c("Heat indicator ×outcome", "β", "CI", "P value", "Model"),
  c("GBDregion", "DCH", "NCH", "CDNCH", "SDI", "Deaths", "Incidence", "DALYs", "Prevalence", "Household air pollution", "Ambient particulate matter pollution", "Fruit intake", "Vegetable intake", "BMI", "LDL"),
  c("Outcome", "age", "sex", "rate difference", "P value"),
  c("Outcome", "age", "sex", "reference group", "rate difference", "P value"),
  c("Outcome", "age", "male−female rate difference", "P value"),
  c("", "Deaths", "Incidence", "DALYS", "Prevalence"),
  c("Sensitivity analysis", "Tests", "P<0.05", "BH q<0.05", "Observations per model", "Countries", "Outcome years")
)
presentation_rows <- c(224L, 21L, 11L, 16L, 12L, 3L, 4L)
for (i in seq_along(presentation_tables)) {
  assert_presentation_table(
    presentation_tables[[i]], presentation_schemas[[i]], presentation_rows[[i]],
    names(presentation_tables)[[i]]
  )
}

markdown_exact_table <- function(x) {
  displayed <- as.data.frame(x, check.names = FALSE)
  c(
    paste0("| ", paste(names(displayed), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(displayed)), collapse = "|"), "|"),
    apply(displayed, 1L, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
  )
}

if (
  nrow(attributable_input_bounds_summary) != 4L ||
    attributable_input_bounds_summary[, any(
      attributable_burden_lower_input > attributable_burden_point |
        attributable_burden_point > attributable_burden_upper_input |
        !is.finite(attributable_burden_lower_input) |
        !is.finite(attributable_burden_point) |
        !is.finite(attributable_burden_upper_input)
    )]
) {
  stop("Attributable input-boundary sensitivity output is invalid.", call. = FALSE)
}

attributable_input_bounds_report <- c(
  "# GBD lower/upper input-boundary sensitivity for attributable burden",
  "",
  "The fitted point-estimate core-M5 4-df spline exposure–response function and its attributable fraction are held fixed. For each country–year, the GBD `lower`, `val`, and `upper` burden inputs are substituted to recalculate attributable burden.",
  "",
  markdown_exact_table(attributable_input_bounds_table),
  "",
  "These bounds isolate sensitivity to the row-level GBD burden inputs. They are deterministic all-lower/all-upper scenarios, not a total 95% uncertainty interval, because the spline coefficients were not refitted to GBD draws and the dependence structure of the GBD uncertainty draws is unavailable.",
  "",
  "Numbers are intentionally rounded to approximately 0.1 million. The global-median CDNCH counterfactual and harmful-only truncation rule remain unchanged in this analysis."
)

if (
  nrow(country_year_twfe_results) != 32L ||
    nrow(covid_exclude_2020_results) != 32L ||
    nrow(reviewer_robustness_comparison) != 32L
) {
  stop("Reviewer-requested robustness outputs must each contain 32 outcome–heat estimates.", call. = FALSE)
}

reviewer_specification_summary <- data.table(
  Specification = c(
    "Original core M5", "Country and year two-way fixed effects",
    "Core M5 excluding 2020"
  ),
  `Country FE` = c("No", "Yes", "No"),
  `Year control` = c(
    "Superregion-specific 3 df natural spline", "Year fixed effects",
    "Superregion-specific 3 df natural spline"
  ),
  `Outcome years` = c("1990–2020", "1990–2020", "1990–2019"),
  `Observations/model` = c("6,138", "6,138", "5,940"),
  `Raw P<0.05` = c(
    sum(results[model_id == core_model_id]$p_value < 0.05),
    sum(country_year_twfe_results$p_value < 0.05),
    sum(covid_exclude_2020_results$p_value < 0.05)
  ),
  `BH q<0.05` = c(
    sum(results[model_id == core_model_id]$p_adjust_bh < 0.05),
    sum(country_year_twfe_results$p_adjust_bh < 0.05),
    sum(covid_exclude_2020_results$p_adjust_bh < 0.05)
  )
)

reviewer_robustness_report <- c(
  "# 国家固定效应与 2020 年 COVID-19 冲击敏感性报告",
  "",
  "## 中文结论摘要",
  "",
  paste0(
    "原 M5 的 32 项检验中有 ",
    sum(results[model_id == core_model_id]$p_value < 0.05),
    " 项原始 P<0.05、",
    sum(results[model_id == core_model_id]$p_adjust_bh < 0.05),
    " 项 BH q<0.05。"
  ),
  paste0(
    "国家 + 年份双向固定效应模型中相应数量降至 ",
    sum(country_year_twfe_results$p_value < 0.05), " 和 ",
    sum(country_year_twfe_results$p_adjust_bh < 0.05),
    "；与原 M5 的系数方向仅 ",
    sum(reviewer_robustness_comparison$m5_twfe_same_sign), "/32 一致。"
  ),
  paste0(
    "剔除 2020 年后仍有 ",
    sum(covid_exclude_2020_results$p_value < 0.05), " 项原始 P<0.05、",
    sum(covid_exclude_2020_results$p_adjust_bh < 0.05),
    " 项 BH q<0.05，且 32/32 项系数方向与原 M5 一致。"
  ),
  "因此，2020 年单年冲击不是原 M5 结果的主要驱动因素；但国家固定效应显著改变了整体结果，应将双向固定效应结果作为重要稳健性边界共同报告。",
  "",
  "## Model specifications",
  "",
  paste0(
    "The original core M5 is retained unchanged: `Heat + 7 covariates + ",
    "ns(year, df=3):factor(superregion)`, with a quasi-Poisson log link and ",
    "`ln(population)` offset. M5 does not contain country fixed effects."
  ),
  "",
  paste0(
    "The additional two-way fixed-effects model is `Heat + 7 covariates + ",
    "factor(country) + factor(year)`, with the same outcome definition, heat ",
    "centering/IQR scaling, quasi-Poisson log link, population offset, and ",
    "country-clustered HC1 inference. It is reported separately because its ",
    "identifying variation differs from M5."
  ),
  "",
  paste0(
    "The COVID-19 sensitivity analysis re-estimates the original M5 after ",
    "excluding outcome year 2020. Heat remains scaled by the full primary-sample ",
    "IQR, so the coefficient is directly comparable with the original M5. This ",
    "analysis is reviewer-requested rather than part of the three original ",
    "prespecified sensitivity analyses."
  ),
  "",
  "All three specifications use country-clustered HC1 standard errors and t inference with 197 degrees of freedom.",
  "",
  "## Overall comparison",
  "",
  markdown_exact_table(reviewer_specification_summary),
  "",
  paste0(
    "Coefficient signs agree with the original M5 for ",
    sum(reviewer_robustness_comparison$m5_twfe_same_sign), "/32 estimates in the country-and-year fixed-effects model and ",
    sum(reviewer_robustness_comparison$m5_no2020_same_sign), "/32 estimates after excluding 2020."
  ),
  "",
  "## Estimate-by-estimate comparison",
  "",
  markdown_exact_table(reviewer_robustness_table),
  "",
  "The complete machine-readable estimates, standard errors, rate ratios, adjusted P values, sample sizes, and diagnostics are provided in the accompanying CSV files."
)

if (
  nrow(cdnch_rch_model_results_long) != 16L ||
    nrow(cdnch_rch_mutual_adjustment_comparison) != 4L ||
    any(!is.finite(
      cdnch_rch_mutual_adjustment_comparison$cdnch_attenuation_percent
    )) ||
    any(!is.finite(
      cdnch_rch_mutual_adjustment_comparison$rch_attenuation_percent
    ))
) {
  stop("CDNCH/RCH mutual-adjustment output is incomplete or invalid.", call. = FALSE)
}

cdnch_rch_mutual_adjustment_report <- c(
  "# CDNCH 与 RCH 相互调整分析报告",
  "",
  "## 模型结构",
  "",
  "三组模型均沿用原核心 M5：准 Poisson 对数链接、`ln(population)` offset、7 项预设协变量、`ns(year, df=3):factor(superregion)` 时间趋势，以及国家聚类 HC1 推断。每个结局均包含 6,138 个国家–年份和 198 个国家。",
  "",
  "1. **CDNCH-only：** `CDNCH + covariates + regional spline trends`。",
  "2. **RCH-only：** `RCH + covariates + regional spline trends`。",
  "3. **CDNCH+RCH：** `CDNCH + RCH + covariates + regional spline trends`。",
  "",
  "CDNCH 与 RCH 均分别按主分析完整样本自身的 IQR 缩放。联合模型中的 CDNCH 系数表示控制 RCH 后每 1 个 CDNCH-IQR 的条件关联；RCH 系数反之。",
  paste0(
    "在 6,138 个主分析国家–年份中，CDNCH 与 RCH 的 Pearson 相关系数为 ",
    formatC(mutual_indicator_correlation, format = "f", digits = 3),
    "。仅由两指标相关性得到的配对 VIF 为 ",
    formatC(mutual_pairwise_vif, format = "f", digits = 2),
    "；在包含 7 项协变量和区域特异时间趋势的完整联合 M5 设计矩阵中，",
    "CDNCH 和 RCH 的 VIF 分别为 ",
    formatC(cdnch_vif$vif, format = "f", digits = 2), " 和 ",
    formatC(rch_vif$vif, format = "f", digits = 2),
    "。联合模型的置信区间需结合这种共线性解释。"
  ),
  "",
  "## 互调前后比较",
  "",
  markdown_exact_table(cdnch_rch_mutual_adjustment_table),
  "",
  "`Change in |β|` 定义为 `100 × (|β_after| − |β_before|) / |β_before|`。负值表示加入另一指标后系数绝对值下降，即部分变异被另一指标吸收；正值表示绝对效应增强。该统计量描述效应幅度的相对变化，不等同于中介比例。方向反转情况另保留在机器可读比较文件中。",
  "",
  paste0(
    "CDNCH 加入 RCH 后有 ",
    sum(cdnch_rch_mutual_adjustment_comparison$cdnch_sign_changed),
    "/4 个结局发生系数方向反转；RCH 加入 CDNCH 后有 ",
    sum(cdnch_rch_mutual_adjustment_comparison$rch_sign_changed),
    "/4 个结局发生方向反转。"
  ),
  "",
  "联合模型中，控制 RCH 后 CDNCH 与 DALYs 和发病的 95% CI 仍不跨 0，而死亡和患病的 95% CI 跨 0；控制 CDNCH 后，RCH 对四类结局的 95% CI 均跨 0。该结果表明两项高度相关热指标的单独模型效应不能解释为彼此独立的效应。",
  "",
  "完整 β、聚类标准误、95% CI、RR、P 值、BH q 值和样本量见 `cdnch_rch_model_results_long.csv`；`Change in |β|`、原衰减比例及方向变化见 `cdnch_rch_mutual_adjustment_comparison.csv`；相关系数和 VIF 见 `cdnch_rch_collinearity_diagnostics.csv`。"
)

fwrite(results, file.path(output_dir, "did_results_full.csv"), bom = TRUE)
fwrite(results, file.path(output_dir, "stage_a_did_results.csv"), bom = TRUE)
fwrite(
  country_year_twfe_results,
  file.path(output_dir, "country_year_twfe_results.csv"), bom = TRUE
)
fwrite(
  covid_exclude_2020_results,
  file.path(output_dir, "m5_excluding_2020_results.csv"), bom = TRUE
)
fwrite(
  reviewer_robustness_comparison,
  file.path(output_dir, "reviewer_robustness_comparison.csv"), bom = TRUE
)
fwrite(
  reviewer_specification_summary,
  file.path(output_dir, "reviewer_robustness_summary.csv"), bom = TRUE
)
fwrite(
  cdnch_m5_twfe_no2020_table,
  file.path(output_dir, "cdnch_m5_twfe_no2020_comparison.csv"), bom = TRUE
)
write_markdown_utf8(
  markdown_exact_table(cdnch_m5_twfe_no2020_table),
  file.path(output_dir, "cdnch_m5_twfe_no2020_comparison.md")
)
fwrite(
  cdnch_rch_model_results_long,
  file.path(output_dir, "cdnch_rch_model_results_long.csv"), bom = TRUE
)
fwrite(
  cdnch_rch_mutual_adjustment_comparison,
  file.path(output_dir, "cdnch_rch_mutual_adjustment_comparison.csv"),
  bom = TRUE
)
fwrite(
  cdnch_rch_collinearity_diagnostics,
  file.path(output_dir, "cdnch_rch_collinearity_diagnostics.csv"),
  bom = TRUE
)
write_markdown_utf8(
  markdown_exact_table(cdnch_rch_mutual_adjustment_table),
  file.path(output_dir, "cdnch_rch_mutual_adjustment_comparison.md")
)
write_markdown_utf8(
  cdnch_rch_mutual_adjustment_report,
  file.path(output_dir, "cdnch_rch_mutual_adjustment_report.md")
)
write_markdown_utf8(
  markdown_exact_table(country_year_twfe_table),
  file.path(output_dir, "country_year_twfe_results.md")
)
write_markdown_utf8(
  markdown_exact_table(m5_excluding_2020_table),
  file.path(output_dir, "m5_excluding_2020_results.md")
)
write_markdown_utf8(
  reviewer_robustness_report,
  file.path(output_dir, "reviewer_robustness_report.md")
)
fwrite(diagnostics, file.path(output_dir, "model_diagnostics.csv"), bom = TRUE)
fwrite(model_definitions, file.path(output_dir, "model_definitions.csv"), bom = TRUE)
fwrite(exposure_scaling, file.path(output_dir, "exposure_scaling.csv"), bom = TRUE)
fwrite(subgroup_exposure_scaling, file.path(output_dir, "subgroup_exposure_scaling.csv"), bom = TRUE)
fwrite(subgroup_m5_results, file.path(output_dir, "subgroup_m5_results.csv"), bom = TRUE)
fwrite(subgroup_m5_summary, file.path(output_dir, "subgroup_m5_summary.csv"), bom = TRUE)
fwrite(subgroup_m5_diagnostics, file.path(output_dir, "subgroup_m5_model_diagnostics.csv"), bom = TRUE)
fwrite(
  subgroup_age_derivation_audit,
  file.path(output_dir, "subgroup_age_derivation_20_70.csv"),
  bom = TRUE
)
fwrite(
  subgroup_interaction_tests,
  file.path(output_dir, "subgroup_interaction_tests.csv"),
  bom = TRUE
)
fwrite(
  subgroup_interaction_diagnostics,
  file.path(output_dir, "subgroup_interaction_model_diagnostics.csv"),
  bom = TRUE
)
fwrite(
  subgroup_difference_contrasts,
  file.path(output_dir, "subgroup_difference_contrasts.csv"),
  bom = TRUE
)
fwrite(
  subgroup_older_joint_tests,
  file.path(output_dir, "subgroup_older_vs_others_joint_tests.csv"),
  bom = TRUE
)
fwrite(sample_flow, file.path(output_dir, "sample_flow.csv"), bom = TRUE)
fwrite(excluded_locations, file.path(output_dir, "excluded_locations.csv"), bom = TRUE)
fwrite(sensitivity_results, file.path(output_dir, "sensitivity_results_m5.csv"), bom = TRUE)
fwrite(sensitivity_summary, file.path(output_dir, "sensitivity_summary_m5.csv"), bom = TRUE)
fwrite(nonlinear_curves, file.path(output_dir, "compound_nonlinear_curves.csv"), bom = TRUE)
fwrite(nonlinear_tests, file.path(output_dir, "compound_nonlinear_tests.csv"), bom = TRUE)
fwrite(attributable_burden, file.path(output_dir, "attributable_cvd_burden_m5.csv"), bom = TRUE)
fwrite(
  attributable_burden_summary,
  file.path(output_dir, "attributable_cvd_burden_m5_annual_summary.csv"),
  bom = TRUE
)
fwrite(
  attributable_burden_overall_summary,
  file.path(output_dir, "attributable_burden_summary_m5.csv"),
  bom = TRUE
)
fwrite(
  attributable_input_bounds_summary,
  file.path(output_dir, "attributable_burden_input_bounds_m5.csv"),
  bom = TRUE
)
fwrite(
  attributable_model_uncertainty_summary,
  file.path(output_dir, "attributable_burden_model_uncertainty_m5.csv"),
  bom = TRUE
)
write_markdown_utf8(
  markdown_exact_table(attributable_model_uncertainty_table),
  file.path(output_dir, "attributable_burden_model_uncertainty_m5.md")
)
write_markdown_utf8(
  c(
    attributable_model_uncertainty_report, "",
    markdown_exact_table(attributable_model_uncertainty_table)
  ),
  file.path(output_dir, "attributable_burden_model_uncertainty_report.md")
)
write_markdown_utf8(
  markdown_exact_table(attributable_input_bounds_table),
  file.path(output_dir, "attributable_burden_input_bounds_m5.md")
)
write_markdown_utf8(
  attributable_input_bounds_report,
  file.path(output_dir, "attributable_burden_input_bounds_report.md")
)
fwrite(attributable_deaths, file.path(output_dir, "attributable_cvd_deaths_m5.csv"), bom = TRUE)
fwrite(
  attributable_deaths_annual_summary,
  file.path(output_dir, "attributable_cvd_deaths_m5_annual_summary.csv"),
  bom = TRUE
)
fwrite(outcome_map_data, file.path(output_dir, "map_outcomes_1990_2020_mean.csv"), bom = TRUE)
fwrite(indicator_map_data, file.path(output_dir, "map_heat_indicators_1990_2020_mean.csv"), bom = TRUE)
fwrite(descriptive_superregion_medians, file.path(output_dir, "descriptive_superregion_medians.csv"), bom = TRUE)
for (table_name in names(presentation_tables)) {
  table_data <- presentation_tables[[table_name]]
  fwrite(table_data, file.path(output_dir, paste0(table_name, ".csv")), bom = TRUE)
  write_markdown_utf8(
    markdown_exact_table(table_data),
    file.path(output_dir, paste0(table_name, ".md"))
  )
}

# Compatibility names point to the matching prespecified presentation tables.
write_markdown_utf8(
  markdown_exact_table(table_0_did_results),
  file.path(output_dir, "did_results_table.md")
)
write_markdown_utf8(
  markdown_exact_table(supplement_table_1_superregion),
  file.path(output_dir, "descriptive_superregion_medians.md")
)
write_markdown_utf8(
  markdown_exact_table(supplement_table_4_sensitivity),
  file.path(output_dir, "sensitivity_summary_m5.md")
)
write_markdown_utf8(
  markdown_exact_table(table_2_attributable_summary),
  file.path(output_dir, "attributable_burden_summary_m5.md")
)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"), useBytes = TRUE)

message(
  "Done. Main rows: ", nrow(analysis_data), " (", uniqueN(analysis_data$location_id),
  " countries); linear results: ", nrow(results),
  "; subgroup rate results: ", nrow(subgroup_m5_results),
  "; subgroup interaction tests: ", nrow(subgroup_interaction_tests),
  "; subgroup pairwise difference tests: ", nrow(subgroup_difference_contrasts),
  "; older-vs-others joint tests: ", nrow(subgroup_older_joint_tests),
  "; country-year TWFE results: ", nrow(country_year_twfe_results),
  "; CDNCH/RCH model-result rows: ", nrow(cdnch_rch_model_results_long),
  "; sensitivity results: ", nrow(sensitivity_results),
  "; nonlinear fits: ", nrow(nonlinear_tests),
  "; attributable burden rows: ", nrow(attributable_burden),
  "; superregion descriptive rows: ", nrow(descriptive_superregion_medians),
  "; 1990–2020 mean outcome-map rows: ", nrow(outcome_map_data),
  "; 1990–2020 mean heat-indicator map rows: ", nrow(indicator_map_data), "."
)
