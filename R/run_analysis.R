required_packages <- c(
  "dplyr", "tidyr", "ggplot2", "broom", "sandwich", "lmtest",
  "splines", "lavaan"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Install required packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(broom)
  library(sandwich)
  library(lmtest)
  library(splines)
  library(lavaan)
})

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

data_path <- "data/synthetic_analysis_data.csv"
if (!file.exists(data_path)) {
  source("R/generate_synthetic_data.R")
}

analysis_data <- read.csv(data_path, stringsAsFactors = FALSE, check.names = FALSE) %>%
  mutate(
    exposure_flag = as.integer(exposure_flag),
    moderator_group = factor(moderator_group),
    covariate_02 = factor(covariate_02),
    covariate_04 = factor(covariate_04),
    covariate_06 = factor(covariate_06)
  )

outcome_names <- sprintf("outcome_%02d", seq_len(5))
indicator_names <- sprintf("indicator_%02d", seq_len(6))
continuous_covariates <- c("covariate_01", "covariate_03", "covariate_05")
categorical_covariates <- c("covariate_02", "covariate_04", "covariate_06")

# Quality-control checks -------------------------------------------------------
qc_summary <- data.frame(
  check = c(
    "record_count", "duplicate_record_ids", "invalid_exposure_values",
    "nonpositive_outcome_values", "missing_outcome_cells", "missing_indicator_cells"
  ),
  value = c(
    nrow(analysis_data),
    sum(duplicated(analysis_data$record_id)),
    sum(!analysis_data$exposure_flag %in% c(0L, 1L)),
    sum(as.matrix(analysis_data[outcome_names]) <= 0, na.rm = TRUE),
    sum(is.na(analysis_data[outcome_names])),
    sum(is.na(analysis_data[indicator_names]))
  )
)
write.csv(qc_summary, "outputs/tables/quality_control_summary.csv", row.names = FALSE)

# Descriptive summaries -------------------------------------------------------
summarise_numeric <- function(data, variables, group_label) {
  bind_rows(lapply(variables, function(variable_name) {
    values <- data[[variable_name]]
    data.frame(
      group = group_label,
      variable = variable_name,
      n_observed = sum(!is.na(values)),
      missing_n = sum(is.na(values)),
      mean = mean(values, na.rm = TRUE),
      sd = sd(values, na.rm = TRUE),
      median = median(values, na.rm = TRUE),
      q1 = unname(quantile(values, 0.25, na.rm = TRUE)),
      q3 = unname(quantile(values, 0.75, na.rm = TRUE))
    )
  }))
}

numeric_variables <- c(continuous_covariates, outcome_names)
descriptive_summary <- bind_rows(
  summarise_numeric(analysis_data, numeric_variables, "Overall"),
  summarise_numeric(filter(analysis_data, exposure_flag == 0), numeric_variables, "Exposure_0"),
  summarise_numeric(filter(analysis_data, exposure_flag == 1), numeric_variables, "Exposure_1")
)
write.csv(descriptive_summary, "outputs/tables/descriptive_summary.csv", row.names = FALSE)

missingness_summary <- data.frame(
  variable = c(outcome_names, indicator_names),
  missing_n = vapply(analysis_data[c(outcome_names, indicator_names)], function(x) sum(is.na(x)), numeric(1)),
  missing_percent = 100 * vapply(
    analysis_data[c(outcome_names, indicator_names)],
    function(x) mean(is.na(x)),
    numeric(1)
  )
)
write.csv(missingness_summary, "outputs/tables/missingness_summary.csv", row.names = FALSE)

# Primary and sensitivity models ---------------------------------------------
adjustment_terms <- paste(
  c("exposure_flag", continuous_covariates, categorical_covariates),
  collapse = " + "
)

fit_primary_model <- function(outcome_name) {
  model_formula <- as.formula(paste0("log(", outcome_name, ") ~ ", adjustment_terms))
  fit <- lm(model_formula, data = analysis_data, na.action = na.omit)
  robust_test <- lmtest::coeftest(fit, vcov. = sandwich::vcovHC(fit, type = "HC3"))
  robust_row <- robust_test["exposure_flag", ]
  estimate <- unname(robust_row[1])
  standard_error <- unname(robust_row[2])
  data.frame(
    outcome = outcome_name,
    model = "Log-linear with HC3 SE",
    n = stats::nobs(fit),
    log_estimate = estimate,
    standard_error = standard_error,
    percent_difference = 100 * (exp(estimate) - 1),
    lower_95 = 100 * (exp(estimate - 1.96 * standard_error) - 1),
    upper_95 = 100 * (exp(estimate + 1.96 * standard_error) - 1),
    p_value = unname(robust_row[4])
  )
}

primary_results <- bind_rows(lapply(outcome_names, fit_primary_model))
write.csv(primary_results, "outputs/tables/model_results.csv", row.names = FALSE)

fit_sensitivity_models <- function(outcome_name) {
  gamma_formula <- as.formula(paste0(outcome_name, " ~ ", adjustment_terms))
  gamma_fit <- glm(
    gamma_formula,
    data = analysis_data,
    family = Gamma(link = "log"),
    na.action = na.omit
  )
  gamma_row <- summary(gamma_fit)$coefficients["exposure_flag", ]

  spline_terms <- paste(
    c(
      "exposure_flag", "splines::ns(covariate_01, df = 3)",
      "covariate_03", "covariate_05", categorical_covariates
    ),
    collapse = " + "
  )
  spline_formula <- as.formula(paste0("log(", outcome_name, ") ~ ", spline_terms))
  spline_fit <- lm(spline_formula, data = analysis_data, na.action = na.omit)
  spline_test <- lmtest::coeftest(spline_fit, vcov. = sandwich::vcovHC(spline_fit, type = "HC3"))
  spline_row <- spline_test["exposure_flag", ]

  bind_rows(
    data.frame(
      outcome = outcome_name,
      model = "Gamma GLM with log link",
      n = stats::nobs(gamma_fit),
      percent_difference = 100 * (exp(gamma_row[1]) - 1),
      lower_95 = 100 * (exp(gamma_row[1] - 1.96 * gamma_row[2]) - 1),
      upper_95 = 100 * (exp(gamma_row[1] + 1.96 * gamma_row[2]) - 1),
      p_value = gamma_row[4]
    ),
    data.frame(
      outcome = outcome_name,
      model = "Log-linear with nonlinear covariate_01",
      n = stats::nobs(spline_fit),
      percent_difference = 100 * (exp(spline_row[1]) - 1),
      lower_95 = 100 * (exp(spline_row[1] - 1.96 * spline_row[2]) - 1),
      upper_95 = 100 * (exp(spline_row[1] + 1.96 * spline_row[2]) - 1),
      p_value = spline_row[4]
    )
  )
}

sensitivity_results <- bind_rows(lapply(outcome_names, fit_sensitivity_models))
write.csv(sensitivity_results, "outputs/tables/sensitivity_results.csv", row.names = FALSE)

# Interaction estimates -------------------------------------------------------
fit_interaction_model <- function(outcome_name) {
  interaction_formula <- as.formula(paste0(
    "log(", outcome_name, ") ~ exposure_flag * moderator_group + ",
    paste(c(continuous_covariates, categorical_covariates), collapse = " + ")
  ))
  fit <- lm(interaction_formula, data = analysis_data, na.action = na.omit)
  robust_vcov <- sandwich::vcovHC(fit, type = "HC3")
  coefficients <- coef(fit)

  group_effect <- function(group_name) {
    base_name <- "exposure_flag"
    interaction_name <- paste0("exposure_flag:moderator_group", group_name)
    contrast <- rep(0, length(coefficients))
    names(contrast) <- names(coefficients)
    contrast[base_name] <- 1
    if (group_name != "Group_A" && interaction_name %in% names(contrast)) {
      contrast[interaction_name] <- 1
    }
    estimate <- sum(contrast * coefficients)
    standard_error <- sqrt(as.numeric(t(contrast) %*% robust_vcov %*% contrast))
    data.frame(
      outcome = outcome_name,
      moderator_group = group_name,
      n = stats::nobs(fit),
      percent_difference = 100 * (exp(estimate) - 1),
      lower_95 = 100 * (exp(estimate - 1.96 * standard_error) - 1),
      upper_95 = 100 * (exp(estimate + 1.96 * standard_error) - 1)
    )
  }

  bind_rows(lapply(levels(analysis_data$moderator_group), group_effect))
}

interaction_results <- bind_rows(lapply(outcome_names, fit_interaction_model))
write.csv(interaction_results, "outputs/tables/interaction_results.csv", row.names = FALSE)

# Correlations ----------------------------------------------------------------
correlation_variables <- c(continuous_covariates, outcome_names)
correlation_matrix <- cor(
  analysis_data[correlation_variables],
  use = "pairwise.complete.obs",
  method = "spearman"
)
write.csv(correlation_matrix, "outputs/tables/spearman_correlations.csv")

# CFA and exploratory SEM -----------------------------------------------------
latent_data <- analysis_data %>%
  mutate(
    indicator_05 = -indicator_05,
    indicator_06 = -indicator_06
  ) %>%
  mutate(across(all_of(indicator_names), ~ as.numeric(scale(.x))))

cfa_model <- paste0("factor_1 =~ ", paste(indicator_names, collapse = " + "))
cfa_fit <- lavaan::cfa(
  cfa_model,
  data = latent_data,
  estimator = "MLR",
  missing = "fiml"
)

cfa_fit_table <- data.frame(
  metric = c("cfi", "tli", "rmsea", "srmr"),
  value = as.numeric(lavaan::fitMeasures(cfa_fit, c("cfi", "tli", "rmsea", "srmr")))
)
write.csv(cfa_fit_table, "outputs/tables/cfa_fit.csv", row.names = FALSE)

cfa_loadings <- lavaan::standardizedSolution(cfa_fit) %>%
  filter(op == "=~") %>%
  transmute(
    factor = lhs,
    indicator = rhs,
    standardized_loading = est.std,
    standard_error = se,
    p_value = pvalue
  )
write.csv(cfa_loadings, "outputs/tables/cfa_loadings.csv", row.names = FALSE)

sem_data <- latent_data %>%
  mutate(
    mediator_z = as.numeric(scale(log(outcome_04))),
    covariate_01_z = as.numeric(scale(covariate_01)),
    covariate_03_z = as.numeric(scale(covariate_03)),
    covariate_05_z = as.numeric(scale(covariate_05))
  )

sem_model <- paste0(
  "factor_1 =~ ", paste(indicator_names, collapse = " + "), "\n",
  "mediator_z ~ a*exposure_flag + d1*covariate_01_z + d2*covariate_03_z + d3*covariate_05_z\n",
  "factor_1 ~ b*mediator_z + cprime*exposure_flag + e1*covariate_01_z + e2*covariate_03_z + e3*covariate_05_z\n",
  "indirect := a*b\n",
  "total := cprime + (a*b)"
)

sem_complete <- sem_data %>%
  select(
    all_of(indicator_names), mediator_z, exposure_flag,
    covariate_01_z, covariate_03_z, covariate_05_z
  ) %>%
  drop_na()

sem_fit <- lavaan::sem(
  sem_model,
  data = sem_complete,
  se = "bootstrap",
  bootstrap = 500
)

sem_paths <- lavaan::parameterEstimates(
  sem_fit,
  standardized = TRUE,
  boot.ci.type = "perc"
) %>%
  filter(label %in% c("a", "b", "cprime", "indirect", "total")) %>%
  transmute(
    path = label,
    estimate = est,
    lower_95 = ci.lower,
    upper_95 = ci.upper,
    standardized_estimate = std.all,
    p_value = pvalue
  )
write.csv(sem_paths, "outputs/tables/sem_paths.csv", row.names = FALSE)

# Figures ---------------------------------------------------------------------
outcome_long <- analysis_data %>%
  select(exposure_flag, all_of(outcome_names)) %>%
  pivot_longer(all_of(outcome_names), names_to = "outcome", values_to = "value")

distribution_plot <- ggplot(outcome_long, aes(x = value, fill = factor(exposure_flag))) +
  geom_histogram(
    aes(y = after_stat(density)), bins = 30, alpha = 0.42,
    position = "identity", na.rm = TRUE
  ) +
  geom_density(alpha = 0.15, na.rm = TRUE) +
  facet_wrap(~ outcome, scales = "free", ncol = 2) +
  scale_fill_manual(values = c("#64748b", "#0f766e"), name = "Exposure") +
  labs(title = "Synthetic outcome distributions", x = "Observed value", y = "Density") +
  theme_minimal(base_size = 11)
ggsave("outputs/figures/outcome_distributions.png", distribution_plot, width = 10, height = 8, dpi = 180)

missingness_plot <- ggplot(missingness_summary, aes(x = reorder(variable, missing_percent), y = missing_percent)) +
  geom_col(fill = "#2563eb") +
  coord_flip() +
  labs(title = "Missingness by variable", x = NULL, y = "Missing (%)") +
  theme_minimal(base_size = 11)
ggsave("outputs/figures/missingness.png", missingness_plot, width = 8, height = 5, dpi = 180)

correlation_long <- as.data.frame(as.table(correlation_matrix))
correlation_plot <- ggplot(correlation_long, aes(Var1, Var2, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#2563eb", mid = "white", high = "#dc2626", midpoint = 0, limits = c(-1, 1)) +
  coord_equal() +
  labs(title = "Spearman correlation matrix", x = NULL, y = NULL, fill = "rho") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("outputs/figures/correlation_matrix.png", correlation_plot, width = 8, height = 7, dpi = 180)

forest_plot <- ggplot(primary_results, aes(y = outcome, x = percent_difference)) +
  geom_vline(xintercept = 0, linetype = 2, color = "#64748b") +
  geom_errorbar(
    aes(xmin = lower_95, xmax = upper_95),
    width = 0.18, color = "#334155", orientation = "y"
  ) +
  geom_point(size = 2.8, color = "#0f766e") +
  labs(
    title = "Adjusted exposure associations",
    subtitle = "Log-linear models with HC3 robust 95% confidence intervals",
    x = "Percent difference", y = NULL
  ) +
  theme_minimal(base_size = 11)
ggsave("outputs/figures/adjusted_forest_plot.png", forest_plot, width = 8, height = 5, dpi = 180)

diagnostic_fit <- lm(
  log(outcome_01) ~ exposure_flag + covariate_01 + covariate_02 + covariate_03 +
    covariate_04 + covariate_05 + covariate_06,
  data = analysis_data,
  na.action = na.omit
)
diagnostic_data <- data.frame(
  fitted = fitted(diagnostic_fit),
  standardized_residual = rstandard(diagnostic_fit)
)
diagnostic_plot <- ggplot(diagnostic_data, aes(fitted, standardized_residual)) +
  geom_point(alpha = 0.35, color = "#334155") +
  geom_smooth(method = "loess", se = FALSE, color = "#dc2626") +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(title = "Primary-model diagnostic example", x = "Fitted log outcome", y = "Standardized residual") +
  theme_minimal(base_size = 11)
ggsave("outputs/figures/model_diagnostic.png", diagnostic_plot, width = 7, height = 5, dpi = 180)

cfa_plot <- ggplot(cfa_loadings, aes(x = standardized_loading, y = reorder(indicator, standardized_loading))) +
  geom_vline(xintercept = 0, linetype = 2, color = "#94a3b8") +
  geom_point(size = 3, color = "#7c3aed") +
  labs(title = "Standardized CFA loadings", x = "Standardized loading", y = NULL) +
  theme_minimal(base_size = 11)
ggsave("outputs/figures/cfa_loadings.png", cfa_plot, width = 7, height = 5, dpi = 180)

sem_plot <- ggplot(sem_paths, aes(x = estimate, y = reorder(path, estimate))) +
  geom_vline(xintercept = 0, linetype = 2, color = "#94a3b8") +
  geom_errorbar(
    aes(xmin = lower_95, xmax = upper_95),
    width = 0.18, orientation = "y"
  ) +
  geom_point(size = 3, color = "#b45309") +
  labs(title = "Exploratory SEM path estimates", x = "Estimate with bootstrap 95% interval", y = NULL) +
  theme_minimal(base_size = 11)
ggsave("outputs/figures/sem_paths.png", sem_plot, width = 7, height = 5, dpi = 180)

message("Analysis complete. Tables and figures are in outputs/.")
