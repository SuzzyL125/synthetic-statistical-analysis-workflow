set.seed(481927)

n <- 840L

moderator_group <- sample(
  c("Group_A", "Group_B", "Group_C"),
  size = n,
  replace = TRUE,
  prob = c(0.44, 0.33, 0.23)
)

covariate_01 <- round(pmin(pmax(rnorm(n, 50, 13), 18), 84), 1)
covariate_02 <- rbinom(n, 1, plogis(-0.15 + 0.012 * (covariate_01 - 50)))
covariate_03 <- round(rnorm(n, 0.2 * covariate_02, 1), 3)
covariate_04 <- sample(1:3, n, replace = TRUE, prob = c(0.42, 0.36, 0.22))
covariate_05 <- round(rgamma(n, shape = 3.2, rate = 0.8), 3)
covariate_06 <- rbinom(n, 1, 0.31)

exposure_probability <- plogis(
  -0.55 +
    0.018 * (covariate_01 - 50) +
    0.38 * covariate_02 -
    0.14 * covariate_03 +
    0.25 * (moderator_group == "Group_C")
)
exposure_flag <- rbinom(n, 1, exposure_probability)

group_shift <- c(Group_A = 0, Group_B = 0.05, Group_C = -0.04)[moderator_group]
base_levels <- log(c(18, 42, 9, 27, 64))
exposure_effects <- c(0.07, -0.04, 0.11, 0.06, -0.08)
error_scales <- c(0.34, 0.29, 0.48, 0.37, 0.31)

outcomes <- vapply(seq_len(5), function(j) {
  linear_predictor <-
    base_levels[j] +
    exposure_effects[j] * exposure_flag +
    0.006 * (covariate_01 - 50) +
    0.08 * covariate_02 -
    0.035 * covariate_03 +
    0.025 * (covariate_04 - 1) +
    0.012 * covariate_05 +
    0.04 * covariate_06 +
    group_shift
  exp(linear_predictor + rnorm(n, 0, error_scales[j]))
}, numeric(n))

colnames(outcomes) <- sprintf("outcome_%02d", seq_len(5))

latent_construct <-
  0.18 * exposure_flag -
  0.012 * (covariate_01 - 50) +
  0.14 * covariate_03 -
  0.10 * covariate_06 +
  0.12 * scale(log(outcomes[, 4]))[, 1] +
  rnorm(n, 0, 0.8)

loading_values <- c(0.82, 0.76, 0.70, 0.65, -0.74, -0.68)
indicator_errors <- c(0.48, 0.55, 0.62, 0.69, 0.57, 0.66)
indicators <- vapply(seq_len(6), function(j) {
  loading_values[j] * latent_construct + rnorm(n, 0, indicator_errors[j])
}, numeric(n))
colnames(indicators) <- sprintf("indicator_%02d", seq_len(6))

synthetic_data <- data.frame(
  record_id = sprintf("SYN%04d", seq_len(n)),
  exposure_flag = exposure_flag,
  moderator_group = moderator_group,
  covariate_01 = covariate_01,
  covariate_02 = covariate_02,
  covariate_03 = covariate_03,
  covariate_04 = covariate_04,
  covariate_05 = covariate_05,
  covariate_06 = covariate_06,
  outcomes,
  indicators,
  check.names = FALSE
)

for (j in seq_len(5)) {
  missing_probability <-
    0.025 +
    0.025 * (moderator_group == "Group_C") +
    0.018 * covariate_02 +
    0.008 * j
  missing_index <- runif(n) < missing_probability
  synthetic_data[missing_index, sprintf("outcome_%02d", j)] <- NA_real_
}

for (j in seq_len(6)) {
  missing_probability <- 0.02 + 0.012 * (moderator_group == "Group_B") + 0.003 * j
  missing_index <- runif(n) < missing_probability
  synthetic_data[missing_index, sprintf("indicator_%02d", j)] <- NA_real_
}

dir.create("data", recursive = TRUE, showWarnings = FALSE)
write.csv(
  synthetic_data,
  file = "data/synthetic_analysis_data.csv",
  row.names = FALSE,
  na = ""
)

message("Wrote ", nrow(synthetic_data), " synthetic records.")

