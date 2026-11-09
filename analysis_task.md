# Analysis Task

## Scenario

You receive a synthetic cross-sectional dataset containing one masked binary exposure, five continuous outcomes, six candidate covariates, a three-level grouping variable, and six indicators intended to measure a common latent construct.

The assignment is to build a transparent analysis that another analyst can reproduce and review.

## Questions

1. What are the main data-quality issues, distributional features, and missing-data patterns?
2. How do participant characteristics differ by `exposure_flag`?
3. Is `exposure_flag` associated with each continuous outcome after covariate adjustment?
4. Are the estimates robust to a Gamma model and nonlinear adjustment for `covariate_01`?
5. Does the exposure association differ across `moderator_group`?
6. Do the six generic indicators support a one-factor measurement model?
7. In an exploratory SEM, is there evidence consistent with an indirect pathway through `outcome_04`?

## Required deliverables

- a reproducible analysis script;
- a data dictionary and explicit quality-control checks;
- descriptive tables and distribution plots;
- adjusted effect estimates with uncertainty intervals;
- regression diagnostics and sensitivity analyses;
- CFA and exploratory SEM results;
- a concise report that separates results from interpretation and limitations.

## Analytical expectations

- Pre-specify transformations before interpreting model results.
- Treat extreme values as observations to investigate, not records to delete automatically.
- Use robust uncertainty estimates for the primary log-linear models.
- Describe interaction estimates on a scale that readers can interpret.
- Present latent-variable analyses as measurement and exploratory pathway models, not as proof of causation.
