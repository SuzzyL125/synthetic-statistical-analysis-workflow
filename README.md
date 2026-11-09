# Synthetic Statistical Analysis Workflow

This repository is a self-contained portfolio example of an end-to-end statistical analysis in R. It starts with reproducible data generation, moves through quality control and exploratory analysis, and finishes with regression modeling, sensitivity analyses, latent-variable modeling, and a rendered report.

**Privacy statement:** All records, variable names, data-generating mechanisms, estimates, figures, and conclusions in this repository are synthetic. Nothing is derived from a confidential employer, interview, clinical, academic, or research dataset.

## What the project demonstrates

- a documented synthetic-data generator with a fixed random seed;
- programmatic checks for duplicate identifiers, invalid values, and missingness;
- descriptive summaries overall and by a masked binary exposure;
- log-linear regression with HC3 robust standard errors;
- Gamma regression and nonlinear covariate adjustment as sensitivity analyses;
- effect-modification estimates across masked groups;
- confirmatory factor analysis (CFA) for a generic latent construct;
- exploratory structural equation modeling (SEM);
- publication-ready tables and figures; and
- a reproducible HTML report.

## Repository structure

```text
.
├── R/
│   ├── generate_synthetic_data.R
│   └── run_analysis.R
├── data/
│   └── synthetic_analysis_data.csv
├── outputs/
│   ├── figures/
│   └── tables/
├── report/
│   ├── analysis_report.Rmd
│   └── analysis_report.html
├── analysis_task.md
├── data_dictionary.md
└── render_report.R
```

## Reproduce the project

From the repository root, run:

```r
source("R/generate_synthetic_data.R")
source("R/run_analysis.R")
source("render_report.R")
```

Required packages are `dplyr`, `tidyr`, `ggplot2`, `broom`, `sandwich`, `lmtest`, `splines`, `lavaan`, `rmarkdown`, and `knitr`.

## Interpretation boundary

The analysis is designed to show workflow and reasoning, not to support a substantive claim. The synthetic associations are deliberately generated for demonstration. They should not be interpreted as clinical, policy, or causal evidence.

