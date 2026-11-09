# Synthetic Data Dictionary

| Variable | Type | Description |
|---|---|---|
| `record_id` | character | Synthetic unique record identifier |
| `exposure_flag` | integer | Masked binary exposure, coded 0/1 |
| `moderator_group` | factor | Masked grouping variable with levels Group_A, Group_B, and Group_C |
| `covariate_01` | numeric | Continuous adjustment variable |
| `covariate_02` | integer | Binary adjustment variable, coded 0/1 |
| `covariate_03` | numeric | Continuous adjustment variable |
| `covariate_04` | integer | Three-level adjustment variable, coded 1/2/3 |
| `covariate_05` | numeric | Continuous adjustment variable |
| `covariate_06` | integer | Binary adjustment variable, coded 0/1 |
| `outcome_01`–`outcome_05` | numeric | Positive, right-skewed continuous outcomes |
| `indicator_01`–`indicator_06` | numeric | Generic continuous indicators of one simulated latent construct |

## Missingness

Missing values are introduced synthetically and vary modestly across outcomes and indicators. The mechanism is intended to make quality-control and complete-case decisions visible; it does not represent any real data-generating system.

