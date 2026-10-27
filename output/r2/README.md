# IPSM R2 refit: sampling month and repeated measures

Produced by `R/10_repeated_measures_refit.R`. Nothing in the R1 pipeline (`R/01`–`R/09`, `output/tables`, `output/figures*`) is modified; the R1 state is tagged `ipsm-R1`.

Primary specification (AE request), per species x response:

- null: `y ~ plot + month + (1 | station)`
- disease: `y ~ plot + month + z_ani + (1 | station)`
- Length and root:shoot: Gaussian LMM (lme4/lmerTest), AICc from ML fits, slope t-test from REML with Satterthwaite df.
- Germination: beta-binomial GLMM (glmmTMB) on seed counts, replacing quasibinomial/QAICc.

Checks: B = `plot * month` fixed (slope from within plot-date spatial gradient only); C = germination as binomial GLMM with observation-level random effect.

| File | Contents |
|---|---|
| `regression_r2_primary.csv` | Primary fits (new Table 1) |
| `regression_r2_checks.csv` | Checks B and C |
| `comparison_R1_vs_R2.csv` | R1 vs R2 slope, p, ΔAICc side by side |
| `sensitivity_r2.csv` / `_summary.csv` | 18-kernel sweep under the primary specification (216 fits) |
| `index_confounding.csv` | Share of index variance explained by plot, month, plot + month, plot x month |
| `index_by_plot_month.csv` | Index mean and SD per plot-date cell |
| `sessionInfo.txt` | R and package versions used |

`singular = TRUE` means the station variance was estimated at or near zero; the fit then reduces to the corresponding fixed-effects model with month.

Committed outputs were generated under R 4.5.3 (lme4 2.0.6, lmerTest 3.2.1, glmmTMB 1.1.14; full list in `sessionInfo.txt`). A run under R 4.3.3 gave the same primary results to <1e-11; the largest difference anywhere was 0.0013 in one check-model p value.
