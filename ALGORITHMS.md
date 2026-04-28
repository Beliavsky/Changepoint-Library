# Algorithms And Fortran Modules

This document maps the main changepoint algorithms to the Fortran modules that implement them. The R and Python scripts in this repository are primarily benchmarks and reference checks; the Fortran modules and `xsim_*` drivers are the standalone implementation side.

For a benchmark-case inventory, see [`COVERAGE.md`](COVERAGE.md). For project setup, see [`README.md`](README.md).

## Conventions

- Most routines use 1-based Fortran indexing.
- Most comparison programs print terminal endpoints in the style used by the reference package being replicated. Some APIs return internal changepoints only, while some drivers print the final endpoint as part of the segmentation.
- File-reading examples usually live in `xsim_*_file.f90` programs and use `compare_io.f90`.
- Many solver names begin with `solve_...`; these are the best entry points when using modules directly.
- This is a research library. Subroutine signatures are usable, but not yet presented as a stable packaged API.

## Core Offline Changepoint Methods

| Algorithm family | Main module | Public entry points | Example drivers |
| --- | --- | --- | --- |
| PELT, AMOC, window, bottom-up, dynamic programming, mean/variance costs | `changepoint.f90` | broad internal changepoint toolkit used by many drivers | `xsim_mean_pelt_file.f90`, `xsim_mean_amoc_file.f90`, `xsim_mean_bottomup_file.f90`, `xsim_mean_dynp_file.f90`, `xsim_window_file.f90` |
| Binary segmentation by generic cost or mean-shift cost | `changepoint_binseg.f90` | `solve_binseg_cost`, `solve_binseg_mean_shift` | `xsim_mean_binseg_file.f90` |
| Mean-variance cost comparison and metric support | `changepoint_metrics.f90` | `precision_recall_metric`, `hausdorff_metric`, `randindex_metric`, `print_metric_block` | `xmetrics_file.f90` |

The central `changepoint.f90` module is used by many standalone executables for univariate mean changes, variance changes, mean-variance changes, regression AMOC, PELT-style optimization, and related cost functions. For direct use, start from the corresponding `xsim_*_file.f90` driver and follow the `use changepoint_mod` calls.

## High-Dimensional And Regression Dynamic Programming

| Algorithm family | Main module | Public entry points | Example drivers |
| --- | --- | --- | --- |
| Regression dynamic programming and local refinement | `changepoints_pkg.f90` | module-level solvers used by DP, CV-DP, and local-refinement drivers | `xsim_dp_regression_file.f90`, `xsim_cv_dp_regression_file.f90`, `xsim_local_refine_regression_file.f90` |
| Dynamic programming with decreased update / DpDu variants | `changepoints_pkg_dpdu.f90` | module-level solvers used by DpDu drivers | `xsim_dpdu_regression_file.f90`, `xsim_dpdu2_regression_file.f90`, `xsim_cv_dpdu_regression_file.f90`, `xsim_local_refine_dpdu_regression_file.f90` |
| VAR(1)-style segmentation | `changepoints_pkg_var1.f90` | `solve_dp_var1_1d`, `solve_cv_dp_var1_1d`, `solve_local_refine_cv_var1_1d` | `xsim_dp_var1_file.f90`, `xsim_cv_dp_var1_file.f90`, `xsim_local_refine_cv_var1_file.f90` |

These modules cover more specialized dynamic-programming experiments: polynomial or regression costs, cross-validation variants, local refinement, and VAR-like one-dimensional examples.

## Sequential CPM Methods

| Algorithm family | Main module | Public entry points | Example drivers |
| --- | --- | --- | --- |
| Parametric and nonparametric CPM detectors | `cpm_pkg.f90` | `solve_cpm_detect_*_1d`, `solve_cpm_process_*_1d`, `solve_cpm_batch_*_1d` | `xsim_cpm_student_file.f90`, `xsim_cpm_joint_file.f90`, `xsim_cpm_ks_file.f90`, `xsim_cpm_cvm_file.f90` |

`cpm_pkg.f90` implements three usage modes for most statistics:

| Mode | Meaning |
| --- | --- |
| `detect` | online-style first alarm with change point and detection time |
| `process` | repeated detection over a stream, resetting after detections |
| `batch` | retrospective best split and statistic vector |

Implemented CPM statistics include Student, Bartlett, Joint, JointAdjusted, JointHawkins, Exponential, ExponentialAdjusted, Poisson, Mann-Whitney, Mood, Lepage, FET, KS, and CVM.

## Energy And Nonparametric ECP Methods

| Algorithm family | Main module | Public entry points | Example drivers |
| --- | --- | --- | --- |
| Energy divisive segmentation | `ecp_pkg.f90` | `solve_e_divisive_fixedk_2d`, `solve_e_divisive_permtest_2d` | `xsim_ecp_edivisive_file.f90`, `xsim_ecp_edivisive_full_file.f90` |
| Energy agglomerative segmentation | `ecp_pkg.f90` | `solve_e_agglo_default_2d` | `xsim_ecp_eagglo_file.f90` |
| CP3O exact / approximate energy segmentation | `ecp_pkg.f90` | `solve_e_cp3o_2d`, `solve_e_cp3o_delta_2d` | `xsim_ecp_cp3o_file.f90`, `xsim_ecp_cp3o_delta_file.f90` |

The ECP routines are matrix-oriented and generally expect observations by row. They use pairwise energy distances and return estimated segment endpoints or changepoint locations in the style of the corresponding reference method.

## Nonparametric Empirical-Distribution PELT

| Algorithm family | Main module | Public entry points | Example drivers |
| --- | --- | --- | --- |
| `changepoint.np`-style empirical-distribution PELT | `changepoint_np_pkg.f90` | `solve_cpt_np_pelt_1d`, `changepoint_np_penalty_value` | `xsim_changepointnp_file.f90` |

This module implements the empirical CDF quantile-grid cost path used by `changepoint.np::cpt.np` comparisons. It returns changepoints plus diagnostic arrays such as last-changepoint and last-likelihood paths in the driver.

## Bayesian Change Point Models

| Algorithm family | Main module | Public entry points | Example drivers |
| --- | --- | --- | --- |
| Bayesian univariate mean changepoints | `bcp_pkg.f90` | `solve_bcp_series_1d` | `xsim_bcp_file.f90` |
| Bayesian multivariate changepoints | `bcp_pkg.f90` | `solve_bcp_series_2d` | `xsim_bcp_multivar_file.f90` |
| Bayesian regression changepoints | `bcp_pkg.f90` | `solve_bcp_regression_1x`, `solve_bcp_regression_1x_probe` | `xsim_bcp_reg_file.f90`, `xsim_probe_bcp_reg_file.f90` |

The BCP routines use deterministic Wichmann-Hill seed handling in the comparison drivers so that posterior means and posterior probabilities can be compared against R package output.

## Bayesian Segmented Regression Draw Evaluation

| Algorithm family | Main module | Public entry points | Example drivers |
| --- | --- | --- | --- |
| `mcp` demo draw evaluation | `mcp_pkg.f90` | `eval_mcp_demo_draws` | `xsim_mcp_demo_file.f90` |
| Mean and variance segmented models | `mcp_pkg.f90` | `eval_mcp_sigma_draws` | `xsim_mcp_sigma_file.f90` |
| Autoregressive segmented models | `mcp_pkg.f90` | `eval_mcp_ar1_draws`, `eval_mcp_ar_change_draws` | `xsim_mcp_ar_file.f90`, `xsim_mcp_ar_change_file.f90` |
| AR plus changing variance models | `mcp_pkg.f90` | `eval_mcp_arsigma_draws` | `xsim_mcp_arsigma_file.f90` |

`mcp_pkg.f90` does not run JAGS. It evaluates posterior draw files produced by the R `mcp` workflow and computes fitted means, sigma paths, changepoint means, and parameter summaries in Fortran.

## Structural Breaks And Segmented Regression

| Algorithm family | Main module | Public entry points | Example drivers |
| --- | --- | --- | --- |
| `strucchange`-style breakpoints, F statistics, empirical fluctuation processes | `strucchange_pkg.f90` | module-level solvers used by structural-break drivers | `xsim_strucchange_breakpoints_file.f90`, `xsim_strucchange_fstats_file.f90`, `xsim_strucchange_efp_file.f90` |
| Recursive CUSUM / MOSUM and MEFP variants | `strucchange_pkg.f90` | module-level solvers used by monitoring drivers | `xsim_strucchange_rec_cusum_file.f90`, `xsim_strucchange_ols_mosum_file.f90`, `xsim_strucchange_mefp_file.f90` |
| Linear segmented regression and step changes | `segmented_pkg.f90` | `eval_segmented_lm1_draws`, `eval_segmented_lm2_draws`, `eval_stepmented_lm1_draws` | `xsim_segmented_file.f90`, `xsim_segmented_multi_file.f90`, `xsim_stepmented_file.f90` |

`strucchange_pkg.f90` depends on `pca_jacobi.f90` for symmetric eigen decompositions used in some structural-break calculations.

## Online, Subspace, Density-Ratio, And Classifier Detectors

| Algorithm family | Main module | Public entry points | Example drivers |
| --- | --- | --- | --- |
| SST and ESST subspace methods | `changepoynt.f90` | `solve_sst_naive_1d`, `solve_esst_exact_1d` | `xsim_sst_file.f90`, `xsim_esst_file.f90` |
| Gaussian RuLSIF and BOCPD mean detector | `changepoynt.f90` | `solve_rulsif_gaussian_1d`, `solve_bocpd_gaussian_mean_1d` | `xsim_rulsif_file.f90`, `xsim_bocpd_changepoynt_file.f90` |
| BOCPD variants built around shared changepoint code | `changepoint.f90` | used by BOCPD drivers | `xsim_bocpd_file.f90`, `xsim_bocpd_poisson_file.f90`, `xsim_bocpd_beta_bernoulli_file.f90`, `xsim_bocd_file.f90` |
| Classifier and density-ratio sliding-window detectors | `roerich.f90` | `solve_sliding_windows_energy_1d`, `solve_cpdc_qda_klsym_1d`, `solve_cpdc_qda_klsym_cv_1d`, `solve_rulsif_linear_pesym_1d` | `xsim_roerich_windows_file.f90`, `xsim_roerich_classifier_file.f90`, `xsim_roerich_classifier_cv_file.f90`, `xsim_roerich_rulsif_file.f90` |
| ClaSP / CLASP-like profile methods | `claspy.f90` | module-level solvers used by ClaSP drivers | `xsim_clasp_file.f90`, `xsim_clasp_alt_file.f90`, `xsim_clasp_ensemble_file.f90`, `xsim_streaming_clasp_file.f90` |

These modules are useful for time-series subsequence and distribution-shift settings where a classical mean-shift cost is not the target model.

## BEAST-Style Piecewise Trend And Seasonal Models

| Algorithm family | Main module | Public entry points | Example drivers |
| --- | --- | --- | --- |
| Trend-only BEAST-style segmentation | `rbeast.f90` | `solve_beast_trend_only_1d` | `xsim_beast_file.f90` |
| Seasonal harmonic BEAST-style segmentation | `rbeast.f90` | `solve_beast_harmonic_1d` | `xsim_beast_seasonal_file.f90` |
| Irregular-time trend segmentation | `rbeast.f90` | `solve_beast_irreg_trend_only_1d` | `xsim_beast_irreg_file.f90` |
| Multivariate regular BEAST123-style example | `rbeast.f90` | `solve_beast123_regular_mv` | `xsim_beast123_file.f90` |

These routines use dynamic programming over piecewise linear trend fits, optional harmonic components, and posterior-like changepoint probability summaries in the comparison drivers.

## Utility Modules

| Module | Purpose |
| --- | --- |
| `kind.f90` | shared numeric kind definitions |
| `util.f90` | display, timing, string, cumulative-sum, and file helpers |
| `compare_io.f90` | simple series and matrix readers for comparison executables |
| `compare_sim.f90` | simulation-comparison helpers |
| `basic_stats.f90` | mean, variance, standard deviation, kurtosis, correlations, covariance summaries |
| `median.f90` | median helper |
| `arma.f90` | AR/MA/ARMA fitting and prediction helpers used by some experiments |
| `pca_jacobi.f90` | Jacobi symmetric eigensolver and PCA helper |
| `random.f90` | random-number utilities |
| `dataframe_index_date.f90`, `df_index_date_ops_mod.f90`, `date.f90`, `io_utils.f90` | date-indexed data-frame and CSV utilities used by return-series examples |

## How To Use A Module Directly

The fastest route is to start from a matching driver:

1. Pick an algorithm family in this file.
2. Open the listed `xsim_*_file.f90` driver.
3. Copy the `use ..._mod` line and the solver call.
4. Replace the `compare_io` file-reading block with your own array construction.
5. Preserve the same indexing convention when interpreting changepoints.

For example, to use the CPM KS batch detector, inspect `xsim_cpm_ks_file.f90`, which calls the KS routines in `cpm_pkg.f90`. To use energy divisive segmentation on a matrix, inspect `xsim_ecp_edivisive_file.f90`, which calls `solve_e_divisive_fixedk_2d` from `ecp_pkg.f90`.

## Build Notes

Most modules can be built through the top-level `Makefile`. Build all registered executables with:

```powershell
make
```

Build one example driver with:

```powershell
make xsim_cpm_ks_file
```

If you compile a standalone program manually, include every module object needed by the driver. For example, `strucchange_pkg.o` requires `pca_jacobi.o`, and most comparison drivers also require `kind.o`, `util.o`, and either `compare_io.o` or `compare_sim.o`.
