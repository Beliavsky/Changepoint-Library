# Coverage

This file lists the comparison cases registered in `xrun_compare.py`. Each case usually has a deterministic data generator, an R or Python reference implementation, and a Fortran executable that prints comparable summaries or checksums.

Build the full executable set with:

```powershell
make
```

Build one executable directly with:

```powershell
make xsim_changepointnp_file
```

Run one comparison case with:

```powershell
python xrun_compare.py xchangepointnp
```

## Classical Offline Changepoint Cases

| Case | Reference scripts | Fortran executable |
| --- | --- | --- |
| `xbinseg` | `xbinseg_file.py`, `xbinseg_file.R` | `xsim_mean_binseg_file.exe` |
| `xbottomup` | `xbottomup_file.py`, `xbottomup_file.R` | `xsim_mean_bottomup_file.exe` |
| `xdynp` | `xdynp_file.py`, `xdynp_file.R` | `xsim_mean_dynp_file.exe` |
| `xmean_shift` | `xmean_shift_file.py`, `xmean_shift_file.R` | `xsim_mean_pelt_file.exe` |
| `xmultivar` | `xmultivar_file.py`, `xmultivar_file.R` | `xsim_multivar_pelt_file.exe` |
| `xwindow` | `xwindow_file.py`, `xwindow_file.R` | `xsim_window_file.exe` |
| `xamoc` | `xamoc_file.R` | `xsim_mean_amoc_file.exe` |
| `xreg` | `xreg_file.R` | `xsim_reg_amoc_file.exe` |
| `xcss` | `xcss_file.R` | `xsim_var_css_amoc_file.exe` |
| `xcusum` | `xcusum_file.R` | `xsim_mean_cusum_amoc_file.exe` |
| `xpoisson` | `xpoisson_file.R` | `xsim_meanvar_poisson_amoc_file.exe` |
| `xexp` | `xexp_file.R` | `xsim_meanvar_exp_amoc_file.exe` |
| `xgamma` | `xgamma_file.R` | `xsim_meanvar_gamma_amoc_file.exe` |
| `xtext_like` | `xtext_like_file.py`, `xtext_like_file.R` | `xsim_text_like_pelt_file.exe` |
| `xvariance_regime` | `xvariance_regime_file.py`, `xvariance_regime_file.R` | `xsim_variance_regime_pelt_file.exe` |
| `xkernel` | `xkernel_file.py`, `xkernel_file_kerSeg.R` | `xsim_kernel_file.exe` |
| `xcosts_mean_variance` | `xcosts_mean_variance_file.py`, `xcosts_mean_variance_file.R` | `xcosts_mean_variance_file.exe` |
| `xmetrics` | `xmetrics_file.py`, `xmetrics_file.R` | `xmetrics_file.exe` |

## Wild Binary Segmentation And Related Cases

| Case | Reference scripts | Fortran executable |
| --- | --- | --- |
| `xwbs` | `xwbs_file.R` | `xsim_wbs_univar_file.exe` |
| `xwbs_rob` | `xwbs_rob_file.R` | `xsim_wbs_univar_rob_file.exe` |
| `xarc` | `xarc_file.R` | `xsim_arc_univar_file.exe` |
| `xaarc` | `xaarc_file.R` | `xsim_aarc_univar_file.exe` |
| `xonline_univar` | `xonline_univar_file.R` | `xsim_online_univar_file.exe` |
| `xbs_cov` | `xbs_cov_file.R` | `xsim_bs_cov_file.exe` |
| `xwbsip_cov` | `xwbsip_cov_file.R` | `xsim_wbsip_cov_file.exe` |
| `xseeded_binseg` | `xseeded_binseg_file.py` | `xsim_seeded_binseg_file.exe` |

## Dynamic Programming, Regression, And Refinement Cases

| Case | Reference scripts | Fortran executable |
| --- | --- | --- |
| `xdp_regression` | `xdp_regression_file.R` | `xsim_dp_regression_file.exe` |
| `xdpdu_regression` | `xdpdu_regression_file.R` | `xsim_dpdu_regression_file.exe` |
| `xdpdu2_regression` | `xdpdu2_regression_file.R` | `xsim_dpdu2_regression_file.exe` |
| `xdp_poly` | `xdp_poly_file.R` | `xsim_dp_poly_file.exe` |
| `xcv_dp_poly` | `xcv_dp_poly_file.R` | `xsim_cv_dp_poly_file.exe` |
| `xlocal_refine_poly` | `xlocal_refine_poly_file.R` | `xsim_local_refine_poly_file.exe` |
| `xcv_dpdu_regression` | `xcv_dpdu_regression_file.R` | `xsim_cv_dpdu_regression_file.exe` |
| `xlocal_refine_dpdu_regression` | `xlocal_refine_dpdu_regression_file.R` | `xsim_local_refine_dpdu_regression_file.exe` |
| `xdp_var1` | `xdp_var1_file.R` | `xsim_dp_var1_file.exe` |
| `xcv_dp_var1` | `xcv_dp_var1_file.R` | `xsim_cv_dp_var1_file.exe` |
| `xlocal_refine_cv_var1` | `xlocal_refine_cv_var1_file.R` | `xsim_local_refine_cv_var1_file.exe` |
| `xci_regression` | `xci_regression_file.R` | `xsim_ci_regression_file.exe` |
| `xcv_dp_regression` | `xcv_dp_regression_file.R` | `xsim_cv_dp_regression_file.exe` |
| `xlocal_refine_regression` | `xlocal_refine_regression_file.R` | `xsim_local_refine_regression_file.exe` |
| `xlocal_refine_univar` | `xlocal_refine_univar_file.R` | `xsim_local_refine_univar_file.exe` |
| `xdp_univar` | `xdp_univar_file.R` | `xsim_dp_univar_file.exe` |
| `xcv_dp_univar` | `xcv_dp_univar_file.R` | `xsim_cv_dp_univar_file.exe` |

## `strucchange`-Style Structural Break Cases

| Case | Reference script | Fortran executable |
| --- | --- | --- |
| `xstrucchange_breakpoints` | `xstrucchange_breakpoints_file.R` | `xsim_strucchange_breakpoints_file.exe` |
| `xstrucchange_fstats` | `xstrucchange_fstats_file.R` | `xsim_strucchange_fstats_file.exe` |
| `xstrucchange_efp` | `xstrucchange_efp_file.R` | `xsim_strucchange_efp_file.exe` |
| `xstrucchange_rec_cusum` | `xstrucchange_rec_cusum_file.R` | `xsim_strucchange_rec_cusum_file.exe` |
| `xstrucchange_ols_mosum` | `xstrucchange_ols_mosum_file.R` | `xsim_strucchange_ols_mosum_file.exe` |
| `xstrucchange_rec_mosum` | `xstrucchange_rec_mosum_file.R` | `xsim_strucchange_rec_mosum_file.exe` |
| `xstrucchange_mefp` | `xstrucchange_mefp_file.R` | `xsim_strucchange_mefp_file.exe` |
| `xstrucchange_mefp_ols_mosum` | `xstrucchange_mefp_ols_mosum_file.R` | `xsim_strucchange_mefp_ols_mosum_file.exe` |
| `xstrucchange_re` | `xstrucchange_re_file.R` | `xsim_strucchange_re_file.exe` |
| `xstrucchange_mefp_re` | `xstrucchange_mefp_re_file.R` | `xsim_strucchange_mefp_re_file.exe` |
| `xstrucchange_me` | `xstrucchange_me_file.R` | `xsim_strucchange_me_file.exe` |
| `xstrucchange_mefp_me` | `xstrucchange_mefp_me_file.R` | `xsim_strucchange_mefp_me_file.exe` |

## CPM Cases

| Case | Reference script | Fortran executable |
| --- | --- | --- |
| `xcpm_student` | `xcpm_student_file.R` | `xsim_cpm_student_file.exe` |
| `xcpm_bartlett` | `xcpm_bartlett_file.R` | `xsim_cpm_bartlett_file.exe` |
| `xcpm_joint` | `xcpm_joint_file.R` | `xsim_cpm_joint_file.exe` |
| `xcpm_joint_adjusted` | `xcpm_joint_adjusted_file.R` | `xsim_cpm_joint_adjusted_file.exe` |
| `xcpm_joint_hawkins` | `xcpm_joint_hawkins_file.R` | `xsim_cpm_joint_hawkins_file.exe` |
| `xcpm_exponential` | `xcpm_exponential_file.R` | `xsim_cpm_exponential_file.exe` |
| `xcpm_exponential_adjusted` | `xcpm_exponential_adjusted_file.R` | `xsim_cpm_exponential_adjusted_file.exe` |
| `xcpm_poisson` | `xcpm_poisson_file.R` | `xsim_cpm_poisson_file.exe` |
| `xcpm_mw` | `xcpm_mw_file.R` | `xsim_cpm_mw_file.exe` |
| `xcpm_mood` | `xcpm_mood_file.R` | `xsim_cpm_mood_file.exe` |
| `xcpm_lepage` | `xcpm_lepage_file.R` | `xsim_cpm_lepage_file.exe` |
| `xcpm_fet` | `xcpm_fet_file.R` | `xsim_cpm_fet_file.exe` |
| `xcpm_ks` | `xcpm_ks_file.R` | `xsim_cpm_ks_file.exe` |
| `xcpm_cvm` | `xcpm_cvm_file.R` | `xsim_cpm_cvm_file.exe` |

## `ecp` Cases

| Case | Reference script | Fortran executable |
| --- | --- | --- |
| `xecp_edivisive` | `xecp_edivisive_file.R` | `xsim_ecp_edivisive_file.exe` |
| `xecp_edivisive_full` | `xecp_edivisive_full_file.R` | `xsim_ecp_edivisive_full_file.exe` |
| `xecp_eagglo` | `xecp_eagglo_file.R` | `xsim_ecp_eagglo_file.exe` |
| `xecp_cp3o` | `xecp_cp3o_file.R` | `xsim_ecp_cp3o_file.exe` |
| `xecp_cp3o_delta` | `xecp_cp3o_delta_file.R` | `xsim_ecp_cp3o_delta_file.exe` |

## Bayesian And Nonparametric R Package Cases

| Case | Reference script | Fortran executable |
| --- | --- | --- |
| `xchangepointnp` | `xchangepointnp_file.R` | `xsim_changepointnp_file.exe` |
| `xbcp` | `xbcp_file.R` | `xsim_bcp_file.exe` |
| `xbcp_multivar` | `xbcp_multivar_file.R` | `xsim_bcp_multivar_file.exe` |
| `xbcp_reg` | `xbcp_reg_file.R` | `xsim_bcp_reg_file.exe` |

## `mcp` Cases

| Case | Reference script | Fortran executable |
| --- | --- | --- |
| `xmcp_demo` | `xmcp_demo_file.R` | `xsim_mcp_demo_file.exe` |
| `xmcp_fit` | `xmcp_fit_file.R` | `xsim_mcp_fit_file.exe` |
| `xmcp_sigma` | `xmcp_sigma_file.R` | `xsim_mcp_sigma_file.exe` |
| `xmcp_ar` | `xmcp_ar_file.R` | `xsim_mcp_ar_file.exe` |
| `xmcp_ar_change` | `xmcp_ar_change_file.R` | `xsim_mcp_ar_change_file.exe` |
| `xmcp_arsigma` | `xmcp_arsigma_file.R` | `xsim_mcp_arsigma_file.exe` |

## `segmented` Cases

| Case | Reference script | Fortran executable |
| --- | --- | --- |
| `xsegmented` | `xsegmented_file.R` | `xsim_segmented_file.exe` |
| `xsegmented_multi` | `xsegmented_multi_file.R` | `xsim_segmented_multi_file.exe` |
| `xstepmented` | `xstepmented_file.R` | `xsim_stepmented_file.exe` |

## Online, Bayesian Online, And Density-Ratio Cases

| Case | Reference script | Fortran executable |
| --- | --- | --- |
| `xbocpd` | `xbocpd_file.py` | `xsim_bocpd_file.exe` |
| `xbocpd_poisson` | `xbocpd_poisson_file.py` | `xsim_bocpd_poisson_file.exe` |
| `xbocpd_beta_bernoulli` | `xbocpd_beta_bernoulli_file.py` | `xsim_bocpd_beta_bernoulli_file.exe` |
| `xbocd` | `xbocd_file.py` | `xsim_bocd_file.exe` |
| `xbocpd_changepoynt` | `xbocpd_changepoynt_file.py` | `xsim_bocpd_changepoynt_file.exe` |
| `xsst` | `xsst_file.py` | `xsim_sst_file.exe` |
| `xesst` | `xesst_file.py` | `xsim_esst_file.exe` |
| `xrulsif` | `xrulsif_file.py` | `xsim_rulsif_file.exe` |
| `xulsif` | `xulsif_file.py` | `xsim_ulsif_file.exe` |
| `xargpcpd` | `xargpcpd_file.py` | `xsim_argpcpd_gp_file.exe` |
| `xargpcpd_rust` | `xargpcpd_file.py` | `xsim_argpcpd_gp_file.exe` |

## Anomaly, Ensemble, ClaSP, Roerich, And BEAST Cases

| Case | Reference script | Fortran executable |
| --- | --- | --- |
| `xcrops` | `xcrops_file.py` | `xsim_crops_file.exe` |
| `xcapa` | `xcapa_file.py` | `xsim_capa_file.exe` |
| `xclass` | `xclass_file.py` | `xsim_class_file.exe` |
| `xclap` | `xclap_file.py` | `xsim_clap_file.exe` |
| `xclasp_ensemble` | `xclasp_ensemble_file.py` | `xsim_clasp_ensemble_file.exe` |
| `xagglomerative_clap` | `xagglomerative_clap_file.py` | `xsim_agglomerative_clap_file.exe` |
| `xroerich_windows` | `xroerich_windows_file.py` | `xsim_roerich_windows_file.exe` |
| `xroerich_energy` | `xroerich_energy_file.py` | `xsim_roerich_energy_file.exe` |
| `xroerich_classifier` | `xroerich_classifier_file.py` | `xsim_roerich_classifier_file.exe` |
| `xroerich_classifier_cv` | `xroerich_classifier_cv_file.py` | `xsim_roerich_classifier_cv_file.exe` |
| `xroerich_rulsif` | `xroerich_rulsif_file.py` | `xsim_roerich_rulsif_file.exe` |
| `xbeast` | `xbeast_file.py` | `xsim_beast_file.exe` |
| `xbeast_seasonal` | `xbeast_seasonal_file.py` | `xsim_beast_seasonal_file.exe` |
| `xbeast_irreg` | `xbeast_irreg_file.py` | `xsim_beast_irreg_file.exe` |
| `xbeast123` | `xbeast123_file.py` | `xsim_beast123_file.exe` |
| `xstreaming_clasp` | `xstreaming_clasp_file.py` | `xsim_streaming_clasp_file.exe` |
| `xclasp` | `xclasp_file.py` | `xsim_clasp_file.exe` |
| `xclasp_changepoynt` | `xclasp_changepoynt_file.py` | `xsim_clasp_changepoynt_file.exe` |
| `xclasp_alt` | `xclasp_alt_file.py` | `xsim_clasp_alt_file.exe` |

## Notes

- The comparison runner may also create deterministic input files such as `*_data.txt`, `*_draws.txt`, or package-specific permutation files.
- Some cases require optional R packages or Python packages that are not needed for the Fortran build itself.
- `xrun_compare.py --all` runs the registered cases in order, but for practical use it is often better to run one family at a time.
