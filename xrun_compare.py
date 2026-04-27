#!/usr/bin/env python
"""
Generic comparison driver for workflows that study the same data file in multiple languages.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path


CASES = {
    "xbinseg": {
        "data_file": "xbinseg_data.txt",
        "steps": [
            {
                "label": "xbinseg_make_data.py",
                "cmd": ["python", "xbinseg_make_data.py", "{data_file}"],
            },
            {
                "label": "xbinseg_file.py",
                "cmd": ["python", "xbinseg_file.py", "{data_file}"],
            },
            {
                "label": "xbinseg_file.R",
                "cmd": ["Rscript", "xbinseg_file.R", "{data_file}"],
            },
            {
                "label": "xsim_mean_binseg_file.exe",
                "cmd": ["xsim_mean_binseg_file.exe", "{data_file}"],
                "requires": ["xsim_mean_binseg_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_mean_binseg_file "
                    "kind.o util.o changepoint.o changepoint_binseg.o xsim_mean_binseg_file.o"
                ),
            },
        ],
    },
    "xbottomup": {
        "data_file": "xbottomup_data.txt",
        "steps": [
            {
                "label": "xbottomup_make_data.py",
                "cmd": ["python", "xbottomup_make_data.py", "{data_file}"],
            },
            {
                "label": "xbottomup_file.py",
                "cmd": ["python", "xbottomup_file.py", "{data_file}"],
            },
            {
                "label": "xbottomup_file.R",
                "cmd": ["Rscript", "xbottomup_file.R", "{data_file}"],
            },
            {
                "label": "xsim_mean_bottomup_file.exe",
                "cmd": ["xsim_mean_bottomup_file.exe", "{data_file}"],
                "requires": ["xsim_mean_bottomup_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_mean_bottomup_file "
                    "kind.o util.o changepoint.o xsim_mean_bottomup_file.o"
                ),
            },
        ],
    },
    "xdynp": {
        "data_file": "xdynp_data.txt",
        "steps": [
            {
                "label": "xdynp_make_data.py",
                "cmd": ["python", "xdynp_make_data.py", "{data_file}"],
            },
            {
                "label": "xdynp_file.py",
                "cmd": ["python", "xdynp_file.py", "{data_file}"],
            },
            {
                "label": "xdynp_file.R",
                "cmd": ["Rscript", "xdynp_file.R", "{data_file}"],
            },
            {
                "label": "xsim_mean_dynp_file.exe",
                "cmd": ["xsim_mean_dynp_file.exe", "{data_file}"],
                "requires": ["xsim_mean_dynp_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_mean_dynp_file "
                    "kind.o util.o changepoint.o xsim_mean_dynp_file.o"
                ),
            },
        ],
    },
    "xmean_shift": {
        "data_file": "xmean_shift_data.txt",
        "steps": [
            {
                "label": "xmean_shift_make_data.py",
                "cmd": ["python", "xmean_shift_make_data.py", "{data_file}"],
            },
            {
                "label": "xmean_shift_file.py",
                "cmd": ["python", "xmean_shift_file.py", "{data_file}"],
            },
            {
                "label": "xmean_shift_file.R",
                "cmd": ["Rscript", "xmean_shift_file.R", "{data_file}"],
            },
            {
                "label": "xsim_mean_pelt_file.exe",
                "cmd": ["xsim_mean_pelt_file.exe", "{data_file}"],
                "requires": ["xsim_mean_pelt_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_mean_pelt_file "
                    "kind.o util.o changepoint.o xsim_mean_pelt_file.o"
                ),
            },
        ],
    },
    "xmultivar": {
        "data_file": "xmultivar_data.txt",
        "steps": [
            {
                "label": "xmultivar_make_data.py",
                "cmd": ["python", "xmultivar_make_data.py", "{data_file}"],
            },
            {
                "label": "xmultivar_file.py",
                "cmd": ["python", "xmultivar_file.py", "{data_file}"],
            },
            {
                "label": "xmultivar_file.R",
                "cmd": ["Rscript", "xmultivar_file.R", "{data_file}"],
            },
            {
                "label": "xsim_multivar_pelt_file.exe",
                "cmd": ["xsim_multivar_pelt_file.exe", "{data_file}"],
                "requires": ["xsim_multivar_pelt_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_multivar_pelt_file "
                    "kind.o util.o changepoint.o xsim_multivar_pelt_file.o"
                ),
            },
        ],
    },
    "xwindow": {
        "data_file": "xwindow_data.txt",
        "steps": [
            {
                "label": "xwindow_make_data.py",
                "cmd": ["python", "xwindow_make_data.py", "{data_file}"],
            },
            {
                "label": "xwindow_file.py",
                "cmd": ["python", "xwindow_file.py", "{data_file}"],
            },
            {
                "label": "xwindow_file.R",
                "cmd": ["Rscript", "xwindow_file.R", "{data_file}"],
            },
            {
                "label": "xsim_window_file.exe",
                "cmd": ["xsim_window_file.exe", "{data_file}"],
                "requires": ["xsim_window_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_window_file "
                    "kind.o util.o changepoint.o xsim_window_file.o"
                ),
            },
        ],
    },
    "xamoc": {
        "data_file": "xamoc_data.txt",
        "steps": [
            {
                "label": "xamoc_make_data.py",
                "cmd": ["python", "xamoc_make_data.py", "{data_file}"],
            },
            {
                "label": "xamoc_file.R",
                "cmd": ["Rscript", "xamoc_file.R", "{data_file}"],
            },
            {
                "label": "xsim_mean_amoc_file.exe",
                "cmd": ["xsim_mean_amoc_file.exe", "{data_file}"],
                "requires": ["xsim_mean_amoc_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_mean_amoc_file "
                    "kind.o util.o median.o arma.o compare_io.o changepoint.o xsim_mean_amoc_file.o"
                ),
            },
        ],
    },
    "xreg": {
        "data_file": "xreg_data.txt",
        "steps": [
            {
                "label": "xreg_make_data.py",
                "cmd": ["python", "xreg_make_data.py", "{data_file}"],
            },
            {
                "label": "xreg_file.R",
                "cmd": ["Rscript", "xreg_file.R", "{data_file}"],
            },
            {
                "label": "xsim_reg_amoc_file.exe",
                "cmd": ["xsim_reg_amoc_file.exe", "{data_file}"],
                "requires": ["xsim_reg_amoc_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_reg_amoc_file "
                    "kind.o util.o median.o arma.o compare_io.o changepoint.o xsim_reg_amoc_file.o"
                ),
            },
        ],
    },
    "xcss": {
        "data_file": "xcss_data.txt",
        "steps": [
            {
                "label": "xcss_make_data.py",
                "cmd": ["python", "xcss_make_data.py", "{data_file}"],
            },
            {
                "label": "xcss_file.R",
                "cmd": ["Rscript", "xcss_file.R", "{data_file}"],
            },
            {
                "label": "xsim_var_css_amoc_file.exe",
                "cmd": ["xsim_var_css_amoc_file.exe", "{data_file}"],
                "requires": ["xsim_var_css_amoc_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_var_css_amoc_file "
                    "kind.o util.o median.o arma.o compare_io.o changepoint.o xsim_var_css_amoc_file.o"
                ),
            },
        ],
    },
    "xcusum": {
        "data_file": "xcusum_data.txt",
        "steps": [
            {
                "label": "xcusum_make_data.py",
                "cmd": ["python", "xcusum_make_data.py", "{data_file}"],
            },
            {
                "label": "xcusum_file.R",
                "cmd": ["Rscript", "xcusum_file.R", "{data_file}"],
            },
            {
                "label": "xsim_mean_cusum_amoc_file.exe",
                "cmd": ["xsim_mean_cusum_amoc_file.exe", "{data_file}"],
                "requires": ["xsim_mean_cusum_amoc_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_mean_cusum_amoc_file "
                    "kind.o util.o median.o arma.o compare_io.o changepoint.o xsim_mean_cusum_amoc_file.o"
                ),
            },
        ],
    },
    "xpoisson": {
        "data_file": "xpoisson_data.txt",
        "steps": [
            {
                "label": "xpoisson_make_data.py",
                "cmd": ["python", "xpoisson_make_data.py", "{data_file}"],
            },
            {
                "label": "xpoisson_file.R",
                "cmd": ["Rscript", "xpoisson_file.R", "{data_file}"],
            },
            {
                "label": "xsim_meanvar_poisson_amoc_file.exe",
                "cmd": ["xsim_meanvar_poisson_amoc_file.exe", "{data_file}"],
                "requires": ["xsim_meanvar_poisson_amoc_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_meanvar_poisson_amoc_file "
                    "kind.o util.o median.o arma.o compare_io.o changepoint.o xsim_meanvar_poisson_amoc_file.o"
                ),
            },
        ],
    },
    "xexp": {
        "data_file": "xexp_data.txt",
        "steps": [
            {
                "label": "xexp_make_data.py",
                "cmd": ["python", "xexp_make_data.py", "{data_file}"],
            },
            {
                "label": "xexp_file.R",
                "cmd": ["Rscript", "xexp_file.R", "{data_file}"],
            },
            {
                "label": "xsim_meanvar_exp_amoc_file.exe",
                "cmd": ["xsim_meanvar_exp_amoc_file.exe", "{data_file}"],
                "requires": ["xsim_meanvar_exp_amoc_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_meanvar_exp_amoc_file "
                    "kind.o util.o median.o arma.o compare_io.o changepoint.o xsim_meanvar_exp_amoc_file.o"
                ),
            },
        ],
    },
    "xgamma": {
        "data_file": "xgamma_data.txt",
        "steps": [
            {
                "label": "xgamma_make_data.py",
                "cmd": ["python", "xgamma_make_data.py", "{data_file}"],
            },
            {
                "label": "xgamma_file.R",
                "cmd": ["Rscript", "xgamma_file.R", "{data_file}"],
            },
            {
                "label": "xsim_meanvar_gamma_amoc_file.exe",
                "cmd": ["xsim_meanvar_gamma_amoc_file.exe", "{data_file}"],
                "requires": ["xsim_meanvar_gamma_amoc_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_meanvar_gamma_amoc_file "
                    "kind.o util.o median.o arma.o compare_io.o changepoint.o xsim_meanvar_gamma_amoc_file.o"
                ),
            },
        ],
    },
    "xwbs": {
        "data_file": "xwbs_data.txt",
        "steps": [
            {
                "label": "xwbs_make_data.py",
                "cmd": ["python", "xwbs_make_data.py", "{data_file}"],
            },
            {
                "label": "xwbs_file.R",
                "cmd": ["Rscript", "xwbs_file.R", "{data_file}"],
            },
            {
                "label": "xsim_wbs_univar_file.exe",
                "cmd": ["xsim_wbs_univar_file.exe", "{data_file}"],
                "requires": ["xsim_wbs_univar_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_wbs_univar_file "
                    "kind.o util.o median.o arma.o compare_io.o changepoints_pkg.o xsim_wbs_univar_file.o"
                ),
            },
        ],
    },
    "xwbs_rob": {
        "data_file": "xwbs_rob_data.txt",
        "steps": [
            {
                "label": "xwbs_rob_make_data.py",
                "cmd": ["python", "xwbs_rob_make_data.py", "{data_file}"],
            },
            {
                "label": "xwbs_rob_file.R",
                "cmd": ["Rscript", "xwbs_rob_file.R", "{data_file}"],
            },
            {
                "label": "xsim_wbs_univar_rob_file.exe",
                "cmd": ["xsim_wbs_univar_rob_file.exe", "{data_file}"],
                "requires": ["xsim_wbs_univar_rob_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_wbs_univar_rob_file "
                    "kind.o util.o median.o arma.o compare_io.o changepoints_pkg.o xsim_wbs_univar_rob_file.o"
                ),
            },
        ],
    },
    "xarc": {
        "data_file": "xarc_data.txt",
        "steps": [
            {
                "label": "xarc_make_data.py",
                "cmd": ["python", "xarc_make_data.py", "{data_file}"],
            },
            {
                "label": "xarc_file.R",
                "cmd": ["Rscript", "xarc_file.R", "{data_file}"],
            },
            {
                "label": "xsim_arc_univar_file.exe",
                "cmd": ["xsim_arc_univar_file.exe", "{data_file}"],
                "requires": ["xsim_arc_univar_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_arc_univar_file "
                    "kind.o util.o median.o arma.o compare_io.o changepoints_pkg.o xsim_arc_univar_file.o"
                ),
            },
        ],
    },
    "xaarc": {
        "data_file": "xaarc_data.txt",
        "steps": [
            {
                "label": "xaarc_make_data.py",
                "cmd": ["python", "xaarc_make_data.py", "{data_file}"],
            },
            {
                "label": "xaarc_file.R",
                "cmd": ["Rscript", "xaarc_file.R", "{data_file}"],
            },
            {
                "label": "xsim_aarc_univar_file.exe",
                "cmd": ["xsim_aarc_univar_file.exe", "{data_file}"],
                "requires": ["xsim_aarc_univar_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_aarc_univar_file "
                    "kind.o util.o median.o arma.o compare_io.o changepoints_pkg.o xsim_aarc_univar_file.o"
                ),
            },
        ],
    },
    "xonline_univar": {
        "data_file": "xonline_univar_data.txt",
        "steps": [
            {
                "label": "xonline_univar_make_data.py",
                "cmd": ["python", "xonline_univar_make_data.py", "{data_file}"],
            },
            {
                "label": "xonline_univar_file.R",
                "cmd": ["Rscript", "xonline_univar_file.R", "{data_file}"],
            },
            {
                "label": "xsim_online_univar_file.exe",
                "cmd": ["xsim_online_univar_file.exe", "{data_file}"],
                "requires": ["xsim_online_univar_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_online_univar_file "
                    "kind.o util.o median.o arma.o compare_io.o changepoints_pkg.o xsim_online_univar_file.o"
                ),
            },
        ],
    },
    "xbs_cov": {
        "data_file": "xbs_cov_data.txt",
        "steps": [
            {
                "label": "xbs_cov_make_data.py",
                "cmd": ["python", "xbs_cov_make_data.py", "{data_file}"],
            },
            {
                "label": "xbs_cov_file.R",
                "cmd": ["Rscript", "xbs_cov_file.R", "{data_file}"],
            },
            {
                "label": "xsim_bs_cov_file.exe",
                "cmd": ["xsim_bs_cov_file.exe", "{data_file}"],
                "requires": ["xsim_bs_cov_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_bs_cov_file "
                    "kind.o util.o median.o arma.o compare_io.o changepoints_pkg.o xsim_bs_cov_file.o"
                ),
            },
        ],
    },
    "xwbsip_cov": {
        "data_file": "xwbsip_cov_data.txt",
        "steps": [
            {
                "label": "xwbsip_cov_make_data.py",
                "cmd": ["python", "xwbsip_cov_make_data.py", "{data_file}"],
            },
            {
                "label": "xwbsip_cov_file.R",
                "cmd": ["Rscript", "xwbsip_cov_file.R", "{data_file}"],
            },
            {
                "label": "xsim_wbsip_cov_file.exe",
                "cmd": ["xsim_wbsip_cov_file.exe", "{data_file}"],
                "requires": ["xsim_wbsip_cov_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_wbsip_cov_file "
                    "kind.o util.o median.o arma.o compare_io.o changepoints_pkg.o xsim_wbsip_cov_file.o"
                ),
            },
        ],
    },
    "xdp_regression": {
        "data_file": "xdp_regression_data.txt",
        "steps": [
            {
                "label": "xdp_regression_make_data.py",
                "cmd": ["python", "xdp_regression_make_data.py", "{data_file}"],
            },
            {
                "label": "xdp_regression_file.R",
                "cmd": ["Rscript", "xdp_regression_file.R", "{data_file}"],
            },
            {
                "label": "xsim_dp_regression_file.exe",
                "cmd": ["xsim_dp_regression_file.exe", "{data_file}"],
                "requires": ["xsim_dp_regression_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_dp_regression_file "
                    "kind.o util.o median.o arma.o compare_io.o changepoints_pkg.o xsim_dp_regression_file.o"
                ),
            },
        ],
    },
    "xdpdu_regression": {
        "data_file": "xdpdu_regression_data.txt",
        "steps": [
            {
                "label": "xdpdu_regression_make_data.py",
                "cmd": ["python", "xdpdu_regression_make_data.py", "{data_file}"],
            },
            {
                "label": "xdpdu_regression_file.R",
                "cmd": ["Rscript", "xdpdu_regression_file.R", "{data_file}"],
            },
            {
                "label": "xsim_dpdu_regression_file.exe",
                "cmd": ["xsim_dpdu_regression_file.exe", "{data_file}"],
                "requires": ["xsim_dpdu_regression_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_dpdu_regression_file "
                    "kind.o util.o median.o arma.o compare_io.o changepoints_pkg.o xsim_dpdu_regression_file.o"
                ),
            },
        ],
    },
    "xdpdu2_regression": {
        "data_file": "xdpdu2_regression_data.txt",
        "steps": [
            {
                "label": "xdpdu2_regression_make_data.py",
                "cmd": ["python", "xdpdu2_regression_make_data.py", "{data_file}"],
            },
            {
                "label": "xdpdu2_regression_file.R",
                "cmd": ["Rscript", "xdpdu2_regression_file.R", "{data_file}"],
            },
            {
                "label": "xsim_dpdu2_regression_file.exe",
                "cmd": ["xsim_dpdu2_regression_file.exe", "{data_file}"],
                "requires": ["xsim_dpdu2_regression_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_dpdu2_regression_file "
                    "kind.o util.o median.o arma.o compare_io.o changepoints_pkg_dpdu.o xsim_dpdu2_regression_file.o"
                ),
            },
        ],
    },
    "xdp_poly": {
        "data_file": "xdp_poly_data.txt",
        "steps": [
            {
                "label": "xdp_poly_make_data.py",
                "cmd": ["python", "xdp_poly_make_data.py", "{data_file}"],
            },
            {
                "label": "xdp_poly_file.R",
                "cmd": ["Rscript", "xdp_poly_file.R", "{data_file}"],
            },
            {
                "label": "xsim_dp_poly_file.exe",
                "cmd": ["xsim_dp_poly_file.exe", "{data_file}"],
                "requires": ["xsim_dp_poly_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_dp_poly_file "
                    "kind.o util.o median.o arma.o compare_io.o changepoints_pkg.o xsim_dp_poly_file.o"
                ),
            },
        ],
    },
    "xcv_dp_poly": {
        "data_file": "xcv_dp_poly_data.txt",
        "steps": [
            {
                "label": "xcv_dp_poly_make_data.py",
                "cmd": ["python", "xcv_dp_poly_make_data.py", "{data_file}"],
            },
            {
                "label": "xcv_dp_poly_file.R",
                "cmd": ["Rscript", "xcv_dp_poly_file.R", "{data_file}"],
            },
            {
                "label": "xsim_cv_dp_poly_file.exe",
                "cmd": ["xsim_cv_dp_poly_file.exe", "{data_file}"],
                "requires": ["xsim_cv_dp_poly_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_cv_dp_poly_file "
                    "kind.o util.o median.o arma.o compare_io.o changepoints_pkg.o xsim_cv_dp_poly_file.o"
                ),
            },
        ],
    },
    "xlocal_refine_poly": {
        "data_file": "xcv_dp_poly_data.txt",
        "steps": [
            {
                "label": "xcv_dp_poly_make_data.py",
                "cmd": ["python", "xcv_dp_poly_make_data.py", "{data_file}"],
            },
            {
                "label": "xlocal_refine_poly_file.R",
                "cmd": ["Rscript", "xlocal_refine_poly_file.R", "{data_file}"],
            },
            {
                "label": "xsim_local_refine_poly_file.exe",
                "cmd": ["xsim_local_refine_poly_file.exe", "{data_file}"],
                "requires": ["xsim_local_refine_poly_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_local_refine_poly_file "
                    "kind.o util.o median.o arma.o compare_io.o changepoints_pkg.o xsim_local_refine_poly_file.o"
                ),
            },
        ],
    },
    "xcv_dpdu_regression": {
        "data_file": "xcv_dpdu_regression_data.txt",
        "steps": [
            {
                "label": "xcv_dpdu_regression_make_data.py",
                "cmd": ["python", "xcv_dpdu_regression_make_data.py", "{data_file}"],
            },
            {
                "label": "xcv_dpdu_regression_file.R",
                "cmd": ["Rscript", "xcv_dpdu_regression_file.R", "{data_file}"],
            },
            {
                "label": "xsim_cv_dpdu_regression_file.exe",
                "cmd": ["xsim_cv_dpdu_regression_file.exe", "{data_file}"],
                "requires": ["xsim_cv_dpdu_regression_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_cv_dpdu_regression_file "
                    "kind.o util.o median.o arma.o compare_io.o changepoints_pkg_dpdu.o xsim_cv_dpdu_regression_file.o"
                ),
            },
        ],
    },
    "xlocal_refine_dpdu_regression": {
        "data_file": "xdpdu_regression_data.txt",
        "steps": [
            {
                "label": "xdpdu_regression_make_data.py",
                "cmd": ["python", "xdpdu_regression_make_data.py", "{data_file}"],
            },
            {
                "label": "xlocal_refine_dpdu_regression_file.R",
                "cmd": ["Rscript", "xlocal_refine_dpdu_regression_file.R", "{data_file}"],
            },
            {
                "label": "xsim_local_refine_dpdu_regression_file.exe",
                "cmd": ["xsim_local_refine_dpdu_regression_file.exe", "{data_file}"],
                "requires": ["xsim_local_refine_dpdu_regression_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_local_refine_dpdu_regression_file "
                    "kind.o util.o median.o arma.o compare_io.o changepoints_pkg_dpdu.o xsim_local_refine_dpdu_regression_file.o"
                ),
            },
        ],
    },
    "xdp_var1": {
        "data_file": "xdp_var1_data.txt",
        "steps": [
            {
                "label": "xdp_var1_make_data.py",
                "cmd": ["python", "xdp_var1_make_data.py", "{data_file}"],
            },
            {
                "label": "xdp_var1_file.R",
                "cmd": ["Rscript", "xdp_var1_file.R", "{data_file}"],
            },
            {
                "label": "xsim_dp_var1_file.exe",
                "cmd": ["xsim_dp_var1_file.exe", "{data_file}"],
                "requires": ["xsim_dp_var1_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_dp_var1_file "
                    "kind.o util.o median.o arma.o compare_io.o changepoints_pkg_var1.o xsim_dp_var1_file.o"
                ),
            },
        ],
    },
    "xcv_dp_var1": {
        "data_file": "xcv_dp_var1_data.txt",
        "steps": [
            {
                "label": "xcv_dp_var1_make_data.py",
                "cmd": ["python", "xcv_dp_var1_make_data.py", "{data_file}"],
            },
            {
                "label": "xcv_dp_var1_file.R",
                "cmd": ["Rscript", "xcv_dp_var1_file.R", "{data_file}"],
            },
            {
                "label": "xsim_cv_dp_var1_file.exe",
                "cmd": ["xsim_cv_dp_var1_file.exe", "{data_file}"],
                "requires": ["xsim_cv_dp_var1_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_cv_dp_var1_file "
                    "kind.o util.o median.o arma.o compare_io.o changepoints_pkg_var1.o xsim_cv_dp_var1_file.o"
                ),
            },
        ],
    },
    "xlocal_refine_cv_var1": {
        "data_file": "xcv_dp_var1_data.txt",
        "steps": [
            {
                "label": "xcv_dp_var1_make_data.py",
                "cmd": ["python", "xcv_dp_var1_make_data.py", "{data_file}"],
            },
            {
                "label": "xlocal_refine_cv_var1_file.R",
                "cmd": ["Rscript", "xlocal_refine_cv_var1_file.R", "{data_file}"],
            },
            {
                "label": "xsim_local_refine_cv_var1_file.exe",
                "cmd": ["xsim_local_refine_cv_var1_file.exe", "{data_file}"],
                "requires": ["xsim_local_refine_cv_var1_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_local_refine_cv_var1_file "
                    "kind.o util.o median.o arma.o compare_io.o changepoints_pkg_var1.o xsim_local_refine_cv_var1_file.o"
                ),
            },
        ],
    },
    "xci_regression": {
        "data_file": "xdpdu_regression_data.txt",
        "steps": [
            {
                "label": "xdpdu_regression_make_data.py",
                "cmd": ["python", "xdpdu_regression_make_data.py", "{data_file}"],
            },
            {
                "label": "xci_regression_file.R",
                "cmd": ["Rscript", "xci_regression_file.R", "{data_file}"],
            },
            {
                "label": "xsim_ci_regression_file.exe",
                "cmd": ["xsim_ci_regression_file.exe", "{data_file}"],
                "requires": ["xsim_ci_regression_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_ci_regression_file "
                    "kind.o util.o median.o arma.o compare_io.o changepoints_pkg_dpdu.o xsim_ci_regression_file.o"
                ),
            },
        ],
    },
    "xstrucchange_breakpoints": {
        "data_file": "xstrucchange_breakpoints_data.txt",
        "steps": [
            {
                "label": "xstrucchange_breakpoints_make_data.py",
                "cmd": ["python", "xstrucchange_breakpoints_make_data.py", "{data_file}"],
            },
            {
                "label": "xstrucchange_breakpoints_file.R",
                "cmd": ["Rscript", "xstrucchange_breakpoints_file.R", "{data_file}"],
            },
            {
                "label": "xsim_strucchange_breakpoints_file.exe",
                "cmd": ["xsim_strucchange_breakpoints_file.exe", "{data_file}"],
                "requires": ["xsim_strucchange_breakpoints_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_strucchange_breakpoints_file"
                ),
            },
        ],
    },
    "xstrucchange_fstats": {
        "data_file": "xstrucchange_fstats_data.txt",
        "steps": [
            {
                "label": "xstrucchange_fstats_make_data.py",
                "cmd": ["python", "xstrucchange_fstats_make_data.py", "{data_file}"],
            },
            {
                "label": "xstrucchange_fstats_file.R",
                "cmd": ["Rscript", "xstrucchange_fstats_file.R", "{data_file}"],
            },
            {
                "label": "xsim_strucchange_fstats_file.exe",
                "cmd": ["xsim_strucchange_fstats_file.exe", "{data_file}"],
                "requires": ["xsim_strucchange_fstats_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_strucchange_fstats_file"
                ),
            },
        ],
    },
    "xstrucchange_efp": {
        "data_file": "xstrucchange_efp_data.txt",
        "steps": [
            {
                "label": "xstrucchange_efp_make_data.py",
                "cmd": ["python", "xstrucchange_efp_make_data.py", "{data_file}"],
            },
            {
                "label": "xstrucchange_efp_file.R",
                "cmd": ["Rscript", "xstrucchange_efp_file.R", "{data_file}"],
            },
            {
                "label": "xsim_strucchange_efp_file.exe",
                "cmd": ["xsim_strucchange_efp_file.exe", "{data_file}"],
                "requires": ["xsim_strucchange_efp_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_strucchange_efp_file"
                ),
            },
        ],
    },
    "xstrucchange_rec_cusum": {
        "data_file": "xstrucchange_rec_cusum_data.txt",
        "steps": [
            {
                "label": "xstrucchange_rec_cusum_make_data.py",
                "cmd": ["python", "xstrucchange_rec_cusum_make_data.py", "{data_file}"],
            },
            {
                "label": "xstrucchange_rec_cusum_file.R",
                "cmd": ["Rscript", "xstrucchange_rec_cusum_file.R", "{data_file}"],
            },
            {
                "label": "xsim_strucchange_rec_cusum_file.exe",
                "cmd": ["xsim_strucchange_rec_cusum_file.exe", "{data_file}"],
                "requires": ["xsim_strucchange_rec_cusum_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_strucchange_rec_cusum_file"
                ),
            },
        ],
    },
    "xstrucchange_ols_mosum": {
        "data_file": "xstrucchange_ols_mosum_data.txt",
        "steps": [
            {
                "label": "xstrucchange_ols_mosum_make_data.py",
                "cmd": ["python", "xstrucchange_ols_mosum_make_data.py", "{data_file}"],
            },
            {
                "label": "xstrucchange_ols_mosum_file.R",
                "cmd": ["Rscript", "xstrucchange_ols_mosum_file.R", "{data_file}"],
            },
            {
                "label": "xsim_strucchange_ols_mosum_file.exe",
                "cmd": ["xsim_strucchange_ols_mosum_file.exe", "{data_file}"],
                "requires": ["xsim_strucchange_ols_mosum_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_strucchange_ols_mosum_file"
                ),
            },
        ],
    },
    "xstrucchange_rec_mosum": {
        "data_file": "xstrucchange_rec_mosum_data.txt",
        "steps": [
            {
                "label": "xstrucchange_rec_mosum_make_data.py",
                "cmd": ["python", "xstrucchange_rec_mosum_make_data.py", "{data_file}"],
            },
            {
                "label": "xstrucchange_rec_mosum_file.R",
                "cmd": ["Rscript", "xstrucchange_rec_mosum_file.R", "{data_file}"],
            },
            {
                "label": "xsim_strucchange_rec_mosum_file.exe",
                "cmd": ["xsim_strucchange_rec_mosum_file.exe", "{data_file}"],
                "requires": ["xsim_strucchange_rec_mosum_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_strucchange_rec_mosum_file"
                ),
            },
        ],
    },
    "xstrucchange_mefp": {
        "data_file": "xstrucchange_mefp_data.txt",
        "steps": [
            {
                "label": "xstrucchange_mefp_make_data.py",
                "cmd": ["python", "xstrucchange_mefp_make_data.py", "{data_file}"],
            },
            {
                "label": "xstrucchange_mefp_file.R",
                "cmd": ["Rscript", "xstrucchange_mefp_file.R", "{data_file}"],
            },
            {
                "label": "xsim_strucchange_mefp_file.exe",
                "cmd": ["xsim_strucchange_mefp_file.exe", "{data_file}"],
                "requires": ["xsim_strucchange_mefp_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_strucchange_mefp_file"
                ),
            },
        ],
    },
    "xstrucchange_mefp_ols_mosum": {
        "data_file": "xstrucchange_mefp_ols_mosum_data.txt",
        "steps": [
            {
                "label": "xstrucchange_mefp_ols_mosum_make_data.py",
                "cmd": ["python", "xstrucchange_mefp_ols_mosum_make_data.py", "{data_file}"],
            },
            {
                "label": "xstrucchange_mefp_ols_mosum_file.R",
                "cmd": ["Rscript", "xstrucchange_mefp_ols_mosum_file.R", "{data_file}"],
            },
            {
                "label": "xsim_strucchange_mefp_ols_mosum_file.exe",
                "cmd": ["xsim_strucchange_mefp_ols_mosum_file.exe", "{data_file}"],
                "requires": ["xsim_strucchange_mefp_ols_mosum_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_strucchange_mefp_ols_mosum_file"
                ),
            },
        ],
    },
    "xstrucchange_re": {
        "data_file": "xstrucchange_re_data.txt",
        "steps": [
            {
                "label": "xstrucchange_re_make_data.py",
                "cmd": ["python", "xstrucchange_re_make_data.py", "{data_file}"],
            },
            {
                "label": "xstrucchange_re_file.R",
                "cmd": ["Rscript", "xstrucchange_re_file.R", "{data_file}"],
            },
            {
                "label": "xsim_strucchange_re_file.exe",
                "cmd": ["xsim_strucchange_re_file.exe", "{data_file}"],
                "requires": ["xsim_strucchange_re_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_strucchange_re_file"
                ),
            },
        ],
    },
    "xstrucchange_mefp_re": {
        "data_file": "xstrucchange_mefp_re_data.txt",
        "steps": [
            {
                "label": "xstrucchange_mefp_re_make_data.py",
                "cmd": ["python", "xstrucchange_mefp_re_make_data.py", "{data_file}"],
            },
            {
                "label": "xstrucchange_mefp_re_file.R",
                "cmd": ["Rscript", "xstrucchange_mefp_re_file.R", "{data_file}"],
            },
            {
                "label": "xsim_strucchange_mefp_re_file.exe",
                "cmd": ["xsim_strucchange_mefp_re_file.exe", "{data_file}"],
                "requires": ["xsim_strucchange_mefp_re_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_strucchange_mefp_re_file"
                ),
            },
        ],
    },
    "xstrucchange_me": {
        "data_file": "xstrucchange_me_data.txt",
        "steps": [
            {
                "label": "xstrucchange_me_make_data.py",
                "cmd": ["python", "xstrucchange_me_make_data.py", "{data_file}"],
            },
            {
                "label": "xstrucchange_me_file.R",
                "cmd": ["Rscript", "xstrucchange_me_file.R", "{data_file}"],
            },
            {
                "label": "xsim_strucchange_me_file.exe",
                "cmd": ["xsim_strucchange_me_file.exe", "{data_file}"],
                "requires": ["xsim_strucchange_me_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_strucchange_me_file"
                ),
            },
        ],
    },
    "xstrucchange_mefp_me": {
        "data_file": "xstrucchange_mefp_me_data.txt",
        "steps": [
            {
                "label": "xstrucchange_mefp_me_make_data.py",
                "cmd": ["python", "xstrucchange_mefp_me_make_data.py", "{data_file}"],
            },
            {
                "label": "xstrucchange_mefp_me_file.R",
                "cmd": ["Rscript", "xstrucchange_mefp_me_file.R", "{data_file}"],
            },
            {
                "label": "xsim_strucchange_mefp_me_file.exe",
                "cmd": ["xsim_strucchange_mefp_me_file.exe", "{data_file}"],
                "requires": ["xsim_strucchange_mefp_me_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_strucchange_mefp_me_file"
                ),
            },
        ],
    },
    "xcpm_student": {
        "data_file": "xcpm_student_data.txt",
        "steps": [
            {
                "label": "xcpm_student_make_data.py",
                "cmd": ["python", "xcpm_student_make_data.py", "{data_file}"],
            },
            {
                "label": "xcpm_student_file.R",
                "cmd": ["Rscript", "xcpm_student_file.R", "{data_file}"],
            },
            {
                "label": "xsim_cpm_student_file.exe",
                "cmd": ["xsim_cpm_student_file.exe", "{data_file}"],
                "requires": ["xsim_cpm_student_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_cpm_student_file"
                ),
            },
        ],
    },
    "xcpm_bartlett": {
        "data_file": "xcpm_bartlett_data.txt",
        "steps": [
            {
                "label": "xcpm_bartlett_make_data.py",
                "cmd": ["python", "xcpm_bartlett_make_data.py", "{data_file}"],
            },
            {
                "label": "xcpm_bartlett_file.R",
                "cmd": ["Rscript", "xcpm_bartlett_file.R", "{data_file}"],
            },
            {
                "label": "xsim_cpm_bartlett_file.exe",
                "cmd": ["xsim_cpm_bartlett_file.exe", "{data_file}"],
                "requires": ["xsim_cpm_bartlett_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_cpm_bartlett_file"
                ),
            },
        ],
    },
    "xcpm_joint": {
        "data_file": "xcpm_joint_data.txt",
        "steps": [
            {
                "label": "xcpm_joint_make_data.py",
                "cmd": ["python", "xcpm_joint_make_data.py", "{data_file}"],
            },
            {
                "label": "xcpm_joint_file.R",
                "cmd": ["Rscript", "xcpm_joint_file.R", "{data_file}"],
            },
            {
                "label": "xsim_cpm_joint_file.exe",
                "cmd": ["xsim_cpm_joint_file.exe", "{data_file}"],
                "requires": ["xsim_cpm_joint_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_cpm_joint_file"
                ),
            },
        ],
    },
    "xcpm_joint_adjusted": {
        "data_file": "xcpm_joint_adjusted_data.txt",
        "steps": [
            {
                "label": "xcpm_joint_adjusted_make_data.py",
                "cmd": ["python", "xcpm_joint_adjusted_make_data.py", "{data_file}"],
            },
            {
                "label": "xcpm_joint_adjusted_file.R",
                "cmd": ["Rscript", "xcpm_joint_adjusted_file.R", "{data_file}"],
            },
            {
                "label": "xsim_cpm_joint_adjusted_file.exe",
                "cmd": ["xsim_cpm_joint_adjusted_file.exe", "{data_file}"],
                "requires": ["xsim_cpm_joint_adjusted_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_cpm_joint_adjusted_file"
                ),
            },
        ],
    },
    "xcpm_joint_hawkins": {
        "data_file": "xcpm_joint_hawkins_data.txt",
        "steps": [
            {
                "label": "xcpm_joint_hawkins_make_data.py",
                "cmd": ["python", "xcpm_joint_hawkins_make_data.py", "{data_file}"],
            },
            {
                "label": "xcpm_joint_hawkins_file.R",
                "cmd": ["Rscript", "xcpm_joint_hawkins_file.R", "{data_file}"],
            },
            {
                "label": "xsim_cpm_joint_hawkins_file.exe",
                "cmd": ["xsim_cpm_joint_hawkins_file.exe", "{data_file}"],
                "requires": ["xsim_cpm_joint_hawkins_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_cpm_joint_hawkins_file"
                ),
            },
        ],
    },
    "xcpm_exponential": {
        "data_file": "xcpm_exponential_data.txt",
        "steps": [
            {
                "label": "xcpm_exponential_make_data.py",
                "cmd": ["python", "xcpm_exponential_make_data.py", "{data_file}"],
            },
            {
                "label": "xcpm_exponential_file.R",
                "cmd": ["Rscript", "xcpm_exponential_file.R", "{data_file}"],
            },
            {
                "label": "xsim_cpm_exponential_file.exe",
                "cmd": ["xsim_cpm_exponential_file.exe", "{data_file}"],
                "requires": ["xsim_cpm_exponential_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_cpm_exponential_file"
                ),
            },
        ],
    },
    "xcpm_exponential_adjusted": {
        "data_file": "xcpm_exponential_adjusted_data.txt",
        "steps": [
            {
                "label": "xcpm_exponential_adjusted_make_data.py",
                "cmd": ["python", "xcpm_exponential_adjusted_make_data.py", "{data_file}"],
            },
            {
                "label": "xcpm_exponential_adjusted_file.R",
                "cmd": ["Rscript", "xcpm_exponential_adjusted_file.R", "{data_file}"],
            },
            {
                "label": "xsim_cpm_exponential_adjusted_file.exe",
                "cmd": ["xsim_cpm_exponential_adjusted_file.exe", "{data_file}"],
                "requires": ["xsim_cpm_exponential_adjusted_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_cpm_exponential_adjusted_file"
                ),
            },
        ],
    },
    "xcpm_poisson": {
        "data_file": "xcpm_poisson_data.txt",
        "steps": [
            {
                "label": "xcpm_poisson_make_data.py",
                "cmd": ["python", "xcpm_poisson_make_data.py", "{data_file}"],
            },
            {
                "label": "xcpm_poisson_file.R",
                "cmd": ["Rscript", "xcpm_poisson_file.R", "{data_file}"],
            },
            {
                "label": "xsim_cpm_poisson_file.exe",
                "cmd": ["xsim_cpm_poisson_file.exe", "{data_file}"],
                "requires": ["xsim_cpm_poisson_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_cpm_poisson_file"
                ),
            },
        ],
    },
    "xcpm_mw": {
        "data_file": "xcpm_mw_data.txt",
        "steps": [
            {
                "label": "xcpm_mw_make_data.py",
                "cmd": ["python", "xcpm_mw_make_data.py", "{data_file}"],
            },
            {
                "label": "xcpm_mw_file.R",
                "cmd": ["Rscript", "xcpm_mw_file.R", "{data_file}"],
            },
            {
                "label": "xsim_cpm_mw_file.exe",
                "cmd": ["xsim_cpm_mw_file.exe", "{data_file}"],
                "requires": ["xsim_cpm_mw_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_cpm_mw_file"
                ),
            },
        ],
    },
    "xcpm_mood": {
        "data_file": "xcpm_mood_data.txt",
        "steps": [
            {
                "label": "xcpm_mood_make_data.py",
                "cmd": ["python", "xcpm_mood_make_data.py", "{data_file}"],
            },
            {
                "label": "xcpm_mood_file.R",
                "cmd": ["Rscript", "xcpm_mood_file.R", "{data_file}"],
            },
            {
                "label": "xsim_cpm_mood_file.exe",
                "cmd": ["xsim_cpm_mood_file.exe", "{data_file}"],
                "requires": ["xsim_cpm_mood_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_cpm_mood_file"
                ),
            },
        ],
    },
    "xcpm_lepage": {
        "data_file": "xcpm_lepage_data.txt",
        "steps": [
            {
                "label": "xcpm_lepage_make_data.py",
                "cmd": ["python", "xcpm_lepage_make_data.py", "{data_file}"],
            },
            {
                "label": "xcpm_lepage_file.R",
                "cmd": ["Rscript", "xcpm_lepage_file.R", "{data_file}"],
            },
            {
                "label": "xsim_cpm_lepage_file.exe",
                "cmd": ["xsim_cpm_lepage_file.exe", "{data_file}"],
                "requires": ["xsim_cpm_lepage_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_cpm_lepage_file"
                ),
            },
        ],
    },
    "xcpm_fet": {
        "data_file": "xcpm_fet_data.txt",
        "steps": [
            {
                "label": "xcpm_fet_make_data.py",
                "cmd": ["python", "xcpm_fet_make_data.py", "{data_file}"],
            },
            {
                "label": "xcpm_fet_file.R",
                "cmd": ["Rscript", "xcpm_fet_file.R", "{data_file}"],
            },
            {
                "label": "xsim_cpm_fet_file.exe",
                "cmd": ["xsim_cpm_fet_file.exe", "{data_file}"],
                "requires": ["xsim_cpm_fet_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_cpm_fet_file"
                ),
            },
        ],
    },
    "xcpm_ks": {
        "data_file": "xcpm_ks_data.txt",
        "steps": [
            {
                "label": "xcpm_ks_make_data.py",
                "cmd": ["python", "xcpm_ks_make_data.py", "{data_file}"],
            },
            {
                "label": "xcpm_ks_file.R",
                "cmd": ["Rscript", "xcpm_ks_file.R", "{data_file}"],
            },
            {
                "label": "xsim_cpm_ks_file.exe",
                "cmd": ["xsim_cpm_ks_file.exe", "{data_file}"],
                "requires": ["xsim_cpm_ks_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_cpm_ks_file"
                ),
            },
        ],
    },
    "xcpm_cvm": {
        "data_file": "xcpm_cvm_data.txt",
        "steps": [
            {
                "label": "xcpm_cvm_make_data.py",
                "cmd": ["python", "xcpm_cvm_make_data.py", "{data_file}"],
            },
            {
                "label": "xcpm_cvm_file.R",
                "cmd": ["Rscript", "xcpm_cvm_file.R", "{data_file}"],
            },
            {
                "label": "xsim_cpm_cvm_file.exe",
                "cmd": ["xsim_cpm_cvm_file.exe", "{data_file}"],
                "requires": ["xsim_cpm_cvm_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_cpm_cvm_file"
                ),
            },
        ],
    },
    "xecp_edivisive": {
        "data_file": "xecp_edivisive_data.txt",
        "steps": [
            {
                "label": "xecp_edivisive_make_data.py",
                "cmd": ["python", "xecp_edivisive_make_data.py", "{data_file}"],
            },
            {
                "label": "xecp_edivisive_file.R",
                "cmd": ["Rscript", "xecp_edivisive_file.R", "{data_file}"],
            },
            {
                "label": "xsim_ecp_edivisive_file.exe",
                "cmd": ["xsim_ecp_edivisive_file.exe", "{data_file}"],
                "requires": ["xsim_ecp_edivisive_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_ecp_edivisive_file"
                ),
            },
        ],
    },
    "xecp_edivisive_full": {
        "data_file": "xecp_edivisive_full_data.txt",
        "steps": [
            {
                "label": "xecp_edivisive_full_make_data.py",
                "cmd": ["python", "xecp_edivisive_full_make_data.py", "{data_file}"],
            },
            {
                "label": "xecp_edivisive_full_file.R",
                "cmd": ["Rscript", "xecp_edivisive_full_file.R", "{data_file}"],
            },
            {
                "label": "xsim_ecp_edivisive_full_file.exe",
                "cmd": ["xsim_ecp_edivisive_full_file.exe", "{data_file}"],
                "requires": ["xsim_ecp_edivisive_full_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_ecp_edivisive_full_file"
                ),
            },
        ],
    },
    "xecp_eagglo": {
        "data_file": "xecp_eagglo_data.txt",
        "steps": [
            {
                "label": "xecp_eagglo_make_data.py",
                "cmd": ["python", "xecp_eagglo_make_data.py", "{data_file}"],
            },
            {
                "label": "xecp_eagglo_file.R",
                "cmd": ["Rscript", "xecp_eagglo_file.R", "{data_file}"],
            },
            {
                "label": "xsim_ecp_eagglo_file.exe",
                "cmd": ["xsim_ecp_eagglo_file.exe", "{data_file}"],
                "requires": ["xsim_ecp_eagglo_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_ecp_eagglo_file"
                ),
            },
        ],
    },
    "xecp_cp3o": {
        "data_file": "xecp_cp3o_data.txt",
        "steps": [
            {
                "label": "xecp_cp3o_make_data.py",
                "cmd": ["python", "xecp_cp3o_make_data.py", "{data_file}"],
            },
            {
                "label": "xecp_cp3o_file.R",
                "cmd": ["Rscript", "xecp_cp3o_file.R", "{data_file}"],
            },
            {
                "label": "xsim_ecp_cp3o_file.exe",
                "cmd": ["xsim_ecp_cp3o_file.exe", "{data_file}"],
                "requires": ["xsim_ecp_cp3o_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_ecp_cp3o_file"
                ),
            },
        ],
    },
    "xecp_cp3o_delta": {
        "data_file": "xecp_cp3o_delta_data.txt",
        "steps": [
            {
                "label": "xecp_cp3o_delta_make_data.py",
                "cmd": ["python", "xecp_cp3o_delta_make_data.py", "{data_file}"],
            },
            {
                "label": "xecp_cp3o_delta_file.R",
                "cmd": ["Rscript", "xecp_cp3o_delta_file.R", "{data_file}"],
            },
            {
                "label": "xsim_ecp_cp3o_delta_file.exe",
                "cmd": ["xsim_ecp_cp3o_delta_file.exe", "{data_file}"],
                "requires": ["xsim_ecp_cp3o_delta_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_ecp_cp3o_delta_file"
                ),
            },
        ],
    },
    "xcv_dp_regression": {
        "data_file": "xcv_dp_regression_data.txt",
        "steps": [
            {
                "label": "xcv_dp_regression_make_data.py",
                "cmd": ["python", "xcv_dp_regression_make_data.py", "{data_file}"],
            },
            {
                "label": "xcv_dp_regression_file.R",
                "cmd": ["Rscript", "xcv_dp_regression_file.R", "{data_file}"],
            },
            {
                "label": "xsim_cv_dp_regression_file.exe",
                "cmd": ["xsim_cv_dp_regression_file.exe", "{data_file}"],
                "requires": ["xsim_cv_dp_regression_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_cv_dp_regression_file "
                    "kind.o util.o median.o arma.o compare_io.o changepoints_pkg.o xsim_cv_dp_regression_file.o"
                ),
            },
        ],
    },
    "xlocal_refine_regression": {
        "data_file": "xdp_regression_data.txt",
        "steps": [
            {
                "label": "xdp_regression_make_data.py",
                "cmd": ["python", "xdp_regression_make_data.py", "{data_file}"],
            },
            {
                "label": "xlocal_refine_regression_file.R",
                "cmd": ["Rscript", "xlocal_refine_regression_file.R", "{data_file}"],
            },
            {
                "label": "xsim_local_refine_regression_file.exe",
                "cmd": ["xsim_local_refine_regression_file.exe", "{data_file}"],
                "requires": ["xsim_local_refine_regression_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_local_refine_regression_file "
                    "kind.o util.o median.o arma.o compare_io.o changepoints_pkg.o xsim_local_refine_regression_file.o"
                ),
            },
        ],
    },
    "xlocal_refine_univar": {
        "data_file": "xwbs_data.txt",
        "steps": [
            {
                "label": "xwbs_make_data.py",
                "cmd": ["python", "xwbs_make_data.py", "{data_file}"],
            },
            {
                "label": "xlocal_refine_univar_file.R",
                "cmd": ["Rscript", "xlocal_refine_univar_file.R", "{data_file}"],
            },
            {
                "label": "xsim_local_refine_univar_file.exe",
                "cmd": ["xsim_local_refine_univar_file.exe", "{data_file}"],
                "requires": ["xsim_local_refine_univar_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_local_refine_univar_file "
                    "kind.o util.o median.o arma.o compare_io.o changepoints_pkg.o xsim_local_refine_univar_file.o"
                ),
            },
        ],
    },
    "xdp_univar": {
        "data_file": "xdp_univar_data.txt",
        "steps": [
            {
                "label": "xdp_univar_make_data.py",
                "cmd": ["python", "xdp_univar_make_data.py", "{data_file}"],
            },
            {
                "label": "xdp_univar_file.R",
                "cmd": ["Rscript", "xdp_univar_file.R", "{data_file}"],
            },
            {
                "label": "xsim_dp_univar_file.exe",
                "cmd": ["xsim_dp_univar_file.exe", "{data_file}"],
                "requires": ["xsim_dp_univar_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_dp_univar_file "
                    "kind.o util.o median.o arma.o compare_io.o changepoints_pkg.o xsim_dp_univar_file.o"
                ),
            },
        ],
    },
    "xcv_dp_univar": {
        "data_file": "xcv_dp_univar_data.txt",
        "steps": [
            {
                "label": "xcv_dp_univar_make_data.py",
                "cmd": ["python", "xcv_dp_univar_make_data.py", "{data_file}"],
            },
            {
                "label": "xcv_dp_univar_file.R",
                "cmd": ["Rscript", "xcv_dp_univar_file.R", "{data_file}"],
            },
            {
                "label": "xsim_cv_dp_univar_file.exe",
                "cmd": ["xsim_cv_dp_univar_file.exe", "{data_file}"],
                "requires": ["xsim_cv_dp_univar_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_cv_dp_univar_file "
                    "kind.o util.o median.o arma.o compare_io.o changepoints_pkg.o xsim_cv_dp_univar_file.o"
                ),
            },
        ],
    },
    "xtext_like": {
        "data_file": "xtext_like_data.txt",
        "steps": [
            {
                "label": "xtext_like_make_data.py",
                "cmd": ["python", "xtext_like_make_data.py", "{data_file}"],
            },
            {
                "label": "xtext_like_file.py",
                "cmd": ["python", "xtext_like_file.py", "{data_file}"],
            },
            {
                "label": "xtext_like_file.R",
                "cmd": ["Rscript", "xtext_like_file.R", "{data_file}"],
            },
            {
                "label": "xsim_text_like_pelt_file.exe",
                "cmd": ["xsim_text_like_pelt_file.exe", "{data_file}"],
                "requires": ["xsim_text_like_pelt_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_text_like_pelt_file "
                    "kind.o util.o changepoint.o xsim_text_like_pelt_file.o"
                ),
            },
        ],
    },
    "xvariance_regime": {
        "data_file": "xvariance_regime_data.txt",
        "steps": [
            {
                "label": "xvariance_regime_make_data.py",
                "cmd": ["python", "xvariance_regime_make_data.py", "{data_file}"],
            },
            {
                "label": "xvariance_regime_file.py",
                "cmd": ["python", "xvariance_regime_file.py", "{data_file}"],
            },
            {
                "label": "xvariance_regime_file.R",
                "cmd": ["Rscript", "xvariance_regime_file.R", "{data_file}"],
            },
            {
                "label": "xsim_variance_regime_pelt_file.exe",
                "cmd": ["xsim_variance_regime_pelt_file.exe", "{data_file}"],
                "requires": ["xsim_variance_regime_pelt_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_variance_regime_pelt_file "
                    "kind.o util.o changepoint.o xsim_variance_regime_pelt_file.o"
                ),
            },
        ],
    },
    "xkernel": {
        "data_file": "xkernel_data.txt",
        "steps": [
            {
                "label": "xkernel_make_data.py",
                "cmd": ["python", "xkernel_make_data.py", "{data_file}"],
            },
            {
                "label": "xkernel_file.py",
                "cmd": ["python", "xkernel_file.py", "{data_file}"],
            },
            {
                "label": "xkernel_file_kerSeg.R",
                "cmd": ["Rscript", "xkernel_file_kerSeg.R", "{data_file}"],
            },
            {
                "label": "xsim_kernel_file.exe",
                "cmd": ["xsim_kernel_file.exe", "{data_file}"],
                "requires": ["xsim_kernel_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_kernel_file "
                    "kind.o util.o changepoint.o xsim_kernel_file.o"
                ),
            },
        ],
    },
    "xcosts_mean_variance": {
        "data_file": "xcosts_mean_variance_data.txt",
        "steps": [
            {
                "label": "xcosts_mean_variance_make_data.py",
                "cmd": ["python", "xcosts_mean_variance_make_data.py", "{data_file}"],
            },
            {
                "label": "xcosts_mean_variance_file.py",
                "cmd": ["python", "xcosts_mean_variance_file.py", "{data_file}"],
            },
            {
                "label": "xcosts_mean_variance_file.R",
                "cmd": ["Rscript", "xcosts_mean_variance_file.R", "{data_file}"],
            },
            {
                "label": "xcosts_mean_variance_file.exe",
                "cmd": ["xcosts_mean_variance_file.exe", "{data_file}"],
                "requires": ["xcosts_mean_variance_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xcosts_mean_variance_file "
                    "kind.o util.o median.o changepoint.o xcosts_mean_variance_file.o"
                ),
            },
        ],
    },
    "xmetrics": {
        "data_file": "xmetrics_data.txt",
        "steps": [
            {
                "label": "xmetrics_make_data.py",
                "cmd": ["python", "xmetrics_make_data.py", "{data_file}"],
            },
            {
                "label": "xmetrics_file.py",
                "cmd": ["python", "xmetrics_file.py", "{data_file}"],
            },
            {
                "label": "xmetrics_file.R",
                "cmd": ["Rscript", "xmetrics_file.R", "{data_file}"],
            },
            {
                "label": "xmetrics_file.exe",
                "cmd": ["xmetrics_file.exe", "{data_file}"],
                "requires": ["xmetrics_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xmetrics_file "
                    "kind.o util.o median.o changepoint.o changepoint_binseg.o xmetrics_file.o"
                ),
            },
        ],
    },
    "xbocpd": {
        "data_file": "xbocpd_data.txt",
        "steps": [
            {
                "label": "xbocpd_make_data.py",
                "cmd": ["python", "xbocpd_make_data.py", "{data_file}"],
            },
            {
                "label": "xbocpd_file.py",
                "cmd": ["python", "xbocpd_file.py", "{data_file}"],
            },
            {
                "label": "xsim_bocpd_file.exe",
                "cmd": ["xsim_bocpd_file.exe", "{data_file}"],
                "requires": ["xsim_bocpd_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_bocpd_file"
                ),
            },
        ],
    },
    "xbocpd_poisson": {
        "data_file": "xbocpd_poisson_data.txt",
        "steps": [
            {
                "label": "xbocpd_poisson_make_data.py",
                "cmd": ["python", "xbocpd_poisson_make_data.py", "{data_file}"],
            },
            {
                "label": "xbocpd_poisson_file.py",
                "cmd": ["python", "xbocpd_poisson_file.py", "{data_file}"],
            },
            {
                "label": "xsim_bocpd_poisson_file.exe",
                "cmd": ["xsim_bocpd_poisson_file.exe", "{data_file}"],
                "requires": ["xsim_bocpd_poisson_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_bocpd_poisson_file"
                ),
            },
        ],
    },
    "xbocpd_beta_bernoulli": {
        "data_file": "xbocpd_beta_bernoulli_data.txt",
        "steps": [
            {
                "label": "xbocpd_beta_bernoulli_make_data.py",
                "cmd": ["python", "xbocpd_beta_bernoulli_make_data.py", "{data_file}"],
            },
            {
                "label": "xbocpd_beta_bernoulli_file.py",
                "cmd": ["python", "xbocpd_beta_bernoulli_file.py", "{data_file}"],
            },
            {
                "label": "xsim_bocpd_beta_bernoulli_file.exe",
                "cmd": ["xsim_bocpd_beta_bernoulli_file.exe", "{data_file}"],
                "requires": ["xsim_bocpd_beta_bernoulli_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_bocpd_beta_bernoulli_file"
                ),
            },
        ],
    },
    "xbocd": {
        "data_file": "xbocd_data.txt",
        "steps": [
            {
                "label": "xbocd_make_data.py",
                "cmd": ["python", "xbocd_make_data.py", "{data_file}"],
            },
            {
                "label": "xbocd_file.py",
                "cmd": ["python", "xbocd_file.py", "{data_file}"],
            },
            {
                "label": "xsim_bocd_file.exe",
                "cmd": ["xsim_bocd_file.exe", "{data_file}"],
                "requires": ["xsim_bocd_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_bocd_file"
                ),
            },
        ],
    },
    "xbocpd_changepoynt": {
        "data_file": "xbocpd_changepoynt_data.txt",
        "steps": [
            {
                "label": "xbocpd_changepoynt_make_data.py",
                "cmd": ["python", "xbocpd_changepoynt_make_data.py", "{data_file}"],
            },
            {
                "label": "xbocpd_changepoynt_file.py",
                "cmd": ["python", "xbocpd_changepoynt_file.py", "{data_file}"],
            },
            {
                "label": "xsim_bocpd_changepoynt_file.exe",
                "cmd": ["xsim_bocpd_changepoynt_file.exe", "{data_file}"],
                "requires": ["xsim_bocpd_changepoynt_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_bocpd_changepoynt_file"
                ),
            },
        ],
    },
    "xsst": {
        "data_file": "xsst_data.txt",
        "steps": [
            {
                "label": "xsst_make_data.py",
                "cmd": ["python", "xsst_make_data.py", "{data_file}"],
            },
            {
                "label": "xsst_file.py",
                "cmd": ["python", "xsst_file.py", "{data_file}"],
            },
            {
                "label": "xsim_sst_file.exe",
                "cmd": ["xsim_sst_file.exe", "{data_file}"],
                "requires": ["xsim_sst_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_sst_file"
                ),
            },
        ],
    },
    "xesst": {
        "data_file": "xesst_data.txt",
        "steps": [
            {
                "label": "xesst_make_data.py",
                "cmd": ["python", "xesst_make_data.py", "{data_file}"],
            },
            {
                "label": "xesst_file.py",
                "cmd": ["python", "xesst_file.py", "{data_file}"],
            },
            {
                "label": "xsim_esst_file.exe",
                "cmd": ["xsim_esst_file.exe", "{data_file}"],
                "requires": ["xsim_esst_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_esst_file"
                ),
            },
        ],
    },
    "xrulsif": {
        "data_file": "xrulsif_data.txt",
        "steps": [
            {
                "label": "xrulsif_make_data.py",
                "cmd": ["python", "xrulsif_make_data.py", "{data_file}"],
            },
            {
                "label": "xrulsif_file.py",
                "cmd": ["python", "xrulsif_file.py", "{data_file}"],
            },
            {
                "label": "xsim_rulsif_file.exe",
                "cmd": ["xsim_rulsif_file.exe", "{data_file}"],
                "requires": ["xsim_rulsif_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_rulsif_file"
                ),
            },
        ],
    },
    "xulsif": {
        "data_file": "xulsif_data.txt",
        "steps": [
            {
                "label": "xulsif_make_data.py",
                "cmd": ["python", "xulsif_make_data.py", "{data_file}"],
            },
            {
                "label": "xulsif_file.py",
                "cmd": ["python", "xulsif_file.py", "{data_file}"],
            },
            {
                "label": "xsim_ulsif_file.exe",
                "cmd": ["xsim_ulsif_file.exe", "{data_file}"],
                "requires": ["xsim_ulsif_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_ulsif_file"
                ),
            },
        ],
    },
    "xargpcpd": {
        "data_file": "xargpcpd_data.txt",
        "steps": [
            {
                "label": "xargpcpd_make_data.py",
                "cmd": ["python", "xargpcpd_make_data.py", "{data_file}"],
            },
            {
                "label": "xargpcpd_file.py",
                "cmd": ["python", "xargpcpd_file.py", "{data_file}"],
            },
            {
                "label": "xsim_argpcpd_gp_file.exe",
                "cmd": ["xsim_argpcpd_gp_file.exe", "{data_file}"],
                "requires": ["xsim_argpcpd_gp_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_argpcpd_gp_file"
                ),
            },
        ],
    },
    "xargpcpd_rust": {
        "data_file": "xargpcpd_data.txt",
        "steps": [
            {
                "label": "xargpcpd_make_data.py",
                "cmd": ["python", "xargpcpd_make_data.py", "{data_file}"],
            },
            {
                "label": "xargpcpd_file.py",
                "cmd": ["python", "xargpcpd_file.py", "{data_file}"],
            },
            {
                "label": "xsim_argpcpd_gp_file.exe",
                "cmd": ["xsim_argpcpd_gp_file.exe", "{data_file}", "rust"],
                "requires": ["xsim_argpcpd_gp_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_argpcpd_gp_file"
                ),
            },
        ],
    },
    "xseeded_binseg": {
        "data_file": "xseeded_binseg_data.txt",
        "steps": [
            {
                "label": "xseeded_binseg_make_data.py",
                "cmd": ["python", "xseeded_binseg_make_data.py", "{data_file}"],
            },
            {
                "label": "xseeded_binseg_file.py",
                "cmd": ["python", "xseeded_binseg_file.py", "{data_file}"],
            },
            {
                "label": "xsim_seeded_binseg_file.exe",
                "cmd": ["xsim_seeded_binseg_file.exe", "{data_file}"],
                "requires": ["xsim_seeded_binseg_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_seeded_binseg_file"
                ),
            },
        ],
    },
    "xcrops": {
        "data_file": "xcrops_data.txt",
        "steps": [
            {
                "label": "xcrops_make_data.py",
                "cmd": ["python", "xcrops_make_data.py", "{data_file}"],
            },
            {
                "label": "xcrops_file.py",
                "cmd": ["python", "xcrops_file.py", "{data_file}"],
            },
            {
                "label": "xsim_crops_file.exe",
                "cmd": ["xsim_crops_file.exe", "{data_file}"],
                "requires": ["xsim_crops_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_crops_file"
                ),
            },
        ],
    },
    "xcapa": {
        "data_file": "xcapa_data.txt",
        "steps": [
            {
                "label": "xcapa_make_data.py",
                "cmd": ["python", "xcapa_make_data.py", "{data_file}"],
            },
            {
                "label": "xcapa_file.py",
                "cmd": ["python", "xcapa_file.py", "{data_file}"],
            },
            {
                "label": "xsim_capa_file.exe",
                "cmd": ["xsim_capa_file.exe", "{data_file}"],
                "requires": ["xsim_capa_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_capa_file"
                ),
            },
        ],
    },
    "xclass": {
        "data_file": "xclass_data.txt",
        "steps": [
            {
                "label": "xclass_make_data.py",
                "cmd": ["python", "xclass_make_data.py", "{data_file}"],
            },
            {
                "label": "xclass_file.py",
                "cmd": ["python", "xclass_file.py", "{data_file}"],
            },
            {
                "label": "xsim_class_file.exe",
                "cmd": ["xsim_class_file.exe", "{data_file}"],
                "requires": ["xsim_class_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_class_file"
                ),
            },
        ],
    },
    "xclap": {
        "data_file": "xclap_data.txt",
        "steps": [
            {
                "label": "xclap_make_data.py",
                "cmd": ["python", "xclap_make_data.py", "{data_file}"],
            },
            {
                "label": "xclap_file.py",
                "cmd": ["python", "xclap_file.py", "{data_file}"],
            },
            {
                "label": "xsim_clap_file.exe",
                "cmd": ["xsim_clap_file.exe", "{data_file}"],
                "requires": ["xsim_clap_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_clap_file"
                ),
            },
        ],
    },
    "xclasp_ensemble": {
        "data_file": "xclasp_ensemble_data.txt",
        "steps": [
            {
                "label": "xclasp_ensemble_make_data.py",
                "cmd": ["python", "xclasp_ensemble_make_data.py", "{data_file}"],
            },
            {
                "label": "xclasp_ensemble_file.py",
                "cmd": ["python", "xclasp_ensemble_file.py", "{data_file}"],
            },
            {
                "label": "xsim_clasp_ensemble_file.exe",
                "cmd": ["xsim_clasp_ensemble_file.exe", "{data_file}"],
                "requires": ["xsim_clasp_ensemble_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_clasp_ensemble_file"
                ),
            },
        ],
    },
    "xagglomerative_clap": {
        "data_file": "xagglomerative_clap_data.txt",
        "steps": [
            {
                "label": "xagglomerative_clap_make_data.py",
                "cmd": ["python", "xagglomerative_clap_make_data.py", "{data_file}"],
            },
            {
                "label": "xagglomerative_clap_file.py",
                "cmd": ["python", "xagglomerative_clap_file.py", "{data_file}"],
            },
            {
                "label": "xsim_agglomerative_clap_file.exe",
                "cmd": ["xsim_agglomerative_clap_file.exe", "{data_file}"],
                "requires": ["xsim_agglomerative_clap_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_agglomerative_clap_file"
                ),
            },
        ],
    },
    "xroerich_windows": {
        "data_file": "xroerich_windows_data.txt",
        "steps": [
            {
                "label": "xroerich_windows_make_data.py",
                "cmd": ["python", "xroerich_windows_make_data.py", "{data_file}"],
            },
            {
                "label": "xroerich_windows_file.py",
                "cmd": ["python", "xroerich_windows_file.py", "{data_file}"],
            },
            {
                "label": "xsim_roerich_windows_file.exe",
                "cmd": ["xsim_roerich_windows_file.exe", "{data_file}"],
                "requires": ["xsim_roerich_windows_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_roerich_windows_file"
                ),
            },
        ],
    },
    "xroerich_energy": {
        "data_file": "xroerich_energy_data.txt",
        "steps": [
            {
                "label": "xroerich_energy_make_data.py",
                "cmd": ["python", "xroerich_energy_make_data.py", "{data_file}"],
            },
            {
                "label": "xroerich_energy_file.py",
                "cmd": ["python", "xroerich_energy_file.py", "{data_file}"],
            },
            {
                "label": "xsim_roerich_energy_file.exe",
                "cmd": ["xsim_roerich_energy_file.exe", "{data_file}"],
                "requires": ["xsim_roerich_energy_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_roerich_energy_file"
                ),
            },
        ],
    },
    "xroerich_classifier": {
        "data_file": "xroerich_classifier_data.txt",
        "steps": [
            {
                "label": "xroerich_classifier_make_data.py",
                "cmd": ["python", "xroerich_classifier_make_data.py", "{data_file}"],
            },
            {
                "label": "xroerich_classifier_file.py",
                "cmd": ["python", "xroerich_classifier_file.py", "{data_file}"],
            },
            {
                "label": "xsim_roerich_classifier_file.exe",
                "cmd": ["xsim_roerich_classifier_file.exe", "{data_file}"],
                "requires": ["xsim_roerich_classifier_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_roerich_classifier_file"
                ),
            },
        ],
    },
    "xroerich_classifier_cv": {
        "data_file": "xroerich_classifier_cv_data.txt",
        "steps": [
            {
                "label": "xroerich_classifier_cv_make_data.py",
                "cmd": ["python", "xroerich_classifier_cv_make_data.py", "{data_file}"],
            },
            {
                "label": "xroerich_classifier_cv_file.py",
                "cmd": ["python", "xroerich_classifier_cv_file.py", "{data_file}"],
            },
            {
                "label": "xsim_roerich_classifier_cv_file.exe",
                "cmd": ["xsim_roerich_classifier_cv_file.exe", "{data_file}"],
                "requires": ["xsim_roerich_classifier_cv_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_roerich_classifier_cv_file"
                ),
            },
        ],
    },
    "xroerich_rulsif": {
        "data_file": "xroerich_rulsif_data.txt",
        "steps": [
            {
                "label": "xroerich_rulsif_make_data.py",
                "cmd": ["python", "xroerich_rulsif_make_data.py", "{data_file}"],
            },
            {
                "label": "xroerich_rulsif_file.py",
                "cmd": ["python", "xroerich_rulsif_file.py", "{data_file}"],
            },
            {
                "label": "xsim_roerich_rulsif_file.exe",
                "cmd": ["xsim_roerich_rulsif_file.exe", "{data_file}"],
                "requires": ["xsim_roerich_rulsif_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_roerich_rulsif_file"
                ),
            },
        ],
    },
    "xbeast": {
        "data_file": "xbeast_data.txt",
        "steps": [
            {
                "label": "xbeast_make_data.py",
                "cmd": ["python", "xbeast_make_data.py", "{data_file}"],
            },
            {
                "label": "xbeast_file.py",
                "cmd": ["python", "xbeast_file.py", "{data_file}"],
            },
            {
                "label": "xsim_beast_file.exe",
                "cmd": ["xsim_beast_file.exe", "{data_file}"],
                "requires": ["xsim_beast_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_beast_file"
                ),
            },
        ],
    },
    "xbeast_seasonal": {
        "data_file": "xbeast_seasonal_data.txt",
        "steps": [
            {
                "label": "xbeast_seasonal_make_data.py",
                "cmd": ["python", "xbeast_seasonal_make_data.py", "{data_file}"],
            },
            {
                "label": "xbeast_seasonal_file.py",
                "cmd": ["python", "xbeast_seasonal_file.py", "{data_file}"],
            },
            {
                "label": "xsim_beast_seasonal_file.exe",
                "cmd": ["xsim_beast_seasonal_file.exe", "{data_file}"],
                "requires": ["xsim_beast_seasonal_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_beast_seasonal_file"
                ),
            },
        ],
    },
    "xbeast_irreg": {
        "data_file": "xbeast_irreg_data.txt",
        "steps": [
            {
                "label": "xbeast_irreg_make_data.py",
                "cmd": ["python", "xbeast_irreg_make_data.py", "{data_file}"],
            },
            {
                "label": "xbeast_irreg_file.py",
                "cmd": ["python", "xbeast_irreg_file.py", "{data_file}"],
            },
            {
                "label": "xsim_beast_irreg_file.exe",
                "cmd": ["xsim_beast_irreg_file.exe", "{data_file}"],
                "requires": ["xsim_beast_irreg_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_beast_irreg_file"
                ),
            },
        ],
    },
    "xbeast123": {
        "data_file": "xbeast123_data.txt",
        "steps": [
            {
                "label": "xbeast123_make_data.py",
                "cmd": ["python", "xbeast123_make_data.py", "{data_file}"],
            },
            {
                "label": "xbeast123_file.py",
                "cmd": ["python", "xbeast123_file.py", "{data_file}"],
            },
            {
                "label": "xsim_beast123_file.exe",
                "cmd": ["xsim_beast123_file.exe", "{data_file}"],
                "requires": ["xsim_beast123_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_beast123_file"
                ),
            },
        ],
    },
    "xstreaming_clasp": {
        "data_file": "xstreaming_clasp_data.txt",
        "steps": [
            {
                "label": "xstreaming_clasp_make_data.py",
                "cmd": ["python", "xstreaming_clasp_make_data.py", "{data_file}"],
            },
            {
                "label": "xstreaming_clasp_file.py",
                "cmd": ["python", "xstreaming_clasp_file.py", "{data_file}"],
            },
            {
                "label": "xsim_streaming_clasp_file.exe",
                "cmd": ["xsim_streaming_clasp_file.exe", "{data_file}"],
                "requires": ["xsim_streaming_clasp_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_streaming_clasp_file"
                ),
            },
        ],
    },
    "xclasp": {
        "data_file": "xclasp_data.txt",
        "steps": [
            {
                "label": "xclasp_make_data.py",
                "cmd": ["python", "xclasp_make_data.py", "{data_file}"],
            },
            {
                "label": "xclasp_file.py",
                "cmd": ["python", "xclasp_file.py", "{data_file}"],
            },
            {
                "label": "xsim_clasp_file.exe",
                "cmd": ["xsim_clasp_file.exe", "{data_file}"],
                "requires": ["xsim_clasp_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_clasp_file"
                ),
            },
        ],
    },
    "xclasp_changepoynt": {
        "data_file": "xclasp_changepoynt_data.txt",
        "steps": [
            {
                "label": "xclasp_changepoynt_make_data.py",
                "cmd": ["python", "xclasp_changepoynt_make_data.py", "{data_file}"],
            },
            {
                "label": "xclasp_changepoynt_file.py",
                "cmd": ["python", "xclasp_changepoynt_file.py", "{data_file}"],
            },
            {
                "label": "xsim_clasp_changepoynt_file.exe",
                "cmd": ["xsim_clasp_changepoynt_file.exe", "{data_file}"],
                "requires": ["xsim_clasp_changepoynt_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  gfortran -O3 -Wall -static -g -o xsim_clasp_changepoynt_file "
                    "kind.o util.o median.o arma.o compare_io.o claspy.o xsim_clasp_changepoynt_file.o"
                ),
            },
        ],
    },
    "xclasp_alt": {
        "data_file": "xclasp_alt_data.txt",
        "steps": [
            {
                "label": "xclasp_alt_make_data.py",
                "cmd": ["python", "xclasp_alt_make_data.py", "{data_file}"],
            },
            {
                "label": "xclasp_alt_file.py",
                "cmd": ["python", "xclasp_alt_file.py", "{data_file}"],
            },
            {
                "label": "xsim_clasp_alt_file.exe",
                "cmd": ["xsim_clasp_alt_file.exe", "{data_file}"],
                "requires": ["xsim_clasp_alt_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_clasp_alt_file"
                ),
            },
        ],
    },
    "xchangepointnp": {
        "data_file": "xchangepointnp_data.txt",
        "steps": [
            {
                "label": "xchangepointnp_make_data.py",
                "cmd": ["python", "xchangepointnp_make_data.py", "{data_file}"],
            },
            {
                "label": "xchangepointnp_file.R",
                "cmd": ["Rscript", "xchangepointnp_file.R", "{data_file}"],
            },
            {
                "label": "xsim_changepointnp_file.exe",
                "cmd": ["xsim_changepointnp_file.exe", "{data_file}"],
                "requires": ["xsim_changepointnp_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_changepointnp_file"
                ),
            },
        ],
    },
    "xbcp": {
        "data_file": "xbcp_data.txt",
        "steps": [
            {
                "label": "xbcp_make_data.py",
                "cmd": ["python", "xbcp_make_data.py", "{data_file}"],
            },
            {
                "label": "xbcp_file.R",
                "cmd": ["Rscript", "xbcp_file.R", "{data_file}"],
            },
            {
                "label": "xsim_bcp_file.exe",
                "cmd": ["xsim_bcp_file.exe", "{data_file}"],
                "requires": ["xsim_bcp_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_bcp_file"
                ),
            },
        ],
    },
    "xbcp_multivar": {
        "data_file": "xbcp_multivar_data.txt",
        "steps": [
            {
                "label": "xbcp_multivar_make_data.py",
                "cmd": ["python", "xbcp_multivar_make_data.py", "{data_file}"],
            },
            {
                "label": "xbcp_multivar_file.R",
                "cmd": ["Rscript", "xbcp_multivar_file.R", "{data_file}"],
            },
            {
                "label": "xsim_bcp_multivar_file.exe",
                "cmd": ["xsim_bcp_multivar_file.exe", "{data_file}"],
                "requires": ["xsim_bcp_multivar_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_bcp_multivar_file"
                ),
            },
        ],
    },
    "xbcp_reg": {
        "data_file": "xbcp_reg_data.txt",
        "steps": [
            {
                "label": "xbcp_reg_make_data.py",
                "cmd": ["python", "xbcp_reg_make_data.py", "{data_file}"],
            },
            {
                "label": "xbcp_reg_file.R",
                "cmd": ["Rscript", "xbcp_reg_file.R", "{data_file}"],
            },
            {
                "label": "xsim_bcp_reg_file.exe",
                "cmd": ["xsim_bcp_reg_file.exe", "{data_file}"],
                "requires": ["xsim_bcp_reg_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_bcp_reg_file"
                ),
            },
        ],
    },
    "xmcp_demo": {
        "data_file": "xmcp_demo_data.txt",
        "steps": [
            {
                "label": "xmcp_demo_make_data.py",
                "cmd": ["python", "xmcp_demo_make_data.py", "{data_file}"],
            },
            {
                "label": "xmcp_demo_file.R",
                "cmd": ["Rscript", "xmcp_demo_file.R", "{data_file}"],
            },
            {
                "label": "xsim_mcp_demo_file.exe",
                "cmd": ["xsim_mcp_demo_file.exe", "{data_file}"],
                "requires": ["xsim_mcp_demo_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_mcp_demo_file"
                ),
            },
        ],
    },
    "xmcp_fit": {
        "data_file": "xmcp_fit_data.txt",
        "steps": [
            {
                "label": "xmcp_fit_make_data.py",
                "cmd": ["python", "xmcp_fit_make_data.py", "{data_file}"],
            },
            {
                "label": "xmcp_fit_file.R",
                "cmd": ["Rscript", "xmcp_fit_file.R", "{data_file}"],
            },
            {
                "label": "xsim_mcp_fit_file.exe",
                "cmd": ["xsim_mcp_fit_file.exe", "{data_file}"],
                "requires": ["xsim_mcp_fit_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_mcp_fit_file"
                ),
            },
        ],
    },
    "xmcp_sigma": {
        "data_file": "xmcp_sigma_data.txt",
        "steps": [
            {
                "label": "xmcp_sigma_make_data.py",
                "cmd": ["python", "xmcp_sigma_make_data.py", "{data_file}"],
            },
            {
                "label": "xmcp_sigma_file.R",
                "cmd": ["Rscript", "xmcp_sigma_file.R", "{data_file}"],
            },
            {
                "label": "xsim_mcp_sigma_file.exe",
                "cmd": ["xsim_mcp_sigma_file.exe", "{data_file}"],
                "requires": ["xsim_mcp_sigma_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_mcp_sigma_file"
                ),
            },
        ],
    },
    "xmcp_ar": {
        "data_file": "xmcp_ar_data.txt",
        "steps": [
            {
                "label": "xmcp_ar_make_data.py",
                "cmd": ["python", "xmcp_ar_make_data.py", "{data_file}"],
            },
            {
                "label": "xmcp_ar_file.R",
                "cmd": ["Rscript", "xmcp_ar_file.R", "{data_file}"],
            },
            {
                "label": "xsim_mcp_ar_file.exe",
                "cmd": ["xsim_mcp_ar_file.exe", "{data_file}"],
                "requires": ["xsim_mcp_ar_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_mcp_ar_file"
                ),
            },
        ],
    },
    "xmcp_ar_change": {
        "data_file": "xmcp_ar_change_data.txt",
        "steps": [
            {
                "label": "xmcp_ar_change_make_data.py",
                "cmd": ["python", "xmcp_ar_change_make_data.py", "{data_file}"],
            },
            {
                "label": "xmcp_ar_change_file.R",
                "cmd": ["Rscript", "xmcp_ar_change_file.R", "{data_file}"],
            },
            {
                "label": "xsim_mcp_ar_change_file.exe",
                "cmd": ["xsim_mcp_ar_change_file.exe", "{data_file}"],
                "requires": ["xsim_mcp_ar_change_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_mcp_ar_change_file"
                ),
            },
        ],
    },
    "xmcp_arsigma": {
        "data_file": "xmcp_arsigma_data.txt",
        "steps": [
            {
                "label": "xmcp_arsigma_make_data.py",
                "cmd": ["python", "xmcp_arsigma_make_data.py", "{data_file}"],
            },
            {
                "label": "xmcp_arsigma_file.R",
                "cmd": ["Rscript", "xmcp_arsigma_file.R", "{data_file}"],
            },
            {
                "label": "xsim_mcp_arsigma_file.exe",
                "cmd": ["xsim_mcp_arsigma_file.exe", "{data_file}"],
                "requires": ["xsim_mcp_arsigma_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_mcp_arsigma_file"
                ),
            },
        ],
    },
    "xsegmented": {
        "data_file": "xsegmented_data.txt",
        "steps": [
            {
                "label": "xsegmented_make_data.py",
                "cmd": ["python", "xsegmented_make_data.py", "{data_file}"],
            },
            {
                "label": "xsegmented_file.R",
                "cmd": ["Rscript", "xsegmented_file.R", "{data_file}"],
            },
            {
                "label": "xsim_segmented_file.exe",
                "cmd": ["xsim_segmented_file.exe", "{data_file}"],
                "requires": ["xsim_segmented_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_segmented_file"
                ),
            },
        ],
    },
    "xsegmented_multi": {
        "data_file": "xsegmented_multi_data.txt",
        "steps": [
            {
                "label": "xsegmented_multi_make_data.py",
                "cmd": ["python", "xsegmented_multi_make_data.py", "{data_file}"],
            },
            {
                "label": "xsegmented_multi_file.R",
                "cmd": ["Rscript", "xsegmented_multi_file.R", "{data_file}"],
            },
            {
                "label": "xsim_segmented_multi_file.exe",
                "cmd": ["xsim_segmented_multi_file.exe", "{data_file}"],
                "requires": ["xsim_segmented_multi_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_segmented_multi_file"
                ),
            },
        ],
    },
    "xstepmented": {
        "data_file": "xstepmented_data.txt",
        "steps": [
            {
                "label": "xstepmented_make_data.py",
                "cmd": ["python", "xstepmented_make_data.py", "{data_file}"],
            },
            {
                "label": "xstepmented_file.R",
                "cmd": ["Rscript", "xstepmented_file.R", "{data_file}"],
            },
            {
                "label": "xsim_stepmented_file.exe",
                "cmd": ["xsim_stepmented_file.exe", "{data_file}"],
                "requires": ["xsim_stepmented_file.exe"],
                "build_hint": (
                    "Build it first with:\n"
                    "  make -f Makefile.xcorr xsim_stepmented_file"
                ),
            },
        ],
    },
}


def run_step(index: int, total: int, step: dict[str, object], data_file: str) -> float:
    label = str(step["label"])
    print(f"[{index}/{total}] Running {label}", flush=True)

    for required in step.get("requires", []):
        if not Path(required).exists():
            print(f"{required} not found.")
            hint = step.get("build_hint")
            if hint:
                print(hint)
            raise SystemExit(1)

    cmd = [str(part).format(data_file=data_file) for part in step["cmd"]]
    t0 = time.perf_counter()
    completed = subprocess.run(cmd, check=False)
    elapsed = time.perf_counter() - t0
    if completed.returncode != 0:
        raise SystemExit(completed.returncode)
    return elapsed


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("cases", nargs="*", choices=sorted(CASES))
    parser.add_argument("--all", action="store_true", dest="run_all")
    parser.add_argument("--data-file", dest="data_file", default=None)
    args = parser.parse_args()

    if args.run_all and args.cases:
        parser.error("use either case names or --all, not both")
    if not args.run_all and not args.cases:
        parser.error("specify one or more case names, or use --all")

    case_names = list(CASES) if args.run_all else args.cases

    ncases = len(case_names)
    for icase, case_name in enumerate(case_names, start=1):
        case = CASES[case_name]
        data_file = args.data_file or case["data_file"]
        steps = case["steps"]

        if ncases > 1:
            print(f"=== {case_name} ({icase}/{ncases}) ===", flush=True)

        timings: list[tuple[str, float]] = []
        total = len(steps)

        try:
            for i, step in enumerate(steps, start=1):
                elapsed = run_step(i, total, step, data_file)
                timings.append((str(step["label"]), elapsed))
                if i < total:
                    print(flush=True)
        except SystemExit as exc:
            print()
            print("Comparison run failed.")
            return int(exc.code) if isinstance(exc.code, int) else 1

        print()
        labels = [label for label, _ in timings]
        values = [elapsed for _, elapsed in timings]
        widths = [max(len(label), len(f"{elapsed:.3f}")) for label, elapsed in timings]
        row_label_width = len("absolute")

        print("time elapsed (s)")
        label_line = " " * (row_label_width + 2) + "  ".join(
            f"{label:>{width}}" for label, width in zip(labels, widths)
        )
        value_line = "  ".join(f"{elapsed:>{width}.3f}" for elapsed, width in zip(values, widths))
        print(label_line)
        print(f"absolute  {value_line}")

        baseline = values[1] if len(values) >= 2 else values[0]
        rel_line = "  ".join(
            f"{(elapsed / baseline):>{width}.3f}" for elapsed, width in zip(values, widths)
        )
        print(f"relative  {rel_line}")
        print("\nCompleted successfully.")
        if icase < ncases:
            print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
