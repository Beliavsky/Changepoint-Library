# Changepoint-Library

`Changepoint-Library` is a collection of changepoint, breakpoint, and segmentation methods implemented primarily in Fortran, with Python and R scripts used to generate data and compare against established reference packages.

The project is comparison-driven: most workflows generate deterministic data, run a reference implementation, run the Fortran implementation, and print matching summaries or checksums.

## Scope

The repository covers a broad range of changepoint methods:

| Area | Reference packages / methods |
| --- | --- |
| Classical changepoints | mean, variance, mean-variance, AMOC, PELT, binary segmentation, dynamic programming |
| Structural breaks | `strucchange`-style breakpoints, F statistics, empirical fluctuation processes |
| Sequential CPM methods | `cpm` families including Student, Bartlett, Joint, Exponential, Poisson, MW, Mood, Lepage, FET, KS, CVM |
| Energy and nonparametric methods | `ecp` methods including `e.divisive`, `e.agglo`, `e.cp3o`, and `e.cp3o_delta` |
| Nonparametric PELT | `changepoint.np::cpt.np` empirical-distribution path |
| Bayesian changepoints | `bcp` univariate, multivariate, and regression paths |
| Bayesian segmented regression | `mcp` Gaussian, `sigma(...)`, `ar(...)`, and combined AR plus variance models |
| Classical segmented regression | `segmented()` and `stepmented()` |
| Other detectors | BOCPD, BEAST, ClaSP, RuLSIF/uLSIF, Roerich-style methods, and related experiments |

## Repository Layout

The files follow a regular naming convention:

| Pattern | Meaning |
| --- | --- |
| `*_pkg.f90` | reusable Fortran implementation modules |
| `xsim_*_file.f90` | Fortran comparison executables that read data files |
| `x*_make_data.py` | deterministic data generators |
| `x*_file.R` | R reference implementations |
| `x*_file.py` | Python reference implementations |
| `xrun_compare.py` | comparison runner for registered cases |
| `Makefile.xcorr` | main build file for Fortran executables |

There are hundreds of source files. The best entry point is usually the comparison case name in `xrun_compare.py`, then the matching generator, reference script, and Fortran driver.

## Requirements

Typical workflows need:

- `gfortran`
- `make`
- Python 3
- NumPy for many Python data generators
- R for R reference comparisons
- selected R packages depending on the comparison case, such as `changepoint`, `strucchange`, `cpm`, `ecp`, `changepoint.np`, `bcp`, `mcp`, and `segmented`

Some optional comparisons require additional package-specific dependencies, for example JAGS for `mcp` sampling through `rjags`.

## Build

Build individual Fortran comparison executables with `Makefile.xcorr`:

```powershell
make -f Makefile.xcorr xsim_changepointnp_file
```

For another example:

```powershell
make -f Makefile.xcorr xsim_mcp_arsigma_file
```

The makefile contains targets for the comparison executables used by `xrun_compare.py`.

## Running Comparisons

Use `xrun_compare.py` with a registered case name:

```powershell
python xrun_compare.py xchangepointnp
```

A comparison usually performs three steps:

1. Generate a deterministic data file.
2. Run the R or Python reference implementation.
3. Run the Fortran executable and print comparable summaries.

Representative cases:

```powershell
python xrun_compare.py xcpm_ks
python xrun_compare.py xecp_edivisive_full
python xrun_compare.py xchangepointnp
python xrun_compare.py xbcp_reg
python xrun_compare.py xmcp_arsigma
python xrun_compare.py xstepmented
```

Build the corresponding executable first when a case requires one:

```powershell
make -f Makefile.xcorr xsim_cpm_ks_file
python xrun_compare.py xcpm_ks
```

## Validation Style

The project emphasizes numerical agreement with reference implementations. Most comparison scripts print:

- detected changepoints or breakpoints
- fitted parameter summaries
- posterior means or probabilities for Bayesian methods
- checksums for vectors, matrices, statistics, or fitted values

Tiny floating-point differences are expected across languages and compilers.

## Notes

This repository is organized as a research and replication library rather than a single public API. Many files are standalone experiments or comparison drivers. Stable reuse points are the Fortran modules such as `changepoint.f90`, `cpm_pkg.f90`, `ecp_pkg.f90`, `bcp_pkg.f90`, `mcp_pkg.f90`, `segmented_pkg.f90`, and related package modules.
