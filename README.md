# Changepoint-Library

`Changepoint-Library` is a collection of changepoint, breakpoint, and segmentation methods implemented primarily in Fortran, with Python and R scripts used to generate data and compare against established reference packages.

The project is comparison-driven: most workflows generate deterministic data, run a reference implementation, run the Fortran implementation, and print matching summaries or checksums.

See [`COVERAGE.md`](COVERAGE.md) for the full list of registered comparison cases, and [`ALGORITHMS.md`](ALGORITHMS.md) for a Fortran module guide.

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
| `Makefile` | main build file for Fortran executables |

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

Build individual Fortran comparison executables with `Makefile`:

```powershell
make xsim_changepointnp_file
```

For another example:

```powershell
make xsim_mcp_arsigma_file
```

The makefile contains targets for the comparison executables used by `xrun_compare.py`.

## Running Comparisons

Use `xrun_compare.py` with a registered case name:

```powershell
python xrun_compare.py xchangepointnp
```

Run all registered cases, or a limited prefix of them:

```powershell
python xrun_compare.py --all
python xrun_compare.py --limit 10
```

By default, a failed comparison is reported and the runner continues with later cases. Use `--fail-fast` to stop at the first failed case.

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
make xsim_cpm_ks_file
python xrun_compare.py xcpm_ks
```

## Timing Results

The uploaded `results.txt` file is an example full run of:

```powershell
python xrun_compare.py --all
```

At the end of the run, `xrun_compare.py` prints a detector-time summary by implementation language. Python data-generation scripts named `*_make_data.py` are excluded from the language comparison because they create shared input files used by all implementations. Their cost is reported separately as setup time.

From the uploaded `results.txt`:

```text
summary detector time by language (s)
language  steps    total   mean  median  geomean  share  avg_rank
Python       42  347.769  8.280   1.987    2.482  0.635     2.071
R            86  163.418  1.900   0.609    0.740  0.298     2.035
Fortran     117   36.629  0.313   0.171    0.172  0.067     1.043

setup time (s)
steps  total   mean  median  geomean
  117  84.329  0.721   0.231    0.345
```

The timing columns mean:

- `steps`: number of detector/reference steps run for that language.
- `total`: total wall-clock seconds spent in that language across the run.
- `mean`: arithmetic average per detector step.
- `median`: typical detector step time, robust to one very slow case.
- `geomean`: geometric mean per detector step, useful for multiplicative speed comparisons.
- `share`: fraction of detector runtime spent in that language, excluding setup.
- `avg_rank`: average per-case speed rank, with `1` fastest. Ranks compare language totals within each case, then average over cases where the language appears.

In this run, the Fortran implementations covered all 117 cases, had the smallest total detector time, and had an average speed rank near 1. The Python and R rows are benchmark/reference timings only; they do not include the shared data-generation setup.

## Validation Style

The project emphasizes numerical agreement with reference implementations. Most comparison scripts print:

- detected changepoints or breakpoints
- fitted parameter summaries
- posterior means or probabilities for Bayesian methods
- checksums for vectors, matrices, statistics, or fitted values

Tiny floating-point differences are expected across languages and compilers.

## Notes

This repository is organized as a research and replication library rather than a single public API. Many files are standalone experiments or comparison drivers. Stable reuse points are the Fortran modules such as `changepoint.f90`, `cpm_pkg.f90`, `ecp_pkg.f90`, `bcp_pkg.f90`, `mcp_pkg.f90`, `segmented_pkg.f90`, and related package modules.
