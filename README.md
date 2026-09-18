# par2circadian

An R implementation of AR(2) eigenvalue persistence analysis for gene
expression time series — the R port of the Python package
[`par2-circadian`](https://github.com/mickwh2764/par2discovery)
(RRID:SCR_028837).

Circadian tools ask *does this gene cycle?* `par2circadian` asks *how long does
its expression remember?* It fits

```
x[t] = phi1 * x[t-1] + phi2 * x[t-2] + e[t]
```

per gene and returns `|lambda|`, the modulus of the largest root of
`r^2 - phi1*r - phi2 = 0`: one number per gene measuring how strongly recent
expression constrains future expression.

## Install

```r
# install.packages("par2circadian")            # once on CRAN
# remotes::install_github("mickwh2764/par2circadian-R")
```

Depends on base R and `stats` only.

## Use

```r
library(par2circadian)

fit <- fit_ar2(expression_vector)          # one gene
fit$eigenvalue; fit$root_type

fits <- fit_ar2_batch(expression_matrix)   # genes x timepoints
discover_hierarchy(fits)$gearbox_gap       # Clock - Target median |lambda|

fit_ar2(short_series, n_bootstrap = 2000, seed = 1)$eigenvalue_ci
```

`vignette("par2circadian")` works through a mouse liver series end to end.

## What it is not

`|lambda|` is **not** a rhythmicity test, and should not be reported as one.
Red noise with no rhythm scores high — correctly, because it is persistent — and
a clean fixed-period cosine scores low. On noise-driven oscillators at 24
timepoints it detected 92% of true oscillators against 40–62% for JTK_CYCLE,
ARSER and RAIN; on a clean 24 h cosine those methods reached 100% and
`|lambda|` reached 12%. It answers a different question: use it alongside a
rhythm-detection method, not instead of one.

Two further limits worth respecting: single-gene estimates from short series are
imprecise (roughly a 0.4-wide 95% interval at 24 timepoints, biased upward
below that), so report `n_bootstrap` intervals; and a high `|lambda|` with real
roots is persistent but *not* oscillatory.

## Agreement with the Python implementation

The test suite refits a 300-gene subset of GSE11923 (mouse liver, hourly for
48 h) and compares every gene against what `par2-circadian` returns. All 300
agree exactly on `|lambda|`, `phi1`, `phi2`, `R^2`, half-life, eigenperiod and
root type at the 6 decimal places both implementations report — the worst-gene
difference is 0, not merely small. Comparing the worst gene rather than a mean
or a correlation is the point: a wrong estimator still correlates at 0.999 with
the right one. The fixtures and the script that regenerates them are in
`inst/extdata/` and `data-raw/`. Bootstrap intervals are the exception — R and
NumPy draw from different generators, so widths are comparable but endpoints
are not.

## Citation

```r
citation("par2circadian")
```

Cite the method preprint (doi:10.21203/rs.3.rs-9283100/v1) alongside the
software.

## Licence

Apache-2.0. Test fixtures derive from GEO series GSE11923 (Hughes et al. 2009),
a public dataset, subset for testing only.
