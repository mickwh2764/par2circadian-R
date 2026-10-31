# par2circadian 0.2.0

* New `timescale()` converts `|lambda|` to the decay time constant
  tau = -sampling_interval / ln|lambda| in real time units, which, unlike
  `|lambda|`, can be compared across sampling designs.
* `fit_ar2()` and `fit_ar2_batch()` gain `sampling_interval`. When supplied, the
  fit carries `timescale`; with `n_bootstrap` too, `timescale_ci` and
  `timescale_imprecise`, which is `TRUE`, with a warning, when the interval is
  unbounded or spans more than tenfold. Matches par2-circadian 1.3.0.
* Existing outputs are unchanged when `sampling_interval` is not supplied.

# par2circadian 0.1.1

* `bootstrap_ar2()` and `fit_ar2()` no longer restore `.Random.seed` in the
  global environment after a seeded bootstrap. Supplying `seed` now simply
  calls `set.seed(seed)`; leave it `NULL` to draw from the caller's stream.
* Documented the return value of `CORE_CLOCK_GENES` and `KNOWN_TARGET_GENES`.
* `fit_ar2()` no longer caps `|lambda|` at 2. An unstable fit now reports the
  spectral radius of its fitted coefficients, matching par2-circadian 1.2.1.
  Stable fits (`|lambda|` < 1) are unchanged.

# par2circadian 0.1.0

First release: the R implementation of the AR(2) eigenvalue persistence method,
ported from the Python package par2-circadian 1.2.0.

* `fit_ar2()` for one series and `fit_ar2_batch()` for a genes-by-timepoints
  matrix, returning `|lambda|`, the AR(2) coefficients, R^2, root type,
  half-life and eigenperiod.
* `bootstrap_ar2()` residual-bootstrap confidence intervals, reachable through
  `fit_ar2(n_bootstrap =)`, for the short series circadian designs produce.
* `discover_hierarchy()`, `gearbox_gap()` and `classify_gene_layer()` for the
  Clock > Target > Background stratification; `classify_dynamics()` for
  root-type-aware labels.
* `half_life()` and `eigenperiod()` convert to real time units.
* Every gene of a 300-gene mouse liver dataset (GSE11923) is tested against the
  Python implementation and matches exactly at reported precision. Bootstrap
  intervals are not expected to match endpoint-for-endpoint, because R and NumPy
  use different random generators.
