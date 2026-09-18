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
