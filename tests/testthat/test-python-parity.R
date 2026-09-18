read_fixture <- function(name) {
  path <- system.file("extdata", name, package = "par2circadian")
  skip_if(path == "", paste("fixture", name, "not installed"))
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

# The bar is bit-level, not "close": the fit can be got subtly wrong in ways a
# correlation or a mean error hides (Yule-Walker instead of OLS, biased variance
# divisor, no mean-centring, a different lag alignment), each of which leaves
# r > 0.999 while moving individual genes by enough to cross a persistence
# threshold. So compare the worst gene, not the average one.
#
# Both implementations round their reported values to 6 decimal places, so the
# observed worst-gene difference over the fixture is exactly zero. The tolerance
# here is slack for the last reported digit only -- anything a wrong estimator
# would produce is orders of magnitude larger.
TOLERANCE <- 1e-12

test_that("every gene matches the Python implementation at reported precision", {
  series <- read_fixture("gse11923_subset.csv")
  reference <- read_fixture("python_reference.csv")

  values <- as.matrix(series[, -1L, drop = FALSE])
  fitted <- fit_ar2_batch(values, gene_names = series[[1L]])

  expect_equal(nrow(fitted), nrow(reference))

  fitted <- fitted[match(reference$gene, fitted$gene), , drop = FALSE]
  expect_identical(fitted$gene, reference$gene)
  expect_identical(fitted$root_type, reference$root_type)

  worst <- function(column) {
    a <- fitted[[column]]
    b <- reference[[column]]
    expect_identical(is.na(a), is.na(b), info = column)
    keep <- !is.na(a)
    if (!any(keep)) {
      return(0)
    }
    max(abs(a[keep] - b[keep]))
  }

  for (column in c("eigenvalue", "phi1", "phi2", "r2", "half_life", "eigenperiod")) {
    expect_lt(worst(column), TOLERANCE, label = paste("max abs difference in", column))
  }

  # And at the precision both implementations report, the agreement is exact.
  expect_identical(fitted$eigenvalue, reference$eigenvalue)
  expect_identical(fitted$phi1, reference$phi1)
  expect_identical(fitted$phi2, reference$phi2)
})

test_that("the fixture exercises both root types and a spread of persistence", {
  reference <- read_fixture("python_reference.csv")
  expect_true(all(c("Real", "Complex") %in% reference$root_type))
  expect_gt(diff(range(reference$eigenvalue)), 0.3)
})
