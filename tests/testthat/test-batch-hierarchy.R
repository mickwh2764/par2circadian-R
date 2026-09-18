simulate_matrix <- function(persistence, n = 24, seed = 1) {
  set.seed(seed)
  t(vapply(persistence, function(p) {
    as.numeric(stats::arima.sim(list(ar = c(0, p^2)), n = n))
  }, numeric(n)))
}

test_that("fit_ar2_batch labels, orders and returns one row per fitted gene", {
  m <- simulate_matrix(c(0.7, 0.6, 0.5))
  rownames(m) <- c("Per2", "Myc", "Actb")
  out <- fit_ar2_batch(m)

  expect_s3_class(out, "par2_batch")
  expect_identical(nrow(out), 3L)
  expect_identical(sort(out$gene), sort(rownames(m)))
  expect_true(all(diff(out$eigenvalue) <= 0))
})

test_that("fit_ar2_batch drops missing values and skips genes left too short", {
  m <- simulate_matrix(rep(0.6, 3), n = 12)
  m[2L, 1:8] <- NA
  rownames(m) <- c("Keep1", "TooShort", "Keep2")
  out <- fit_ar2_batch(m)
  expect_identical(sort(out$gene), c("Keep1", "Keep2"))
})

test_that("fit_ar2_batch names genes when the matrix has no rownames", {
  out <- fit_ar2_batch(simulate_matrix(c(0.6, 0.5)))
  expect_identical(sort(out$gene), c("Gene_1", "Gene_2"))
  expect_error(fit_ar2_batch(simulate_matrix(c(0.6, 0.5)), gene_names = "only-one"),
    "one entry per row"
  )
})

test_that("fit_ar2_batch returns an empty frame rather than failing on unfittable input", {
  out <- fit_ar2_batch(matrix(1:6, nrow = 2))
  expect_identical(nrow(out), 0L)
  expect_true(all(c("gene", "eigenvalue", "root_type") %in% names(out)))
})

test_that("batch fits are identical to fitting each gene alone", {
  m <- simulate_matrix(c(0.7, 0.5), seed = 8)
  out <- fit_ar2_batch(m, gene_names = c("A", "B"))
  expect_identical(out$eigenvalue[out$gene == "A"], fit_ar2(m[1L, ])$eigenvalue)
  expect_identical(out$eigenvalue[out$gene == "B"], fit_ar2(m[2L, ])$eigenvalue)
})

test_that("classify_gene_layer uses the clock and target sets", {
  expect_identical(classify_gene_layer("Per2"), "Clock")
  expect_identical(classify_gene_layer("ARNTL"), "Clock")
  expect_identical(classify_gene_layer("Myc"), "Target")
  expect_identical(classify_gene_layer("Actb"), "Background")
})

test_that("discover_hierarchy recovers an imposed Clock > Target > Background order", {
  results <- data.frame(
    gene = c("Per2", "Arntl", "Myc", "Wee1", "Actb", "Gapdh"),
    eigenvalue = c(0.75, 0.71, 0.58, 0.56, 0.47, 0.45),
    stringsAsFactors = FALSE
  )
  h <- discover_hierarchy(results)

  expect_identical(h$clock_median, 0.73)
  expect_identical(h$target_median, 0.57)
  expect_identical(h$background_median, 0.46)
  expect_identical(h$gearbox_gap, 0.16)
  expect_true(h$hierarchy_preserved)
  expect_identical(h$health_grade, "A")
  expect_identical(c(h$n_clock, h$n_target, h$n_background), c(2L, 2L, 2L))
  expect_identical(h$clock_genes$gene, c("Per2", "Arntl"))
})

test_that("a collapsed gap grades F and is not reported as preserved", {
  results <- data.frame(
    gene = c("Per2", "Myc", "Actb"),
    eigenvalue = c(0.50, 0.50, 0.60),
    stringsAsFactors = FALSE
  )
  h <- discover_hierarchy(results)
  expect_identical(h$health_grade, "F")
  expect_false(h$hierarchy_preserved)
})

test_that("discover_hierarchy accepts custom gene sets and rejects bad input", {
  results <- data.frame(
    gene = c("kaiA", "kaiB", "other"),
    eigenvalue = c(0.8, 0.6, 0.4),
    stringsAsFactors = FALSE
  )
  h <- discover_hierarchy(results, clock_genes = c("kaiA"), target_genes = c("kaiB"))
  expect_identical(h$n_clock, 1L)
  expect_identical(h$gearbox_gap, 0.2)
  expect_error(discover_hierarchy(data.frame(x = 1)), "'gene' and 'eigenvalue'")
})

test_that("gearbox_gap matches the difference of medians and is 0 when empty", {
  expect_equal(gearbox_gap(c(0.72, 0.70), c(0.61, 0.59)), 0.11)
  expect_identical(gearbox_gap(numeric(0), c(0.5)), 0)
})
