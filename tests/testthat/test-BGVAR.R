# Small, reproducible models keep input tests independent of bundled datasets.
bgvar_inputs <- function() {
  set.seed(1209)
  Data <- list(AA = matrix(rnorm(160), 80, 2),
               BB = matrix(rnorm(160), 80, 2))
  Data <- lapply(Data, function(x) { colnames(x) <- c("y", "p"); x })
  W <- matrix(c(0, 1, 1, 0), 2, dimnames = list(names(Data), names(Data)))
  Ex <- matrix(rnorm(160), 80, 2,
               dimnames = list(NULL, c("AA.z", "BB.z")))
  list(Data = Data, W = W, Ex = Ex)
}

bgvar_fit <- function(...) {
  inputs <- bgvar_inputs()
  args <- list(Data = inputs$Data, W = inputs$W, draws = 20, burnin = 10,
               prior = "MN", SV = FALSE, eigen = FALSE, verbose = FALSE)
  overrides <- list(...)
  args[names(overrides)] <- overrides
  set.seed(921)
  # bgvar currently writes some progress output even with verbose = FALSE.
  output <- capture.output(result <- do.call(bgvar, args))
  attr(result, "test_output") <- output
  result
}

test_that("Data and W accept their documented representations", {
  inputs <- bgvar_inputs()
  baseline <- bgvar_fit()
  expect_s3_class(baseline, "bgvar")
  expect_equal(dim(baseline$xglobal), c(80L, 4L))
  matrix_data <- do.call(cbind, inputs$Data)
  colnames(matrix_data) <- c("AA.y", "AA.p", "BB.y", "BB.p")
  expect_equal(bgvar_fit(Data = matrix_data)$stacked.results,
               baseline$stacked.results)
  expect_equal(bgvar_fit(W = list(W = inputs$W))$stacked.results,
               baseline$stacked.results)
  reordered <- inputs$W[c("BB", "AA"), c("BB", "AA")]
  expect_equal(bgvar_fit(W = reordered)$stacked.results,
               baseline$stacked.results)
  time_data <- lapply(inputs$Data, ts, start = c(2000, 1), frequency = 4)
  expect_equal(bgvar_fit(Data = time_data)$xglobal, baseline$xglobal)
})

test_that("Data and W reject missing values, incompatible sizes and names", {
  inputs <- bgvar_inputs()
  bad_data <- inputs$Data
  bad_data$AA[1, 1] <- NA_real_
  expect_error(bgvar_fit(Data = bad_data), "contains NAs")
  bad_data <- inputs$Data
  bad_data$BB <- bad_data$BB[-1, ]
  expect_error(bgvar_fit(Data = bad_data), "same sample size")
  names(bad_data) <- c("AAA", "BB")
  expect_error(bgvar_fit(Data = bad_data), "exactly two characters")
  expect_error(bgvar_fit(W = 1), "argument 'W'")
  bad_w <- inputs$W
  bad_w[1, 2] <- NA_real_
  expect_error(bgvar_fit(W = bad_w), "weight matrix you have provided contains NAs")
  expect_error(bgvar_fit(W = matrix(0, 3, 3)), "same dimension")
  dimnames(bad_w) <- list(c("AA", "CC"), c("AA", "CC"))
  bad_w[1, 2] <- 1
  expect_error(bgvar_fit(W = bad_w), "same country names")
})

test_that("plag supports common and separate lag orders", {
  common <- bgvar_fit(plag = 2)
  separate <- bgvar_fit(plag = c(2, 1))
  expect_equal(common$args$plag, 2)
  expect_equal(separate$args$plag, c(2, 1))
  expect_true(any(grepl("_lag2", rownames(separate$cc.results$coeffs$AA))))
  expect_error(bgvar_fit(plag = "one"), "lags as numeric")
  expect_error(bgvar_fit(plag = NA_real_), "number of lags")
  expect_error(bgvar_fit(plag = c(1, 1, 1)), "One lag length")
})

test_that("draws, burnin and thin control posterior storage", {
  fit <- bgvar_fit(draws = 24, burnin = 12, thin = 3,
                   expert = list(save.country.store = TRUE))
  expect_equal(fit$args$draws, 24)
  expect_equal(fit$args$burnin, 12)
  expect_equal(fit$args$thindraws, 8)
  expect_equal(dim(fit$stacked.results$A_large)[3], 8L)
  expect_equal(bgvar_fit(thin = 0.5)$args$thin, 2)
  adjusted <- bgvar_fit(draws = 20, thin = 3)
  expect_equal(20 %% adjusted$args$thin, 0)
  for (argument in c("draws", "burnin")) {
    for (value in list("ten", -1, c(10, 20))) {
      expect_error(do.call(bgvar_fit, setNames(list(value), argument)),
                   "draws and burnin")
    }
  }
})

test_that("all priors and hyperparameter overrides can be estimated", {
  for (prior in c("MN", "SSVS", "NG", "HS")) {
    fit <- bgvar_fit(prior = prior, expert = list(save.shrink.store = TRUE))
    expect_s3_class(fit, "bgvar")
    expect_equal(fit$args$prior, prior)
    expect_true(all(is.finite(fit$stacked.results$A_large)))
  }
  fit <- bgvar_fit(prior = "SSVS", hyperpara = list(p_i = 0.3, q_ij = 0.7),
                   expert = list(save.shrink.store = TRUE))
  expect_equal(fit$args$hyperpara, list(p_i = 0.3, q_ij = 0.7))
  expect_false(is.null(fit$cc.results$PIP))
  expect_error(bgvar_fit(prior = "unknown"), "prior options")
  expect_warning(bgvar_fit(hyperpara = list(unknown_parameter = 1)),
                 "no valid hyperparameter")
})

test_that("SV, hold.out and trend alter the requested model", {
  fit <- bgvar_fit(SV = TRUE, hold.out = 5, trend = TRUE,
                   expert = list(save.country.store = TRUE,
                                 save.vola.store = TRUE))
  expect_true(fit$args$SV)
  expect_equal(nrow(fit$xglobal), 75L)
  expect_equal(nrow(fit$args$yfull), 80L)
  expect_true("trend" %in% rownames(fit$cc.results$coeffs$AA))
  expect_false(is.null(fit$cc.results$store$AA$pars_store))
  baseline <- bgvar_fit()
  expect_false(baseline$args$SV)
  expect_false("trend" %in% rownames(baseline$cc.results$coeffs$AA))
})

test_that("eigen accepts logical and numeric trimming settings", {
  for (threshold in list(TRUE, 1.05)) {
    fit <- bgvar_fit(eigen = threshold)
    expect_true(length(fit$stacked.results$F.eigen) > 0)
    expect_true(all(fit$stacked.results$F.eigen < 1.05))
  }
  expect_equal(bgvar_fit(eigen = FALSE)$args$thindraws, 20)
})

test_that("Ex accepts matrix and list inputs and checks their contents", {
  inputs <- bgvar_inputs()
  matrix_fit <- bgvar_fit(Ex = inputs$Ex)
  ex_list <- list(AA = inputs$Ex[, 1, drop = FALSE],
                  BB = inputs$Ex[, 2, drop = FALSE])
  ex_list <- lapply(ex_list, function(x) { colnames(x) <- "z"; x })
  expect_equal(bgvar_fit(Ex = ex_list)$stacked.results,
               matrix_fit$stacked.results)
  expect_true("z" %in% rownames(matrix_fit$cc.results$coeffs$AA))
  expect_error(bgvar_fit(Ex = 1), "argument 'Ex'")
  bad_ex <- inputs$Ex
  bad_ex[1, 1] <- NA_real_
  expect_error(bgvar_fit(Ex = bad_ex), "data for exogenous variables you have submitted contains NAs")
  expect_error(bgvar_fit(Ex = inputs$Ex[-1, ]), "not equally long")
  bad_ex <- inputs$Ex
  colnames(bad_ex) <- NULL
  expect_error(bgvar_fit(Ex = bad_ex), "column names")
  colnames(bad_ex) <- c("z1", "z2")
  expect_error(bgvar_fit(Ex = bad_ex), "with a point")
  colnames(bad_ex) <- c("CC.z", "BB.z")
  expect_error(bgvar_fit(Ex = bad_ex), "country names")
  duplicate <- inputs$Data$AA[, 1, drop = FALSE]
  colnames(duplicate) <- "AA.z"
  expect_error(bgvar_fit(Ex = duplicate), "also contained")
})

test_that("expert selects R estimation, storage and a custom apply function", {
  calls <- 0L
  custom_apply <- function(X, FUN, ...) {
    calls <<- calls + 1L
    lapply(X, FUN, ...)
  }
  fit <- bgvar_fit(expert = list(use_R = TRUE, applyfun = custom_apply,
                                 save.country.store = TRUE))
  expect_s3_class(fit, "bgvar")
  expect_equal(calls, 1L)
  expect_equal(names(fit$cc.results$store), c("AA", "BB"))
  expect_s3_class(bgvar_fit(expert = list(cores = 1)), "bgvar")
  expect_error(bgvar_fit(expert = list(cores = "invalid")), "argument 'cores'")
})

test_that("verbose changes progress output without changing estimates", {
  quiet <- bgvar_fit(verbose = FALSE)
  loud <- bgvar_fit(verbose = TRUE)
  expect_equal(loud$stacked.results, quiet$stacked.results)
  expect_true(any(grepl("Stacking of global model", attr(loud, "test_output"))))
  expect_false(any(grepl("Stacking of global model", attr(quiet, "test_output"))))
})
