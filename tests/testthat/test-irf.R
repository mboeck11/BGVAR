test_that("Cholesky and generalized IRFs have finite, ordered posterior quantiles", {
  for (ident in c("chol", "girf")) {
    fit <- postprocess_irf(ident)
    expect_s3_class(fit, "bgvar.irf")
    expect_equal(fit$ident, ident)
    expect_equal(dim(fit$posterior), c(4L, 5L, 1L, 3L))
    expect_true(all(is.finite(fit$posterior)))
    expect_true(all(fit$posterior[, , , 1] <= fit$posterior[, , , 2]))
    expect_true(all(fit$posterior[, , , 2] <= fit$posterior[, , , 3]))
    expect_equal(as.numeric(fit$posterior["AA.y", 1, 1, ]), rep(1, 3))
    expect_equal(dim(fit$IRF_store), c(4L, 1L, 5L, 30L))
  }
})

test_that("IRF implementations agree and respect shock scaling", {
  for (ident in c("chol", "girf")) {
    compiled <- postprocess_irf(ident)
    fallback <- postprocess_irf(ident, use_R = TRUE)
    expect_equal(compiled$posterior, fallback$posterior, tolerance = 1e-8)
    doubled <- postprocess_irf(ident, scale = 2)
    expect_equal(doubled$posterior, 2 * compiled$posterior, tolerance = 1e-8)
    compact <- postprocess_irf(ident, save = FALSE)
    expect_equal(compact$posterior, compiled$posterior)
    expect_null(compact$IRF_store)
  }
})

test_that("IRFs reject unknown shocks and nonnumeric quantiles", {
  model <- postprocess_model()
  info <- get_shockinfo("chol")
  info$shock <- "XX.unknown"
  expect_error(irf(model, shockinfo = info, verbose = FALSE), "variables available")
  expect_error(irf(model, quantiles = "median", verbose = FALSE), "numeric vector")
})

test_that("sign restrictions produce valid rotations and restricted responses", {
  set.seed(517)
  info <- add_shockinfo(get_shockinfo("sign"), shock = "AA.y",
                       restriction = "AA.y", sign = "<", horizon = 1,
                       prob = 1, scale = -1)
  fit <- irf(postprocess_model(), n.ahead = 4, shockinfo = info,
             expert = list(MaxTries = 1000, save.store = TRUE), verbose = FALSE)
  expect_true(all(is.finite(fit$posterior)))
  expect_true(all(fit$IRF_store["AA.y", 1, 1, ] < 0))
  expect_equal(dim(fit$struc.obj$Rmed), c(4L, 4L))
  expect_equal(unname(crossprod(fit$struc.obj$Rmed)), diag(4), tolerance = 1e-8)
})
