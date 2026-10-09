test_that("Cholesky FEVD shares are finite, nonnegative and sum to one", {
  shocks <- irf(postprocess_model(), n.ahead = 4, verbose = FALSE)
  fit <- fevd(shocks, verbose = FALSE)
  expect_s3_class(fit, "bgvar.fevd")
  expect_true(all(is.finite(fit$FEVD)))
  expect_true(all(fit$FEVD >= -1e-10 & fit$FEVD <= 1 + 1e-10))
  expect_equal(unname(apply(fit$FEVD, c(2, 3), sum)),
               matrix(1, ncol(fit$FEVD), dim(fit$FEVD)[3]), tolerance = 1e-8)
  selected <- fevd(shocks, var.slct = c("AA.y", "BB.p"), verbose = FALSE)
  expect_equal(selected$FEVD, fit$FEVD[, c("Decomp. of AA.y", "Decomp. of BB.p"), , drop = FALSE])
})

test_that("FEVD rejects unknown variables and unsupported identification", {
  expect_error(fevd(postprocess_irf(), var.slct = "XX.y", verbose = FALSE),
               "not contained")
  expect_error(fevd(postprocess_irf("girf"), verbose = FALSE), "only")
  missing <- postprocess_irf()
  missing$ident <- "sign"
  missing$struc.obj$Rmed <- NA_real_
  expect_error(fevd(missing, verbose = FALSE), "No rotation matrix")
})
test_that("sign-restricted FEVD handles shocks in multiple countries", {
  shocks <- postprocess_sign_irf()
  expect_equal(sort(unique(sub("[.].*$", "", shocks$shockinfo$shock))), c("AA", "BB"))
  fit <- fevd(shocks, verbose = FALSE)
  expect_s3_class(fit, "bgvar.fevd")
  expect_equal(fit$rotation.matrix, shocks$struc.obj$Rmed)
  expect_true(all(is.finite(fit$FEVD)))
  expect_true(all(fit$FEVD >= -1e-10 & fit$FEVD <= 1 + 1e-10))
  expect_equal(unname(apply(fit$FEVD, c(2, 3), sum)),
               matrix(1, 4, dim(fit$FEVD)[3]), tolerance = 1e-8)
  selected <- fevd(shocks, var.slct = "BB.y", verbose = FALSE)
  expect_equal(selected$FEVD, fit$FEVD[, "Decomp. of BB.y", , drop = FALSE])
})
