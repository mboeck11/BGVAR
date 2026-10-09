test_that("historical decomposition reconstructs the observed series", {
  shocks <- postprocess_irf()
  fit <- hd(shocks, verbose = FALSE)
  expect_s3_class(fit, "bgvar.hd")
  expect_equal(dim(fit$hd_array), c(79L, 4L, 7L))
  expect_true(all(is.finite(fit$hd_array)))
  expect_equal(unname(apply(fit$hd_array, c(1, 2), sum)),
               unname(fit$xglobal), tolerance = 1e-8)
  expect_equal(dim(fit$struc_shock), c(79L, 4L))
  expect_true(all(is.finite(fit$struc_shock)))
  # Independently recover reduced-form innovations from structural shocks.
  model <- postprocess_model()
  X <- cbind(model$xglobal[-80, ], 1)
  residuals <- model$xglobal[-1, ] - X %*% t(shocks$struc.obj$A)
  sigma <- shocks$struc.obj$Ginv %*% shocks$struc.obj$Smat %*%
    t(shocks$struc.obj$Ginv)
  expect_equal(unname(fit$struc_shock %*% chol(sigma)),
               unname(residuals), tolerance = 1e-8)
})
test_that("sign-restricted HD uses the median rotation by default", {
  shocks <- postprocess_sign_irf()
  fit <- hd(shocks, verbose = FALSE)
  expect_s3_class(fit, "bgvar.hd")
  expect_true(all(is.finite(fit$hd_array)))
  expect_equal(unname(apply(fit$hd_array, c(1, 2), sum)),
               unname(fit$xglobal), tolerance = 1e-8)
  sigma <- shocks$struc.obj$Ginv %*% shocks$struc.obj$Smat %*% t(shocks$struc.obj$Ginv)
  impact <- t(chol(sigma)) %*% shocks$struc.obj$Rmed
  expect_equal(unname(fit$struc_shock %*% t(impact)),
               unname(postprocess_residuals(shocks)), tolerance = 1e-8)
  explicit <- hd(shocks, rotation.matrix = shocks$struc.obj$Rmed, verbose = FALSE)
  expect_equal(fit$hd_array, explicit$hd_array)
  expect_equal(fit$struc_shock, explicit$struc_shock)
  # An explicitly supplied rotation must be preserved, rather than overwritten.
  custom_rotation <- shocks$struc.obj$Rmed
  custom_rotation[, 1] <- -custom_rotation[, 1]
  custom <- hd(shocks, rotation.matrix = custom_rotation, verbose = FALSE)
  custom_impact <- t(chol(sigma)) %*% custom_rotation
  expect_equal(unname(custom$struc_shock %*% t(custom_impact)),
               unname(postprocess_residuals(shocks)), tolerance = 1e-8)
  expect_equal(custom$struc_shock[, 1], -fit$struc_shock[, 1])
  expect_equal(custom$struc_shock[, -1], fit$struc_shock[, -1])
})

test_that("sign-restricted HD rejects a missing median rotation", {
  shocks <- postprocess_sign_irf()
  for (missing in list(NULL, NA_real_)) {
    shocks$struc.obj$Rmed <- missing
    expect_error(hd(shocks, verbose = FALSE), "No rotation matrix available")
  }
})
