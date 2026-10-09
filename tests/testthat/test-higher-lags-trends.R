for (settings in list(list(plag = 2L, trend = FALSE),
                      list(plag = 1L, trend = TRUE),
                      list(plag = 2L, trend = TRUE))) {
  label <- paste("lags", settings$plag, "trend", settings$trend)

  test_that(paste(label, "IRF implementations agree and FEVD shares sum to one"), {
    model <- do.call(postprocess_model, settings)
    compiled <- irf(model, n.ahead = 4, expert = list(save.store = TRUE), verbose = FALSE)
    fallback <- irf(model, n.ahead = 4, expert = list(use_R = TRUE), verbose = FALSE)
    expect_equal(dim(compiled$posterior), c(4L, 5L, 4L, 7L))
    expect_true(all(is.finite(compiled$posterior)))
    expect_equal(compiled$posterior, fallback$posterior, tolerance = 1e-8)
    expect_equal(dim(compiled$struc.obj$Fmat), c(4L, 4L, settings$plag))
    decomposition <- fevd(compiled, verbose = FALSE)
    expect_true(all(is.finite(decomposition$FEVD)))
    expect_equal(unname(apply(decomposition$FEVD, c(2, 3), sum)),
                 matrix(1, 4, dim(decomposition$FEVD)[3]), tolerance = 1e-8)
  })

  test_that(paste(label, "HD reconstructs observations with the correct components"), {
    shocks <- irf(do.call(postprocess_model, settings), n.ahead = 4, verbose = FALSE)
    fit <- hd(shocks, verbose = FALSE)
    expect_equal(dim(fit$hd_array), c(80L - settings$plag, 4L,
                                    7L + as.integer(settings$trend)))
    expect_true(all(is.finite(fit$hd_array)))
    expect_equal(unname(apply(fit$hd_array, c(1, 2), sum)),
                 unname(fit$xglobal), tolerance = 1e-8)
    expect_equal("trend" %in% dimnames(fit$hd_array)[[3]], settings$trend)
    sigma <- shocks$struc.obj$Ginv %*% shocks$struc.obj$Smat %*% t(shocks$struc.obj$Ginv)
    expect_equal(unname(fit$struc_shock %*% chol(sigma)),
                 unname(postprocess_residuals(shocks)), tolerance = 1e-8)
  })

  test_that(paste(label, "forecasts retain finite posterior draws and quantiles"), {
    model <- do.call(postprocess_model, settings)
    set.seed(653)
    fit <- predict(model, n.ahead = 3, save.store = TRUE, verbose = FALSE)
    expect_equal(dim(fit$fcast), c(4L, 3L, 7L))
    expect_equal(dim(fit$pred_store), c(30L, 4L, 3L))
    expect_true(all(is.finite(fit$pred_store)))
    expect_equal(unname(fit$fcast[, , "Q50"]),
                 unname(apply(fit$pred_store, c(2, 3), median)))
  })
}
