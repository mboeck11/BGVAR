test_that("forecast quantiles summarize the stored posterior draws", {
  set.seed(651)
  fit <- predict(postprocess_model(), n.ahead = 3,
                 quantiles = c(0.16, 0.5, 0.84), save.store = TRUE, verbose = FALSE)
  expect_s3_class(fit, "bgvar.pred")
  expect_equal(dim(fit$fcast), c(4L, 3L, 3L))
  expect_true(all(is.finite(fit$fcast)))
  expect_true(all(fit$fcast[, , 1] <= fit$fcast[, , 2]))
  expect_true(all(fit$fcast[, , 2] <= fit$fcast[, , 3]))
  for (i in 1:3) {
    expect_equal(unname(fit$fcast[, , i]),
                 unname(apply(fit$pred_store, c(2, 3), quantile,
                              probs = c(0.16, 0.5, 0.84)[i])))
  }
  set.seed(651)
  compact <- predict(postprocess_model(), n.ahead = 3,
                     quantiles = c(0.16, 0.5, 0.84), verbose = FALSE)
  expect_equal(compact$fcast, fit$fcast)
  expect_null(compact$pred_store)
})

test_that("hard conditional forecasts honor the supplied path", {
  model <- postprocess_model()
  path <- matrix(NA_real_, 3, 4, dimnames = list(NULL, colnames(model$xglobal)))
  path[, "AA.y"] <- c(0.2, 0.4, 0.6)
  set.seed(61)
  fit <- predict(model, n.ahead = 3, constr = path,
                 quantiles = c(0.16, 0.5, 0.84), verbose = FALSE)
  expect_equal(unname(fit$fcast["AA.y", , ]),
               matrix(path[, "AA.y"], 3, 3), tolerance = 1e-8)
  expect_true(all(is.finite(fit$fcast)))
})

test_that("forecasting validates quantiles and conditional matrix dimensions", {
  model <- postprocess_model()
  expect_error(predict(model, quantiles = "median", verbose = FALSE), "numeric vector")
  expect_error(predict(model, n.ahead = 3, constr = matrix(NA, 2, 4), verbose = FALSE),
               "dimensions of 'constr'")
  expect_error(predict(model, n.ahead = 3, constr = matrix(NA, 3, 4),
                       constr_sd = matrix(0, 2, 4), verbose = FALSE),
               "dimensions of 'constr_sd'")
})
