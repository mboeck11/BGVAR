test_that("forecast evaluation matches independently computed means, errors and scores", {
  model <- postprocess_model(hold.out = 3L)
  set.seed(361)
  forecast <- predict(model, n.ahead = 3, save.store = TRUE, verbose = FALSE)
  observed <- model$args$yfull[78:80, , drop = FALSE]
  means <- t(apply(forecast$pred_store, c(2, 3), mean))
  sds <- t(apply(forecast$pred_store, c(2, 3), sd))
  expect_equal(forecast$hold.out.sample, observed)
  expect_equal(unname(t(forecast$lps.stats[, "mean", ])), unname(means))
  expect_equal(unname(t(forecast$lps.stats[, "sd", ])), unname(sds))
  errors <- rmse(forecast)
  scores <- lps(forecast)
  expect_s3_class(errors, "bgvar.rmse")
  expect_s3_class(scores, "bgvar.lps")
  expect_equal(dim(errors), c(3L, 4L))
  expect_equal(dim(scores), c(3L, 4L))
  # The API returns a separate error at each horizon, not one aggregate RMSE.
  expect_equal(as.numeric(errors), as.numeric(abs(observed - means)))
  expect_equal(as.numeric(scores),
               dnorm(as.numeric(observed), as.numeric(means), as.numeric(sds), log = TRUE))
  expect_true(all(is.finite(scores)))
})

test_that("evaluation without hold-out data is rejected", {
  forecast <- predict(postprocess_model(), n.ahead = 3, verbose = FALSE)
  expect_error(rmse(forecast), "hold out sample")
  expect_error(lps(forecast), "hold out sample")
})

test_that("a single forecast horizon retains its evaluation dimensions", {
  forecast <- predict(postprocess_model(hold.out = 1L), n.ahead = 1,
                       save.store = TRUE, verbose = FALSE)
  expect_equal(dim(forecast$lps.stats), c(4L, 2L, 1L))
  expect_equal(dim(rmse(forecast)), c(1L, 4L))
  expect_equal(dim(lps(forecast)), c(1L, 4L))
})

test_that("a forecast shorter than the hold-out sample evaluates the first held-out observations", {
  model <- postprocess_model(hold.out = 5L)
  forecast <- predict(model, n.ahead = 2, verbose = FALSE)
  expect_equal(forecast$hold.out.sample, model$args$yfull[76:77, , drop = FALSE])
})
