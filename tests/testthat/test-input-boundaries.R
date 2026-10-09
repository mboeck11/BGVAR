# Stop before any country sampler is called, so malformed sizes cannot reach
# compiled code or request large allocations. Invalid inputs must fail earlier.
expect_bgvar_boundary_error <- function(argument, value, hyperparameter = FALSE) {
  reached_estimation <- FALSE
  guard <- function(X, FUN, ...) {
    reached_estimation <<- TRUE
    stop("Input reached estimation without validation")
  }
  model <- postprocess_model()
  args <- list(Data = model$args$Data, W = model$args$W,
               prior = "SSVS", draws = 20, burnin = 10, SV = FALSE,
               eigen = FALSE, verbose = FALSE, expert = list(applyfun = guard))
  if (hyperparameter) {
    args$hyperpara <- setNames(list(value), argument)
  } else {
    args[argument] <- list(value)
  }
  error <- NULL
  invisible(capture.output(tryCatch(suppressWarnings(do.call(bgvar, args)),
                                   error = function(e) error <<- e)))
  label <- paste(argument, paste(value, collapse = ","))
  expect_true(inherits(error, "error"), info = label)
  expect_false(reached_estimation, info = label)
}

test_that("valid numerical endpoints remain supported", {
  model <- postprocess_model()
  set.seed(438)
  invisible(capture.output(fit <- bgvar(
    Data = model$args$Data, W = model$args$W, prior = "MN",
    draws = 20, burnin = 0, thin = 0.5, hold.out = 0,
    SV = FALSE, eigen = FALSE, verbose = FALSE
  )))
  expect_s3_class(fit, "bgvar")
  expect_equal(fit$args$burnin, 0)
  expect_equal(fit$args$thin, 2)
  expect_equal(fit$args$thindraws, 10)
  for (fun in list(predict, irf)) {
    result <- fun(model, n.ahead = 1, quantiles = c(0, 1), verbose = FALSE)
    values <- if (inherits(result, "bgvar.pred")) result$fcast else result$posterior
    expect_true(all(is.finite(values)))
  }
})

test_that("invalid lag counts are rejected before estimation", {
  for (value in list(0, -1, 1.5, NA_real_, Inf, NaN, numeric(0), c(1, 2, 3))) {
    expect_bgvar_boundary_error("plag", value)
  }
})

test_that("invalid draws and burn-in counts are rejected before estimation", {
  for (value in list(0, -1, 1.5, NA_real_, Inf, NaN, numeric(0), c(10, 20))) {
    expect_bgvar_boundary_error("draws", value)
  }
  # Zero burn-in is valid; only negative or noninteger/nonfinite counts fail.
  for (value in list(-1, 1.5, NA_real_, Inf, NaN, numeric(0), c(10, 20))) {
    expect_bgvar_boundary_error("burnin", value)
  }
})

test_that("invalid thinning and hold-out counts are rejected before estimation", {
  # Positive reciprocal thinning (e.g. 0.5) is intentionally supported.
  for (value in list(0, -1, NA_real_, Inf, NaN, numeric(0), c(1, 2))) {
    expect_bgvar_boundary_error("thin", value)
  }
  for (value in list(-1, 0.5, 80, 81, NA_real_, Inf, NaN, numeric(0), c(1, 2))) {
    expect_bgvar_boundary_error("hold.out", value)
  }
})

test_that("SSVS probabilities outside the unit interval are rejected", {
  for (argument in c("p_i", "q_ij")) {
    for (value in list(-0.1, 1.1, NA_real_, Inf, NaN, numeric(0), c(0.3, 0.7))) {
      expect_bgvar_boundary_error(argument, value, hyperparameter = TRUE)
    }
  }
})

test_that("prediction and IRFs reject invalid horizons and quantiles", {
  model <- postprocess_model()
  for (name in c("predict", "irf")) {
    fun <- get(name, mode = "function")
    for (horizon in list(0, -1, 1.5, NA_real_, Inf, numeric(0), c(1, 2))) {
      error <- tryCatch(suppressWarnings(fun(model, n.ahead = horizon, verbose = FALSE)),
                        error = identity)
      expect_true(inherits(error, "error"),
                  info = paste(name, "n.ahead", paste(horizon, collapse = ",")))
    }
    for (quantiles in list(-0.1, 1.1, NA_real_, Inf, numeric(0))) {
      error <- tryCatch(fun(model, n.ahead = 2, quantiles = quantiles, verbose = FALSE),
                        error = identity)
      expect_true(inherits(error, "error"),
                  info = paste(name, "quantiles", paste(quantiles, collapse = ",")))
    }
  }
})
