# Call each implementation directly: the public wrapper can silently fall back
# from C++ to R, which would otherwise hide compiled-sampler failures.
prior_sample <- function(implementation, prior, hyperpara = list(),
                         save = TRUE, thin = 1L) {
  set.seed(2718)
  defaults <- list(
    Mstar = 1L, crit_eig = 1, prmean = 0, a_1 = 3, b_1 = 0.3,
    Bsigma = 1, a0 = 25, b0 = 1.5, bmu = 0, Bmu = 100^2,
    lambda1 = 0.1, lambda2 = 0.2, lambda3 = 0.1, lambda4 = 100,
    tau0 = 0.1, tau1 = 3, kappa0 = 0.1, kappa1 = 7,
    p_i = 0.5, q_ij = 0.5, d_lambda = 0.01, e_lambda = 0.01,
    tau_theta = 0.7, sample_tau = TRUE, tau_log = TRUE
  )
  defaults[names(hyperpara)] <- hyperpara
  store <- list(shrink_MN = FALSE, shrink_SSVS = FALSE,
                shrink_NG = FALSE, shrink_HS = FALSE, vola_pars = FALSE)
  store[[paste0("shrink_", prior)]] <- save
  args <- list(
    Yraw = matrix(rnorm(160), 80, 2),
    Wraw = matrix(rnorm(80), 80, 1), Exraw = matrix(0, 1, 1),
    lags = c(1L, 1L), draws = 60L, burnin = 20L, thin = thin,
    cons = TRUE, trend = FALSE, sv = FALSE,
    prior = unname(c(MN = 1L, NG = 3L, HS = 4L)[prior]),
    setting_store = store
  )
  if (implementation == "C++") {
    sampler <- getFromNamespace("BVAR_linear", "BGVAR")
    args$hyperparam <- defaults
  } else {
    sampler <- getFromNamespace(".BVAR_linear_R", "BGVAR")
    args$hyperpara <- defaults
    args$verbose <- FALSE
  }
  do.call(sampler, args)
}

expect_positive_draws <- function(x) {
  expect_true(length(x) > 0)
  expect_true(all(is.finite(x)))
  expect_true(all(x > 0))
}

expect_prior_storage <- function(implementation, prior) {
  saved <- prior_sample(implementation, prior)
  unsaved <- prior_sample(implementation, prior, save = FALSE)
  # Saving diagnostics must not alter posterior simulation.
  expect_equal(saved$A_store, unsaved$A_store)
  expect_equal(saved$L_store, unsaved$L_store)
  expect_true(all(vapply(unsaved[[prior]], length, integer(1)) == 0L))
  expect_equal(dim(saved$A_store), c(5L, 2L, 60L))
  expect_true(all(is.finite(saved$A_store)))
  thinned <- prior_sample(implementation, prior, thin = 3L)
  expect_equal(dim(thinned$A_store), c(5L, 2L, 20L))
  # Both samplers retain the first post-burn-in draw, then every third draw.
  # Compare within each implementation; their random sequences need not agree.
  idx <- seq(1L, 60L, 3L)
  expect_equal(thinned$A_store, saved$A_store[, , idx, drop = FALSE])
  for (name in names(saved[[prior]])) {
    full <- saved[[prior]][[name]]
    selected <- if (length(dim(full)) == 3L) full[, , idx, drop = FALSE] else
      full[, idx, drop = FALSE]
    expect_equal(thinned[[prior]][[name]], selected)
  }
}
