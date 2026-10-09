# Exercise the country samplers directly so a C++ failure cannot silently
# fall back to R and leave the compiled implementation untested.
ssvs_sample <- function(implementation, p_i, q_ij, equal_scales = FALSE,
                        draws = 100L, burnin = 20L, hyperpara_override = list()) {
  set.seed(2718)
  Yraw <- matrix(rnorm(160), 80, 2)
  Wraw <- matrix(rnorm(80), 80, 1)
  hyperpara <- list(
    Mstar = 1L, crit_eig = 1, prmean = 0, a_1 = 3, b_1 = 0.3,
    Bsigma = 1, a0 = 25, b0 = 1.5, bmu = 0, Bmu = 100^2,
    lambda1 = 0.1, lambda2 = 0.2, lambda3 = 0.1, lambda4 = 100,
    tau0 = 0.1, tau1 = 3, kappa0 = 0.1, kappa1 = 7,
    p_i = p_i, q_ij = q_ij, d_lambda = 0.01, e_lambda = 0.01,
    tau_theta = 0.7, sample_tau = TRUE, tau_log = TRUE
  )
  if (equal_scales) {
    hyperpara$tau0 <- hyperpara$tau1 <- 1
    hyperpara$kappa0 <- hyperpara$kappa1 <- 1
  }
  hyperpara[names(hyperpara_override)] <- hyperpara_override
  args <- list(
    Yraw = Yraw, Wraw = Wraw, Exraw = matrix(0, 1, 1),
    lags = c(1L, 1L), draws = draws, burnin = burnin, thin = 1L,
    cons = TRUE, trend = FALSE, sv = FALSE, prior = 2L,
    setting_store = list(shrink_MN = FALSE, shrink_SSVS = TRUE,
                         shrink_NG = FALSE, shrink_HS = FALSE,
                         vola_pars = FALSE)
  )
  if (implementation == "C++") {
    sampler <- getFromNamespace("BVAR_linear", "BGVAR")
    args$hyperparam <- hyperpara
  } else {
    sampler <- getFromNamespace(".BVAR_linear_R", "BGVAR")
    args$hyperpara <- hyperpara
    args$verbose <- FALSE
  }
  result <- do.call(sampler, args)
  list(gamma = result$SSVS$gamma_store,
       # With two equations, only (2, 1) is a sampled covariance indicator.
       omega = result$SSVS$omega_store[2, 1, ])
}

for (implementation in c("C++", "R")) {
  test_that(paste(implementation, "SSVS endpoints select the documented component"), {
    for (p in c(0, 1)) {
      result <- ssvs_sample(implementation, p_i = p, q_ij = 1 - p)
      expect_true(length(result$gamma) > 0)
      expect_true(length(result$omega) > 0)
      expect_true(all(is.finite(result$gamma)))
      expect_true(all(is.finite(result$omega)))
      expect_true(all(result$gamma == p))
      expect_true(all(result$omega == 1 - p))
    }
  })

  test_that(paste(implementation, "SSVS uses inclusion weights for asymmetric priors"), {
    # Equal component densities remove the data from the indicator probability.
    # Only in this controlled case must the inclusion rate equal the prior.
    for (p in c(0.3, 0.5, 0.7)) {
      result <- ssvs_sample(implementation, p_i = p, q_ij = 1 - p,
                            equal_scales = TRUE, draws = 2000L)
      expect_true(all(result$gamma %in% c(0, 1)))
      expect_true(all(result$omega %in% c(0, 1)))
      # More than six binomial standard errors for the smaller omega sample.
      expect_lt(abs(mean(result$gamma) - p), 0.08)
      expect_lt(abs(mean(result$omega) - (1 - p)), 0.08)
    }
  })
}

for (implementation in c("C++", "R")) {
  test_that(paste(implementation, "SSVS covariance probabilities survive underflow"), {
    # At the first sampled covariance coefficient both ordinary densities are
    # zero. Log odds still select the wider component with probability ~1.
    result <- ssvs_sample(implementation, p_i=.3, q_ij=.7,
                          draws=1L, burnin=0L,
                          hyperpara_override=list(kappa0=1e-8, kappa1=1e-7))
    expect_equal(as.numeric(result$omega), 1)
  })
}
