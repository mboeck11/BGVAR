for (implementation in c("C++", "R")) {
  test_that(paste(implementation, "Normal-Gamma samples positive local and global scales"), {
    fit <- prior_sample(implementation, "NG")
    expect_equal(dim(fit$NG$theta_store), c(5L, 2L, 60L))
    expect_positive_draws(fit$NG$theta_store)
    expect_equal(dim(fit$NG$lambda2_store), c(2L, 3L, 60L))
    # Endogenous lag zero and unused covariance rows are structural placeholders.
    expect_positive_draws(fit$NG$lambda2_store[2, 1, ])
    expect_positive_draws(fit$NG$lambda2_store[, 2, ])
    expect_positive_draws(fit$NG$lambda2_store[1, 3, ])
    expect_positive_draws(fit$NG$tau_store[2, 1, ])
    expect_positive_draws(fit$NG$tau_store[, 2, ])
    expect_positive_draws(fit$NG$tau_store[1, 3, ])
  })

  test_that(paste(implementation, "Normal-Gamma honors fixed shape parameters"), {
    for (shape in c(0.4, 1.2)) {
      fit <- prior_sample(implementation, "NG",
                          hyperpara = list(tau_theta = shape, sample_tau = FALSE))
      expect_equal(as.numeric(fit$NG$tau_store[2, 1, ]), rep(shape, 60))
      expect_equal(as.numeric(fit$NG$tau_store[, 2, ]), rep(shape, 120))
      expect_equal(as.numeric(fit$NG$tau_store[1, 3, ]), rep(shape, 60))
    }
    sampled <- prior_sample(implementation, "NG")
    expect_gt(length(unique(sampled$NG$tau_store[2, 1, ])), 1L)
  })

  test_that(paste(implementation, "Normal-Gamma storage and thinning preserve the chain"), {
    expect_prior_storage(implementation, "NG")
  })
}
