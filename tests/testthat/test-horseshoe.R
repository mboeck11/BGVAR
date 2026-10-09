for (implementation in c("C++", "R")) {
  test_that(paste(implementation, "Horseshoe samples local, global and auxiliary scales"), {
    fit <- prior_sample(implementation, "HS")
    # Two equations: four domestic lag coefficients, four foreign coefficients
    # (current and lagged), and one lower-triangular covariance coefficient.
    for (name in names(fit$HS)) {
      expected_rows <- if (grepl("^(lambda|nu)_A_", name)) 4L else 1L
      expect_equal(dim(fit$HS[[name]]), c(expected_rows, 60L))
      expect_positive_draws(fit$HS[[name]])
      expect_gt(length(unique(as.numeric(fit$HS[[name]]))), 1L)
    }
    endo <- fit$HS$lambda_A_endo_store *
      matrix(fit$HS$tau_A_endo_store, 4, 60, byrow = TRUE)
    exo <- fit$HS$lambda_A_exo_store *
      matrix(fit$HS$tau_A_exo_store, 4, 60, byrow = TRUE)
    expect_positive_draws(endo)
    expect_positive_draws(exo)
    expect_false(isTRUE(all.equal(endo, exo)))
  })

  test_that(paste(implementation, "Horseshoe does not use other priors' hyperparameters"), {
    baseline <- prior_sample(implementation, "HS")
    changed <- prior_sample(implementation, "HS",
                            hyperpara = list(p_i = 0.2, q_ij = 0.8,
                                             d_lambda = 2, e_lambda = 3))
    expect_equal(baseline$A_store, changed$A_store)
    expect_equal(baseline$HS, changed$HS)
  })

  test_that(paste(implementation, "Horseshoe storage and thinning preserve the chain"), {
    expect_prior_storage(implementation, "HS")
  })
}
