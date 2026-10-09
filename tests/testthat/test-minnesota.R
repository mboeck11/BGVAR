for (implementation in c("C++", "R")) {
  test_that(paste(implementation, "Minnesota samples positive, changing shrinkage parameters"), {
    fit <- prior_sample(implementation, "MN")
    expect_equal(dim(fit$MN$lambda_store), c(3L, 1L, 60L))
    expect_positive_draws(fit$MN$lambda_store)
    for (i in 1:3) {
      expect_gt(length(unique(fit$MN$lambda_store[i, 1, ])), 1L)
    }
  })

  test_that(paste(implementation, "Minnesota honors its own-lag prior mean"), {
    zero <- prior_sample(implementation, "MN", hyperpara = list(prmean = 0))
    unit <- prior_sample(implementation, "MN", hyperpara = list(prmean = 1))
    expect_false(isTRUE(all.equal(zero$A_store, unit$A_store)))
    # SSVS mixture probabilities must not influence Minnesota estimation.
    other <- prior_sample(implementation, "MN", hyperpara = list(p_i = 0.2, q_ij = 0.8))
    expect_equal(zero$A_store, other$A_store)
  })

  test_that(paste(implementation, "Minnesota storage and thinning preserve the chain"), {
    expect_prior_storage(implementation, "MN")
  })
}
