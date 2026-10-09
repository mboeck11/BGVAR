test_that("Minnesota variances preserve every foreign row and deterministic boundary", {
  variance <- function(ndet, wexo = TRUE, lambda3 = .3, lambda4 = 100) {
    .get_V(k = 2 + if(wexo) 2 + ndet else ndet, M = 2,
           Mstar = 1, Mex = 0, plag = 1, plagstar = 1,
           lambda1 = .1, lambda2 = .2, lambda3 = lambda3,
           lambda4 = lambda4, sigma_sq = c(1, 2), sigma_wex = 1,
           wexo = wexo)
  }
  for (ndet in 0:2) {
    v <- variance(ndet)
    expect_equal(v[3:4, ], outer(c(.03^2, (.03/2)^2), c(1, 2)))
    proposed <- variance(ndet, lambda3 = .6)
    expect_equal(proposed[1:2, ], v[1:2, ])
    expect_equal(proposed[3:4, ], 4*v[3:4, ])
    if(ndet > 0) {
      expect_equal(v[seq.int(5, 4+ndet), , drop=FALSE],
                   matrix(rep(c(100, 200), each=ndet), ndet, 2))
      expect_equal(proposed[seq.int(5, 4+ndet), ], v[seq.int(5, 4+ndet), ])
    }
    noforeign <- variance(ndet, wexo=FALSE)
    expect_equal(noforeign[1:2, ], v[1:2, ])
  }
})

test_that("auxiliary AR variances use all requested lags and residual degrees of freedom", {
  set.seed(905)
  y <- as.numeric(arima.sim(list(ar=c(.5, -.3, .2)), n=120))
  for(lag in 1:3) {
    embedded <- embed(y, lag+1)
    design <- cbind(embedded[, -1, drop=FALSE], seq_len(nrow(embedded)))
    reference <- lm.fit(design, embedded[, 1])
    expected <- sum(reference$residuals^2)/(nrow(design)-reference$rank)
    expect_equal(.ar_residual_variance(y, lag), expected, tolerance=1e-12)
  }
})

test_that("Minnesota implementations agree with unequal lag orders", {
  for(lags in list(c(1L, 1L), c(2L, 1L), c(1L, 2L))) {
    cpp <- prior_sample("C++", "MN", lags=lags)
    r <- prior_sample("R", "MN", lags=lags)
    expect_equal(cpp$MN$lambda_store, r$MN$lambda_store, tolerance=1e-10)
    expect_equal(unname(cpp$A_store), unname(r$A_store), tolerance=1e-10)
    expect_equal(cpp$Sv_store, r$Sv_store, tolerance=1e-10)
    # This fixture has innovations around unit variance, so some log variances
    # must be negative. Storing raw variances makes this assertion fail.
    expect_true(any(cpp$Sv_store < 0))
  }
})

test_that("SSVS probabilities match Bayes' rule and survive density underflow", {
  for(value in c(-1, 0, .2, 1)) {
    spike <- dnorm(value, 0, .1)*.7
    slab <- dnorm(value, 0, 3)*.3
    expect_equal(.ssvs_spike_probability(value, 0, .1, 3, .3),
                 spike/(spike+slab))
  }
  for(value in c(-200, 200, -1e200, 1e200)) {
    expect_equal(.ssvs_spike_probability(value, 0, .1, 3, .3), 0)
    expect_equal(.ssvs_spike_probability(value, 0, .1, 3, 0), 1)
    expect_equal(.ssvs_spike_probability(value, 0, .1, 3, 1), 0)
    expect_equal(.ssvs_spike_probability(value, 0, 1, 1, .3), .7)
  }
})

for(implementation in c("C++", "R")) {
  test_that(paste(implementation, "Horseshoe conditionals respect auxiliary scales and prior means"), {
    fit <- prior_sample(implementation, "HS", hyperpara=list(prmean=1), draws=1000L)
    hs <- fit$HS
    for(block in c("endo", "exo")) {
      tau <- as.numeric(hs[[paste0("tau_A_", block, "_store")]])
      zeta <- as.numeric(hs[[paste0("zeta_A_", block, "_store")]])
      # Under the auxiliary IG(1, 1+1/tau) conditional, this is Exp(1).
      expect_lt(abs(mean((1+1/tau)/zeta)-1), .2)
      rows <- if(block == "endo") 1:2 else 3:4
      coefficients <- fit$A_store[rows, , , drop=FALSE]
      if(block == "endo") {
        coefficients[1, 1, ] <- coefficients[1, 1, ]-1
        coefficients[2, 2, ] <- coefficients[2, 2, ]-1
      }
      coefficients <- matrix(coefficients, 4, 1000)
      local <- hs[[paste0("lambda_A_", block, "_store")]]
      rate <- 1/zeta[-1000] + .5*colSums(coefficients[, -1]^2/local[, -1])
      # The transformed global draw is Gamma((4+1)/2, rate=1).
      expect_lt(abs(mean(rate/tau[-1])-2.5), .3)
    }
  })
}
