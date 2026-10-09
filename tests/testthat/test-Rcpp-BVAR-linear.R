test_that("Rcpp BVAR_linear samples covariance rows with three and four equations", {
  # Call the compiled entry point directly. The public country wrapper catches
  # C++ errors and falls back to R, which would hide this regression.
  cpp <- getFromNamespace("BVAR_linear", "BGVAR")
  for(M in c(3L, 4L)) {
    set.seed(318)
    hyperparam <- list(
      Mstar=1L, crit_eig=1, prmean=0, a_1=3, b_1=.3,
      Bsigma=1, a0=25, b0=1.5, bmu=0, Bmu=100^2,
      lambda1=.1, lambda2=.2, lambda3=.1, lambda4=100,
      tau0=.1, tau1=3, kappa0=.1, kappa1=7,
      p_i=.5, q_ij=.5, d_lambda=.01, e_lambda=.01,
      tau_theta=.7, sample_tau=FALSE, tau_log=FALSE
    )
    args <- list(
      Yraw=matrix(rnorm(100*M), 100, M),
      Wraw=matrix(rnorm(100), 100, 1), Exraw=matrix(0,1,1),
      lags=c(1L,1L), draws=20L, burnin=10L, thin=1L,
      cons=TRUE, trend=FALSE, sv=FALSE,
      hyperparam=hyperparam,
      setting_store=list(shrink_MN=FALSE, shrink_SSVS=FALSE,
                         shrink_NG=FALSE, shrink_HS=FALSE, vola_pars=FALSE)
    )
    for(prior in 1:4) {
      args$prior <- prior
      args$setting_store$shrink_SSVS <- prior == 2L
      set.seed(921)
      fit <- do.call(cpp, args)
      expect_equal(dim(fit$L_store), c(M,M,20L))
      expect_equal(dim(fit$A_store), c(M+3L,M,20L))
      expect_true(all(is.finite(fit$A_store)))
      expect_true(all(is.finite(fit$L_store)))
      # Disabled NG diagnostics must remain empty even when SSVS is stored.
      expect_true(all(vapply(fit$NG, length, integer(1)) == 0L))
      for(draw in 1:20) {
        L <- fit$L_store[,,draw]
        expect_equal(diag(L), rep(1,M))
        expect_true(all(L[upper.tri(L)] == 0))
        expect_equal(fit$res_store[,,draw],
                     fit$Y-fit$X %*% fit$A_store[,,draw], tolerance=1e-12)
      }
      # The final covariance row contains a vector of length >1: the old
      # column-to-row assignment fails here before returning any draws.
      expect_true(all(apply(fit$L_store[M,seq_len(M-1L),], 1, sd) > 0))
      if(prior == 1L) {
        reference_args <- args
        reference_args$hyperparam <- NULL
        reference_args$hyperpara <- hyperparam
        reference_args$verbose <- FALSE
        set.seed(921)
        reference <- do.call(getFromNamespace(".BVAR_linear_R", "BGVAR"), reference_args)
        expect_equal(fit$L_store, reference$L_store, tolerance=1e-10)
        expect_equal(unname(fit$A_store), unname(reference$A_store), tolerance=1e-10)
      }
    }
  }
})
