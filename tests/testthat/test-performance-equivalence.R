test_that("blocked covariance summaries reproduce the original full-array calculation", {
  set.seed(618)
  for(M in c(1L,3L)) for(T in c(1L,7L)) for(draws in c(1L,4L,5L)) {
    L <- array(0,c(M,M,draws))
    for(draw in seq_len(draws)) {
      matrix <- matrix(rnorm(M*M),M,M)
      matrix[upper.tri(matrix)] <- 0
      diag(matrix) <- 1
      L[,,draw] <- matrix
    }
    for(constant in c(FALSE,TRUE)) {
      log_variances <- array(rnorm(T*M*draws),c(T,M,draws))
      if(constant) for(t in seq_len(T)) log_variances[t,,] <- log_variances[1,,]
      names <- paste0("v",seq_len(M))
      original <- array(0,c(T,M,M,draws),dimnames=list(NULL,names,names,NULL))
      for(draw in seq_len(draws)) for(t in seq_len(T)) {
        triangular <- matrix(L[,,draw],M,M)
        original[t,,,draw] <- triangular %*%
          diag(as.numeric(exp(log_variances[t,,draw])),nrow=M) %*% t(triangular)
      }
      fast <- .covariance_summaries(L,log_variances,constant,names)
      expect_equal(fast$draw_medians, apply(original,c(2,3,4),median), tolerance=1e-12)
      expect_equal(fast$posterior_medians, apply(original,c(1,2,3),median), tolerance=1e-12)
    }
  }
})
