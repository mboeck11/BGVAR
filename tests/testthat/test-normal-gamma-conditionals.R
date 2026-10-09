test_that("Normal-Gamma factor conditionals reproduce the complete joint density", {
  factors <- c(1.2, 1.8, .7)
  shapes <- c(.6, .9, 1.1)
  variances <- list(c(.2,.4,.5,.8), c(.3,.6,.8,1,1.2,1.4), c(.4,.7,.9,1.1,1.5))
  names(factors) <- names(shapes) <- paste0("lag.", 1:3)
  grid <- c(.3,.8,1.4,2.5)
  for(index in seq_along(factors)) {
    conditional <- .ng_factor_conditional(index, factors, shapes, variances, .01, .01)
    joint <- vapply(grid, function(value) {
      proposed <- factors
      proposed[index] <- value
      cumulative <- cumprod(proposed)
      sum(vapply(seq_along(variances), function(block)
        sum(dgamma(variances[[block]], shape=shapes[block],
                   rate=shapes[block]*cumulative[block]/2, log=TRUE)), numeric(1))) +
        sum(dgamma(proposed, shape=.01, rate=.01, log=TRUE))
    }, numeric(1))
    gamma <- dgamma(grid, shape=conditional["shape"], rate=conditional["rate"], log=TRUE)
    expect_equal(joint-joint[1], gamma-gamma[1], tolerance=1e-12)
    # The old value of the sampled factor cannot enter its conditional.
    changed <- factors
    changed[index] <- 20
    expect_equal(.ng_factor_conditional(index, changed, shapes, variances, .01, .01), conditional)
  }
  # Two blocks also cover contemporaneous and lagged foreign coefficients.
  expected <- c(shape=.01+4*.6+6*.9,
                rate=.01+.6/2*sum(variances[[1]])+.9/2*1.8*sum(variances[[2]]))
  expect_equal(.ng_factor_conditional(1, factors[1:2], shapes[1:2], variances[1:2], .01, .01), expected)
})

for(implementation in c("C++", "R")) {
  test_that(paste(implementation, "Normal-Gamma sampled factors follow full lag conditionals"), {
    for(lags in list(c(1L,1L), c(2L,1L), c(1L,2L))) {
      # M=2, Mstar=3 deliberately makes foreign blocks rectangular.
      fit <- prior_sample(implementation, "NG", lags=lags, Mstar=3L,
                          draws=1000L, hyperpara=list(sample_tau=FALSE, tau_theta=.7))
      expect_true(all(is.finite(fit$A_store)))
      for(foreign in c(FALSE, TRUE)) {
        block_count <- if(foreign) lags[2]+1L else lags[1]
        width <- if(foreign) 3L else 2L
        offset <- if(foreign) lags[1]*2L else 0L
        column <- if(foreign) 2L else 1L
        factor_rows <- if(foreign) seq_len(block_count) else seq.int(2L, block_count+1L)
        factors <- matrix(fit$NG$lambda2_store[factor_rows,column,], block_count,1000)
        for(index in seq_len(block_count)) {
          rates <- numeric(999)
          # Later blocks have not yet been updated in this sweep; earlier
          # factors already have their current-sweep values.
          for(draw in 2:1000) {
            rate <- .01
            for(block in seq.int(index, block_count)) {
              product <- 1
              for(factor in seq_len(block)) {
                if(factor != index)
                  product <- product*factors[factor,if(factor < index) draw else draw-1L]
              }
              rows <- offset + seq.int((block-1L)*width+1L, block*width)
              rate <- rate + .7/2*product*sum(fit$NG$theta_store[rows,,draw-1L])
            }
            rates[draw-1L] <- rate
          }
          shape <- .01 + (block_count-index+1L)*width*2*.7
          transformed <- rates*factors[index,-1L]
          # Rate*delta has Gamma(shape, rate=1). Fixed seeds and generous
          # six-standard-error bounds keep this reference check reproducible.
          expect_lt(abs(mean(transformed)-shape), 6*sqrt(shape/999))
          expect_lt(abs(mean(transformed^2)-shape*(shape+1)),
                    6*sqrt(2*shape*(shape+1)*(2*shape+3)/999))
        }
      }
    }
  })
}
