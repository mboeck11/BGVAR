# Replay the fixture's RNG state immediately before estimation, so trimming
# compares subsets of the very same posterior chain rather than different fits.
stability_fit <- function(plag, eigen) {
  baseline <- postprocess_model(plag = plag)
  set.seed(821)
  invisible(rnorm(320))
  invisible(capture.output(fit <- bgvar(
    Data = baseline$args$Data, W = baseline$args$W,
    draws = 30, burnin = 20, prior = "MN", plag = plag,
    SV = FALSE, eigen = eigen, verbose = FALSE
  )))
  fit
}

for (plag in c(1L, 2L)) {
  test_that(paste("lag order", plag, "stability is measured on the full companion matrix"), {
    fit <- postprocess_model(plag = plag)
    F <- fit$stacked.results$F_large
    radius <- vapply(seq_len(dim(F)[4]), function(draw) {
      top <- do.call(cbind, lapply(seq_len(plag), function(lag) F[, , lag, draw]))
      companion <- if (plag == 1L) top else
        rbind(top, cbind(diag(4 * (plag - 1L)), matrix(0, 4 * (plag - 1L), 4)))
      max(Mod(eigen(companion, only.values = TRUE)$values))
    }, numeric(1))
    expect_equal(as.numeric(fit$stacked.results$F.eigen), radius, tolerance = 1e-8)
    expect_equal(fit$args$thindraws, 30)
  })

  test_that(paste("lag order", plag, "trimming keeps all posterior arrays aligned"), {
    baseline <- postprocess_model(plag = plag)
    radius <- as.numeric(baseline$stacked.results$F.eigen)
    sorted <- sort(radius)
    cutoff <- mean(sorted[15:16])
    keep <- which(radius < cutoff)
    expect_equal(length(keep), 15L)
    trimmed <- stability_fit(plag, cutoff)
    expect_equal(trimmed$args$thindraws, length(keep))
    expect_equal(trimmed$stacked.results$F.eigen,
                 baseline$stacked.results$F.eigen[keep])
    expect_equal(trimmed$stacked.results$F_large,
                 baseline$stacked.results$F_large[, , , keep, drop = FALSE])
    for (name in c("A_large", "S_large", "Ginv_large")) {
      expect_equal(trimmed$stacked.results[[name]],
                   baseline$stacked.results[[name]][, , keep, drop = FALSE])
    }
    expect_true(all(trimmed$stacked.results$F.eigen < cutoff))
  })

  test_that(paste("lag order", plag, "trimming enforces the ten-draw minimum"), {
    radius <- sort(as.numeric(postprocess_model(plag = plag)$stacked.results$F.eigen))
    ten <- stability_fit(plag, mean(radius[10:11]))
    expect_equal(ten$args$thindraws, 10L)
    invisible(capture.output(expect_error(
      stability_fit(plag, mean(radius[9:10])), "Less than 10 stable draws"
    )))
    # Reject an estimate with no surviving draws before postprocessing is used.
    invisible(capture.output(expect_error(
      stability_fit(plag, min(radius) / 2), "Less than 10 stable draws"
    )))
  })
}

test_that("default trimming is equivalent to the explicit 1.05 threshold", {
  implicit <- stability_fit(1L, TRUE)
  explicit <- stability_fit(1L, 1.05)
  expect_equal(implicit$stacked.results, explicit$stacked.results)
  expect_true(all(implicit$stacked.results$F.eigen < 1.05))
})
