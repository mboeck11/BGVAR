# A shared deterministic fixture; fitting is deferred until a test needs it.
postprocess_model <- local({
  cache <- list()
  function(plag = 1L, trend = FALSE, hold.out = 0L) {
    key <- paste(plag, trend, hold.out, sep = "-")
    if (is.null(cache[[key]])) {
      set.seed(821)
      Data <- list(AA = matrix(rnorm(160), 80, 2),
                   BB = matrix(rnorm(160), 80, 2))
      Data <- lapply(Data, function(x) { colnames(x) <- c("y", "p"); x })
      W <- matrix(c(0, 1, 1, 0), 2, dimnames = list(names(Data), names(Data)))
      invisible(capture.output(cache[[key]] <<- bgvar(
        Data = Data, W = W, draws = 30, burnin = 20, prior = "MN",
        plag = plag, trend = trend, hold.out = hold.out,
        SV = FALSE, eigen = FALSE, verbose = FALSE
      )))
    }
    cache[[key]]
  }
})

postprocess_irf <- function(ident = "chol", scale = 1, use_R = FALSE,
                            save = TRUE) {
  info <- get_shockinfo(ident)
  info$shock <- "AA.y"
  info$scale <- scale
  irf(postprocess_model(), n.ahead = 4, shockinfo = info,
      quantiles = c(0.16, 0.5, 0.84),
      expert = list(use_R = use_R, save.store = save), verbose = FALSE)
}

postprocess_sign_irf <- function() {
  set.seed(517)
  info <- get_shockinfo("sign")
  for (country in c("AA", "BB")) {
    info <- add_shockinfo(info, shock = paste0(country, ".y"),
                          restriction = paste0(country, ".y"), sign = "<",
                          horizon = 1, prob = 1, scale = -1)
  }
  irf(postprocess_model(), n.ahead = 4, shockinfo = info,
      expert = list(MaxTries = 1000), verbose = FALSE)
}

# Independent design matrix for checking innovations and decomposition output.
postprocess_residuals <- function(shocks) {
  data <- shocks$model.obj$xglobal
  p <- max(shocks$model.obj$lags)
  rows <- seq.int(p + 1L, nrow(data))
  X <- do.call(cbind, lapply(seq_len(p), function(lag) data[rows - lag, , drop = FALSE]))
  X <- cbind(X, 1)
  if ("trend" %in% colnames(shocks$struc.obj$A)) X <- cbind(X, seq_along(rows))
  data[rows, , drop = FALSE] - X %*% t(shocks$struc.obj$A)
}
