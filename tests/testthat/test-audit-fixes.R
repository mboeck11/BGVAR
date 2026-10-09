test_that("covariance quantiles transform each posterior draw", {
  model <- postprocess_model()
  S <- model$stacked.results$S_large
  G <- model$stacked.results$Ginv_large
  expected <- S
  for (i in seq_len(dim(S)[3])) expected[,,i] <- G[,,i] %*% S[,,i] %*% t(G[,,i])
  expect_equal(vcov(model), apply(expected, c(1,2), quantile, .5))
  expect_equal(vcov(model, quantile=c(.1,.9)), apply(expected, c(1,2), quantile, c(.1,.9)))
})

test_that("likelihood caching retains all quantiles", {
  model <- postprocess_model(trend=TRUE)
  first <- logLik(model)
  expect_equal(logLik(model), first)
  expect_equal(as.numeric(logLik(model, quantile=.1)),
               unname(quantile(model$args$logLik_draws, .1)))
})

test_that("FEVD retains supplied rotations and HD reports identification", {
  shocks <- postprocess_irf()
  rotation <- diag(4)
  rotation[1:2,1:2] <- matrix(c(0,1,-1,0),2)
  expect_equal(fevd(shocks, rotation.matrix=rotation, verbose=FALSE)$rotation.matrix, rotation)
  expect_error(fevd(shocks, rotation.matrix=matrix(1,4,4), verbose=FALSE), "orthogonal")
  expect_error(fevd(shocks, var.slct=character(), verbose=FALSE), "variables")
  signed <- postprocess_sign_irf()
  expect_output(print(hd(signed, verbose=FALSE)), "Sign-restrictions")
  expect_output(print(hd(shocks, verbose=FALSE)), "Cholesky")
})

test_that("matrix conversion matches country components exactly", {
  x <- matrix(seq_len(24),6,4,dimnames=list(NULL,c("AA.y","AA.p","BB.AA","BB.p")))
  result <- matrix_to_list(x)
  expect_equal(unname(result$AA), unname(x[,1:2]))
  expect_equal(unname(result$BB), unname(x[,3:4]))
})

test_that("shock helpers validate all signs and apply horizon defaults", {
  expect_error(get_shockinfo("invalid"), "arg")
  expect_error(get_shockinfo(nr_rows=0), "positive integer")
  expect_error(add_shockinfo(shock="AA.y", restriction=c("AA.y","AA.p"),
                            sign=c(">","invalid"), horizon=1), "sign")
  info <- suppressWarnings(add_shockinfo(shock="AA.y", restriction="AA.y", sign=">"))
  expect_equal(info$horizon, 1)
})

test_that("plotting supports unequal lags and a single forecast quantile", {
  grDevices::pdf(tempfile(fileext=".pdf"))
  on.exit(grDevices::dev.off())
  expect_error(plot(postprocess_model(plag=c(1,2)), resp="AA.y"), NA)
  prediction <- predict(postprocess_model(), n.ahead=1, quantiles=.5, verbose=FALSE)
  expect_error(plot(prediction, resp="AA.y", quantiles=.5, cut=200), NA)
  expect_error(predict(postprocess_model(), n.ahead=2, constr=rep(NA_real_,8), verbose=FALSE), "dimensions")
})

test_that("parallel core presets resolve to valid counts", {
  normalize <- getFromNamespace(".normalize_cores", "BGVAR")
  expect_null(normalize(NULL))
  expect_gte(normalize("half"), 1)
  expect_gte(normalize("all"), normalize("half"))
  for (invalid in list(0, -1, NA_real_, 1.5, c(1,2), "unknown"))
    expect_error(normalize(invalid), "cores")
})

test_that("Excel sheet exclusions accept integer and double indices", {
  file <- tempfile(fileext=".xlsx")
  file.create(file)
  on.exit(unlink(file))
  testthat::local_mocked_bindings(
    excel_sheets=function(...) c("AA", "BB"),
    excel_format=function(...) "xlsx",
    read_xlsx=function(path, sheet, ...) data.frame(y=1:3, p=4:6),
    .package="BGVAR")
  expect_named(excel_to_list(file, first_column_as_time=FALSE, skipsheet=1), "BB")
  expect_named(excel_to_list(file, first_column_as_time=FALSE, skipsheet=1L), "BB")
  expect_named(excel_to_list(file, first_column_as_time=FALSE, skipsheet="AA"), "BB")
  expect_error(excel_to_list(file, skipsheet=0), "indices")
  expect_error(excel_to_list(file, skipsheet=3), "indices")
})
