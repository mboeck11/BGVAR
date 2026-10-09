# Run separately from testthat: machine-dependent timings are measurements,
# not pass/fail assertions. Example from the package root:
# Rscript scripts/benchmark-runtime.R
# For release timings, install to a temporary library first:
# dir.create("/tmp/bgvar-runtime-library", showWarnings=FALSE)
# R CMD INSTALL --preclean --library=/tmp/bgvar-runtime-library .
# BGVAR_BENCH_LIBRARY=/tmp/bgvar-runtime-library Rscript scripts/benchmark-runtime.R
# When sourcing, call benchmark_bgvar_runtime(draws=1000L, burnin=500L).
benchmark_bgvar_runtime <- function(repetitions=3L, draws=500L, burnin=250L) {
  bench_library <- Sys.getenv("BGVAR_BENCH_LIBRARY")
  if(nzchar(bench_library)) .libPaths(c(bench_library, .libPaths()))
  if(!nzchar(bench_library) && file.exists("DESCRIPTION")) {
    if(!requireNamespace("pkgload", quietly=TRUE)) stop("Install pkgload or specify BGVAR_BENCH_LIBRARY.")
    pkgload::load_all(quiet=TRUE)
    build <- "workspace build (may use debug compiler flags)"
  } else {
    library(BGVAR)
    build <- "installed package; use an optimized local install for release timings"
  }
  ns <- asNamespace("BGVAR")
  cpp <- get("BVAR_linear", ns)
  public <- get("bgvar", ns)
  # Fail immediately on any attempt to use the R fallback. The counter is
  # checked afterward too, because the public wrapper catches sampler errors.
  guard <- new.env(parent=emptyenv())
  guard$calls <- 0L
  tracer <- substitute({ G$calls <- G$calls+1L; stop("C++ failed: R fallback invalidates benchmark.") }, list(G=guard))
  trace(".BVAR_linear_R", where=ns, tracer=tracer, print=FALSE)
  on.exit(untrace(".BVAR_linear_R", where=ns), add=TRUE)
  timing <- function(fun) {
    times <- numeric(repetitions)
    size <- NA_real_
    for(i in seq_len(repetitions)) {
      gc()
      set.seed(410)
      elapsed <- system.time({
        output <- capture.output(result <- fun())
      })[["elapsed"]]
      if(guard$calls > 0L) stop("R fallback occurred; no timing is valid.")
      times[i] <- elapsed
      size <- as.numeric(object.size(result))/2^20
      rm(result, output)
    }
    c(median_seconds=median(times), min_seconds=min(times),
      max_seconds=max(times), returned_MiB=size)
  }
  scenarios <- data.frame(countries=c(4L,8L), variables=c(4L,6L), observations=c(150L,200L))
  results <- list()
  for(scenario in seq_len(nrow(scenarios))) {
    config <- scenarios[scenario,]
    N <- config$countries; M <- config$variables; T <- config$observations
    set.seed(1402)
    names <- sprintf("%02d",seq_len(N))
    Data <- setNames(lapply(seq_len(N), function(i)
      matrix(rnorm(T*M),T,M,dimnames=list(NULL,paste0("v",seq_len(M))))),names)
    W <- matrix(1/(N-1),N,N,dimnames=list(names,names)); diag(W) <- 0
    foreign <- Reduce(`+`, Map(function(x,w) x*w,Data,W[1,]))
    hyper <- list(Mstar=M, crit_eig=1, prmean=0, a_1=3, b_1=.3,
      Bsigma=1,a0=25,b0=1.5,bmu=0,Bmu=100^2,
      lambda1=.1,lambda2=.2,lambda3=.1,lambda4=100,
      tau0=.1,tau1=3,kappa0=.1,kappa1=7,p_i=.5,q_ij=.5,
      d_lambda=.01,e_lambda=.01,tau_theta=.7,sample_tau=FALSE,tau_log=FALSE)
    store <- list(shrink_MN=FALSE,shrink_SSVS=FALSE,shrink_NG=FALSE,shrink_HS=FALSE,vola_pars=FALSE)
    for(SV in c(FALSE,TRUE)) {
      direct <- function() cpp(Data[[1]], foreign, matrix(NA_real_),c(1L,1L),
        draws,burnin,1L,TRUE,FALSE,SV,1L,hyper,store)
      complete <- function() public(Data,W,draws=draws,burnin=burnin,thin=1L,
        prior="MN",SV=SV,eigen=FALSE,verbose=FALSE,
        expert=list(use_R=FALSE,save.country.store=FALSE,save.shrink.store=FALSE))
      # Warm up dispatch and native code outside measured runs.
      set.seed(410); invisible(direct())
      for(target in c("C++ country sampler","Complete bgvar")) {
        measured <- timing(if(target=="C++ country sampler") direct else complete)
        row <- data.frame(target=target,countries=if(target=="Complete bgvar") N else 1L,
          variables=M,observations=T,SV=SV,draws=draws,burnin=burnin,
          repetitions=repetitions,as.list(measured),check.names=FALSE)
        results[[length(results)+1L]] <- row
        print(row,row.names=FALSE)
      }
    }
  }
  results <- do.call(rbind,results)
  write.csv(results,"scripts/benchmark-runtime-results.csv",row.names=FALSE)
  metadata <- c(paste("Build:",build),paste("R:",R.version.string),
    paste("Platform:",R.version$platform),paste("Package:",utils::packageVersion("BGVAR")),
    paste("Library:",getNamespaceInfo(ns,"path")),
    "Synthetic data; Minnesota prior; one domestic/foreign lag; sequential country estimation.",
    "Compilation, data creation, warm-up and explicit garbage collection excluded from timing.",
    "Complete bgvar includes R preprocessing, C++ estimation, summaries and global stacking.",
    "returned_MiB is retained object size, not peak process memory.")
  writeLines(metadata,"scripts/benchmark-runtime-environment.txt")
  invisible(results)
}
if(sys.nframe()==0L) benchmark_bgvar_runtime()
