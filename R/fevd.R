#' @export
"fevd" <- function(x, rotation.matrix=NULL, var.slct=NULL, verbose=TRUE){
  UseMethod("fevd", x)
}

#' @name fevd
#' @title Forecast Error Variance Decomposition
#' @description This function calculates the forecast error variance decomposition (FEVDs) for Cholesky and sign-identified shocks.
#' @usage fevd(x, rotation.matrix=NULL, var.slct=NULL, verbose=TRUE)
#' @details Since the calculations are very time consuming, the FEVDs are based on the posterior median only (as opposed to calculating FEVDs for each MCMC sweep). In case the underlying shock has been identified via sign restrictions, the rotation matrix corresponds to the one that fulfills the sign restrictions at the posterior median of the estimated coefficients. More precisely, the algorithm searches for 50 rotation matrices that fulfill the sign restrictions at the \emph{posterior median} of the coefficients and then singles out the rotation matrix that minimizes the distance to the median of the impulse responses as suggested in Fry and Pagan (2011).
#' @param x an object of class \code{bgvar.irf}.
#' @param rotation.matrix If \code{NULL} and the \code{x} has been fitted via sign restrictions, the rotation matrix is used that minimizes the distance to the median impulse responses at the posterior median.
#' @param var.slct character vector that contains the variables for which forecast error variance decomposition should be performed. If \code{NULL} the FEVD is computed for the whole system, which is very time consuming.
#' @param verbose If set to \code{FALSE} it suppresses printing messages to the console.
#' @return Returns a list with two elements \describe{
#' \item{\code{FEVD}}{  an array of size (K times horizon times N), where K are all variables in the system, horizon is the specified impulse response horizon and N is the size of the decomposed structural variables (if \code{var.slct=NULL} then K=N).}
#' \item{\code{xglobal}}{ used data of the model.}
#' }
#' @author Maximilian Boeck, Martin Feldkircher, Florian Huber
#' @seealso \code{\link{bgvar}}, \code{\link{irf}}
#' @examples
#' \donttest{
#' set.seed(123)
#' library(BGVAR)
#' data(testdata)
#' model.eer<-bgvar(Data=testdata,W=W.test,prior="MN",
#'                  draws=50,burnin=50,plag=1,eigen=TRUE)
#'                       
#' # US monetary policy shock
#' shockinfo <- get_shockinfo("chol")
#' shockinfo$shock <- "US.stir"; shockinfo$scale <- -100
#' irf.chol.us.mp<-irf(model.eer,n.ahead=48,shockinfo=shockinfo)
#' 
#' # calculates FEVD for variables US.Dp and EA.y
#' fevd.us.mp=fevd(irf.chol.us.mp,var.slct=c("US.Dp","EA.y"))
#' 
#' # US monetary policy shock with sign restrictions
#' shockinfo <- get_shockinfo("sign")
#' shockinfo <- add_shockinfo(shockinfo, shock="US.stir", 
#'                            restriction=c("US.y","US.Dp"), 
#'                            sign=c("<","<"), horizon=c(1,1), 1, 100)
#' irf.sign.us.mp<-irf(model.eer,n.ahead=24,shockinfo=shockinfo)
#' 
#' # calculates FEVD for variables US.Dp and EA.y
#' fevd.us.mp=fevd(irf.sign.us.mp,var.slct=c("US.Dp","EA.y"))
#' }
#' @export
fevd.bgvar.irf <- function(x, rotation.matrix=NULL, var.slct=NULL, verbose=TRUE){
  start.fevd <- Sys.time()
  if(verbose) cat("Start computing forecast error variance decomposition of Bayesian Global Vector Autoregression.\n\n")
  #------------------------------ get stuff -------------------------------------------------------#
  xglobal     = x$model.obj$xglobal
  lags        = x$model.obj$lags
  pmax        = max(lags)
  ident       = x$ident
  Traw        = nrow(xglobal)
  bigK        = ncol(xglobal)
  xdat        = xglobal[(pmax+1):Traw,,drop=FALSE]
  bigT        = nrow(x)
  A           = x$struc.obj$A
  Fmat        = x$struc.obj$Fmat
  Ginv        = x$struc.obj$Ginv
  Smat        = x$struc.obj$S
  Rmed        = x$struc.obj$Rmed
  Sigma_u     = Ginv%*%Smat%*%t(Ginv)
  horizon     = dim(x$posterior)[2]
  varNames    = colnames(xglobal)
  shock       = x$shock
  sign.constr = x$sign.constr
  cN          = unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x)x[1]))
  N           = length(cN)
  if(ident=="sign"){
    if(verbose) cat("Identification scheme: Sign-restrictions provided.\n")
    shock.cN <- strsplit(unique(x$shockinfo$shock),".",fixed=TRUE)[[1]][1]
  }else if(ident=="chol"){
    if(verbose) cat("Identification scheme: Short-run restrictions via Cholesky decomposition.\n")
    shock.cN <- unique(x$shockinfo$shock)
  }
  #-------------------- some checks ------------------------------------------------------------------#
  if(!ident%in%c("sign","chol")){
    stop("FEVD implemented for shocks identified via cholesky ordering or sign restrictions only.")
  }
  if(!is.null(var.slct)){
    if(!all(var.slct%in%varNames)){
      stop("One of the variables you want to decompose is not contained in the system. Please re-specify!")
    }
  }
  if((is.null(Rmed) || any(is.na(Rmed))) && ident == "sign"){
    stop("No rotation matrix available. Please supply rotation matrix or re-estimate IRFs with sign-restrictions.")
  }else if(ident=="sign" && !is.null(Rmed)){
    rotation.matrix = Rmed
  }else{
    rotation.matrix = diag(bigK)
  }
  if(is.null(var.slct)){
    if(verbose) cat("FEVD computed for all variables.\n\n")
    var.slct<-varNames
  }else{
    var.print <- var.slct[1]
    if(length(var.slct)>1) for(kk in 2:length(var.slct)) var.print <- paste(var.print,", ",var.slct[kk],sep="")
    if(verbose) cat(paste("FEVD computed for the following variables: ",var.print,".\n",sep=""))
  }
  #----------------------------------------------------------------------------------------------------#
  P0G <- diag(bigK); colnames(P0G) <- rownames(P0G) <- varNames
  gcov <- Smat
  for(cc in 1:N){
    if(cN[cc] %in% shock.cN){
      idx           <- grep(cN[cc],varNames)
      P0G[idx,idx]  <- t(chol(gcov[idx,idx,drop=FALSE])) # calculate local cholesky factor of gcov
      gcov[idx,idx] <- diag(length(idx)) #set vcv matrix to identity for coutnry where shock occurs
    }
  }
  # create dynamic multiplier
  PHIx <- array(0,c(bigK,bigK,pmax+horizon+1)); dimnames(PHIx)[[1]] <- dimnames(PHIx)[[2]] <- varNames
  PHIx[,,pmax+1]  <-  diag(bigK)
  for (ihor in (pmax+2):(pmax+horizon+1)){
    acc = matrix(0,bigK,bigK)
    for (pp in 1:pmax){
      acc  <-  acc + Fmat[,,pp]%*%PHIx[,,ihor-pp]
    }
    PHIx[,,ihor]  <-  acc
  }
  PHI  <-  PHIx[,,(pmax+1):(pmax+horizon+1)]
  #----------------------------------------------------------------------------------------------------#
  if(verbose) cat("Start computing FEVDs...\n")
  vslct      <- diag(bigK)
  invGSigmau <- Ginv%*%gcov # use gcov not SIGMA_u from irf.obj
  invGSinvG  <- invGSigmau%*%t(Ginv)
  scale      <-  1/diag(gcov)
  
  FEVDres  <-  array(0,dim=c(bigK,length(var.slct),horizon+1))
  dimnames(FEVDres) <- list(varNames,paste("Decomp. of",var.slct),0:horizon)
  
  for(zz in 1:length(var.slct)){
    eslct <-matrix(0,bigK,1);rownames(eslct) <- varNames
    eslct[var.slct[zz],1] <- 1
    
    num  <-  matrix(0,bigK,horizon+1)
    den  <-  matrix(0,bigK,horizon+1)
    
    N <- 1
    while (N<=horizon+1){
      for (l in 1:N){
        acc1  <-  t((t(eslct)%*%rotation.matrix%*%PHI[,,l]%*%invGSigmau%*%vslct)^2)
        num[,N]  <-  num[,N] + acc1
        acc2  <-  (t(eslct)%*%rotation.matrix%*%PHI[,,l]%*%invGSinvG%*%t(rotation.matrix%*%PHI[,,l])%*%eslct)
        den[,N]  <-  den[,N] + matrix(1,bigK,1)*as.numeric(acc2)
      }
      FEVDres[,paste("Decomp. of",var.slct[zz]),N]  <-  (scale*num[,N])/den[,N]
      N <- N+1
    }
  }
  #------------------------------------------------------------------------------------------------------
  out <- structure(list(FEVD=FEVDres,
                        xglobal=xglobal,
                        rotation.matrix=rotation.matrix),
                   class="bgvar.fevd", type="fevd")
  if(verbose) cat(paste("\nSize of FEVD object: ", format(object.size(FEVDres),unit="MB")))
  end.fevd <- Sys.time()
  diff.fevd <- difftime(end.fevd,start.fevd,units="mins")
  mins.fevd <- round(diff.fevd,0); secs.fevd <- round((diff.fevd-floor(diff.fevd))*60,0)
  if(verbose) cat(paste("\nNeeded time for computation: ",mins.fevd," ",ifelse(mins.fevd==1,"min","mins")," ",secs.fevd, " ",ifelse(secs.fevd==1,"second.","seconds.\n"),sep=""))
  return(out)
}

#' @export
"gfevd" <- function(x, n.ahead=24, running=TRUE, applyfun=NULL, cores=NULL, verbose=TRUE){
  UseMethod("gfevd", x)
}

#' @name gfevd
#' @title Generalized Forecast Error Variance Decomposition
#' @description This function calculates a complete generalized forecast error variance decomposition (GFEVDs) based on generalized impulse response functions akin to Lanne-Nyberg (2016). The Lanne-Nyberg (2016) corrected GFEVD sum up to unity.
#' @method gfevd bgvar
#' @export
#' @usage gfevd(x, n.ahead=24, running=TRUE, applyfun=NULL, cores=NULL, verbose=TRUE)
#' @param x an object of class \code{bgvar}.
#' @param n.ahead the forecast horizon.
#' @param running Default is set to \code{TRUE} and implies that only a running mean over the posterior draws is calculated. A full analysis including posterior bounds is likely to cause memory issues.
#' @param applyfun Allows for user-specific apply function, which has to have the same interface than \code{lapply}. If \code{cores=NULL} then \code{lapply} is used, if set to a numeric either \code{parallel::parLapply()} is used on Windows platforms and \code{parallel::mclapply()} on non-Windows platforms.
#' @param cores Specifies the number of cores which should be used. Default is set to \code{NULL} and \code{applyfun} is used.
#' @param verbose If set to \code{FALSE} it suppresses printing messages to the console.
#' @return Returns a list with two elements \describe{
#' \item{\code{GFEVD}}{ a three or four-dimensional array, with the first dimension referring to the K time series that are decomposed into contributions of K time series (second dimension) for \code{n.ahead} forecast horizons. In case \code{running=TRUE} only the posterior mean else also its 16\% and 84\% credible intervals is contained in the fourth dimension.}
#' \item{\code{xglobal}}{ used data of the model.}
#' }
#' @author Maximilian Boeck, Martin Feldkircher
#' @references 
#' Lanne, M. and H. Nyberg (2016) \emph{Generalized Forecast Error Variance Decomposition for Linear and Nonlinear Multivariate Models.} Oxford Bulletin of Economics and Statistics, Vol. 78(4), pp. 595-603.
#' @seealso \code{\link{bgvar}}.
#' @examples 
#' \donttest{
#' library(BGVAR)
#' data(testdata)
#' model.eer<-bgvar(Data=testdata, W=W.test, draws=50, burnin=50,
#'                  plag=1, prior="SSVS", eigen=TRUE)
#'                       
#' GFEVD<-gfevd(model.eer, n.ahead=24)
#' }
#' @importFrom abind adrop
gfevd.bgvar<-function(x,n.ahead=24,running=TRUE,applyfun=NULL,cores=NULL,verbose=TRUE){
  start.gfevd <- Sys.time()
  if(verbose) cat("\nStart computing generalized forecast error variance decomposition of Bayesian Global Vector Autoregression.\n\n")
  #------------------------------ get stuff -------------------------------------------------------#
  lags        <- x$args$lags
  pmax        <- max(lags)
  xglobal     <- x$xglobal
  Traw        <- nrow(xglobal)
  bigK        <- ncol(xglobal)
  Kbig        <- pmax*bigK
  bigT        <- Traw-pmax
  A_large     <- x$stacked.results$A_large
  F_large     <- x$stacked.results$F_large
  S_large     <- x$stacked.results$S_large
  Ginv_large  <- x$stacked.results$Ginv_large
  F.eigen     <- x$stacked.results$F.eigen
  x           <- xglobal[(pmax+1):Traw,,drop=FALSE]
  thindraws   <- length(F.eigen) ### prior: draws
  varNames    <- colnames(xglobal)
  #------------------------------ prepare applyfun --------------------------------------------------------#
  if(is.null(applyfun)) {
    applyfun <- if(is.null(cores)) {
      lapply
    } else {
      if(.Platform$OS.type == "windows") {
        cl_cores <- parallel::makeCluster(cores)
        on.exit(parallel::stopCluster(cl_cores))
        function(X, FUN, ...) parallel::parLapply(cl = cl_cores, X, FUN, ...)
      } else {
        function(X, FUN, ...) parallel::mclapply(X, FUN, ..., mc.cores =
                                                   cores)
      }
    }
  }
  if(is.null(cores)) cores <- 1
  #-----------------------------------------------------------------------------------------------------#
  if(running){
    GFEVD_post <- array(0,dim=c(bigK,bigK,n.ahead)); dimnames(GFEVD_post)<-list(varNames, paste("Decomp. of",varNames),0:(n.ahead-1))
    if(verbose) cat(paste("Start computation on ", cores, " cores", " (",thindraws," stable draws in total).",sep=""),"\n")
    
    imp.obj <- applyfun(1:thindraws,function(irep){
      irfa <- .irf.girf.sims(invG=Ginv_large[,,irep],
                             lF=adrop(F_large[,,,irep,drop=FALSE],drop=4),
                             gcov=S_large[,,irep],
                             x,horizon=n.ahead)$impl
      GFEVD <- .mk_fevd.sims(irfa)
      return(list(GFEVD=GFEVD))
    })
    thindraws2 <- 0
    for(irep in 1:thindraws){
      if(!any(is.na(imp.obj[[irep]]$GFEVD))){
        GFEVD_post<-GFEVD_post+imp.obj[[irep]]$GFEVD
        thindraws2 <- thindraws2+1
      }
    }
    GFEVD_post<-GFEVD_post/thindraws2
  }else{ #-------------------------HERE DO FULL CALCULATION INCLUDING BOUNDS- VERY MEMORY INTENSIVE!!!-----
    GFEVD_draws<-array(NA,dim=c(thindraws,bigK,bigK,n.ahead))
    GFEVD_post <-array(NA,dim=c(bigK,bigK,n.ahead,3))
    dimnames(GFEVD_post)[[1]]<-dimnames(GFEVD_post)[[2]]<-varNames
    dimnames(GFEVD_post)[[3]]<-0:(n.ahead-1)
    dimnames(GFEVD_post)[[4]]<-c("low16","median","high84")
    if(verbose) cat(paste("Start computation on ", cores, " cores", " (",thindraws," stable draws in total).",sep=""),"\n")
    imp.obj <- applyfun(1:thindraws,function(irep){
      irfa <- .irf.girf.sims(invG=Ginv_large[,,irep],
                             lF=adrop(F_large[,,,irep,drop=FALSE],drop=1),
                             gcov=S_large[,,irep],
                             x,horizon=n.ahead)$impl
      GFEVD <- .mk_fevd.sims(irfa)
      return(list(GFEVD=GFEVD))
    })
    for(irep in 1:thindraws){
      GFEVD_draws[irep,,,] <- imp.obj[[irep]]$GFEVD
    }
    GFEVD_post[,,,1] <- apply(GFEVD_draws,c(2,3,4),quantile,.16,na.rm=TRUE)
    GFEVD_post[,,,2] <- apply(GFEVD_draws,c(2,3,4),median,na.rm=TRUE)
    GFEVD_post[,,,3] <- apply(GFEVD_draws,c(2,3,4),quantile,.16,na.rm=TRUE)
  }
  #----------------------------------------------------------------------------------------------------------------
  out <- structure(list(FEVD=GFEVD_post,
                        xglobal=xglobal,
                        R=NULL),
                   class="bgvar.fevd", type="gfevd")
  if(!running) out$GFEVD_store <- GFEVD_draws
  if(verbose) cat(paste("Size of IRF object:", format(object.size(out),unit="MB")))
  end.gfevd <- Sys.time()
  diff.gfevd <- difftime(end.gfevd,start.gfevd,units="mins")
  mins.gfevd <- round(diff.gfevd,0); secs.gfevd <- round((diff.gfevd-floor(diff.gfevd))*60,0)
  if(verbose) cat(paste("\nNeeded time for computation: ",mins.gfevd," ",ifelse(mins.gfevd==1,"min","mins")," ",secs.gfevd, " ",ifelse(secs.gfevd==1,"second.","seconds.\n"),sep=""))
  return(out)
}

#' @method print bgvar.fevd
#' @export
print.bgvar.fevd <- function(x, ...){
  cat("---------------------------------------------------------------------------------------")
  cat("\n")
  cat("Object contains forecast error variance decomposition of object estimated with 'bgvar':")
  cat("\n")
  cat(paste0("Size of FEVD containing forecast error variance decompositions: ",dim(x$FEVD)[[1]]," x ",dim(x$FEVD)[[2]]," x ",dim(x$FEVD)[[3]],"."))
  cat("\n")
  if(attributes(x)$type=="fevd"){
    cat("Identification scheme: ")
    if(is.null(x$rotation.matrix)){
      cat("Short-run restrictions via Cholesky decomposition.")
    }else{
      cat("Sign-restrictions.")
    }
  }else if(attributes(x)$type=="gfevd"){
    cat("Identification scheme: Generalized - no identification scheme employed.")
  }
  cat("\n")
  cat(paste0("Size ob object: ",format(object.size(x),unit="MB")))
  cat("\n")
  cat("---------------------------------------------------------------------------------------")
  
  return(invisible(x))
}
