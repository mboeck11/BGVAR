#' @export
"fevd" <- function(x, ...){
  UseMethod("fevd", x)
}

#' @name fevd
#' @title Forecast Error Variance Decomposition
#' @description This function calculates the forecast error variance decomposition (FEVDs) for Cholesky and sign-identified shocks.
#' @usage fevd(x, R=NULL, var.slct=NULL, verbose=TRUE)
#' @details Since the calculations are very time consuming, the FEVDs are based on the posterior median only (as opposed to calculating FEVDs for each MCMC sweep). In case the underlying shock has been identified via sign restrictions, the rotation matrix corresponds to the one that fulfills the sign restrictions at the posterior median of the estimated coefficients. More precisely, the algorithm searches for 50 rotation matrices that fulfill the sign restrictions at the \emph{posterior median} of the coefficients and then singles out the rotation matrix that minimizes the distance to the median of the impulse responses as suggested in Fry and Pagan (2011).
#' @param x an object of class \code{bgvar.irf}.
#' @param R If \code{NULL} and the \code{x} has been fitted via sign restrictions, the rotation matrix is used that minimizes the distance to the median impulse responses at the posterior median.
#' @param var.slct character vector that contains the variables for which forecast error variance decomposition should be performed. If \code{NULL} the FEVD is computed for the whole system, which is very time consuming.
#' @param verbose If set to \code{FALSE} it suppresses printing messages to the console.
#' @return Returns a list with two elements \itemize{
#' \item{\code{FEVD}}{  an array of size (K times horizon times N), where K are all variables in the system, horizon is the specified impulse response horizon and N is the size of the decomposed structural variables (if \code{var.slct=NULL} then K=N).}
#' \item{\code{xglobal}}{ used data of the model.}
#' }
#' @author Maximilian Boeck, Martin Feldkircher, Florian Huber
#' @seealso \code{\link{IRF}}
#' @examples
#' \dontshow{
#' library(BGVAR)
#' data(eerData)
#' cN<-c("EA","US","UK")
#' eerData<-eerData[cN]
#' W.trade0012<-apply(W.trade0012[cN,cN],2,function(x)x/rowSums(W.trade0012[cN,cN]))
#' model.ssvs.eer<-bgvar(Data=eerData,W=W.trade0012,draws=50,burnin=50,plag=1,
#'                       prior="SSVS",thin=1,eigen=TRUE)
#'                       
#' # US monetary policy shock
#' shocks<-list();shocks$var="stir";shocks$cN<-"US";shocks$ident="chol";shocks$scal=-100
#' irf.chol.us.mp<-irf(model.ssvs.eer,shock=shocks,nhor=48)
#' 
#' # calculates FEVD for variables US.Dp and EA.y
#' fevd.us.mp=fevd(irf.chol.us.mp,var.slct=c("US.Dp","EA.y"))
#' 
#' # US monetary policy shock with sign restrictions
#' sign.constr<-list()
#' sign.constr$shock1$shock             <- c("US.stir")
#' sign.constr$shock1$restrictions$res1 <- c("US.y")
#' sign.constr$shock1$restrictions$res2 <- c("US.Dp")
#' sign.constr$shock1$sign              <- c(">","<","<")
#' sign.constr$shock1$rest.horz         <- c(1,1,1)
#' sign.constr$shock1$constr            <- c(1,1,1)
#' sign.constr$shock1$scal              <- +100 
#' sign.constr$MaxTries<-200
#' irf.sign.us.mp<-irf(model.ssvs.eer,sign.constr=sign.constr,nhor=24)
#' 
#' # calculates FEVD for variables US.Dp and EA.y
#' fevd.us.mp=fevd(irf.sign.us.mp,var.slct=c("US.Dp","EA.y"))
#' }
#' \donttest{
#' library(BGVAR)
#' data(eerData)
#' model.ssvs.eer<-bgvar(Data=eerData,W=W.trade0012,draws=100,burnin=100,plag=1,
#'                       prior="SSVS",thin=1,eigen=TRUE)
#'                       
#' # US monetary policy shock
#' shocks<-list();shocks$var="stir";shocks$cN<-"US";shocks$ident="chol";shocks$scal=-100
#' irf.chol.us.mp<-irf(model.ssvs.eer,shock=shocks,nhor=48)
#' 
#' # calculates FEVD for variables US.Dp and EA.y
#' fevd.us.mp=fevd(irf.chol.us.mp,var.slct=c("US.Dp","EA.y"))
#' 
#' # US monetary policy shock with sign restrictions
#' sign.constr<-list()
#' sign.constr$shock1$shock             <- c("US.stir")
#' sign.constr$shock1$restrictions$res1 <- c("US.y")
#' sign.constr$shock1$restrictions$res2 <- c("US.Dp")
#' sign.constr$shock1$sign              <- c(">","<","<")
#' sign.constr$shock1$rest.horz         <- c(1,1,1)
#' sign.constr$shock1$constr            <- c(1,1,1)
#' sign.constr$shock1$scal              <- +100 
#' sign.constr$MaxTries<-200
#' irf.sign.us.mp<-irf(model.ssvs.eer,sign.constr=sign.constr,nhor=24)
#' 
#' # calculates FEVD for variables US.Dp and EA.y
#' fevd.us.mp=fevd(irf.sign.us.mp,var.slct=c("US.Dp","EA.y"))
#' }
#' # NOT RUN - calculates FEVDs for all variables in the system, very time consuming
#' \dontrun{
#  fevd.us.mp=fevd(irf.chol.us.mp,var.slct=NULL)
#' }
#' @export
#' @rdname fevd
fevd.bgvar.irf <- function(x,R=NULL,var.slct=NULL,verbose=TRUE){
  start.fevd <- Sys.time()
  if(verbose) cat("\nStart computing forecast error variance decomposition of Bayesian Global Vector Autoregression.\n\n")
  #------------------------------ get stuff -------------------------------------------------------#
  xglobal <- x$model.obj$xglobal
  plag    <- x$model.obj$plag
  ident   <- x$ident
  Traw    <- nrow(xglobal)
  bigK    <- ncol(xglobal)
  xdat    <- xglobal[(plag+1):Traw,,drop=FALSE]
  bigT    <- nrow(x)
  A       <- x$struc.obj$A
  Fmat    <- x$struc.obj$Fmat
  Ginv    <- x$struc.obj$Ginv
  Smat    <- x$struc.obj$S
  Sigma_u <- Ginv%*%Smat%*%t(Ginv)
  horizon <- dim(x$posterior)[2]
  varNames<- colnames(xglobal)
  shock   <- x$shock
  sign.constr <- x$sign.constr
  cN <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x)x[1]))
  N  <- length(cN)
  if(ident=="sign"){
    if(verbose) cat("Identification scheme: Sign-restrictions provided.\n")
    shock.cN <- sapply(strsplit(unlist(lapply(sign.constr,function(l)l$shock)),".",fixed=TRUE),function(x)x[1])
  }else if(ident=="chol"){
    if(verbose) cat("Identification scheme: Short-run restrictions via Cholesky decomposition.\n")
    shock.cN <- shock$cN
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
  if(is.null(var.slct)){
    if(verbose) cat("FEVD computed for all variables.\n\n")
    var.slct<-varNames
  }else{
    var.print <- var.slct[1]
    if(length(var.slct)>1) for(kk in 2:length(var.slct)) var.print <- paste(var.print,", ",var.slct[kk],sep="")
    if(verbose) cat(paste("FEVD computed for the following variables: ",var.print,".\n",sep=""))
  }
  if(ident=="sign" && is.null(R)){
    R <- x$struc.obj$Rmed
  }else if(ident=="chol"){
    R<-diag(bigK)
  }
  rownames(R) <- colnames(R) <- varNames
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
  PHIx <- array(0,c(bigK,bigK,plag+horizon+1)); dimnames(PHIx)[[1]] <- dimnames(PHIx)[[2]] <- varNames
  PHIx[,,plag+1]  <-  diag(bigK)
  for (ihor in (plag+2):(plag+horizon+1)){
    acc = matrix(0,bigK,bigK)
    for (pp in 1:plag){
      acc  <-  acc + Fmat[,,pp]%*%PHIx[,,ihor-pp]
    }
    PHIx[,,ihor]  <-  acc
  }
  PHI  <-  PHIx[,,(plag+1):(plag+horizon+1)]
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
        acc1  <-  t((t(eslct)%*%R%*%PHI[,,l]%*%invGSigmau%*%vslct)^2)
        num[,N]  <-  num[,N] + acc1
        acc2  <-  (t(eslct)%*%R%*%PHI[,,l]%*%invGSinvG%*%t(R%*%PHI[,,l])%*%eslct)
        den[,N]  <-  den[,N] + matrix(1,bigK,1)*as.numeric(acc2)
      }
      FEVDres[,paste("Decomp. of",var.slct[zz]),N]  <-  (scale*num[,N])/den[,N]
      N <- N+1
    }
  }
  #------------------------------------------------------------------------------------------------------
  out <- structure(list(FEVD=FEVDres,
                        xglobal=xglobal),
                   class="bgvar.fevd", type="fevd")
  if(verbose) cat(paste("\nSize of FEVD object: ", format(object.size(FEVDres),unit="MB")))
  end.fevd <- Sys.time()
  diff.fevd <- difftime(end.fevd,start.fevd,units="mins")
  mins.fevd <- round(diff.fevd,0); secs.fevd <- round((diff.fevd-floor(diff.fevd))*60,0)
  if(verbose) cat(paste("\nNeeded time for computation: ",mins.fevd," ",ifelse(mins.fevd==1,"min","mins")," ",secs.fevd, " ",ifelse(secs.fevd==1,"second.","seconds.\n"),sep=""))
  return(out)
}

#' @export
"gfevd" <- function(x, ...){
  UseMethod("gfevd", x)
}

#' @name gfevd
#' @title Generalized Forecast Error Variance Decomposition
#' @description This function calculates a complete generalized forecast error variance decomposition (GFEVDs) based on generalized impulse response functions akin to Lanne-Nyberg (2016). The Lanne-Nyberg (2016) corrected GFEVD sum up to unity.
#' @method gfevd bgvar
#' @export
#' @usage gfevd(x, nhor=24, running=TRUE, multithread=FALSE, verbose=TRUE)
#' @param x an object of class \code{bgvar}.
#' @param nhor the forecast horizon.
#' @param running Default is set to \code{TRUE} and implies that only a running mean over the posterior draws is calculated. A full analysis including posterior bounds is likely to cause memory issues.
#' @param applyfun Allows for user-specific apply function, which has to have the same interface than \code{lapply}. If \code{cores=NULL} then \code{lapply} is used, if set to a numeric either \code{parallel::parLapply()} is used on Windows platforms and \code{parallel::mclapply()} on non-Windows platforms.
#' @param cores Specifies the number of cores which should be used. Default is set to \code{NULL} and \code{applyfun} is used.
#' @param verbose If set to \code{FALSE} it suppresses printing messages to the console.
#' @return Returns a list with two elements \itemize{
#' \item{\code{GFEVD}}{ a three or four-dimensional array, with the first dimension referring to the K time series that are decomposed into contributions of K time series (second dimension) for \code{nhor} forecast horizons. In case \code{running=TRUE} only the posterior mean else also its 16\% and 84\% credible intervals is contained in the fourth dimension.}
#' \item{\code{xglobal}}{ used data of the model.}
#' }
#' @author Maximilian Boeck, Martin Feldkircher
#' @references 
#' Lanne, M. and H. Nyberg (2016) \emph{Generalized Forecast Error Variance Decomposition for Linear and Nonlinear Multivariate Models.} Oxford Bulletin of Economics and Statistics, Vol. 78(4), pp. 595-603.
#' @seealso \code{\link{bgvar}}.
#' @examples 
#' \dontshow{
#' library(BGVAR)
#' data(eerData)
#' cN<-c("EA","US","UK")
#' eerData<-eerData[cN]
#' W.trade0012<-apply(W.trade0012[cN,cN],2,function(x)x/rowSums(W.trade0012[cN,cN]))
#' 
#' model.ssvs.eer<-bgvar(Data=eerData,W=W.trade0012,draws=100,burnin=100,plag=1,
#'                       prior="SSVS",thin=1,eigen=TRUE)
#'                       
#' GFEVD<-gfevd(model.ssvs.eer,nhor=24,running=TRUE)
#' }
#' \donttest{
#' library(BGVAR)
#' data(eerData)
#' model.ssvs.eer<-bgvar(Data=eerData,W=W.trade0012,draws=100,burnin=100,plag=1,
#'                       prior="SSVS",thin=1,eigen=TRUE)
#'                       
#' # Calculates running mean GFEVDs for all variables in the system 
#' GFEVD<-gfevd(model.ssvs.eer,nhor=24,running=TRUE)
#' }
#' @importFrom abind adrop
#' @importFrom parallel parLapply mclapply
gfevd.bgvar<-function(x,nhor=24,running=TRUE,applyfun=NULL,cores=NULL,verbose=TRUE){
  start.gfevd <- Sys.time()
  if(verbose) cat("\nStart computing generalized forecast error variance decomposition of Bayesian Global Vector Autoregression.\n\n")
  #------------------------------ get stuff -------------------------------------------------------#
  plag        <- x$args$plag
  xglobal     <- x$xglobal
  Traw        <- nrow(xglobal)
  bigK        <- ncol(xglobal)
  Kbig        <- plag*bigK
  bigT        <- Traw-plag
  A_large     <- x$stacked.results$A_large
  F_large     <- x$stacked.results$F_large
  S_large     <- x$stacked.results$S_large
  Ginv_large  <- x$stacked.results$Ginv_large
  F.eigen     <- x$stacked.results$F.eigen
  x           <- xglobal[(plag+1):Traw,,drop=FALSE]
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
  if(is.null(cores)) {cores <- 1}
  #-----------------------------------------------------------------------------------------------------#
  if(running){
    GFEVD_post <- array(0,dim=c(bigK,bigK,nhor)); dimnames(GFEVD_post)<-list(varNames, paste("Decomp. of",varNames),0:(nhor-1))
    if(verbose) cat(paste("Start computation on ", cores, " cores", " (",thindraws," stable draws in total).",sep=""),"\n")
    
    imp.obj <- applyfun(1:thindraws,function(irep){
      irfa <- .irf.girf.sims(invG=Ginv_large[irep,,],lF=adrop(F_large[irep,,,,drop=FALSE],drop=1),gcov=S_large[irep,,],
                             x,horizon=nhor)$impl
      GFEVD <- .mk_fevd.sims(irfa)
      return(list(GFEVD=GFEVD))
    })
    for(irep in 1:thindraws){
      GFEVD_post<-GFEVD_post+imp.obj[[irep]]$GFEVD
    }
    GFEVD_post<-GFEVD_post/thindraws
  }else{ #-------------------------HERE DO FULL CALCULATION INCLUDING BOUNDS- VERY MEMORY INTENSIVE!!!-----
    GFEVD_draws<-array(NA,dim=c(thindraws,bigK,bigK,nhor))
    GFEVD_post <-array(NA,dim=c(bigK,bigK,nhor,3))
    dimnames(GFEVD_post)[[1]]<-dimnames(GFEVD_post)[[2]]<-varNames
    dimnames(GFEVD_post)[[3]]<-0:(nhor-1)
    dimnames(GFEVD_post)[[4]]<-c("low16","median","high84")
    if(verbose) cat(paste("Start computation on ", cores, " cores", " (",thindraws," stable draws in total).",sep=""),"\n")
    imp.obj <- applyfun(1:thindraws,function(irep){
      irfa <- .irf.girf.sims(invG=Ginv_large[irep,,],lF=adrop(F_large[irep,,,,drop=FALSE],drop=1),gcov=S_large[irep,,],
                             x,horizon=nhor)$impl
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
                        xglobal=xglobal),
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
  
}
