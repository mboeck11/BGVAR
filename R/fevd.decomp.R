#' @name fevd.decomp
#' @title Forecast Error Variance Decomposition
#' @description This function calculates the forecast error variance decomposition (FEVDs) for Cholesky and sign-identified shocks.
#' @usage fevd.decomp(obj, R=NULL, var.slct=NULL, verbose=TRUE)
#' @details Since the calculations are very time consuming, the FEVDs are based on the posterior median only (as opposed to calculating FEVDs for each MCMC sweep). In case the underlying shock has been identified via sign restrictions, the rotation matrix corresponds to the one that fulfills the sign restrictions at the posterior median of the estimated coefficients. More precisely, the algorithm searches for 50 rotation matrices that fulfill the sign restrictions at the \emph{posterior median} of the coefficients and then singles out the rotation matrix that minimizes the distance to the median of the impulse responses as suggested in Fry and Pagan (2011).
#' @param obj an object of class \code{bgvar.irf}.
#' @param R If \code{NULL} and the \code{obj} has been fitted via sign restrictions, the rotation matrix is used that minimizes the distance to the median impulse responses at the posterior median.
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
#' model.ssvs.eer<-bgvar(Data=eerData,W=W.trade0012,saves=50,burns=50,plag=1,
#'                       prior="SSVS",thin=1,eigen=TRUE)
#'                       
#' # US monetary policy shock
#' shocks<-list();shocks$var="stir";shocks$cN<-"US";shocks$ident="chol";shocks$scal=-100
#' irf.chol.us.mp<-IRF(obj=model.ssvs.eer,shock=shocks,nhor=48)
#' 
#' # calculates FEVD for variables US.Dp and EA.y
#' fevd.us.mp=fevd.decomp(obj=irf.chol.us.mp,var.slct=c("US.Dp","EA.y"))
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
#' irf.sign.us.mp<-IRF(obj=model.ssvs.eer,sign.constr=sign.constr,nhor=24)
#' 
#' # calculates FEVD for variables US.Dp and EA.y
#' fevd.us.mp=fevd.decomp(obj=irf.sign.us.mp,var.slct=c("US.Dp","EA.y"))
#' }
#' \donttest{
#' library(BGVAR)
#' data(eerData)
#' model.ssvs.eer<-bgvar(Data=eerData,W=W.trade0012,saves=100,burns=100,plag=1,
#'                       prior="SSVS",thin=1,eigen=TRUE)
#'                       
#' # US monetary policy shock
#' shocks<-list();shocks$var="stir";shocks$cN<-"US";shocks$ident="chol";shocks$scal=-100
#' irf.chol.us.mp<-IRF(obj=model.ssvs.eer,shock=shocks,nhor=48)
#' 
#' # calculates FEVD for variables US.Dp and EA.y
#' fevd.us.mp=fevd.decomp(obj=irf.chol.us.mp,var.slct=c("US.Dp","EA.y"))
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
#' irf.sign.us.mp<-IRF(obj=model.ssvs.eer,sign.constr=sign.constr,nhor=24)
#' 
#' # calculates FEVD for variables US.Dp and EA.y
#' fevd.us.mp=fevd.decomp(obj=irf.sign.us.mp,var.slct=c("US.Dp","EA.y"))
#' }
#' # NOT RUN - calculates FEVDs for all variables in the system, very time consuming
#' \dontrun{
#  fevd.us.mp=fevd.decomp(obj=irf.chol.us.mp,var.slct=NULL)
#' }
#' @export
fevd.decomp <- function(obj,R=NULL,var.slct=NULL,verbose=TRUE){
  start.fevd <- Sys.time()
  if(!inherits(obj, "bgvar.irf")) {stop("Please provide a `bgvar.irf` object.")}
  if(verbose) cat("\nStart computing forecast error variance decomposition of Bayesian Global Vector Autoregression.\n\n")
  #------------------------------ get stuff -------------------------------------------------------#
  xglobal <- obj$model.obj$xglobal
  plag    <- obj$model.obj$plag
  ident   <- obj$ident
  Traw    <- nrow(xglobal)
  bigK    <- ncol(xglobal)
  x       <- xglobal[(plag+1):Traw,,drop=FALSE]
  bigT    <- nrow(x)
  A       <- obj$struc.obj$A
  Fmat    <- obj$struc.obj$Fmat
  Ginv    <- obj$struc.obj$Ginv
  Smat    <- obj$struc.obj$S
  Sigma_u <- Ginv%*%Smat%*%t(Ginv)
  horizon <- dim(obj$posterior)[2]
  varNames<- colnames(xglobal)
  shock   <- obj$shock
  sign.constr <- obj$sign.constr
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
    R <- obj$struc.obj$Rmed
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
                   class="bgvar.fevd")
  if(verbose) cat(paste("\nSize of FEVD object: ", format(object.size(FEVDres),unit="MB")))
  end.fevd <- Sys.time()
  diff.fevd <- difftime(end.fevd,start.fevd,units="mins")
  mins.fevd <- round(diff.fevd,0); secs.fevd <- round((diff.fevd-floor(diff.fevd))*60,0)
  if(verbose) cat(paste("\nNeeded time for computation: ",mins.fevd," ",ifelse(mins.fevd==1,"min","mins")," ",secs.fevd, " ",ifelse(secs.fevd==1,"second.","seconds.\n"),sep=""))
  return(out)
}

#' @name gfevd.decomp
#' @title Generalized Forecast Error Variance Decomposition
#' @description This function calculates a complete generalized forecast error variance decomposition (GFEVDs) based on generalized impulse response functions akin to Lanne-Nyberg (2016). The Lanne-Nyberg (2016) corrected GFEVD sum up to unity.
#' @usage gfevd.decomp(obj, nhor=24, running=TRUE, multithread=FALSE, verbose=TRUE)
#' @param obj an object of class \code{bgvar}.
#' @param nhor the forecast horizon.
#' @param running Default is set to \code{TRUE} and implies that only a running mean over the posterior draws is calculated. A full analysis including posterior bounds is likely to cause memory issues.
#' @param multithread If set to \code{TRUE} parallel computing using the packages \code{\link{foreach}} and \code{\link{doParallel}}. Number of cores is set to maximum number of cores in the computer. This option is recommended when working with sign restrictions to speed up computations. Default is set to \code{FALSE} and thus no parallelization.
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
#' model.ssvs.eer<-bgvar(Data=eerData,W=W.trade0012,saves=100,burns=100,plag=1,
#'                       prior="SSVS",thin=1,eigen=TRUE)
#'                       
#' GFEVD<-gfevd.decomp(model.ssvs.eer,nhor=24,running=TRUE)
#' }
#' \donttest{
#' library(BGVAR)
#' data(eerData)
#' model.ssvs.eer<-bgvar(Data=eerData,W=W.trade0012,saves=100,burns=100,plag=1,
#'                       prior="SSVS",thin=1,eigen=TRUE)
#'                       
#' # Calculates running mean GFEVDs for all variables in the system 
#' GFEVD<-gfevd.decomp(model.ssvs.eer,nhor=24,running=TRUE)
#' }
#' @export
#' @importFrom abind adrop
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel detectCores
gfevd.decomp<-function(obj,nhor=24,running=TRUE,multithread=FALSE,verbose=TRUE){
  start.gfevd <- Sys.time()
  if(!inherits(obj, "bgvar")) {stop("Please provide a `bgvar` object.")}
  if(verbose) cat("\nStart computing generalized forecast error variance decomposition of Bayesian Global Vector Autoregression.\n\n")
  #------------------------------ get stuff -------------------------------------------------------#
  plag        <- obj$args$plag
  xglobal     <- obj$xglobal
  Traw        <- nrow(xglobal)
  bigK        <- ncol(xglobal)
  Kbig        <- plag*bigK
  bigT        <- Traw-plag
  A_large     <- obj$stacked.results$A_large
  F_large     <- obj$stacked.results$F_large
  S_large     <- obj$stacked.results$S_large
  Ginv_large  <- obj$stacked.results$Ginv_large
  F.eigen     <- obj$stacked.results$F.eigen
  x           <- xglobal[(plag+1):Traw,,drop=FALSE]
  thinsaves   <- length(F.eigen) ### prior: saves
  varNames    <- colnames(xglobal)
  #-----------------------------------------------------------------------------------------------------#
  if(running){
    GFEVD_post <- array(0,dim=c(bigK,bigK,nhor)); dimnames(GFEVD_post)<-list(varNames, paste("Decomp. of",varNames),0:(nhor-1))
    if(multithread){
      numCores <- detectCores()
      registerDoParallel(cores=numCores)
      if(verbose) cat(paste("Start computation on ", numCores, " cores", " (",thinsaves," stable draws in total).",sep=""),"\n")
      imp.obj <-foreach(irep=1:thinsaves) %dopar% {
        irfa <- .irf.girf.sims(invG=Ginv_large[irep,,],lF=adrop(F_large[irep,,,,drop=FALSE],drop=1),gcov=S_large[irep,,],
                               x,horizon=nhor)$impl
        GFEVD <- .mk_fevd.sims(irfa)
        return(list(GFEVD=GFEVD))
      }
      for(irep in 1:thinsaves){
        GFEVD_post<-GFEVD_post+imp.obj[[irep]]$GFEVD
      }
      GFEVD_post<-GFEVD_post/thinsaves
    }else{
      if(verbose) cat(paste("Start computation on single core", " (",thinsaves," stable draws in total).",sep=""),"\n")
      for(irep in 1:thinsaves){
        irfa <- .irf.girf.sims(invG=Ginv_large[irep,,],lF=adrop(F_large[irep,,,,drop=FALSE],drop=1),gcov=S_large[irep,,],
                               x,horizon=nhor)$impl
        GFEVD_post<-GFEVD_post+.mk_fevd.sims(irfa)
      }
      GFEVD_post<-GFEVD_post/thinsaves
    }
  }else{ #-------------------------HERE DO FULL CALCULATION INCLUDING BOUNDS- VERY MEMORY INTENSIVE!!!-----
    GFEVD_draws<-array(NA,dim=c(thinsaves,bigK,bigK,nhor))
    GFEVD_post <-array(NA,dim=c(bigK,bigK,nhor,3))
    dimnames(GFEVD_post)[[1]]<-dimnames(GFEVD_post)[[2]]<-varNames
    dimnames(GFEVD_post)[[3]]<-0:(nhor-1)
    dimnames(GFEVD_post)[[4]]<-c("low16","median","high84")
    if(multithread){# use parallel computing for IRF analysis
      numCores <- detectCores()
      registerDoParallel(cores=numCores)
      if(verbose) cat(paste("Start computation on ", numCores, " cores", " (",thinsaves," stable draws in total).",sep=""),"\n")
      imp.obj <-foreach(irep=1:thinsaves) %dopar%{
        irfa <- .irf.girf.sims(invG=Ginv_large[irep,,],lF=adrop(F_large[irep,,,,drop=FALSE],drop=1),gcov=S_large[irep,,],
                               x,horizon=nhor)$impl
        GFEVD <- .mk_fevd.sims(irfa)
        return(list(GFEVD=GFEVD))
      }
      for(irep in 1:thinsaves){
        GFEVD_draws[irep,,,] <- imp.obj[[irep]]$GFEVD
      }
    }else{ #carry out impulse response function analysis in a loop
      if(verbose) cat(paste("Start computation on single core", " (",thinsaves," stable draws in total).",sep=""),"\n")
      for(irep in 1:thinsaves){
        irfa <- .irf.girf.sims(invG=Ginv_large[irep,,],lF=adrop(F_large[irep,,,,drop=FALSE],drop=1),gcov=S_large[irep,,],
                               x,horizon=nhor)$impl
        GFEVD_draws[irep,,,]<-.mk_fevd.sims(irfa)
        
      }
    }
    GFEVD_post[,,,1] <- apply(GFEVD_draws,c(2,3,4),quantile,.16,na.rm=TRUE)
    GFEVD_post[,,,2] <- apply(GFEVD_draws,c(2,3,4),mean,na.rm=TRUE)
    GFEVD_post[,,,3] <- apply(GFEVD_draws,c(2,3,4),quantile,.16,na.rm=TRUE)
  }
  #----------------------------------------------------------------------------------------------------------------
  out <- structure(list(GFEVD=GFEVD_post,
                        xglobal=xglobal),
                   class="bgvar.fevd")
  if(!running) out$GFEVD_store <- GFEVD_draws
  if(verbose) cat(paste("Size of IRF object:", format(object.size(out),unit="MB")))
  end.gfevd <- Sys.time()
  diff.gfevd <- difftime(end.gfevd,start.gfevd,units="mins")
  mins.gfevd <- round(diff.gfevd,0); secs.gfevd <- round((diff.gfevd-floor(diff.gfevd))*60,0)
  if(verbose) cat(paste("\nNeeded time for computation: ",mins.gfevd," ",ifelse(mins.gfevd==1,"min","mins")," ",secs.gfevd, " ",ifelse(secs.gfevd==1,"second.","seconds.\n"),sep=""))
  return(out)
}

#' @name plot.bgvar.fevd
#' @title Plotting Function for Forecast Error Variance Decomposition
#' @description  Plots the decomposition of a specific time series into selected structural shocks.
#' @param x an object of class \code{bgvar.fevd}.
#' @param ... additional arguments.
#' @param ts specify the decomposed time series to be plotted.
#' @param k.max plots the k series with the highest for the decomposition of \code{ts}.
#' @return No return value.
#' @author Maximilian Boeck
#' @examples
#' \donttest{
#' library(BGVAR)
#' data(eerData)
#' model.ssvs.eer<-bgvar(Data=eerData,W=W.trade0012,saves=100,burns=100,plag=1,
#'                       prior="SSVS",thin=1,eigen=TRUE)
#'                       
#' # US monetary policy shock
#' shocks<-list();shocks$var="stir";shocks$cN<-"US";shocks$ident="chol";shocks$scal=-100
#' irf.chol.us.mp<-IRF(obj=model.ssvs.eer,shock=shocks,nhor=48)
#' 
#' # calculates FEVD for variables US.Dp and EA.y
#' fevd.us.mp=fevd.decomp(obj=irf.chol.us.mp,var.slct=c("US.Dp","EA.y"))
#' 
#' plot(fevd.us.mp, ts="US.Dp", k.max=10)
#' }
#' @seealso \code{\link{IRF}}
#' @importFrom graphics abline matplot polygon
#' @importFrom MASS Null
#' @importFrom stats rnorm
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
plot.bgvar.fevd<-function(x, ..., ts, k.max=10){
  # restore user par settings on exit
  oldpar <- par(no.readonly=TRUE)
  on.exit(par(oldpar))
  if(!inherits(x, "bgvar.fevd")) {stop("Please provide a `bgvar.fevd` object.")}
  fevd      <- x[[1]]
  xglobal   <- x$xglobal
  varNames  <- colnames(xglobal)
  varAll    <- varNames
  cN        <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x) x[1]))
  vars      <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x) x[2]))
  nhor      <- dim(fevd)[[3]]
  
  ts        <- paste("Decomp. of ", ts,sep="")
  if(length(ts)>1){
    stop("Please provide just one time series in 'ts'.")
  }
  if(!(ts%in%dimnames(fevd)[[2]])){
    stop("Please provide time series present in dataset.")
  }
  if(is.numeric(k.max)){
    mean <- apply(fevd,c(1,2),mean)
    resp <- names(sort(mean[,ts], decreasing=TRUE))[1:k.max]
  }
  varNames <- list()
  for(kk in 1:ceiling(k.max/10)){
    if(kk*10>k.max) kk.max <- k.max else kk.max <- kk*10
    varNames[[kk]] <- resp[((kk-1)*10+1):kk.max]
  }
  for(kk in 1:length(varNames)){
    rows <- length(varNames[[kk]])/2
    if(rows<1) cols <- 1 else cols <- 2
    if(rows%%1!=0) rows <- ceiling(rows)
    if(rows%%1!=0) rows <- ceiling(rows)
    # update par settings
    newpar <- oldpar
    if(prod(oldpar$mfrow)<(rows*cols)) newpar$mfrow <- c(rows,cols)
    newpar$mar <- bgvar.env$mar
    par(newpar)
    for(kkk in 1:length(varNames[[kk]])){
      idx <- grep(varNames[[kk]][kkk],varAll)
      x<-fevd[idx,ts,]
      b<-range(x); b1<-b[1]; b2<-b[2]
      plot.ts(x,col=bgvar.env$plot$col.50,xaxt="n",yaxt="n",lwd=bgvar.env$plot.lwd.line,ylab="",
              main=varAll[idx],cex.main=bgvar.env$plot.cex.main,cex.axis=bgvar.env$plot$cex.axis,
              cex.lab=bgvar.env$plot$cex.lab,lty=1,ylim=c(b1,b2),xlab="")
      
      axis(2, at=seq(b1,b2,length.out=5), labels=format(seq(b1,b2,length.out=5),digits=2,nsmall=1),cex.axis=1.2,las=1)
      axisindex<-seq(1,length(x),by=4)
      axis(side=1, las=1,at=axisindex, labels=c(0:length(x))[axisindex], cex.axis=1.6,tick=FALSE)
      abline(v=axisindex,col=bgvar.env$plot$col.tick,lty=bgvar.env$plot$lty.tick)
    }
    if(kk<length(varNames)) readline(prompt="Press enter for next group...")
  }
}
