#' @export
"hd" <- function(x, rotation.matrix=NULL, verbose=TRUE){
  UseMethod("hd", x)
}

#' @name hd
#' @title Historical Decomposition
#' @description A function that calculates historical decomposition (HD) of the time series and the structural error.
#' @method hd bgvar.irf
#' @export
#' @usage hd(x, rotation.matrix=NULL, verbose=TRUE)
#' @param x an item fitted by \code{irf}.
#' @param rotation.matrix If \code{NULL} and the \code{irf.bgvar} object has been fitted via sign restrictions, the rotation matrix is used that minimizes the distance to the median impulse responses at the posterior median.
#' @param verbose If set to \code{FALSE} it suppresses printing messages to the console.
#' @details To save computational time as well as due to storage limits, both functions are based on the posterior median (as opposed to calculating HDs and the structural error for each draw of the MCMC chain). In case the shock has been identified via sign restrictions, a rotation matrix has to be selected to calculate both statistics. If not specified otherwise (via \code{R}), the algorithm searches for 50 rotation matrices that fulfill the sign restrictions at the \emph{posterior median} of the coefficients and then singles out the rotation matrix that minimizes the distance to the median of the impulse responses as suggested in Fry and Pagan (2011).
#' @return Returns a list with the following objects \describe{
#' \item{\code{hd_array}}{ is a three-dimensional array with the first dimension referring to the K time series, the second to the T observations and the third dimensions containing the contribution of the shocks in explaining historically deviations in the time series from their trend. The third dimension is K+3, since the last three entries contain the contributions of the constant, the initial condition and a residual component that the contributions sum up to the original time series. If a trend i specified in the model the third dimension is K+3 with trend ordered after the constant.}
#' \item{\code{struc.shcok}}{ contains the structural shock.}
#' \item{\code{x}}{ is a matrix object that contains the original time series, which is of dimension K times (T-plag).}
#' }
#' @author Maximilian Boeck, Martin Feldkircher, Florian Huber
#' @seealso \code{\link{bgvar}} and \code{\link{irf}}.
#' @examples
#' \donttest{
#' library(BGVAR)
#' data(testdata)
#' model.eer<-bgvar(Data=testdata, W=W.test, draws=50, burnin=50, 
#'                  plag=1, prior="SSVS", eigen=TRUE)
#'                  
#' # US monetary policy shock
#' shockinfo <- get_shockinfo("chol")
#' shockinfo$shock <- "US.stir"; shockinfo$scale <- -100
#' irf.chol.us.mp<-irf(model.eer,n.ahead=48,shockinfo=shockinfo)
#' 
#' # calculates historical decomposition
#' HD <- hd(irf.chol.us.mp)
#' }
#' @references 
#' Fry, R. and A. Pagan (2011) \emph{Sign restrictions in Structural Vector Autoregressions: A Critical Review}. Journal of Economic Literature, Vol. 49(4), pp. 938-960.
hd.bgvar.irf<-function(x, rotation.matrix=NULL, verbose=TRUE){
  start.hd <- Sys.time()
  if(verbose) cat("Start computing historical decomposition of Bayesian Global Vector Autoregression.\n\n")
  #------------------------------ get stuff -------------------------------------------------------#
  xglobal <- x$model.obj$xglobal
  lags    <- x$model.obj$lags
  pmax    <- max(lags)
  ident   <- x$ident
  Traw    <- nrow(xglobal)
  bigK    <- ncol(xglobal)
  xdat    <- xglobal[(pmax+1):Traw,,drop=FALSE]
  bigT    <- nrow(xdat)
  ALPHA   <- x$struc.obj$A
  Ginv    <- x$struc.obj$Ginv
  Smat    <- x$struc.obj$Smat
  Sigma_u <- Ginv%*%Smat%*%t(Ginv)
  varNames<- colnames(xglobal)
  trend   <- FALSE
  if(!is.null(rotation.matrix)){
    rotation.matrix<-x$struc.obj$Rmed
  }else{
    rotation.matrix<-diag(bigK)
  }
  rownames(rotation.matrix) <- colnames(rotation.matrix) <- varNames
  #------------------------checks-------------------------------------------------------------------#
  if(ident=="girf"){
    message("Historical decomposition of the time series not implemented for GIRFs since cross-correlation is unequal to zero (and hence decompositions do not sum up to original time series).")
    return(list(hd_array=NA,struc.shock=vv,xglobal=xglobal))
  }
  #------ initialize objects -----------------------------------------------------------------------#
  if("trend"%in%dimnames(ALPHA)[[2]]) trend <- TRUE
  struc_post <- array(NA,dim=c(bigT,bigK))
  hd_array   <- array(0,c(bigK,bigT,(bigK+3+trend)))
  if(trend){
    dimnames(hd_array)<-list(rownames(x), NULL, c(paste("contribution of shock to", c(varNames)),"constant","trend","initial cond.","residual"))
  }else{
    dimnames(hd_array)<-list(rownames(x), NULL, c(paste("contribution of shock to", c(varNames)),"constant","initial cond.","residual"))
  }
  #------------------------------------------------------------------------------------------------#
  Rinv       <- solve(rotation.matrix)
  Sigchol_u  <- t(chol(Sigma_u))
  Sigcholinv <- solve(Sigchol_u)
  
  vv <- matrix(0,bigT,bigK,dimnames=list(NULL,varNames))
  YY <- xglobal[(pmax+1):Traw,]
  XX <- cbind(.mlag(xglobal,pmax),1)
  XX <- XX[(pmax+1):nrow(XX),]
  if(trend) XX <- cbind(XX,seq(1,bigT))
  
  strMat <- Rinv%*%Sigcholinv
  
  for (t in 1:bigT){
    Yhat <- strMat%*%YY[t,]
    Xhat <- strMat%*%ALPHA%*%XX[t,]
    
    PHI <- Rinv%*%Sigcholinv%*%ALPHA
    vv[t,] <- Yhat-Xhat
  }
  #Start historical decompositions -------------------------------------------------------------------------------#
  if(verbose) cat("Start computing HDs...\n")
  HDshock_big <- array(0,c(pmax*bigK,bigT,bigK))
  HDconst_big <- matrix(0,pmax*bigK,bigT)
  HDinit_big  <- matrix(0,pmax*bigK,bigT)
  HDshock     <- array(0,c(bigK,bigT,bigK))
  HDinit      <- matrix(0,bigK,bigT)
  HDconst     <- matrix(0,bigK,bigT)
  if(trend){
    HDtrend_big <- matrix(0,pmax*bigK,bigT)
    HDtrend     <- matrix(0,bigK,bigT)
  }
  
  solveA <- (Sigchol_u%*%rotation.matrix) #Depends on identification, if Cholesky then solveA = t(chol(SIGMA)), where SIGMA is the VC of the global model
  eps <- (YY-XX%*%t(ALPHA))%*%t(solve(solveA)) #Atilda is the matrix of autoregressive coefficients of the global model
  Fcomp <-  .get_companion(ALPHA[,1:(bigK*pmax)],varndxv = c(bigK,0,pmax))$MM#Fcomp is the companion matrix (used in the eigenvalue stuff without the constant)
  
  invA_big <- matrix(0,bigK*pmax,bigK)  #M is the number of endogenous variables ; p is the number of lags
  invA_big[1:bigK,] <- solveA
  Icomp <- cbind(diag(bigK),matrix(0,bigK,(pmax-1)*bigK))
  for (nn in 2:bigT){
    for (jj in 1:bigK){
      eps_big <- matrix(0,bigK,1)
      eps_big[jj,] <- eps[nn,jj]
      HDshock_big[,nn,jj] <- (invA_big)%*%eps_big+Fcomp%*%HDshock_big[,nn-1,jj]
      HDshock[,nn,jj] <- Icomp%*%HDshock_big[,nn,jj]
    }
    #Initial value
    HDinit_big[,1] <- XX[1,1:(pmax*bigK)]
    HDinit[,1] <- Icomp%*%HDinit_big[,1]
    HDinit_big[,nn] <- Fcomp%*%HDinit_big[,nn-1]
    HDinit[,nn] <- Icomp%*%HDinit_big[,nn]
    
    #Constant
    CC <- matrix(0,bigK*pmax,1)
    CC[1:bigK] <- t(ALPHA)[(bigK*pmax)+1,]
    HDconst_big[,nn] <- CC+Fcomp%*%HDconst_big[,nn-1]
    HDconst[,nn] <- Icomp%*%HDconst_big[,nn]
    
    # Trend
    if(trend){
      TT <- matrix(0,bigK*pmax,1)
      TT[1:bigK] <- t(ALPHA)[(bigK*pmax)+2,]
      HDtrend_big[,nn] <- TT+Fcomp%*%HDtrend_big[,nn-1]
      HDtrend[,nn] <- Icomp%*%HDtrend_big[,nn]
    }
  }
  hd_array[,,1:bigK]   <- HDshock_big
  hd_array[,,(bigK+1)] <- HDconst_big
  if(trend) hd_array[,,(bigK+1+trend)] <- HDtrend_big
  hd_array[,,(bigK+2+trend)] <- HDinit_big
  hd_array[,,(bigK+3+trend)] <- (t(xdat)-apply(hd_array,c(1,2),sum)) # residual part
  #----------------------------------------------------------------------------------#
  hd_array <- aperm(hd_array,c(2,1,3))
  out      <- structure(list(hd_array=hd_array,struc_shock=vv,xglobal=xdat, R=NULL), class="bgvar.hd")
  if(verbose) cat(paste("Size of object:", format(object.size(out),unit="MB")))
  end.hd <- Sys.time()
  diff.hd <- difftime(end.hd,start.hd,units="mins")
  mins.hd <- round(diff.hd,0); secs.hd <- round((diff.hd-floor(diff.hd))*60,0)
  if(verbose) cat(paste("\nNeeded time for computation: ",mins.hd," ",ifelse(mins.hd==1,"min","mins")," ",secs.hd, " ",ifelse(secs.hd==1,"second.","seconds.\n"),sep=""))
  return(out)
}

#' @method print bgvar.hd
#' @export
print.bgvar.hd <- function(x, ...){
  cat("---------------------------------------------------------------------------------------")
  cat("\n")
  cat("Object contains historical decomposition of object estimated with 'bgvar':")
  cat("\n")
  cat(paste0("Size of hd_array containing historical decompositions: ",dim(x$hd_array)[[1]]," x ",dim(x$hd_array)[[2]]," x ",dim(x$hd_array)[[3]],"."))
  cat("\n")
  cat(paste0("Size of struc_shock containing structural errors: ",dim(x$struc_shock)[[1]]," x ",dim(x$struc_shock)[[2]],"."))
  cat("\n")
  cat("Identification scheme: ")
  if(is.null(x$R)){
    cat("Short-run restrictions via Cholesky decomposition.")
  }else{
    cat("Sign-restrictions.")
  }
  cat("\n")
  cat(paste0("Size ob object: ",format(object.size(x),unit="MB")))
  cat("\n")
  cat("---------------------------------------------------------------------------------------")
  
  return(invisible(x))
}