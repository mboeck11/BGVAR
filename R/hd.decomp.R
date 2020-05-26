#' @name hd.decomp
#' @title Historical Decomposition
#' @description A function that calculates historical decomposition (HD) of the time series and the structural error.
#' @usage hd.decomp(obj, R=NULL)
#' @param obj an item fitted by \code{IRF}.
#' @param R If \code{NULL} and the \code{irf.bgvar} object has been fitted via sign restrictions, the rotation matrix is used that minimizes the distance to the median impulse responses at the posterior median.
#' @details To save computational time as well as due to storage limits, both functions are based on the posterior median (as opposed to calculating HDs and the structural error for each draw of the MCMC chain). In case the shock has been identified via sign restrictions, a rotation matrix has to be selected to calculate both statistics. If not specified otherwise (via \code{R}), the algorithm searches for 50 rotation matrices that fulfill the sign restrictions at the \emph{posterior median} of the coefficients and then singles out the rotation matrix tha minimizes the distance to the median of the impulse responses as suggested in Fry and Pagan (2011).
#' @return Returns a list with the following objects \itemize{
#' \item{\code{hd_array}}{ is a three-dimensional array with the first dimension referring to the K time series, the second to the T observations and the third dimensions containing the contribution of the shocks in explaining historically deviations in the time series from their trend. The third dimension is K+3, since the last three entries contain the contributions of the constant, the initial condition and a residual component that the contributions sum up to the original time series. If a trend i specified in the model the third dimension is K+3 with trend ordered after the constant.}
#' \item{\code{struc.shcok}}{ contains the structural shock.}
#' \item{\code{x}}{ is a matrix object that contains the original time series, which is of dimension K times (T-plag).}
#' }
#' @author Maximilian Boeck, Martin Feldkircher, Florian Huber
#' @seealso \code{\link{bgvar}} and \code{\link{IRF}}.
#' @examples 
#' \donttest{
#' set.seed(571)
#' library(BGVAR)
#' data(eerData)
#' model.ssvs.eer<-bgvar(Data=eerData,W=W.trade0012,saves=100,burns=100,plag=1,
#'                       prior="SSVS",thin=1,eigen=TRUE)
#' # US monetary policy shock
#' shocks<-list();shocks$var="stir";shocks$cN<-"US";shocks$ident="chol";shocks$scal=-100
#' irf.chol.us.mp <- IRF(obj=model.ssvs.eer,shock=shocks,nhor=48)
#' 
#' HD <- hd.decomp(irf.chol.us.mp)
#' # summing them up should get you back the original time series
#' org.ts<-apply(HD$hd_array,c(1,2),sum)
#' matplot(cbind(HD$x[,1],org.ts[,1]),type="l",ylab="")
#' legend("bottomright",c("hd series","original"),col=c("black","red"),lty=c(1,2),bty="n")
#' }
#' @references 
#' Fry, R. and A. Pagan (2011) \emph{Sign restrictions in Structural Vector Autoregressions: A Critical Review}. Journal of Economic Literature, Vol. 49(4), pp. 938-960.
#' @export
hd.decomp<-function(obj, R=NULL){
  start.hd <- Sys.time()
  if(!inherits(obj, "bgvar.irf")) {stop("Please provide a `bgvar.irf` object.")}
  cat("\nStart computing historical decomposition of Bayesian Global Vector Autoregression.\n\n")
  #------------------------------ get stuff -------------------------------------------------------#
  xglobal <- obj$model.obj$xglobal
  plag    <- obj$model.obj$plag
  ident   <- obj$ident
  Traw    <- nrow(xglobal)
  bigK    <- ncol(xglobal)
  x       <- xglobal[(plag+1):Traw,,drop=FALSE]
  bigT    <- nrow(x)
  ALPHA   <- obj$struc.obj$A
  Ginv    <- obj$struc.obj$Ginv
  Smat    <- obj$struc.obj$Smat
  Sigma_u <- Ginv%*%Smat%*%t(Ginv)
  varNames<- colnames(xglobal)
  trend   <- FALSE
  if(!is.null(R)){
    R<-obj$struc.obj$Rmed
  }else{
    R<-diag(bigK)
  }
  rownames(R) <- colnames(R) <- varNames
  #------------------------checks-------------------------------------------------------------------#
  if(ident=="girf"){
    print("Historical decomposition of the time series not implemented for GIRFs since cross-correlation is unequal to zero (and hence decompositions do not sum up to original time series).")
    return(list(hd_array=NA,struc.shock=vv,xglobal=xglobal) )
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
  Rinv       <- solve(R)
  Sigchol_u  <- t(chol(Sigma_u))
  Sigcholinv <- solve(Sigchol_u)
  
  vv <- matrix(0,bigT,bigK,dimnames=list(NULL,varNames))
  YY <- xglobal[(plag+1):Traw,]
  XX <- cbind(.mlag(xglobal,plag),1)
  XX <- XX[(plag+1):nrow(XX),]
  if(trend) XX <- cbind(XX,seq(1,bigT))
  
  strMat <- Rinv%*%Sigcholinv
  
  for (t in 1:bigT){
    Yhat <- strMat%*%YY[t,]
    Xhat <- strMat%*%ALPHA%*%XX[t,]
    
    PHI <- Rinv%*%Sigcholinv%*%ALPHA
    vv[t,] <- Yhat-Xhat
  }
  #Start historical decompositions -------------------------------------------------------------------------------#
  cat("Start computing HDs...\n")
  HDshock_big <- array(0,c(plag*bigK,bigT,bigK))
  HDconst_big <- matrix(0,plag*bigK,bigT)
  HDinit_big  <- matrix(0,plag*bigK,bigT)
  HDshock     <- array(0,c(bigK,bigT,bigK))
  HDinit      <- matrix(0,bigK,bigT)
  HDconst     <- matrix(0,bigK,bigT)
  if(trend){
    HDtrend_big <- matrix(0,plag*bigK,bigT)
    HDtrend     <- matrix(0,bigK,bigT)
  }
  
  #NOTE FOR MARTIN: IF SIGN RESTRICTIONS: R = ROTATION ELSE R = diag(M)
  solveA <- (Sigchol_u%*%R) #Depends on identification, if Cholesky then solveA = t(chol(SIGMA)), where SIGMA is the VC of the global model
  eps <- (YY-XX%*%t(ALPHA))%*%t(solve(solveA)) #Atilda is the matrix of autoregressive coefficients of the global model
  Fcomp <-  .get_companion(ALPHA[,1:(bigK*plag)],varndxv = c(bigK,0,plag))$MM#Fcomp is the companion matrix (used in the eigenvalue stuff without the constant)
  
  invA_big <- matrix(0,bigK*plag,bigK)  #M is the number of endogenous variables ; p is the number of lags
  invA_big[1:bigK,] <- solveA
  Icomp <- cbind(diag(bigK),matrix(0,bigK,(plag-1)*bigK))
  for (nn in 2:bigT){
    for (jj in 1:bigK){
      eps_big <- matrix(0,bigK,1)
      eps_big[jj,] <- eps[nn,jj]
      HDshock_big[,nn,jj] <- (invA_big)%*%eps_big+Fcomp%*%HDshock_big[,nn-1,jj]
      HDshock[,nn,jj] <- Icomp%*%HDshock_big[,nn,jj]
    }
    #Initial value
    HDinit_big[,1] <- XX[1,1:(plag*bigK)]
    HDinit[,1] <- Icomp%*%HDinit_big[,1]
    HDinit_big[,nn] <- Fcomp%*%HDinit_big[,nn-1]
    HDinit[,nn] <- Icomp%*%HDinit_big[,nn]
    
    #Constant
    CC <- matrix(0,bigK*plag,1)
    CC[1:bigK] <- t(ALPHA)[(bigK*plag)+1,]
    HDconst_big[,nn] <- CC+Fcomp%*%HDconst_big[,nn-1]
    HDconst[,nn] <- Icomp%*%HDconst_big[,nn]
    
    # Trend
    if(trend){
      TT <- matrix(0,bigK*plag,1)
      TT[1:bigK] <- t(ALPHA)[(bigK*plag)+2,]
      HDtrend_big[,nn] <- TT+Fcomp%*%HDtrend_big[,nn-1]
      HDtrend[,nn] <- Icomp%*%HDtrend_big[,nn]
    }
  }
  hd_array[,,1:bigK]   <- HDshock_big
  hd_array[,,(bigK+1)] <- HDconst_big
  if(trend) hd_array[,,(bigK+1+trend)] <- HDtrend_big
  hd_array[,,(bigK+2+trend)] <- HDinit_big
  hd_array[,,(bigK+3+trend)] <- (t(x)-apply(hd_array,c(1,2),sum)) # residual part
  #----------------------------------------------------------------------------------#
  hd_array <- aperm(hd_array,c(2,1,3))
  out      <- structure(list(hd_array=hd_array,struc.shock=vv,x=x), class="bgvar.hd")
  cat(paste("Size of object:", format(object.size(out),unit="MB")))
  end.hd <- Sys.time()
  diff.hd <- difftime(end.hd,start.hd,units="mins")
  mins.hd <- round(diff.hd,0); secs.hd <- round((diff.hd-floor(diff.hd))*60,0)
  cat(paste("\nNeeded time for computation: ",mins.hd," ",ifelse(mins.hd==1,"min","mins")," ",secs.hd, " ",ifelse(secs.hd==1,"second.","seconds.\n"),sep=""))
  return(out)
}