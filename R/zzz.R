#' @name bgvar.sim
#' @title Simulating a Global Vector Autoregression
#' @description This function is used to produce simulated realizations which follow a Global Vector Autorgression (GVAR). It will also automatically simulate coefficients. All parameters can also be set by the user.
#' @usage bgvar.sim(len, M, N, plag=1, cons=FALSE, trend=FALSE, SV=FALSE)
#' @details For testing purposes, this function enables to simulate time series processes which can be described by a Global Vector Autoregression. Since stability conditions are not checked, it is only implemented for \code{M=3}.
#' @param len length of the simulated time series.
#' @param M number of endogenous variables.
#' @param N number of countries.
#' @param plag number of lags.
#' @param cons logical indicating whether to include an intercept. Default set to \code{FALSE}.
#' @param trend logical indicating whether to include an intercept. Default set to \code{FALSE}.
#' @param SV logical indicating whether the process should be simulated with or without stochastic volatility. Default set to \code{FALSE}.
#' @return Returns a list with the following elements
#' @author Maximilian Boeck
#' @seealso 
#' \code{\link{bgvar}} for estimation of a \code{bgvar} object.
#' @examples 
#' library(BGVAR)
#' sim <- bgvar.sim(len=200, M=3, N=4, plag=2, cons=TRUE, trend=FALSE, SV=TRUE)
#' Data = sim$obs$xglobal
#' W    = sim$obs$W
#' @importFrom abind adrop
#' @importFrom bayesm rdirichlet
#' @importFrom stats rnorm
#' @importFrom stochvol svsim
#' @noRd
bgvar.sim <- function(len, M, N, plag=1, cons=FALSE, trend=FALSE, SV=FALSE){
  idi.para <- t(matrix(rep(c(-5,0.9,0.01),M*N),3,M*N))
  # -----------------------------------------------------------------------------------------
  # simluate shocks
  shock <- matrix(NA, len, M*N)
  vol_true <- matrix(NA, len, M*N)
  if(SV) {
    for(mm in 1:(M*N)) {
      temp <- svsim(len=len, mu=idi.para[mm,1], phi=idi.para[mm,2], sigma=idi.para[mm,3])
      shock[,mm] <- temp$y
      vol_true[,mm] <- temp$vol
    }
  } else {
    for(mm in 1:(M*N)) {
      shock[,mm]    <- rnorm(len, 0, exp(-5))
      vol_true[,mm] <- exp(-5)
    }
  }
  
  xglobal <- matrix(0,len,M*N)
  ident.country <- t(matrix(seq(1,M*N),M,N))
  cN   <- paste0(letters[1:N],letters[1:N],sep="")
  vars <- c("y","inf","stir")
  colnames(xglobal) <- colnames(shock) <- paste0(rep(cN,each=M),".",vars)
  
  # state 0
  Theta.mean   <- matrix(c(0.7,-0.03,0.08,-0.11,0.4,0,0.1,0.03,0.7),M,M)
  Lambda0.mean <- matrix(c(0.05,0.01,0.07,0,0.03,0,0.01,0,0.04),M,M)
  Lambda.mean  <- matrix(c(0.1,0,-0.05,0.02,0.08,-0.01,0,0,0.12),M,M)
  
  max(abs(Re(eigen(Theta.mean)$values)))
  max(abs(Re(eigen(Lambda0.mean)$values)))
  max(abs(Re(eigen(Lambda.mean)$values)))
  
  V.Theta   <- diag(M)*1e-7
  V.Lambda0 <- diag(M)*1e-7
  V.Lambda  <- diag(M)*1e-7
  
  # -----------------------------------------------------------------------------------------
  # get weights
  # weights <- list()
  # for(cc in 1:N){
  #   temp <- matrix(0,M,N); colnames(temp) <- cN; rownames(temp) <- vars
  #   for(ii in 1:length(vars)){
  #     temp[ii,-cc] <- rdirichlet(rep(1/(N-1),N-1))
  #   }
  #   weights[[cc]] <- temp
  # }
  # names(weights) <- cN
  weights <- matrix(0,N,N); rownames(weights) <- colnames(weights) <- cN
  for(cc in 1:N){
    weights[cc,-cc] <- rdirichlet(rep(1/(N-1),N-1)*100)
  }
  # -----------------------------------------------------------------------------------------
  # create weight matrix
  W <- list()
  for(cc in 1:N){
    Wnew <- matrix(0,2*M,M*N); colnames(Wnew) <- colnames(xglobal)
    rownames(Wnew) <- c(vars,paste0(vars,"*"))
    diag(Wnew[,grep(cN[cc],colnames(Wnew))]) <- 1
    for(ii in 1:M){
      #Wnew[paste0(vars[ii],"*"),grep(vars[ii],colnames(Wnew))] <- weights[[cc]][grep(vars[ii],rownames(weights[[cc]])),]
      Wnew[paste0(vars[ii],"*"),grep(vars[ii],colnames(Wnew))] <- weights[cc,]
    }
    W[[cc]] <- Wnew
  }
  names(W) <- cN
  # -----------------------------------------------------------------------------------------
  # simulate coefficients
  Theta   <- Lambda <- list()
  for(pp in 1:plag){
    Theta[[pp]]  <- array(NA,dim=c(M,M,N))
    Lambda[[pp]] <- array(NA,dim=c(M,M,N))
    dimnames(Theta[[pp]])[[1]] <- vars
    dimnames(Lambda[[pp]])[[1]] <- paste0(vars,"*")
    dimnames(Theta[[pp]])[[2]] <- dimnames(Lambda[[pp]])[[2]] <- vars
    dimnames(Theta[[pp]])[[3]] <- dimnames(Lambda[[pp]])[[3]] <- cN
  }
  Lambda0 <- array(NA,dim=c(M,M,N))
  a0l     <- array(NA,dim=c(1,M,N))
  a1l     <- array(NA,dim=c(1,M,N))
  for(cc in 1:N){
    for(pp in 1:plag){
      if(pp==1){
        Theta[[pp]][,,cc]   <- matrix(as.vector(Theta.mean) + kronecker(diag(M),t(chol(V.Theta))) %*% rnorm(M*M), M, M)
        Lambda[[pp]][,,cc]  <- matrix(as.vector(Lambda.mean) + kronecker(diag(M),t(chol(V.Lambda))) %*% rnorm(M*M), M, M)
      }else{
        Theta[[pp]][,,cc]   <- matrix(kronecker(diag(M),t(chol(V.Theta))) %*% rnorm(M*M), M, M)
        Lambda[[pp]][,,cc]  <- matrix(kronecker(diag(M),t(chol(V.Lambda))) %*% rnorm(M*M), M, M)
      }
    }
    Lambda0[,,cc] <- matrix(as.vector(Lambda0.mean) + kronecker(diag(M),t(chol(V.Lambda0))) %*% rnorm(M*M), M, M)
    a0l[1,,cc] <- runif(M,-.3,.3)
    a1l[1,,cc] <- runif(M,1e-10,5e-4)
  }
  names(Theta) <- names(Lambda) <- paste0("lag.",seq(1,plag))
  # -----------------------------------------------------------------------------------------
  # get global solution
  a0     <- NULL
  a1     <- NULL
  G      <- NULL
  S_post <- list()
  for (cc in 1:N){
    A   <- cbind(diag(M),-t(Lambda0[,,cc]))
    for(pp in 1:plag){
      assign(paste0("B",pp),cbind(t(Theta[[pp]][,,cc]),t(Lambda[[pp]][,,cc])))
      if(cc==1) assign(paste0("H",pp), get(paste0("B",pp))%*%W[[cc]])
      if(cc>1)  assign(paste0("H",pp), rbind(get(paste0("H",pp)),get(paste0("B",pp))%*%W[[cc]]))
    }
    G            <- rbind(G,A%*%W[[cc]])
    a0           <- rbind(a0,t(adrop(a0l[,,cc,drop=FALSE],drop=3)))
    a1           <- rbind(a1,t(adrop(a1l[,,cc,drop=FALSE],drop=3)))
    S_post[[cc]] <- crossprod(shock[,grep(cN[cc],colnames(shock))])
  }
  G.inv   <- solve(G)
  b0      <- G.inv%*%a0
  b1      <- G.inv%*%a1
  F_large <- NULL
  F_sum   <- matrix(0,M*N,M*N)
  for (pp in 1:plag){
    assign(paste("F",pp,sep=""),G.inv%*%get(paste0("H",pp)))
    F_large <- cbind(F_large,get(paste0("F",pp)))
    F_sum   <- F_sum + get(paste0("F",pp))
  }
  Cm      <- rbind(F_large,cbind(diag(M*N*(plag-1)),matrix(0,M*N,M*N)))
  if(cons){
    F_large <- cbind(F_large,b0)
  }else{
    a0l <- NULL; a0 <- NULL; b0 <- matrix(0,M*N,1)
  }
  if(trend){
    F_large <- cbind(F_large,b1)
  }else{
    a1l <- NULL; a1 <- NULL; b1 <- matrix(0,M*N,1)
  }
  mu <- solve(diag(M*N)-F_sum)%*%b0
  #max(abs(Re(eigen(Cm)$values)))

  ##### need global solution
  for(pp in 1:plag) xglobal[pp,] <- mu + pp*b1 + shock[pp,]
  for (tt in (plag+1):len){
    xlag <- NULL
    for(pp in 1:plag) xlag <- cbind(xlag,xglobal[tt-1,,drop=FALSE])
    if(cons) xlag <- cbind(xlag,1)
    if(trend) xlag <- cbind(xlag,tt)
    xglobal[tt,] <- xlag%*% t(F_large) + shock[tt,]
  }
  
  
  true.global <- list(F_large=F_large, G.inv=G.inv, S_post=S_post)
  true.cc     <- list(a0l=a0l,a1l=a1l,Theta=Theta,Lambda0=Lambda0,Lambda=Lambda,vol_true)
  obs         <- list(xglobal=xglobal,W=weights)
  
  return(list(obs=obs,true.global=true.global,true.cc=true.cc))
}

#' @noRd
"irfcf" <- function(x, shockvar, resp, n.ahead=24, save.store=FALSE, verbose=TRUE){
  UseMethod("irfcf", x)
}

#' @name irfcf
#' @title Counterfactual Analysis
#' @description Function to perform counterfactual analysis. It enables to neutralize the response of a specific variable to a given shock.
#' @export
#' @usage irfcf(x, shockvar, resp, n.ahead=24, save.store=FALSE, verbose=TRUE)
#' @param x an object of class \code{bgvar}.
#' @param shockvar structural shock of interest.
#' @param resp response variable to neutralize.
#' @param n.ahead forecasting horizon.
#' @param save.store If set to \code{TRUE} the full posterior is returned. Default is set to \code{FALSE} in order to save storage.
#' @param verbose If set to \code{FALSE} it suppresses printing messages to the console.
#' @return Returns a list of class \code{bgvar.irf} with the following elements: \describe{
#' \item{\code{posterior}}{ is a four-dimensional array (K times K times n.ahead times 7) that contains 7 quantiles of the posterior distribution of the impulse response functions: the 50\% ("low25" and "high75"), the 68\% ("low16" and "high84") and the 90\% ("low05" and "high95") credible sets along with the posterior median ("median").}
#' \item{\code{struc.obj}}{ is a list object that contains posterior quantitites needed when calculating historical decomposition and structural errors via \code{hd.decomp}.\describe{
#' \item{\code{A}}{ median posterior of global coefficient matrix.}
#' \item{\code{Ginv}}{ median posterior of matrix \code{Ginv}, which describes contemporaneous relationships between countries.}
#' \item{\code{S}}{ posterior median of matrix with country variance-covariance matrices on the main diagonal.}
#' }}
#' \item{\code{model.obj}}{ is a list object that contains model-specific information, in particular\describe{
#' \item{\code{xglobal}}{ used data of the model.}
#' \item{\code{plag}}{ used lag specification of the model.}
#' }}
#' \item{\code{IRF_store}}{ is a four-dimensional array (K times n.ahead times nr. of shock times draws) which stores the whole posterior distribution. Exists only if \code{save.irf.store=TRUE}.}
#' }
#' @author Maximilian Boeck, Martin Feldkircher
#' @examples
#' \dontrun{
#' library(BGVAR)
#' data(eerDatasmall)
#' model.ssvs.eer<-bgvar(Data=eerDatasmall,W=W.trade0012.small,draws=100,burnin=100,
#'                       plag=1,prior="SSVS",eigen=TRUE)
#' # very time-consuming
#' irfcf <- irfcf(model.ssvs.eer,shockvar="US.stir",resp="US.rer",n.ahead=24)
#' }
#' @noRd
#' @importFrom stats quantile
irfcf.bgvar.irf <- function(x, shockvar, resp, n.ahead=24, save.store=FALSE, verbose=TRUE){
  start.irf <- Sys.time()
  if(verbose) cat("\nStart counterfactual analysis of Bayesian Global Vector Autoregression.\n\n")
  #----------------get stuff-------------------------------------------------------#
  plag        <- x$args$plag
  xglobal     <- x$xglobal
  bigK        <- ncol(xglobal)
  A_large     <- x$stacked.results$A_large
  F_large     <- x$stacked.results$F_large
  S_large     <- x$stacked.results$S_large
  Ginv_large  <- x$stacked.results$Ginv_large
  F.eigen     <- x$stacked.results$F.eigen
  thindraws   <- length(F.eigen)
  varNames    <- colnames(xglobal)
  #----------------------checks-----------------------------------------------------#
  if(length(shockvar)!=1&&length(resp)!=1){
    stop("Please specify only one shock and one response variable to neutralize.")
  }
  if(!(shockvar%in%varNames)){
    stop("Please respecify shockvar. Variable not contained in dataset.")
  }
  if(!(resp%in%varNames)){
    stop("Please respecify response variable. Variable not contained in dataset.")
  }
  neutR <- which(varNames%in%shockvar)
  neutS <- which(varNames%in%resp)
  if(verbose){
    cat(paste("Shock of interest: ",shockvar,".\n",sep=""))
    cat(paste("Response to neutralize: ",resp,".\n",sep=""))
  }
  #--------------compute-----------------------------------------------------------#
  if(verbose) cat("Start computing...\n")
  IRF_store     <- array(NA, dim=c(thindraws,bigK,bigK,n.ahead))
  dimnames(IRF_store)[[2]] <- dimnames(IRF_store)[[3]] <- colnames(xglobal)
  pb <- txtProgressBar(min = 0, max = thindraws, style = 3)
  for(irep in 1:thindraws){
    Sigma_u <- Ginv_large[irep,,]%*%S_large[irep,,]%*%t(Ginv_large[irep,,])
    irf<-Phi2<- .impulsdtrf(B=adrop(F_large[irep,,,,drop=FALSE],drop=1),
                            smat=t(chol(Sigma_u)),nstep=n.ahead)
    for(h in 1:n.ahead){
      aux<-NULL
      e0<-Phi2[neutR,,h]/irf[neutR,neutS,1] # shocks are vectorized, no loop here
      for(i in 1:bigK){ # loop over varibles / responses
        idx<-c(1:(n.ahead-h+1))
        aux<-rbind(aux,matrix(irf[i,neutS,idx],nrow=bigK,ncol=length(idx),byrow=TRUE)*e0)
      }
      dim(aux)<-c(bigK,bigK,length(idx));aux<-aperm(aux,c(2,1,3))
      Phi2[,,(h:n.ahead)]<-Phi2[,,(h:n.ahead),drop=FALSE]-aux
    }
    IRF_store[irep,,,] <- Phi2
    setTxtProgressBar(pb, irep)
  }
  #---------------------compute posterior----------------------------------------#
  imp_posterior <- array(NA, dim=c(bigK,bigK,n.ahead,7))
  dimnames(imp_posterior)[[1]] <- colnames(xglobal)
  dimnames(imp_posterior)[[2]] <- colnames(xglobal)
  dimnames(imp_posterior)[[3]] <- 1:n.ahead
  dimnames(imp_posterior)[[4]] <- c("low25","low16","low05","median","high75","high84","high95")
  imp_posterior[,,,"low25"] <- apply(IRF_store,c(2,3,4),quantile,.25,na.rm=TRUE)
  imp_posterior[,,,"low16"] <- apply(IRF_store,c(2,3,4),quantile,.16,na.rm=TRUE)
  imp_posterior[,,,"low05"] <- apply(IRF_store,c(2,3,4),quantile,.05,na.rm=TRUE)
  imp_posterior[,,,"median"]<- apply(IRF_store,c(2,3,4),quantile,.50,na.rm=TRUE)
  imp_posterior[,,,"high75"]<- apply(IRF_store,c(2,3,4),quantile,.75,na.rm=TRUE)
  imp_posterior[,,,"high84"]<- apply(IRF_store,c(2,3,4),quantile,.84,na.rm=TRUE)
  imp_posterior[,,,"high95"]<- apply(IRF_store,c(2,3,4),quantile,.95,na.rm=TRUE)
  # other stuff
  A         <- apply(A_large,c(2,3),median)
  Fmat      <- apply(F_large,c(2,3,4),median)
  Ginv      <- apply(Ginv_large,c(2,3),median)
  Smat      <- apply(S_large,c(2,3),median)
  Sigma_u   <- Ginv%*%Smat%*%t(Ginv)
  struc.obj <- list(A=A,Fmat=Fmat,Ginv=Ginv,Smat=Smat)
  model.obj <- list(xglobal=xglobal,plag=plag)
  #--------------------------------- prepare output----------------------------------------------------------------------#
  out <- structure(list("posterior"   = imp_posterior,
                        "struc.obj"   = struc.obj,
                        "model.obj"   = model.obj), 
                   class="bgvar.irf")
  if(save.store){
    out$IRF_store = IRF_store
  }
  if(verbose) cat(paste("\nSize of irf object: ", format(object.size(out),unit="MB")))
  end.irf <- Sys.time()
  diff.irf <- difftime(end.irf,start.irf,units="mins")
  mins.irf <- round(diff.irf,0); secs.irf <- round((diff.irf-floor(diff.irf))*60,0)
  if(verbose) cat(paste("\nNeeded time for impulse response analysis: ",mins.irf," ",ifelse(mins.irf==1,"min","mins")," ",secs.irf, " ",ifelse(secs.irf==1,"second.","seconds.\n"),sep=""))
  return(out)
}

#' @name .divisors
#' @noRd
.divisors <- function (n,div) {
  div <- round(div)
  for(dd in div:1){
    if(n%%div==0) break else div<-div-1
  }
  return(div)
}