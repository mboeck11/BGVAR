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
