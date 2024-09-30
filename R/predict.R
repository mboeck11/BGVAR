#' @name predict
#' @title Predictions
#' @description A function that computes predictions and conditional predictions based on a object of class \code{bgvar}.
#' @details Predictions are performed up to an horizon of \code{n.ahead}. Note that conditional forecasts need a fully identified system. Therefore this function utilizes short-run restrictions via the Cholesky decomposition on the global solution of the variance-covariance matrix of the Bayesian GVAR.
#' @param object An object of class \code{bgvar}.
#' @param ... Additional arguments.
#' @param n.ahead Forecast horizon.
#' @param constr Matrix containing the conditional forecasts of size horizon times K, where horizon corresponds to the forecast horizon specified in \code{pred.obj}, while K is the number of variables in the system. The ordering of the variables have to correspond the ordering of the variables in the system. Rest is just set to NA.
#' @param constr_sd Matrix containing the standard deviations around the conditional forecasts. Must have the same size as \code{constr}.
#' @param quantiles Numeric vector with posterior quantiles. Default is set to compute median along with 68\%/80\%/90\% confidence intervals.
#' @param save.store If set to \code{TRUE} the full distribution is returned. Default is set to \code{FALSE} in order to save storage.
#' @param verbose If set to \code{FALSE} it suppresses printing messages to the console.
#' @return Returns an object of class \code{bgvar.pred} with the following elements \describe{
#' \item{\code{fcast}}{ is a K times n.ahead times Q-dimensional array that contains  Q quantiles of the posterior predictive distribution.}
#' \item{\code{xglobal}}{ is a matrix object of dimension T times N (T # of observations, K # of variables in the system).}
#' \item{\code{n.ahead}}{ specified forecast horizon.}
#' \item{\code{lps.stats}}{ is an array object of dimension K times 2 times n.ahead and contains the mean and standard deviation of the log-predictive scores for each variable and each forecast horizon.}
#' \item{\code{hold.out}}{ if \code{h} is not set to zero, this contains the hold-out sample.}
#' }
#' @examples
#' \donttest{
#' library(BGVAR)
#' data(testdata)
#' model.ssvs <- bgvar(Data=testdata,W=W.test,plag=1,draws=100,burnin=100,
#'                     prior="SSVS")
#' fcast <- predict(model.ssvs, n.ahead=8)
#' 
#' # conditional predictions
#' # et up constraints matrix of dimension n.ahead times K
#' constr <- matrix(NA,nrow=8,ncol=ncol(model.ssvs$xglobal))
#' colnames(constr) <- colnames(model.ssvs$xglobal)
#' constr[1:5,"US.Dp"] <- model.ssvs$xglobal[76,"US.Dp"]
#' 
#' # add uncertainty to conditional forecasts
#' constr_sd <- matrix(NA,nrow=8,ncol=ncol(model.ssvs$xglobal))
#' colnames(constr_sd) <- colnames(model.ssvs$xglobal)
#' constr_sd[1:5,"US.Dp"] <- 0.001
#' 
#' fcast_cond <- predict(model.ssvs, n.ahead=8, constr=constr, constr_sd=constr_sd)
#' }
#' @references 
#' Jarocinski, M. (2010) \emph{Conditional forecasts and uncertainty about forecasts revisions in vector autoregressions.} Economics Letters, Vol. 108(3), pp. 257-259.
#' 
#' Waggoner, D., F. and T. Zha (1999) \emph{Conditional Forecasts in Dynamic Multivariate Models.} Review of Economics and Statistics, Vol. 81(4), pp. 639-561.
#' @importFrom stats rnorm tsp sd
#' @author Maximilian Boeck, Martin Feldkircher, Florian Huber
#' @export
predict.bgvar <- function(object, ..., n.ahead=1, constr=NULL, constr_sd=NULL, quantiles=NULL, save.store=FALSE, verbose=TRUE){
  start.pred <- Sys.time()
  if(!inherits(object, "bgvar")) {stop("Please provide a `bgvar` object.")}
  # check if posterior draws are available
  if(object$args$thindraws == 0){
    cat("Computation of BGVAR has yielded no stable posterior draws!")
    return(invisible(object))
  }
  if(is.null(quantiles)){
    quantiles <- c(.05,.10,.16,.50,.84,.90,.95)
  }
  if(!is.numeric(quantiles)){
    stop("Please provide 'quantiles' as numeric vector.")
  }
  if(verbose) cat("Start computing predictions of Bayesian Global Vector Autoregression.\n\n")
  thindraws  <- object$args$thindraws
  lags       <- object$args$lags
  pmax       <- max(lags)
  xglobal    <- object$xglobal
  S_large    <- object$stacked.results$S_large
  F_large    <- object$stacked.results$F_large
  A_large    <- object$stacked.results$A_large
  Ginv_large <- object$stacked.results$Ginv_large
  F.eigen    <- object$stacked.results$F.eigen
  varNames   <- colnames(xglobal)
  cN         <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x) x[1]))
  vars       <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x) x[2]))
  N          <- length(cN)
  Traw       <- nrow(xglobal)
  bigT       <- Traw-pmax
  bigK       <- ncol(xglobal)
  cons       <- 1
  trend      <- ifelse(object$args$trend,1,0)
  Q          <- length(quantiles)
  flag_cond  <- FALSE
  #---------------------check conditional predictions--------------------------------#
  if(!is.null(constr)){
    if(!all(dim(constr)==c(n.ahead,bigK))){
      stop("Please respecify dimensions of 'constr'.")
    }
    if(!is.null(constr_sd)){
      if(!all(dim(constr_sd)==c(n.ahead,bigK))){
        stop("Please respecify dimensions of 'constr_sd'.")
      }
      constr_sd[is.na(constr_sd)] <- 0
    }else{
      constr_sd <- matrix(0,n.ahead,bigK)
    }
    flag_cond <- TRUE
  }
  #---------------------------------------------------------------------------------#
  varndxv <- c(bigK,cons+trend,pmax)
  nkk     <- (pmax*bigK)+cons+trend
  
  Yn <- xglobal
  Xn <- cbind(.mlag(Yn,pmax),1)
  Xn <- Xn[(pmax+1):Traw,,drop=FALSE]
  Yn <- Yn[(pmax+1):Traw,,drop=FALSE]
  if(trend) Xn <- cbind(Xn,seq(1,bigT))
  
  pred_store <- array(NA,dim=c(thindraws,bigK,n.ahead))
  # start loop here
  if(verbose){
    if(flag_cond)
      cat("Start computing conditional predictions...\n")
    else
      cat("Start computing predictions...\n")
  }
  for(irep in 1:thindraws){
    #Step I: Construct a global VC matrix Omega_t
    Ginv    <- Ginv_large[,,irep]
    Sig_t   <- Ginv%*%(S_large[,,irep])%*%t(Ginv)
    Sig_t   <- as.matrix(Sig_t)
    zt      <- Xn[bigT,]
    z1      <- zt
    Mean00  <- zt
    Sigma00 <- matrix(0,nkk,nkk)
    y2      <- NULL
    
    #gets companion form
    aux   <- .get_companion(A_large[,,irep],varndxv)
    Mm    <- aux$MM
    Jm    <- aux$Jm
    Jsigt <- Jm%*%Sig_t%*%t(Jm)
    # this is the forecast loop
    stop <- FALSE
    for(ih in 1:n.ahead){
      z1      <- Mm%*%z1
      Sigma00 <- Mm%*%Sigma00%*%t(Mm) + Jsigt
      chol_varyt <- try(t(chol(Sigma00[1:bigK,1:bigK])),silent=TRUE)
      if(is(chol_varyt,"matrix")){
        yf <- z1[1:bigK]+chol_varyt%*%rnorm(bigK,0,1)
      }
      if(is(chol_varyt,"try-error")){
        yf <- try(mvrnorm(1,mu=z1[1:bigK],Sigma00[1:bigK,1:bigK]),silent=TRUE)
      }
      if(is(yf,"try-error")){
        stop = TRUE
        break # break inner loop
      }
      y2 <- cbind(y2,yf)
    }
    if(stop){next} # continue outer loop
    pred_store[irep,,] <- y2
  }
  #----------do conditional forecasting -------------------------------------------#
  if(flag_cond){
    cond_store <- array(NA, c(thindraws, bigK, n.ahead))
    dimnames(cond_store)[[2]] <- varNames
    
    if(verbose) pb <- txtProgressBar(min = 0, max = thindraws, style = 3)
    for(irep in 1:thindraws){
      pred    <- pred_store[irep,,]
      Sigma_u <- Ginv_large[,,irep]%*%S_large[,,irep]%*%t(Ginv_large[,,irep])
      chol_varyt <- try(t(chol(Sigma_u)), silent=TRUE)
      if(is(chol_varyt,"try-error")) {next}
      irf     <- .impulsdtrf(B=adrop(F_large[,,,irep,drop=FALSE],drop=4),
                             smat=chol_varyt,nstep=n.ahead)
      
      temp <- as.vector(constr) + rnorm(bigK*n.ahead,0,as.vector(constr_sd))
      constr_use <- matrix(temp,n.ahead,bigK)
      
      v <- sum(!is.na(constr))
      s <- bigK * n.ahead
      r <- c(rep(0, v))
      R <- matrix(0, v, s)
      pos <- 1
      for(i in 1:n.ahead) {
        for(j in 1:bigK) {
          if(is.na(constr_use[i, j])) {next}
          r[pos] <- constr_use[i, j] - pred[j, i]
          for(k in 1:i) {
            R[pos, ((k - 1) * bigK + 1):(k * bigK)] <- irf[j,,(i - k + 1)]
          }
          pos <- pos + 1
        }
      }
      
      R_svd <- svd(R, nu=nrow(R), nv=ncol(R))
      U     <- R_svd[["u"]]
      P_inv <- diag(1/R_svd[["d"]])
      V1    <- R_svd[["v"]][,1:v]
      V2    <- R_svd[["v"]][,(v+1):s]
      eta   <- V1 %*% P_inv %*% t(U) %*% r + V2 %*% rnorm(s-v)
      eta   <- matrix(eta, n.ahead, bigK, byrow=TRUE)
      
      for(ih in 1:n.ahead) {
        temp <- matrix(0, bigK, 1)
        for(k in 1:ih) {
          temp <- temp + irf[, , (ih - k + 1)] %*% t(eta[k , , drop=FALSE])
        }
        cond_store[irep,,ih] <- pred[,ih,drop=FALSE] + temp
      }
      if(verbose) setTxtProgressBar(pb, irep)
    }
  }
  #--------------- compute posteriors ----------------------------------------------#
  imp_posterior<-array(NA,dim=c(bigK,n.ahead,Q), dimnames=list(varNames,seq(1,n.ahead),paste0("Q",quantiles*100)))
  for(qq in 1:Q){
    if(flag_cond){
      imp_posterior[,,qq] <- apply(cond_store,c(2,3),quantile,quantiles[qq],na.rm=TRUE)
    }else{
      imp_posterior[,,qq] <- apply(pred_store,c(2,3),quantile,quantiles[qq],na.rm=TRUE)
    }
  }
  #---------------------------------------------------------------------------------#
  hold.out <- object$args$hold.out
  if(hold.out>n.ahead) hold.out <- n.ahead
  yfull <- object$args$yfull
  if(hold.out>0){
    lps.stats <- array(0,dim=c(bigK,2,hold.out), dimnames=list(colnames(xglobal),c("mean","sd"),seq(1,hold.out)))
    lps.stats[,"mean",] <- apply(pred_store[,,1:hold.out],c(2:3),mean,na.rm=TRUE)
    lps.stats[,"sd",]   <- apply(pred_store[,,1:hold.out],c(2:3),sd,na.rm=TRUE)
    hold.out.sample<-yfull[(nrow(yfull)+1-hold.out):nrow(yfull),,drop=FALSE]
  }else{
    lps.stats<-NULL
    hold.out.sample<-NULL
  }
  #---------------------------------------------------------------------------------#
  rownames(xglobal)<-.timelabel(object$args$time)
  
  out <- structure(list(fcast=imp_posterior,
                        xglobal=xglobal,
                        n.ahead=n.ahead,
                        lps.stats=lps.stats,
                        hold.out.sample=hold.out.sample),
                   class="bgvar.pred")
  if(save.store){
    out$pred_store = pred_store
  }
  if(verbose) cat(paste("\n\nSize of object:", format(object.size(out),unit="MB")))
  end.pred <- Sys.time()
  diff.pred <- difftime(end.pred,start.pred,units="mins")
  mins.pred <- round(diff.pred,0); secs.pred <- round((diff.pred-floor(diff.pred))*60,0)
  if(verbose) cat(paste("\nNeeded time for computation: ",mins.pred," ",ifelse(mins.pred==1,"min","mins")," ",secs.pred, " ",ifelse(secs.pred==1,"second.","seconds.\n"),sep=""))
  return(out)
}

#' @method print bgvar.pred
#' @export
print.bgvar.pred <- function(x, ...){
  cat("---------------------------------------------------------------------------------------")
  cat("\n")
  cat("Object contains predictions of object estimated with 'bgvar':")
  cat("\n")
  cat(paste0("Size of posterior containing predictions: ",dim(x$fcast)[[1]]," x ",dim(x$fcast)[[2]]," x ",dim(x$fcast)[[3]],"."))
  cat("\n")
  cat(paste0("Size ob object: ",format(object.size(x),unit="MB")))
  cat("\n")
  cat("---------------------------------------------------------------------------------------")
  
  return(invisible(x))
}

#' @export
"lps" <- function(object){
  UseMethod("lps", object)
}

#' @name lps
#' @title Compute Log-Predictive Scores
#' @method lps bgvar.pred
#' @description  Computes and prints log-predictive score of an object of class \code{bgvar.predict}.
#' @param object An object of class \code{bgvar.predict}.
#' @param ... Additional arguments.
#' @return Returns an object of class \code{bgvar.lps}, which is a matrix of dimension h times K, whereas h is the forecasting horizon and K is the number of variables in the system.
#' @examples
#' \donttest{
#' library(BGVAR)
#' data(testdata)
#' model.ssvs.eer<-bgvar(Data=testdata,W=W.test,draws=100,burnin=100,
#'                       plag=1,prior="SSVS",eigen=TRUE,hold.out=8)
#' fcast <- predict(model.ssvs.eer,n.ahead=8,save.store=TRUE)
#' lps <- lps(fcast)
#' }
#' @author Maximilian Boeck, Martin Feldkircher
#' @importFrom stats dnorm
#' @export
lps.bgvar.pred <- function(object, ...){
  hold.out <- object$hold.out
  h        <- nrow(hold.out)
  K        <- ncol(hold.out)
  if(is.null(hold.out)){
    stop("Please submit a forecast object that includes a hold out sample for evaluation (set hold.out>0 when estimating the model with bgvar)!")
  }
  lps.stats  <- object$lps.stats
  lps.scores <- matrix(NA,h,K)
  for(i in 1:K){
    lps.scores[,i]<-dnorm(hold.out[,i],mean=lps.stats[i,"mean",],sd=lps.stats[i,"sd",],log=TRUE)
  }
  colnames(lps.scores)<-dimnames(lps.stats)[[1]]
  out <- structure(lps.scores, class="bgvar.lps")
  return(out)
}

#' @export
"rmse" <- function(object){
  UseMethod("rmse", object)
}

#' @name rmse
#' @title Compute Root Mean Squared Errors
#' @method rmse bgvar.pred
#' @description  Computes and prints root mean squared errors (RMSEs) of an object of class \code{bgvar.predict}.
#' @param object An object of class \code{bgvar.predict}.
#' @param ... Additional arguments.
#' @return Returns an object of class \code{bgvar.rmse}, which is a matrix of dimension h times K, whereas h is the forecasting horizon and K is the number of variables in the system.
#' @examples
#' \donttest{
#' library(BGVAR)
#' data(testdata)
#' model.ssvs.eer<-bgvar(Data=testdata,W=W.test,draws=100,burnin=100,
#'                       plag=1,prior="SSVS",eigen=TRUE,hold.out=8)
#' fcast <- predict(model.ssvs.eer,n.ahead=8,save.store=TRUE)
#' rmse <- rmse(fcast)
#' }
#' @author Maximilian Boeck, Martin Feldkircher
#' @importFrom stats dnorm
#' @export
rmse.bgvar.pred <- function(object, ...){
  hold.out <- object$hold.out
  h        <- nrow(hold.out)
  K        <- ncol(hold.out)
  if(is.null(hold.out)){
    stop("Please submit a forecast object that includes a hold out sample for evaluation (set hold.out>0 in fcast)!")
  }
  lps.stats   <- object$lps.stats
  rmse.scores <- matrix(NA,h,K)
  for(i in 1:K){
    rmse.scores[,i]<-sqrt((hold.out[,i]-lps.stats[i,"mean",])^2)
  }
  colnames(rmse.scores)<-dimnames(lps.stats)[[1]]
  out <- structure(rmse.scores, class="bgvar.rmse")
  return(out)
}

#' @method print bgvar.lps
#' @export
print.bgvar.lps<-function(x, ..., resp=NULL){
  h        <- dim(x)[1]
  varNames <- colnames(x)
  cN       <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(y)y[1]))
  vars     <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(y)y[2]))
  cntry    <- round(sapply(cN,function(y)mean(x[grepl(y,colnames(x))])),2)
  Ki       <- unlist(lapply(cN,function(x)length(grep(x,varNames))))
  bigK     <- length(vars)
  
  cat("---------------------------------------------------------------------------")
  cat("\n")
  cat("Log-Predictive Density Scores")
  cat("\n")
  cat(paste0("Available for hold.out times K: ",h," times ",bigK))
  cat("\n")
  cat(paste0("Total: ", round(sum(x),2)))
  cat("\n")
  cat("---------------------------------------------------------------------------")
  
  return(invisible(x))
}

#' @method print bgvar.rmse
#' @export
print.bgvar.rmse<-function(x, ..., resp=NULL){
  h        <- dim(x)[1]
  varNames <- colnames(x)
  cN       <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(y)y[1]))
  vars     <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(y)y[2]))
  cntry    <- round(sapply(cN,function(y)mean(x[grepl(y,colnames(x))])),2)
  Ki       <- unlist(lapply(cN,function(x)length(grep(x,varNames))))
  bigK     <- length(vars)
  
  cat("---------------------------------------------------------------------------")
  cat("\n")
  cat("Root-Mean Squared Error")
  cat("\n")
  cat(paste0("Available for hold.out times K: ",h," times ",bigK))
  cat("\n")
  cat(paste0("Total: ", round(sum(x),2)))
  cat("\n")
  cat("---------------------------------------------------------------------------")
  
  return(invisible(x))
}
