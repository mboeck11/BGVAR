#' @export
"irf" <- function(x, n.ahead=24, shockinfo=NULL, quantiles=NULL, expert=NULL, verbose=TRUE){
  UseMethod("irf", x)
}

#' @name irf
#' @title Impulse Response Function
#' @description This function calculates three alternative ways of dynamic responses, namely generalized impulse response functions (GIRFs) as in Pesaran and Shin (1998), orthogonalized impulse response functions using a Cholesky decomposition and finally impulse response functions given a set of user-specified sign restrictions.
#' @export
#' @usage irf(x, n.ahead=24, shockinfo=NULL, quantiles=NULL, 
#'     expert=NULL, verbose=TRUE)
#' @param x Object of class \code{bgvar}.
#' @param n.ahead Forecasting horizon.
#' @param shockinfo Dataframe with additional information about the nature of shocks. Depending on the \code{ident} argument, the dataframe has to be specified differently. In order to get a dummy version for each identification scheme use \code{\link{get_shockinfo}}.
#' @param quantiles Numeric vector with posterior quantiles. Default is set to compute median along with 68\%/80\%/90\% confidence intervals.
#' @param expert Expert settings, must be provided as list. Default is set to \code{NULL}.\describe{
#' \item{\code{MaxTries}}{ Numeric specifying maximal number of tries for finding a rotation matrix with sign-restrictions. Attention: setting this number very large may results in very long computational times. Default is set to \code{MaxTries=100}.}
#' \item{\code{save.store}}{ If set to \code{TRUE} the full posterior of both, impulses responses and rotation matrices, are returned. Default is set to \code{FALSE} in order to save storage.}
#' \item{\code{use_R}}{ Boolean whether IRF computation should fall back on \code{R} version, otherwise \code{Rcpp} version is used.}
#' \item{\code{applyfun}}{ In case \code{use_R=TRUE}, this allows for user-specific apply function, which has to have the same interface than \code{lapply}. If \code{cores=NULL} then \code{lapply} is used, if set to a numeric either \code{parallel::parLapply()} is used on Windows platforms and \code{parallel::mclapply()} on non-Windows platforms.}
#' \item{\code{cores}}{ Numeric specifying the number of cores which should be used, also \code{all} and \code{half} is possible. By default only one core is used.}
#' }
#' @param verbose If set to \code{FALSE} it suppresses printing messages to the console.
#' @return Returns a list of class \code{bgvar.irf} with the following elements: \describe{
#' \item{\code{posterior}}{ Four-dimensional array (K times n.ahead times number of shocks times Q) that contains Q quantiles of the posterior distribution of the impulse response functions.}
#' \item{\code{shockinfo}}{ Dataframe with details on identification specification.}
#' \item{\code{rot.nr}}{ In case identification is based on sign restrictions (i.e., \code{ident="sign"}), this provides the number of rotation matrices found for the number of posterior draws (save*save_thin).}
#' \item{\code{struc.obj}}{ List object that contains posterior quantitites needed when calculating historical decomposition and structural errors via \code{hd.decomp}.\describe{
#' \item{\code{A}}{ Median posterior of global coefficient matrix.}
#' \item{\code{Ginv}}{ Median posterior of matrix \code{Ginv}, which describes contemporaneous relationships between countries.}
#' \item{\code{S}}{ Posterior median of matrix with country variance-covariance matrices on the main diagonal.}
#' \item{\code{Rmed}}{ Posterior rotation matrix if \code{ident="sign"}.}
#' }}
#' \item{\code{model.obj}}{ List object that contains model-specific information, in particular\describe{
#' \item{\code{xglobal}}{ Data of the model.}
#' \item{\code{lags}}{ Lag specification of the model.}
#' }}
#' \item{\code{IRF_store}}{ Four-dimensional array (K times n.ahead times number of shock times draws) which stores the whole posterior distribution. Exists only if \code{save.store=TRUE}.}
#' \item{\code{R_store}}{ Three-dimensional array (K times K times draws) which stores all rotation matrices. Exists only if \code{save.store=TRUE}.}
#' }
#' @author Maximilian Boeck, Martin Feldkircher, Florian Huber
#' @references 
#' Arias, J.E., Rubio-Ramirez, J.F, and D.F. Waggoner (2018) \emph{Inference Based on SVARs Identified with Sign and Zero Restrictions: Theory and Applications.} Econometrica Vol. 86(2), pp. 685-720.
#' 
#' D'Amico, S. and T. B. King (2017) \emph{What Does Anticipated Monetary Policy Do?} Federal Reserve Bank of Chicago Working paper series, Nr. 2015-10.
#' 
#' Pesaran, H.M. and Y. Shin (1998) \emph{Generalized impulse response analysis in linear multivariate models.} Economics Letters, Volume 58, Issue 1, p. 17-29.
#' @examples
#' \donttest{
#' oldpar <- par(no.readonly = TRUE)
#' # First example, a US monetary policy shock, quarterly data
#' library(BGVAR)
#' data(testdata)
#' # US monetary policy shock
#' model.eer<-bgvar(Data=testdata, W=W.test, draws=50, burnin=50, 
#'                  plag=1, prior="SSVS", eigen=TRUE)
#'
#' # generalized impulse responses
#' shockinfo<-get_shockinfo("girf")
#' shockinfo$shock<-"US.stir"; shockinfo$scale<--100
#' 
#' irf.girf.us.mp<-irf(model.eer, n.ahead=24, shockinfo=shockinfo)
#' 
#' # cholesky identification
#' shockinfo<-get_shockinfo("chol")
#' shockinfo$shock<-"US.stir"; shockinfo$scale<--100
#' 
#' irf.chol.us.mp<-irf(model.eer, n.ahead=24, shockinfo=shockinfo)
#' 
#' # sign restrictions
#' shockinfo <- get_shockinfo("sign")
#' shockinfo <- add_shockinfo(shockinfo, shock="US.stir", restriction=c("US.y","US.Dp"), 
#'                            sign=c("<","<"), horizon=c(1,1), scale=1, prob=1)
#' irf.sign.us.mp<-irf(model.eer, n.ahead=24, shockinfo=shockinfo)
#' 
#' # sign restrictions
#' shockinfo <- get_shockinfo("sign")
#' shockinfo <- add_shockinfo(shockinfo, shock="US.stir", restriction=c("US.y","US.Dp"), 
#' sign=c("<","<"), horizon=c(1,1), scale=1, prob=1)
#' irf.sign.us.mp<-irf(model.eer, n.ahead=24, shockinfo=shockinfo)
#' 
#' #' # sign restrictions with relaxed cross-country restrictions
#' shockinfo <- get_shockinfo("sign")
#' # restriction for other countries holds to 75\%
#' shockinfo <- add_shockinfo(shockinfo, shock="US.stir", restriction=c("US.y","EA.y","UK.y"), 
#'                            sign=c("<","<","<"), horizon=1, scale=1, prob=c(1,0.75,0.75))
#' shockinfo <- add_shockinfo(shockinfo, shock="US.stir", restriction=c("US.Dp","EA.Dp","UK.Dp"),
#'                            sign=c("<","<","<"), horizon=1, scale=1, prob=c(1,0.75,0.75))
#' irf.sign.us.mp<-irf(model.eer, n.ahead=20, shockinfo=shockinfo)
#' }
#' @seealso \code{\link{bgvar}}, \code{\link{get_shockinfo}}, \code{\link{add_shockinfo}}
#' @importFrom abind adrop abind
#' @importFrom stochvol sv_normal sv_beta sv_gamma
#' @importFrom RcppParallel RcppParallelLibs setThreadOptions defaultNumThreads
irf.bgvar <- function(x,n.ahead=24,shockinfo=NULL,quantiles=NULL,expert=NULL,verbose=TRUE){
  start.irf <- Sys.time()
  # get identification
  ident <- attr(shockinfo, "ident")
  if(is.null(ident)){
    ident <- "chol"
  }
  #--------------- checks ------------------------------------------------------------------------------------#
  if(!ident%in%c("chol","girf","sign")){
    stop("Please choose available identification scheme!")
  }
  if(is.null(shockinfo) && ident=="sign"){
    stop("Please provide 'shockinfo' argument.")
  }
  if(is.null(quantiles)){
    quantiles <- c(.05,.10,.16,.50,.84,.90,.95)
  }
  if(!is.numeric(quantiles)){
    stop("Please provide 'quantiles' as numeric vector.")
  }
  if(!is.null(shockinfo)){ # delete double entries
    shockinfo<-shockinfo[!duplicated(shockinfo),]
  }
  #-----------------------------------------------------------------------------------------------------------#
  if(verbose) cat("Start computing impulse response functions of Bayesian Global Vector Autoregression.\n\n")
  #------------------------------ get stuff -------------------------------------------------------#
  lags        <- x$args$lags
  pmax        <- max(lags)
  xglobal     <- x$xglobal
  Traw        <- nrow(xglobal)
  bigK        <- ncol(xglobal)
  bigT        <- Traw-pmax
  A_large     <- x$stacked.results$A_large
  F_large     <- x$stacked.results$F_large
  S_large     <- x$stacked.results$S_large
  Ginv_large  <- x$stacked.results$Ginv_large
  F.eigen     <- x$stacked.results$F.eigen
  thindraws   <- length(F.eigen) ### prior: draws
  Global      <- FALSE
  if(!is.null(shockinfo)) Global <- ifelse(any(shockinfo$global),TRUE,FALSE)
  Rmed        <- NULL
  rot.nr      <- NULL
  xdat        <- xglobal[(pmax+1):Traw,,drop=FALSE]
  varNames    <- colnames(xdat)
  cN          <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x)x[1]))
  vars        <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x)x[2]))
  N           <- length(cN)
  Q           <- length(quantiles)
  # expert settings
  expert.list <- list(MaxTries=100, save.store=FALSE, use_R=FALSE, applyfun=NULL, cores=NULL)
  if(!is.null(expert)){
    if(!(is.null(expert$cores) || is.numeric(expert$cores) || expert$cores%in%c("all","half"))){
      stop("Please provide the expert argument 'cores' in appropriate form. Please recheck.")
    }
    for(n in names(expert))
      expert.list[[n]] <- expert[[n]]
  }
  MaxTries    <- expert.list$MaxTries
  save.store  <- expert.list$save.store
  use_R       <- expert.list$use_R
  applyfun    <- expert.list$applyfun
  cores       <- expert.list$cores
  #---------------------------- identification schemes --------------------------------------------#
  if(ident=="chol")
  {
    if(verbose){
      cat("Identification scheme: Short-run identification via Cholesky decomposition.\n")
    }
    if(is.null(shockinfo)){
      shockinfo <- get_shockinfo("chol", nr_rows = length(varNames))
      shockinfo$shock <- varNames
    }
    if(!all(c("shock","scale")%in%colnames(shockinfo))){
      stop("Please provide appropriate dataframe for argument 'shockinfo'. Respecify.")
    }
    if(!all(shockinfo$shock%in%varNames)){
      stop("Please provide shock of 'shockinfo' only to variables available in the dataset used for estimation. Respecify.")
    }
    irf.fun  <- .irf.chol
    shock.nr <- nrow(shockinfo)
    # shock details
    shocks <- shocknames <- unique(shockinfo$shock)
    select_shocks <- NULL
    for(ss in 1:shock.nr) select_shocks <- c(select_shocks,which(shocks[ss] == varNames))
    scale <- shockinfo$scale[!duplicated(shockinfo$shock)]
    shock.cN  <- sapply(strsplit(shockinfo$shock,".",fixed=TRUE),function(x)x[1])
    shock.var <- sapply(strsplit(shockinfo$shock,".",fixed=TRUE),function(x)x[2])
    shock.idx <- list()
    for(cc in 1:N) shock.idx[[cc]] <- grep(cN[cc],varNames)
    shock.cidx <- cN%in%shock.cN
    if(Global){
      if(length(unique(shock.var[shockinfo$global])) != 1){
        stop("Please indicate global shock for same variables. Respecify.")
      }
      shock.nr <- shock.nr-(sum(shockinfo$global)-1)
      scale.new <- rep(1,shock.nr)
      shocknames <- shocks
      idx_global <- which(shockinfo$global)
      shocknames[idx_global[1]] <- paste0("Global.",unique(shock.var[shockinfo$global]))
      shocknames <- shocknames[-idx_global[c(2:length(idx_global))]]
      shock.global <- list()
      tt <- 1
      for(ss in 1:shock.nr){
        if(shockinfo$global[tt] == TRUE){
          shock.global[[shocknames[ss]]] <- varNames%in%shocks[shockinfo$global]
          scale.new[ss] <- scale[shockinfo$global][1]
          tt<-max(which(shockinfo$global))+1
        }else{
          shock.global[[shocknames[ss]]] <- varNames%in%shocks[tt]
          scale.new[ss] <- scale[tt]
          tt<-tt+1
        }
      }
      scale <- scale.new
    }
    shocklist = list(shock.idx=shock.idx,shock.cidx=shock.cidx,plag=pmax,MaxTries=MaxTries)
  }else if(ident=="girf")
  {
    if(verbose){
      cat("Identification scheme: Generalized impulse responses.\n")
    }
    if(is.null(shockinfo)){
      shockinfo <- get_shockinfo("girf", nr_rows = length(varNames))
      shockinfo$shock <- varNames
    }
    if(!all(c("shock","scale")%in%colnames(shockinfo))){
      stop("Please provide appropriate dataframe for argument 'shockinfo'. Respecify.")
    }
    if(!all(shockinfo$shock%in%varNames)){
      stop("Please provide shock of 'shockinfo' only to variables available in the dataset used for estimation. Respecify.")
    }
    if(!is.null(shockinfo)){
      shocks <- shocknames <- unique(shockinfo$shock)
      scale <- shockinfo$scale[!duplicated(shockinfo$shock)]
    }else{
      shocks <- shocknames <- varNames
      scale <- rep(1,length(shocks))
    }
    irf.fun  <- .irf.girf
    shock.nr <- length(shocks)
    select_shocks <- NULL
    for(ss in 1:shock.nr) select_shocks <- c(select_shocks,which(shocks[ss] == varNames))
    shock.idx <- list()
    for(cc in 1:N) shock.idx[[cc]] <- grep(cN[cc],varNames)
    shock.cidx <- rep(FALSE,N)
    if(Global){
      shock.var <- sapply(strsplit(shockinfo$shock,".",fixed=TRUE),function(x)x[2])
      if(length(unique(shock.var[shockinfo$global])) != 1){
        stop("Please indicate global shock for same variables. Respecify.")
      }
      shock.nr <- shock.nr-(sum(shockinfo$global)-1)
      scale.new <- rep(1,shock.nr)
      shocknames <- shocks
      idx_global <- which(shockinfo$global)
      shocknames[idx_global[1]] <- paste0("Global.",unique(shock.var[shockinfo$global]))
      shocknames <- shocknames[-idx_global[c(2:length(idx_global))]]
      shock.global <- list()
      tt <- 1
      for(ss in 1:shock.nr){
        if(shockinfo$global[tt] == TRUE){
          shock.global[[shocknames[ss]]] <- varNames%in%shocks[shockinfo$global]
          scale.new[ss] <- scale[shockinfo$global][1]
          tt<-max(which(shockinfo$global))+1
        }else{
          shock.global[[shocknames[ss]]] <- varNames%in%shocks[tt]
          scale.new[ss] <- scale[tt]
          tt<-tt+1
        }
      }
      scale <- scale.new
    }
    shocklist = list(shock.idx=shock.idx,shock.cidx=shock.cidx,plag=pmax,MaxTries=MaxTries)
  }else if(ident=="sign")
  {
    # --------------- checks -------------------------------------------------#
    if(!all(c("shock","restriction","sign","horizon","scale","prob")%in%colnames(shockinfo))){
      stop("Please provide columns 'shock', 'restriction', 'sign', 'horizon' and 'scal' in dataframe 'shockinfo'.")
    }
    if(!(all(shockinfo$shock%in%varNames) && all(shockinfo$restriction%in%varNames))){
      stop("Please provide in columns 'shock' and 'restriction' of 'shockinfo' only variable names available in the dataset used for estimation. Respecify.")
    }
    if(!any(shockinfo$sign%in%c(">","<","0","ratio.H","ratio.avg"))){
      stop("Misspecification in 'sign'. Only the following is allowed: <, >, 0, ratio.H, ratio.avg")
    }
    if(verbose){
      cat("Identification scheme: identification via sign-restriction.\n")
    }
    irf.fun<-.irf.sign.zero
    shocks <- shocknames <- unique(shockinfo$shock)
    shock.nr <- length(shocks)
    select_shocks <- NULL
    for(ss in 1:shock.nr) select_shocks <- c(select_shocks,which(shocks[ss] == varNames))
    shock.cN  <- unique(sapply(strsplit(shockinfo$shock,".",fixed=TRUE),function(x)x[1]))
    shock.var <- sapply(strsplit(shockinfo$shock,".",fixed=TRUE),function(x)x[2])
    shock.idx <- list()
    for(cc in 1:N) shock.idx[[cc]] <- grep(cN[cc],varNames)
    shock.cidx <- cN%in%shock.cN
    scale <- shockinfo$scale[!duplicated(shockinfo$shock)]
    if(Global){
      if(length(unique(shock.var[shockinfo$global])) != 1){
        stop("Please indicate global shock for same variables. Respecify.")
      }
      shock.nr <- shock.nr-(sum(shockinfo$global[!duplicated(shockinfo$shock)])-1)
      scale.new <- rep(1,shock.nr)
      shocknames <- shocks
      idx_global <- which(shockinfo$global)
      shocknames[idx_global[1]] <- paste0("Global.",unique(shock.var[shockinfo$global]))
      shocknames <- shocknames[-idx_global[c(2:length(idx_global))]]
      shock.global <- list()
      tt <- 1
      for(ss in 1:shock.nr){
        if(shockinfo$global[tt] == TRUE){
          shock.global[[shocknames[ss]]] <- varNames%in%shocks[shockinfo$global[!duplicated(shockinfo$shock)]]
          scale.new[ss] <- scale[shockinfo$global][1]
          tt<-max(which(shockinfo$global))+1
        }else{
          shock.global[[shocknames[ss]]] <- varNames%in%shockinfo$shock[tt]
          scale.new[ss] <- scale[tt]
          tt<-tt+1
        }
      }
      scale <- scale.new
    }
    # check zero/rationality
    if(any(shockinfo$sign%in%c("0","ratio.H","ratio.avg"))){
      for(ss in 1:length(shocks)){
        idx <- shockinfo$sign[grep(shocks[ss],shockinfo$shock)]%in%c("0","ratio.H","ratio.avg")
        if(!all(sapply(strsplit(shockinfo$restriction[grep(shocks[ss],shockinfo$shock)[idx]],".",fixed=TRUE),function(x)x[1])==shock.cN[ss])){
          stop("Please provide zero and rationality conditions only in same country as the origin of the shock.")
        }
      }
    }
    # adjust for rationality conditions
    if(any(shockinfo$sign=="ratio.H")){
      idx <- which(shockinfo$sign=="ratio.H")
      for(ii in idx){
        Kshock <- nrow(shockinfo)
        Mshock <- as.numeric(shockinfo$horizon[ii])
        shockinfo[(Kshock+1):(Kshock+2),] <- NA
        shockinfo$shock[(Kshock+1):nrow(shockinfo)] <- rep(shockinfo$shock[ii],2)
        shockinfo$restriction[(Kshock+1):nrow(shockinfo)] <- c(shockinfo$restriction[ii], strsplit(shockinfo$restriction[ii],"_")[[1]][1])
        shockinfo$sign[(Kshock+1):nrow(shockinfo)] <- c("0","-1")
        shockinfo$horizon[(Kshock+1):nrow(shockinfo)] <- c(1,Mshock)
        shockinfo$scale[(Kshock+1):nrow(shockinfo)] <- rep(shockinfo$scale[ii],2)
        shockinfo$prob[(Kshock+1):nrow(shockinfo)] <- rep(shockinfo$prob[ii],2)
      }
      shockinfo <- shockinfo[-idx,]
      rownames(shockinfo)<-seq(1,nrow(shockinfo))
    }
    if(any(shockinfo$sign=="ratio.avg")){
      idx <- which(shockinfo$sign=="ratio.avg")
      for(ii in idx){
        Kshock <- nrow(shockinfo)
        Mshock <- as.numeric(shockinfo$horizon[ii])
        shockinfo[(Kshock+1):(Kshock+Mshock),] <- NA
        shockinfo$shock[(Kshock+1):nrow(shockinfo)] <- rep(shockinfo$shock[ii],Mshock)
        shockinfo$restriction[(Kshock+1):nrow(shockinfo)] <- c(shockinfo$restriction[ii],rep(strsplit(shockinfo$restriction[ii],"_")[[1]][1],Mshock-1))
        shockinfo$sign[(Kshock+1):nrow(shockinfo)] <- c("0",rep(-1/(Mshock-1),Mshock-1))
        shockinfo$horizon[(Kshock+1):nrow(shockinfo)] <- seq(1,Mshock)
        shockinfo$scale[(Kshock+1):nrow(shockinfo)] <- rep(shockinfo$scale[ii],Mshock)
        shockinfo$prob[(Kshock+1):nrow(shockinfo)] <- rep(shockinfo$prob[ii],Mshock)
      }
      shockinfo <- shockinfo[-idx,]
      rownames(shockinfo)<-seq(1,nrow(shockinfo))
    }
    #---------------------------------------------------------------------------
    # create Scube and Zcube with signs
    sign.horizon   <- unique(shockinfo$horizon)
    sign.horizon   <- sort(sign.horizon, decreasing=FALSE)
    sign.shockvars <- unique(shockinfo$shock)
    H.restr        <- length(sign.horizon)
    N.restr        <- bigK*H.restr
    S.cube         <- matrix(0, N.restr, bigK) # sign restriction
    P.cube         <- matrix(0, N.restr, bigK) # probability of sign-restriction
    Z.cube         <- array(NA, c(N.restr, N.restr, bigK)) # zero restriction
    dimnames(S.cube)[[1]] <- dimnames(Z.cube)[[1]] <- dimnames(Z.cube)[[2]] <- dimnames(P.cube)[[1]] <- paste(rep(varNames,H.restr),".",rep(sign.horizon,each=bigK),sep="")
    dimnames(S.cube)[[2]] <- dimnames(Z.cube)[[3]] <- dimnames(P.cube)[[2]] <- varNames
    for(vv in 1:length(varNames)){
      Z.temp <- matrix(0, N.restr, N.restr)
      if(varNames[vv]%in%sign.shockvars){
        idx <- which(shockinfo$shock==varNames[vv])
        sign.restr <- shockinfo$restriction[idx]
        sign.signs <- shockinfo$sign[idx]
        sign.horiz <- shockinfo$horizon[idx]
        sign.probs <- shockinfo$prob[idx]
        
        s.point <- which(sign.signs=="<"|sign.signs==">")
        z.point <- seq(1,length(idx))[-s.point]
        
        # own shock: default is positive and for one period
        S.cube[paste(varNames[vv],".1",sep=""),varNames[vv]] <- 1
        P.cube[paste(varNames[vv],".1",sep=""),varNames[vv]] <- 1
        if(length(s.point)>0){
          for(ss in 1:length(s.point)){
            S.cube[paste(sign.restr[s.point[ss]],sign.horiz[s.point[ss]],sep="."),varNames[vv]] <- ifelse(sign.signs[s.point[ss]]=="<",-1,1)
            P.cube[paste(sign.restr[s.point[ss]],sign.horiz[s.point[ss]],sep="."),varNames[vv]] <- sign.probs[s.point[ss]]
          }
        }
        if(length(z.point)>0){
          for(zz in 1:length(z.point)){
            if(sign.signs[z.point[zz]]=="0"){
              grp <- which(sign.horiz[z.point[zz]] == sign.horizon)
              row <- (grp-1)*bigK+which(sign.restr[z.point[zz]]==varNames)
              Z.temp[row,row] <- 1
            }else{ # take row from above
              grp <- which(sign.horiz[z.point[zz]] == sign.horizon)
              col <- (grp-1)*bigK+which(sign.restr[z.point[zz]]==varNames)
              Z.temp[row,col] <- as.numeric(sign.signs[z.point[zz]])
            }
          }
        }
      }
      Z.cube[,,vv] <- Z.temp
    }
    # stuff needed for zero-restriction
    no.zero.restr <- rep(TRUE,N)
    shock.order   <- seq(bigK)
    for(cc in 1:N){
      idx <- shock.idx[[cc]]
      no.zero.restr[cc] <- ifelse(base::sum(abs(Z.cube[,,idx]))>0,FALSE,TRUE)
      shock.names <- names(sort(apply(Z.cube[,,idx], 3, function(x) base::sum(abs(x))), decreasing=TRUE))
      for(kk in 1:length(shock.names)) shock.order[idx[kk]] <- which(shock.names[kk]==varNames)
    }
    #---------------------------------------------------------------------------
    shocklist <- list(shock.idx=shock.idx,shock.cidx=shock.cidx,MaxTries=MaxTries,S.cube=S.cube,Z.cube=Z.cube,P.cube=P.cube,
                      shock.order=shock.order,shock.horz=sign.horizon,plag=pmax,no.zero.restr=no.zero.restr)
    rm(S.cube,Z.cube,P.cube)
  }
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
  #------------------------------ container -------------------------------------------------------#
  # initialize objects to save IRFs, HDs, etc.
  if(ident=="sign"){
    R_store       <- array(NA_real_, dim=c(bigK,bigK,thindraws), dimnames=list(colnames(xglobal),colnames(xglobal),NULL))
  }else{
    R_store <- NULL
  }
  IRF_store     <- array(NA_real_, dim=c(bigK,bigK,n.ahead+1,thindraws), dimnames=list(colnames(xglobal),paste0("shock_",colnames(xglobal)),seq(0,n.ahead),NULL))
  imp_posterior <- array(NA_real_, dim=c(bigK,n.ahead+1,shock.nr,Q))
  dimnames(imp_posterior) <- list(colnames(xglobal),seq(0,n.ahead),shocknames,paste0("Q",quantiles*100))
  #------------------------------ start computing irfs  ---------------------------------------------------#
  start.comp <- Sys.time()
  if(verbose) cat(paste("Start impulse response analysis on ", cores, " core",ifelse(cores>1,"s",""), " (",thindraws," stable draws in total).",sep=""),"\n")
  if(use_R)
  {
    #--------------------------------------------------------------
    # r-version
    counter <- numeric(length=thindraws)
    imp.obj <- applyfun(1:thindraws,function(irep){
      Ginv <- Ginv_large[,,irep]
      Fmat <- adrop(F_large[,,,irep,drop=FALSE],drop=4)
      Smat <- S_large[,,irep]
      imp.obj <- irf.fun(xdat=xdat,plag=pmax,n.ahead=n.ahead,Ginv=Ginv,Fmat=Fmat,Smat=Smat,shocklist=shocklist)
      if(verbose){
        if(ident=="sign"){
          if(!any(is.null(imp.obj$rot))){
            cat("\n",as.character(Sys.time()), "MCMC draw", irep, ": rotation found after ",imp.obj$icounter," tries", "\n")
          }else{
            cat("\n",as.character(Sys.time()), "MCMC draw", irep, ": no rotation found", "\n")
          }
        }
      }
      return(list(impl=imp.obj$impl,rot=imp.obj$rot,icounter=imp.obj$icounter))
    })
    for(irep in 1:thindraws){
      counter[irep] <- imp.obj[[1]]$icounter
      if(imp.obj[[1]]$icounter == MaxTries){
        imp.obj[[1]] <- NULL
      }else{
        IRF_store[,,,irep] <- imp.obj[[1]]$impl
        if(ident=="sign") R_store[,,irep] <- imp.obj[[1]]$rot
        imp.obj[[1]] <- NULL
      }
      if(irep %% 50 == 0) gc() # free up memory
    }
    rm(imp.obj)
  }else{ # cpp-version
    #--------------------------------------------------------------
    # adjust indexes due to different indexation (starting with zero in cpp)
    shocklist$shock.idx<-lapply(shocklist$shock.idx,function(l)l-1)
    shocklist$shock.horz <- shocklist$shock.horz-1
    shocklist$shock.order <- shocklist$shock.order-1
    # type
    type <- ifelse(ident=="chol",1,ifelse(ident=="girf",2,3))
    counter <- numeric(length=thindraws)
    save_rot <- ifelse(ident=="sign",TRUE,FALSE)
    # Rcpp::sourceCpp("./src/irf.cpp")
    # Rcpp::sourceCpp("/users/mboeck/documents/packages/bgvar/src/irf.cpp")
    temp = compute_irf(A_large=A_large,S_large=S_large,Ginv_large=Ginv_large,type=type,nhor=n.ahead+1,thindraws=thindraws,shocklist_in=shocklist,save_rot=save_rot,verbose=verbose)
    for(irep in 1:thindraws){
      counter[irep] <- temp$counter[irep,1]
      if(temp$counter[irep,1] == MaxTries){
        temp$irf[[1]] <- NULL; temp$rot[[1]] <- NULL
      }else{
        IRF_store[,,,irep] <- temp$irf[[1]]
        if(ident=="sign") R_store[,,irep] <- temp$rot[[1]]
        temp$irf[[1]] <- NULL; temp$rot[[1]] <- NULL
      }
      if(irep %% 50 == 0) gc() # free up memory
    }
    rm(temp)
    # transform back to R-version
    shocklist$shock.idx   = lapply(shocklist$shock.idx,function(l)l+1)
    shocklist$shock.horz  = shocklist$shock.horz+1
    shocklist$shock.order = shocklist$shock.order+1
  }
  end.comp <- Sys.time()
  diff.comp <- difftime(end.comp,start.comp,units="mins")
  mins <- round(diff.comp,0); secs <- round((diff.comp-floor(diff.comp))*60,0)
  if(verbose) cat(paste("\nImpulse response analysis took ",mins," ",ifelse(mins==1,"min","mins")," ",secs, " ",ifelse(secs==1,"second.\n","seconds.\n"),sep=""))
  #------------------------------ post processing  ---------------------------------------------------#
  # re-set IRF object in case we have found only a few rotation matrices
  if(ident=="sign")
  {
    idx <- which(counter!=MaxTries)
    rot.nr<-paste("For ", length(idx), " draws out of ", thindraws, " draws, a rotation matrix has been found.")
    if(length(idx)==0){
      stop("No rotation matrix found with imposed sign restrictions. Please respecify.")
    }
    if(verbose) cat(rot.nr)
    # subset posterior draws
    #IRF_store <- IRF_store[,,,idx,drop=FALSE]
    #R_store   <- R_store[,,idx,drop=FALSE]
    Ginv_large<-Ginv_large[,,idx,drop=FALSE]
    A_large   <- A_large[,,idx,drop=FALSE]
    S_large   <- S_large[,,idx,drop=FALSE]
    thindraws <- length(idx)
  }
  # Subset to shocks under consideration
  if(Global){
    impulse <- NULL
    for(ss in 1:shock.nr){
      temp <- apply(IRF_store[,shock.global[[ss]],,,drop=FALSE],c(1,3,4),sum)
      Mean <- temp[which(shock.global[[ss]])[1],1,]
      for(irep in 1:thindraws){
        temp[,,irep]<-(temp[,,irep]/Mean[irep])*scale[ss]
      }
      impulse <- abind(impulse,temp,along=4)
    }
    IRF_store <- aperm(impulse,c(1,4,2,3))
    dimnames(IRF_store)[[2]] <- names(shock.global)
  }else{
    IRF_store <- IRF_store[,select_shocks,,,drop=FALSE]
    for(ss in 1:shock.nr){
      Mean<-IRF_store[select_shocks[ss],ss,1,]
      for(irep in 1:thindraws){
        IRF_store[,ss,,irep]<-(IRF_store[,ss,,irep]/Mean[irep])*scale[ss]
      }
    }
  }
  # Normalization
  for(ss in 1:shock.nr){
    for(qq in 1:Q){
      imp_posterior[,,ss,qq] <- apply(IRF_store[,ss,,],c(1,2),quantile,quantiles[qq],na.rm=TRUE)
    }
  }
  # calculate objects needed for HD and struc shock functions later---------------------------------------------
  # median quantitities
  A       <- apply(A_large,c(1,2),median)
  Fmat    <- apply(F_large,c(1,2,3),median)
  Ginv    <- apply(Ginv_large,c(1,2),median)
  Smat    <- apply(S_large,c(1,2),median)
  Sigma_u <- Ginv%*%Smat%*%t(Ginv)
  if(ident=="sign")
  {
    imp.obj    <- try(irf.fun(xdat=xdat,plag=pmax,n.ahead=n.ahead,Ginv=Ginv,Fmat=Fmat,Smat=Smat,shocklist=shocklist),silent=TRUE)
    if(!is(imp.obj,"try-error")){
      Rmed<-imp.obj$rot
    }else{
      Rmed<-NULL
    }
  }
  struc.obj <- list(A=A,Fmat=Fmat,Ginv=Ginv,Smat=Smat,Rmed=Rmed)
  model.obj <- list(xglobal=xglobal,lags=lags)
  #--------------------------------- prepare output----------------------------------------------------------------------#
  out <- structure(list("posterior"   = imp_posterior,
                        "ident"       = ident,
                        "shockinfo"   = shockinfo,
                        "rot.nr"      = rot.nr,
                        "struc.obj"   = struc.obj,
                        "model.obj"   = model.obj), 
                   class="bgvar.irf")
  if(save.store){
    out$IRF_store = IRF_store
    out$R_store = R_store
  }
  if(verbose) cat(paste("\nSize of irf object: ", format(object.size(out),unit="MB")))
  end.irf <- Sys.time()
  diff.irf <- difftime(end.irf,start.irf,units="mins")
  mins.irf <- round(diff.irf,0); secs.irf <- round((diff.irf-floor(diff.irf))*60,0)
  if(verbose) cat(paste("\nNeeded time for impulse response analysis: ",mins.irf," ",ifelse(mins.irf==1,"min","mins")," ",secs.irf, " ",ifelse(secs.irf==1,"second.","seconds.\n"),sep=""))
  return(out)
}

#' @method print bgvar.irf
#' @export
print.bgvar.irf <- function(x, ...){
  cat("---------------------------------------------------------------------------------------")
  cat("\n")
  cat("Object contains impulse responses of object estimated with 'bgvar':")
  cat("\n")
  cat(paste0("Size of posterior containing impulse responses: ",dim(x$posterior)[[1]]," x ",dim(x$posterior)[[2]]," x ",dim(x$posterior)[[3]]," x ",dim(x$posterior)[[4]],"."))
  cat("\n")
  cat(paste0("Number of shocks identified: ",dim(x$posterior)[[3]],"."))
  cat("\n")
  cat("Identification scheme: ")
  if(x$ident=="chol"){
    cat("Short-run restrictions via Cholesky decomposition.")
  }else if(x$ident=="sign"){
    cat("Sign-restrictions.")
  }else if(x$ident=="girf"){
    cat("Generalized - no identification scheme employed.")
  }
  cat("\n")
  cat(paste0("Size ob object: ",format(object.size(x),unit="MB")))
  cat("\n")
  cat("---------------------------------------------------------------------------------------")
  
  return(invisible(x))
}

#' @name get_shockinfo
#' @title Create \code{shockinfo} argument 
#' @description Creates dummy \code{shockinfo} argument for appropriate use in  \code{irf} function.
#' @param ident Definition of identification scheme, either \code{chol}, \code{girf} or \code{sign}.
#' @param nr_rows Number of rows in the created dataframe.
#' @details Depending on the identification scheme a different \code{shockinfo} argument in the \code{irf} function is needed. To handle this convenient, an appropriate data.frame with is created with this function.
#' @usage get_shockinfo(ident="chol", nr_rows=1)
#' @seealso \code{\link{irf}}
#' @export
get_shockinfo <- function(ident="chol", nr_rows=1){
  if(ident == "chol"){
    df <- data.frame(shock=rep(NA,nr_rows),scale=rep(1,nr_rows),global=rep(FALSE,nr_rows))
    attr(df, "ident") <- "chol"
  }else if(ident == "girf"){
    df <- data.frame(shock=rep(NA,nr_rows),scale=rep(1,nr_rows),global=rep(FALSE,nr_rows))
    attr(df, "ident") <- "girf"
  }else if(ident=="sign"){
    df <- data.frame(shock=rep(NA,nr_rows),restriction=rep(NA,nr_rows),sign=rep(NA,nr_rows),
                     horizon=rep(NA,nr_rows),scale=rep(NA,nr_rows),prob=rep(NA,nr_rows),global=rep(NA,nr_rows))
    attr(df, "ident") <- "sign"
  }
  return(df)
}

#' @name add_shockinfo
#' @title Adding shocks to 'shockinfo' argument
#' @description Adds automatically rows to 'shockinfo' data.frame for appropriate use in \code{irf}.
#' @usage add_shockinfo(shockinfo=NULL, shock=NULL, restriction=NULL, sign=NULL, horizon=NULL, 
#' prob=NULL, scale=NULL, global=NULL, horizon.fillup=TRUE)
#' @param shockinfo Dataframe to append shocks. If \code{shockinfo=NULL} appropriate dataframe for sign-restrictions will be created.
#' @param shock String element. Variable of interest for structural shock. Only possible to add restrictions to one structural shock at a time.
#' @param restriction Character vector with variables that are supposed to be sign restricted.
#' @param sign Character vector with signs.
#' @param horizon Numeric vector with horizons to which restriction should hold. Set \code{horizon.fillup} to \code{FALSE} to just restrict one specific horizon.
#' @param prob Number between zero and one determining the probability with which restriction is supposed to hold.
#' @param scale Scaling parameter.
#' @param global If set to \code{TRUE}, shock is defined as global shock.
#' @param horizon.fillup Default set to \code{TRUE}, horizon specified up to given horizon. Otherwise just one specific horizon is restricted.
#' @details This is only possible for sign restriction, hence if \code{ident="sign"} in \code{get_shockinfo()}.
#' @seealso \code{\link{irf}}
#' @export
add_shockinfo <- function(shockinfo=NULL, shock=NULL, restriction=NULL, sign=NULL, horizon=NULL, prob=NULL, scale=NULL, global=NULL, horizon.fillup=TRUE){
  if(is.null(shockinfo)){
    shockinfo <- get_shockinfo(ident="sign")
  }
  if(is.null(shock)){
    stop("Please specify structural shock. This corresponds to the variable the shock is originating from.")
  }
  if(length(shock)>1){
    stop("Please only specify one structural shock at once.")
  }
  if(is.null(restriction) || is.null(sign)){
    stop("Please specify 'restriction' together with 'sign'.")
  }
  if(length(restriction)!=length(sign)){
    stop("Please provide the arguments 'restriction' and 'sign' with equal length. Please respecify.")
  }
  if(length(restriction)!=length(horizon)){
    if(length(horizon)!=1) stop("Please provide the argument 'horizon' either with length equal to one for all shocks or with an equal length of the restrictions.")
  }
  nr <- length(sign)
  if(!(is.null(restriction) && is.null(sign)) && is.null(horizon)){
    warning("No horizon specified, is set to one, i.e., a shock restriction on impact.")
    horizon <- rep(1,nr)
  }
  if(!any(sign%in%c(">","<","0","ratio.H","ratio.avg"))){
    stop("Misspecification in 'sign'. Only the following is allowed: <, >, 0, ratio.H, ratio.avg")
  }
  if(is.null(scale)){
    warning("Scaling is not specified, set positive.")
    scale <- rep(1,nr)
  }
  if(length(scale)==1) scale <- rep(scale,nr)
  scale <- sign(scale)
  if(length(unique(scale))>1){
    warning("Different scaling supplied. Set to default value: positive.")
    scale <- rep(1,nr)
  }
  if(is.null(prob)){
    warning("Restriction proabilities not specified, set to one.")
    prob <- rep(1,nr)
  }
  if(is.null(global)){
    global <- rep(FALSE,nr)
  }else{
    global <- rep(global,nr)
  }
  if(length(prob)==1) prob <- rep(prob,nr)
  if(length(prob)!=nr || length(scale)!=nr){
    stop("Please specify 'prob' or 'scale' with unit length for all restrictions or equal length than restriction.")
  }
  if(length(horizon)==1 && length(horizon)<nr){
    warning("Only one horizon specified, is used for all horizons.")
    horizon <- rep(horizon,nr)
  }
  for(irep in 1:nr){
    # if horizon is bigger than one
    idx_ratio  <- sign[irep] %in% c("ratio.H","ratio.avg")
    if(horizon[irep]>1 && !idx_ratio && horizon.fillup){
      repetition <- max(horizon[irep])
      # horizon <- c(unlist(sapply(horizon[idx_nr],seq)),horizon[idx_r])
    }else{
      repetition <- 1
    }
    # add to shockinfo
    nt<-ifelse(all(is.na(shockinfo)),0,nrow(shockinfo))
    for(nn in 1:repetition){
      shockinfo[nt+nn,] <- NA
      shockinfo$shock[nt+nn] <- shock
      shockinfo$restriction[nt+nn] <- restriction[irep]
      shockinfo$sign[nt+nn] <- sign[irep]
      if(repetition == 1){
        shockinfo$horizon[nt+nn] <- horizon[irep]
      }else{
        shockinfo$horizon[nt+nn] <- seq(1,horizon[irep])[nn]
      }
      shockinfo$prob[nt+nn] <- prob[irep]
      shockinfo$scale[nt+nn] <- scale[irep]
      shockinfo$global[nt+nn] <- global[irep]
    }
  }
  # delete duplicate lines
  shockinfo<-shockinfo[!duplicated(shockinfo),]
  return(shockinfo)
}
 