#' @name IRF
#' @title Impulse Response Functions
#' @description This function calculates three alternative ways of dynamic responses, namely generalized impulse response functions (GIRFs) as in Pesaran and Shin (1998), orthogonalized impulse response functions using a Cholesky decomposition and finally impulse response functions given a set of user-specified sign restrictions.
#' @usage IRF(obj, nhor=24, shock=NULL, sign.constr=NULL, save.store=FALSE, multithread=FALSE)
#' @param obj an object of class \code{bgvar}.
#' @param nhor forecasting horizon.
#' @param shock This is a list object. It should contain an entry labeled \code{var} that contains the name of the variable to be shocked. Also it should contain a list entry labeled \code{cN} that contains a character (or character vector) of the country (countries) in which the variable should be shocked. Finally it has to contain an entry labeled \code{ident} that is either \code{chol} if the shock is based on a short-run identification scheme done with the Cholesky decomposition or \code{girf} if generalized impulse responses should be calculated. In case impulses should be normalized (e.g., a +100bp increase in the interest rate), add another entry \code{scal} that contains a numeric value of the desired impact normalization.
#' @param sign.constr the user should submit a list containing the following entries \itemize{
#' \item{\code{shock1}}{ is a list object that defines sign restrictions for a particular shock. \itemize{
#' \item{\code{shockvar}}{ is a character vector containing the country and variable to shock separated by a dot. Example, "AT.ltir" (long-term interest rates in Austria).}
#' \item{\code{restrictions}}{ is a list containing the variables to restrict. Can have several sub-list restrictions, e.g., \code{sign.constr$shock1$restricionts$rest1=c("DE.y","AT.y")}, \code{sign.constr$shock1$restricionts$rest2=c("NL.p","AT.p")}, putting restrictions on GDP in Germany and Austria and a second set of restrictions on prices in the Netherlands and Austria.}
#' \item{\code{sign}}{ is a character vector of length of set of restrictions + 1, specifying the signs to impose. Use either \code{>}, \code{<} or \code{0}. The latter implements zero restrictions according to Arias et al. (2019). First entry is for the shock, say \code{AT.ltir} should go up, the following entries refer to the restrictions. \code{sign.constr$shock1$sign=c(">", "<", "<")} would impose \code{AT.ltir} to increase, and variables specified in \code{sign.constr$shock1$restricionts$rest1} and \code{sign.constr$shock1$restricionts$rest2} to decrease.}
#' \item{\code{rest.horz}}{ is a vector with same length as slot \code{sign} above and specifies the length of periods the restrictions are imposed. If \code{rest.horz} is 1, only impact restrictions are considered.}
#' \item{\code{constr}}{ is a vector with same length as slot \code{sign} above with elements lying between \code{0} and \code{1}. It specifies the percentage of countries for which cross-country restrictions have to hold. If no cross-country restrictions are supplied, set all elements of \code{constr} to 1.}
#' \item{\code{scal}}{ optional numeric in case impact normalization is desired.}}
#' #' \item{\code{MaxTries}}{ Optional numeric corresponding to the maximum tries to search for a rotation matrix that fulfills the user-specified restrictions. Default is set to 7500. After \code{MaxTries} unsuccessful tries the algorithm sets the impulse response for that specific posterior draw to \code{NA}.}
#' \item{\code{shock2}}{ define a second list with the same arguments as \code{shock1} to identify a second shock. Can be used iteratively to identify multiple shocks.}}}
#' @param save.store If set to \code{TRUE} the full posterior is returned. Default is set to \code{FALSE} in order to save storage.
#' @param multithread If set to \code{TRUE} parallel computing using the packages \code{\link{foreach}} and \code{\link{doParallel}}. Number of cores is set to maximum number of cores in the computer. This option is recommended when working with sign restrictions to speed up computations. Default is set to \code{FALSE} and thus no parallelization.
#' @return Returns a list of class \code{bgvar.irf} with the following elements: \itemize{
#' \item{\code{posterior}}{ is a four-dimensional array (K times nhor times nr. of shocks times 7) that contains 7 quantiles of the posterior distribution of the impulse response functions: the 50\% ("low25" and "high75"), the 68\% ("low16" and "high84") and the 90\% ("low05" and "high95") credible sets along with the posterior median ("median").}
#' \item{\code{rot.nr}}{ in case identification is based on sign restrictions (i.e., \code{ident="sign"}), this provides the number of rotation matrices found for the number of posterior draws (save*save_thin).}
#' \item{\code{shock}}{ in case identification is based on Cholesky decomposition (i.e. \code{ident="chol"}), this gives back the details of the identification specification.}
#' \item{\code{sign.constr}}{ in case identification is based on sign restrictions (i.e. \code{ident="sign"}), this gives back the set of sign restrictions specified by the user.}
#' \item{\code{ident}}{ character giving back the chosen identification scheme.}
#' \item{\code{struc.obj}}{ is a list object that contains posterior quantitites needed when calculating historical decomposition and structural errors via \code{hd.decomp}.\itemize{
#' \item{\code{A}}{ median posterior of global coefficient matrix.}
#' \item{\code{Ginv}}{ median posterior of matrix \code{Ginv}, which describes contemporaneous relationships between countries.}
#' \item{\code{S}}{ posterior median of matrix with country variance-covariance matrices on the main diagonal.}
#' \item{\code{xglobal}}{ dataset for whole model.}
#' \item{\code{plag}}{ specified lag length.}
#' \item{\code{Rmed}}{ posterior rotation matrix if \code{ident="sign"}.}
#' }}
#' \item{\code{model.obj}}{ is a list object that contains model-specific information, in particular\itemize{
#' \item{\code{xglobal}}{ used data of the model.}
#' \item{\code{plag}}{ used lag specification of the model.}
#' }}
#' \item{\code{IRF_store}}{ is a four-dimensional array (K times nhor times nr. of shock times saves) which stores the whole posterior distribution. Exists only if \code{save.irf.store=TRUE}.}
#' }
#' @author Maximilian Boeck, Martin Feldkircher, Florian Huber
#' @export 
#' @references 
#' Arias, J.E., Rubio-Ramirez, J.F, and D.F. Waggoner (2018) \emph{Inference Based on SVARs Identified with Sign and Zero Restrictions: Theory and Applications.} Econometrica Vol. 86(2), pp. 685-720.
#' 
#' D'Amico, S. and T. B. King (2017) \emph{What Does Anticipated Monetary Policy Do?} Federal Reserve Bank of Chicago Working paper series, Nr. 2015-10.
#' 
#' Pesaran, H.M. and Y. Shin (1998) \emph{Generalized impulse response analysis in linear multivariate models.} Economics Letters, Volume 58, Issue 1, p. 17-29.
#' @examples
#' \donttest{
#' set.seed(571)
#' # First example, a US monetary policy shock, quarterly data
#' library(BGVAR)
#' data(eerData)
#' model.ssvs.eer<-bgvar(Data=eerData,W=W.trade0012,saves=100,burns=100,plag=1,prior="SSVS",
#'                       eigen=TRUE)
#' # US monetary policy shock
#' shocks<-list();shocks$var="stir";shocks$cN<-"US";shocks$ident="chol";shocks$scal=-100
#' irf.chol.us.mp<-IRF(obj=model.ssvs.eer,shock=shocks,nhor=24)
#' # plots an impulse response function
#' plot(irf.chol.us.mp,resp="US")
#' 
#' # calculates generalized impulse response functions for the same shock as above
#' shocks$ident="girf"
#' irf.girf.ssvs<-IRF(obj=model.ssvs.eer,shock=shocks,nhor=24)
#' plot(irf.girf.ssvs,resp="US.y")
#' # Shock to first ordered variable yields same responses of Cholesky and GIRF
#' shocks<-list();shocks$var="y";shocks$cN<-"US";shocks$ident="chol";shocks$scal<--100
#' irf.chol<-IRF(model.ssvs.eer,shock=shocks,nhor=24)
#' shocks$ident<-"girf"
#' irf.girf<-IRF(model.ssvs.eer,shock=shocks,nhor=24)
#' matplot(cbind(irf.chol$posterior["US.y",,1,"median"],
#'               irf.girf$posterior["US.y",,1,"median"]),
#'         type="l",ylab="")
#' matplot(cbind(irf.chol$posterior["US.Dp",,1,"median"],
#'               irf.girf$posterior["US.Dp",,1,"median"]),
#'         type="l",ylab="")
#' matplot(cbind(irf.chol$posterior["EA.y",,1,"median"],
#'               irf.girf$posterior["EA.y",,1,"median"]),
#'         type="l",ylab="")
#' 
#' sign.constr<-list()
#' # the variable to shock, can be imposed for more than one country
#' sign.constr$shock1$shock<-c("US.stir") 
#' # but should be the same variable for all of them 
#' sign.constr$shock1$restrictions$res1<-c("US.y")
#' sign.constr$shock1$restrictions$res2<-c("US.Dp")
#' # first entry is for the shock, following entries for the restrictions 
#' # (ltir should go up, y and p go down)
#' sign.constr$shock1$sign<-c(">","<","<")
#' # nr. of time periods restrictions are imposed, first entry is for the shock, 
#' # following entries for the restrictions 
#' sign.constr$shock1$rest.horz<-c(1,1,1) 
#' # are constraints binding for all (1) countries specified for 
#' # at least 50\% of the countries (0.5), or 75\% (0.75)
#' sign.constr$shock1$constr<-c(1,1,1) 
#' # a minus 100 bp shock to long-term interest rates (on average)
#' sign.constr$shock1$scal=+100 
#' sign.constr$MaxTries<-200
#' irf.sign.us.mp<-IRF(obj=model.ssvs.eer,sign.constr=sign.constr,nhor=24)
#' plot(irf.sign.us.mp,resp=c("US"))
#' 
#' # second example, cross-country restrictions, multiple shocks and ECB country model
#' data(monthlyData);monthlyData$OC<-NULL
#' OE.weights <- list(EB=EA.weights)
#' model.ssvs<-bgvar(Data=monthlyData,W=W,saves=100,burns=100,plag=1,prior="SSVS",
#'                   thin=1,eigen=TRUE,OE.weights=OE.weights)
#' EA_countries <- c("AT", "BE", "DE","ES", "FI","FR", "IE", "IT", "NL", "PT","GR","SK")
#' # A simultaneous Cholesky shock to long-term interest rates in the euro area countries, 
#' # scaled to amount to -100 basis points (on average over the EA countries).
#' # Note that the ordering of the variables influences the response, the ordering is exactly as 
#' # in the country models, to use a different order you have re-estimate the model (by bgvar)
#' shocks<-list();shocks$var="ltir";shocks$cN<-EA_countries;shocks$ident="chol";shocks$scal=-100
#' irf.chol.ssvs<-IRF(obj=model.ssvs,shock=shocks,nhor=48)
#' plot(irf.chol.ssvs,resp=c("AT"))
#' # imposes sign restrictions on the cross-section and for a global shock (long-term interest rates)
#' sign.constr<-list()
#' # the variable to shock, can be imposed for more than one country
#' sign.constr$shock1$shock<-c(paste(EA_countries[-c(3,12)],".ltir",sep="")) 
#' # but should be the same variable for all of them
#' # restrictions (industrial production should decrease for selected countries)
#' sign.constr$shock1$restrictions$res1<-paste(EA_countries,".y",sep="")
#' # another set of restrictions (inflation  should decrease for selected countries) 
#' sign.constr$shock1$restrictions$res2<-paste(EA_countries,".p",sep="")
#' # first entry is for the shock, following entries for the restrictions 
#' # (ltir should go up, y and p go down) 
#' sign.constr$shock1$sign<-c(">","<","<") 
#' # nr. of time periods restrictions are imposed, first entry is for the shock, 
#' # following entries for the restrictions
#' sign.constr$shock1$rest.horz<-c(1,1,1) 
#' # are constraints binding for all (1) countries specified or for 
#' # at least 50\% of the countries (0.5), or 75\% (0.75)
#' sign.constr$shock1$constr<-c(1,0.5,0.5) 
#' # a minus 100 bp shock to long-term interest rates (on average)
#' sign.constr$shock1$scal=-100 
#' sign.constr$MaxTries<-200
#' irf.sign.ssvs<-IRF(obj=model.ssvs,nhor=48,sign.constr=sign.constr)
#' plot(irf.sign.ssvs,resp=c("AT"))
#'  
#' # Same example but using a local (German) shock and cross-country restrictions.
#' # Note that the ordering of the variables influences the response, 
#' # the ordering is exactly as in the country models, to use a different order you have re-estimate
#' # the model (by bgvar)
#' shocks<-list();shocks$var="ltir";shocks$cN<-EA_countries;shocks$ident="chol";shocks$scal=-100
#' irf.chol.ssvs<-IRF(obj=model.ssvs,shock=shocks,nhor=24)
#'  
#' # imposes sign restrictions on the cross-section and for a global shock (long-term interest rates)
#' sign.constr<-list()
#' sign.constr$shock1$shock<-c("DE.ltir") # the variable to shock, 
#' # can be imposed for more than one country
#' #but should be the same variable for all of them
#' # restrictions (industrial production should decrease for selected countries)
#' sign.constr$shock1$restrictions$res1<-paste(EA_countries,".y",sep="") 
#' # another set of restrictions (inflation  should decrease for selected countries)
#' sign.constr$shock1$restrictions$res2<-paste(EA_countries,".p",sep="")
#' # first entry is for the shock, following entries for the restrictions 
#' # (ltir should go up, y and p go down)
#' sign.constr$shock1$sign<-c(">","<","<") 
#' # nr. of time periods restrictions are imposed, 
#' # first entry is for the shock, following entries for the restrictions
#' sign.constr$shock1$rest.horz<-c(2,2,1) 
#' # are constraints binding for all (1) countries specified or for 
#' # at least 50\% of the countries (0.5), or 75\% (0.75)
#' sign.constr$shock1$constr<-c(1,0.5,0.5) 
#' # a minus 100 bp shock to long-term interest rates (on average)
#' sign.constr$shock1$scal=-100
#' sign.constr$MaxTries<-200
#' irf.sign.ssvs<-IRF(obj=model.ssvs,nhor=24,sign.constr=sign.constr)
#'  
#' # Example with zero restriction (Arias et al., 2018) and 
#' # rationality conditions (D'Amico and King, 2017).
#' data("eerDataspf")
#' model.ssvs.eer<-bgvar(Data=eerDataspf,W=W.trade0012.spf,saves=300,burns=300,
#'                       plag=1,prior="SSVS",thin=1,eigen=TRUE)
#' sign.constr<-list()
#' sign.constr$shock1$shock<-"US.stir_t+4"
#' sign.constr$shock1$restrictions$res1<-"US.Dp_t+4"
#' sign.constr$shock1$restrictions$res2<-"US.stir"
#' sign.constr$shock1$restrictions$res3<-"US.y_t+4"
#' # rationality condition: US.stir_t+4 on impact is equal to average of
#' # IRF of US.stir between horizon 1 and 4 (defined with rest.horz, but as period 5!)
#' sign.constr$shock1$restrictions$res4<-"US.stir_t+4"
#' # rationality condition: US.Dp_t+4 on impact is equal to H-step ahead IRF of US.Dp in 
#' # horizon 4 (defined with rest.horz, but as period 5!)
#' sign.constr$shock1$restrictions$res5<-"US.Dp_t+4"
#' # rationality condition: US.y_t+4 on impact is equal to H-step ahead IRF of US.y in 
#' # horizon 4 (defined with rest.horz, but as period 5!)
#' sign.constr$shock1$restrictions$res6<-"US.y_t+4"
#' sign.constr$shock1$sign<-c(">","<","0","<","ratio.avg","ratio.H","ratio.H")
#' sign.constr$shock1$rest.horz<-c(1,1,1,1,5,5,5)
#' sign.constr$shock1$constr<-c(1,1,1,1,1,1,1)
#' sign.constr$shock1$scal=0.1
#' sign.constr$MaxTries<-100
#' irf.sign<-IRF(obj=model.ssvs.eer,nhor=20,sign.constr=sign.constr)
#' par(mfrow=c(4,1),mar=c(5.1,4.1,4.1,2.1))
#' # rationality condition: US.stir_t+4 on impact is equal to average of IRF of 
#' # US.stir between horizon 1 and 4
#' matplot(cbind(irf.sign$posterior["US.stir_t+4",,1,"median"],
#'               irf.sign$posterior["US.stir",,1,"median"]),
#'         type="l",ylab="",main="stir")
#' abline(h=mean(irf.sign$posterior["US.stir",2:5,1,"median"]))
#' abline(v=c(1,5),lty=3,col="grey")
#' # rationality condition: US.y_t+4 on impact is equal to H-step ahead IRF of US.y in horizon 4
#' matplot(cbind(irf.sign$posterior["US.y_t+4",,1,"median"],
#'               irf.sign$posterior["US.y",,1,"median"]),
#'         type="l",ylab="",main="y")
#' abline(h=irf.sign$posterior["US.y_t+4",1,1,"median"])
#' abline(v=5,lty=3,col="grey")
#' # rationality condition: US.Dp_t+4 on impact is equal to H-step ahead IRF of US.Dp in horizon 4
#' matplot(cbind(irf.sign$posterior["US.Dp_t+4",,1,"median"],
#'               irf.sign$posterior["US.Dp",,1,"median"]),
#'         type="l",ylab="",main="Dp")
#' abline(h=irf.sign$posterior["US.Dp_t+4",1,1,"median"])
#' abline(v=5,lty=3,col="grey")
#' par(mar=rep(0,4))
#' plot("1",type="n",axes=FALSE)
#' legend("center",c("expectation","actual"),lty=1:2,col=c("black","red"),bty="n",ncol=2)
#' }
#' @importFrom abind adrop
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel detectCores
IRF <- function(obj,nhor=24,shock=NULL,sign.constr=NULL,save.store=FALSE,multithread=FALSE){
  start.irf <- Sys.time()
  if(!inherits(obj, "bgvar")) {stop("Please provide a `bgvar` object.")}
  cat("\nStart computing impulse response functions of Bayesian Global Vector Autoregression.\n\n")
  #------------------------------ get stuff -------------------------------------------------------#
  plag        <- obj$args$plag
  xglobal     <- obj$xglobal
  Traw        <- nrow(xglobal)
  bigK        <- ncol(xglobal)
  bigT        <- Traw-plag
  A_large     <- obj$stacked.results$A_large
  F_large     <- obj$stacked.results$F_large
  S_large     <- obj$stacked.results$S_large
  Ginv_large  <- obj$stacked.results$Ginv_large
  F.eigen     <- obj$stacked.results$F.eigen
  thinsaves   <- length(F.eigen) ### prior: saves
  Global      <- FALSE
  x           <- xglobal[(plag+1):Traw,,drop=FALSE]
  #------------------------------ user checks  ---------------------------------------------------#
  # checks general
  if(is.null(shock)&&is.null(sign.constr)){
    stop("Please provide either the list objects shock OR sign.constr to run GIRF / Cholesky or sign restriction based impulse response analysis.")
  }
  if(!is.null(shock)&&!is.null(sign.constr)){
    stop("Please provide either the list objects shock OR sign.constr to run GIRF / Cholesky or sign restriction based impulse response analysis.")
  }
  # checks cholesky/girf
  if(!is.null(shock)){
    if(!all(c("var","ident","cN")%in%names(shock))){
      stop("Please provide list entries shock, ccode and ident for the object shock. See the help files for more details.")
    }
    if(!shock$ident%in%c("girf","chol")){
      stop("Please provide identification scheme.")
    }
    shockvar <- shock$var
    cN       <- shock$cN
    ss       <- paste(cN,shockvar,sep=".")
    if(!all(ss%in%colnames(xglobal))){
      stop("Not all variables you have selected to shock are contained in the model. Please re-specify.")
    }
    if(length(ss)>1){
      Global=TRUE
    }
    ident    <- shock$ident
    scal     <- shock$scal
    if(is.null(scal)) scal <- 1
    ss<-list(ss)
    Rmed <- NULL
  }
  # checks identification via sign restrictions
  if(!is.null(sign.constr)){
    ident<-"sign"
    # check MaxTries, if no MaxTries, set it to 7500
    MaxTries <- 7500
    shock.nr <- length(sign.constr)
    if(!is.null(sign.constr$MaxTries)){
      MaxTries <- sign.constr$MaxTries
      sign.constr$MaxTries <- NULL
    }
    shock.nr<-length(sign.constr)
    # check whether sign restriction list is correctly specified, for each shock have to specify
    tt<-all(sapply(lapply(sign.constr,names),function(x) all(c("shock","sign","constr","restrictions")%in%x)))
    if(!tt){
      stop("For each shock (i.e., first layer in the list), please specify lists named shock, sign, constr and restrictions. See the details and examples in the manual.")
    }
    # check scaling, if no scaling, set it to 1
    for(kk in 1:shock.nr){
      sign.constr[[kk]]$scal<-ifelse(is.null(sign.constr[[kk]]$scal),1,sign.constr[[kk]]$scal)
    }
    scal <- unlist(sapply(sign.constr,function(x) x$scal))
    ss   <- lapply(sign.constr,function(x) x$shock)
    if(length(sign.constr[[1]]$shock)>1){ # in case you have specified a global shock
      Global=TRUE
    }else{
      Global=FALSE
    }
  }else{# in case we do not have sign restrictions, we can only identify one shock meaningfully
    shock.nr<-1
  }
  #------------------------------ assign irf function  ---------------------------------------------------#
  if(ident=="sign"){
    irf<-.irf.sign.zero
    # first check whether all elements of sign restriction list have been specified
    res_len<-unlist(lapply(sign.constr, function(l){
      if(is.list(l$restrictions)){return(length(l$restrictions))}else{return(0)}
    }))
    cc   <- lapply(ss,function(x)strsplit(x,".",fixed=TRUE))
    varNames <- colnames(xglobal)
    vars <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x)x[2]))
    cN   <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x)x[1]))
    cc.list <- list()
    for(ii in 1:length(cc)){
      temp <- c()
      for(iii in 1:length(cc[[ii]])){
        temp <- c(temp,cc[[ii]][[iii]][1])
      }
      cc.list[[ii]] <- temp
    }
    ccode <- unique(unlist(cc.list))
    sign.constr.new <- list()
    sign.constr.new.counter <- 1
    if(!any(ccode%in%cN)){
      stop("One of the shocks you specified is for a country not available in the data.")
    }
    for(m in 1:shock.nr){
      aux<-sign.constr[[m]]
      if(!all(c("shock","restrictions","sign","rest.horz","constr")%in%names(aux))){
        stop("Please specify for the object 'sign.contr' slots labeled 'shock', 'restrictions', 'sign', 'rest.horz' and 'constr'. A slot 'scal' is optional. See the help files for more information.")
      }
      mm<-1+length(aux$restrictions) # one for the shock, rest for the restrictions
      if(length(aux$sign)!=mm || length(aux$rest.horz)!=mm || length(aux$constr)!=mm){
        stop("Please specify M +1 signs, rest.horz and constr with M denoting the nr. of restrictions.")
      }
      if(!any(unlist(aux$restrictions)%in%varNames)){
        stop("Please restrict variables available in the dataset. Respecify.")
      }
      #-------------------------------------------------------
      if(res_len[m]>0){
        crosscountrylength <- unlist(lapply(sign.constr[[m]]$restrictions,length))
        sign.constr[[m]]$restrictions <- as.character(unlist(sign.constr[[m]]$restrictions))
        var<-unlist(lapply(strsplit(sign.constr[[m]]$restrictions,".",fixed=TRUE),function(y) y[[2]]))
        sign.constr[[m]]$sign<-c(sign.constr[[m]]$sign[1],
                                 rep(sign.constr[[m]]$sign[-1],
                                     times=crosscountrylength))
        sign.constr[[m]]$rest.horz<-c(sign.constr[[m]]$rest.horz[1],
                                      rep(sign.constr[[m]]$rest.horz[-1],
                                          times=crosscountrylength))
        sign.constr[[m]]$constr<-c(sign.constr[[m]]$constr[1],
                                   rep(sign.constr[[m]]$constr[-1],
                                       times=crosscountrylength))
      }
      if(length(cc.list[[m]])>1){
        sign.constr[[m]]$restrictions <- paste(ccode,".",rep(sign.constr[[m]]$restrictions,
                                                              each=length(cc.list[[m]])),sep="")
        sign.constr[[m]]$sign <- paste(ccode,".",rep(sign.constr[[m]]$sign,each=length(cc.list[[m]])),sep="")
        sign.constr[[m]]$rest.horz <- paste(ccode,".",rep(sign.constr[[m]]$rest.horz,
                                                           each=length(cc.list[[m]])),sep="")
        sign.constr[[m]]$constr <- paste(ccode,".",rep(sign.constr[[m]]$constr,
                                                        each=length(cc.list[[m]])),sep="")
        ### make single shocks
        for(cc in 1:length(cc.list[[m]])){
          temp              <- list()
          cntry             <- cc.list[[m]][cc]
          head              <- grep(cntry,sign.constr[[m]]$shock)
          body_rest         <- grep(paste("^",cntry,sep=""),sign.constr[[m]]$restrictions)
          body_else         <- grep(paste("^",cntry,sep=""),sign.constr[[m]]$sign)
          temp$shock        <- sign.constr[[m]]$shock[head]
          temp$restrictions <- gsub(paste("^",cntry,"\\.",sep=""),"",
                                    sign.constr[[m]]$restrictions[body_rest])
          temp$sign         <- gsub(paste("^",cntry,"\\.",sep=""),"",sign.constr[[m]]$sign[body_else])
          temp$rest.horz    <- as.numeric(gsub(paste("^",cntry,"\\.",sep=""),"",
                                               sign.constr[[m]]$rest.horz[body_else]))
          temp$constr       <- as.numeric(gsub(paste("^",cntry,"\\.",sep=""),"",
                                               sign.constr[[m]]$constr[body_else]))
          sign.constr.new[[sign.constr.new.counter]] <- temp
          sign.constr.new.counter<-sign.constr.new.counter+1
        }
      }else{
        sign.constr.new[[m]] <- sign.constr[[m]]
      }
    }
    sign.constr <- sign.constr.new
    names(sign.constr) <- paste("shock",seq(1,length(sign.constr)),sep="")
    for(m in 1:length(sign.constr)){
      if(any(sign.constr[[m]]$sign=="ratio.H"|sign.constr[[m]]$sign=="ratio.avg")){
        ratios <- which(sign.constr[[m]]$sign=="ratio.H"|sign.constr[[m]]$sign=="ratio.avg")
        temp.restr <- sign.constr[[m]]$restrictions[ratios-1]
        temp.signs <- sign.constr[[m]]$sign[ratios]
        temp.horiz <- sign.constr[[m]]$rest.horz[ratios]
        temp.constr <- sign.constr[[m]]$constr[ratios]
        temp.list  <- list(restrictions=c(),
                           sign=c(),
                           rest.horz=c(),
                           constr=c())
        for(ii in 1:length(ratios)){
          var_exp  <- temp.restr[ii]
          var_real <- strsplit(temp.restr[ii],"_")[[1]][1]
          if(temp.signs[ii]=="ratio.avg"){
            horizons <- seq(2,temp.horiz[ii]) # not on impact
          }else if(temp.signs[ii]=="ratio.H"){
            horizons <- temp.horiz[ii]
          }
          #### 0 restriction on expectations
          temp.list$restrictions <- c(temp.list$restrictions,var_exp)
          temp.list$sign         <- c(temp.list$sign,"0")
          temp.list$rest.horz    <- c(temp.list$rest.horz,1)
          temp.list$constr       <- c(temp.list$constr,temp.constr[ii])
          #### -1 restriction indicates sum
          if(length(horizons)==1){
            temp.list$restrictions <- c(temp.list$restrictions,var_real)
            temp.list$sign         <- c(temp.list$sign,"-1")
            temp.list$rest.horz    <- c(temp.list$rest.horz,horizons)
            temp.list$constr       <- c(temp.list$constr,temp.constr[ii])
          } else if(length(horizons)>1){ ### -1/H restriction indicates average
            temp.list$restrictions <- c(temp.list$restr,rep(var_real,length(horizons)))
            share <- as.character(-1/length(horizons))
            temp.list$sign         <- c(temp.list$sign,rep(share,length(horizons)))
            temp.list$rest.horz    <- c(temp.list$rest.horz,horizons)
            temp.list$constr       <- c(temp.list$constr,rep(temp.constr[ii],length(horizons)))
          }
        }
        # kick old stuff
        sign.constr[[m]]$restrictions <- sign.constr[[m]]$restrictions[-c(ratios-1)]
        sign.constr[[m]]$sign <- sign.constr[[m]]$sign[-ratios]
        sign.constr[[m]]$rest.horz <- sign.constr[[m]]$rest.horz[-ratios]
        sign.constr[[m]]$constr <- sign.constr[[m]]$constr[-ratios]
        # get new stuff in
        sign.constr[[m]]$restrictions <- c(sign.constr[[m]]$restrictions,temp.list$restrictions)
        sign.constr[[m]]$sign <- c(sign.constr[[m]]$sign,temp.list$sign)
        sign.constr[[m]]$rest.horz <- c(sign.constr[[m]]$rest.horz,temp.list$rest.horz)
        sign.constr[[m]]$constr <- c(sign.constr[[m]]$constr,temp.list$constr)
      }else if(any(sign.constr[[m]]$sign=="0")){
        zeros <- which(sign.constr[[m]]$sign=="0")
        temp.restr  <- sign.constr[[m]]$restrictions[zeros-1]
        temp.signs  <- sign.constr[[m]]$sign[zeros]
        temp.horiz  <- sign.constr[[m]]$rest.horz[zeros]
        temp.constr <- sign.constr[[m]]$constr[zeros]
        temp.list   <- list(restrictions=c(),
                            sign=c(),
                            rest.horz=c(),
                            constr=c())
        for(ii in 1:length(zeros)){
          var_zero <- temp.restr[ii]
          horizons <- seq(1,temp.horiz[ii])
          #### 0 restriction on expectations
          temp.list$restrictions <- c(temp.list$restrictions,rep(var_zero,length(horizons)))
          temp.list$sign         <- c(temp.list$sign,rep("0",length(horizons)))
          temp.list$rest.horz    <- c(temp.list$rest.horz,horizons)
          temp.list$constr       <- c(temp.list$constr,rep(temp.constr[ii],length(horizons)))
        }
        # kick old stuff
        sign.constr[[m]]$restrictions <- sign.constr[[m]]$restrictions[-c(zeros-1)]
        sign.constr[[m]]$sign <- sign.constr[[m]]$sign[-zeros]
        sign.constr[[m]]$rest.horz <- sign.constr[[m]]$rest.horz[-zeros]
        sign.constr[[m]]$constr <- sign.constr[[m]]$constr[-zeros]
        # get new stuff in
        sign.constr[[m]]$restrictions <- c(sign.constr[[m]]$restrictions,temp.list$restrictions)
        sign.constr[[m]]$sign <- c(sign.constr[[m]]$sign,temp.list$sign)
        sign.constr[[m]]$rest.horz <- c(sign.constr[[m]]$rest.horz,temp.list$rest.horz)
        sign.constr[[m]]$constr <- c(sign.constr[[m]]$constr,temp.list$constr)
      }else if(any(sign.constr[[m]]$sign==">"|sign.constr[[m]]$sign=="<")){
        signs <- which(sign.constr[[m]]$sign==">"|sign.constr[[m]]$sign=="<")
        if(any(signs==1)){
          temp.restr <- c(sign.constr[[m]]$shock,sign.constr[[m]]$restrictions[signs-1])
          first <- TRUE
        }else{
          temp.restr <- sign.constr[[m]]$restrictions[signs-1]
        }
        temp.signs <- sign.constr[[m]]$sign[signs]
        temp.horiz <- sign.constr[[m]]$rest.horz[signs]
        temp.constr<- sign.constr[[m]]$constr[signs]
        temp.list  <- list(restrictions=c(),
                           sign=c(),
                           rest.horz=c(),
                           constr=c())
        for(ii in 1:length(signs)){
          var_horz <- seq(1,temp.horiz[ii])
          temp.list$restrictions <- c(temp.list$restrictions,
                                      rep(temp.restr[ii],length(var_horz)))
          temp.list$sign         <- c(temp.list$sign,
                                      rep(temp.signs[ii],length(var_horz)))
          temp.list$rest.horz    <- c(temp.list$rest.horz,var_horz)
          temp.list$constr       <- c(temp.list$constr,
                                      rep(temp.constr[ii],length(var_horz)))
        }
        # kick old stuff
        sign.constr[[m]]$restrictions <- sign.constr[[m]]$restrictions[-c(signs-1)]
        sign.constr[[m]]$sign         <- sign.constr[[m]]$sign[-signs]
        sign.constr[[m]]$rest.horz    <- sign.constr[[m]]$rest.horz[-signs]
        sign.constr[[m]]$constr       <- sign.constr[[m]]$constr[-signs]
        # get new stuff in
        if(first){
          sign.constr[[m]]$restrictions <- c(sign.constr[[m]]$restrictions,
                                             temp.list$restrictions[-1])
        }else{
          sign.constr[[m]]$restrictions <- c(sign.constr[[m]]$restrictions,
                                             temp.list$restrictions)
        }
        sign.constr[[m]]$sign         <- c(sign.constr[[m]]$sign,
                                           temp.list$sign)
        sign.constr[[m]]$rest.horz    <- c(sign.constr[[m]]$rest.horz,
                                           temp.list$rest.horz)
        sign.constr[[m]]$constr       <- c(sign.constr[[m]]$constr,
                                           temp.list$constr)
      }
    }
    strg.list <- NULL
  }else if(ident=="chol"){
    cat("Identification scheme: Short-run identification via Cholesky decomposition.\n")
    cat(paste("Structural shock of interest: ", ss[[1]],".\n",sep=""))
    irf <- .irf.chol
    MaxTries<-str<-sign.constr<-Srots<-rot.nr<-NULL
  }else{
    cat("Identification scheme: Generalized impulse responses.\n")
    cat(paste("Structural shock of interest: ", ss[[1]],".\n",sep=""))
    irf <- .irf.girf
    MaxTries<-str<-sign.constr<-rot.nr<-NULL
  }
  
  # initialize objects to save IRFs, HDs, etc.
  R_store       <- array(NA, dim=c(thinsaves,bigK,bigK))
  IRF_store     <- array(NA, dim=c(thinsaves,bigK,nhor+1,shock.nr));dimnames(IRF_store)[[2]] <- colnames(xglobal)
  imp_posterior <- array(NA, dim=c(bigK,nhor+1,shock.nr,7))
  dimnames(imp_posterior)[[1]] <- colnames(xglobal)
  dimnames(imp_posterior)[[2]] <- 1:(nhor+1)
  dimnames(imp_posterior)[[3]] <- paste("shock",1:shock.nr,sep="_")
  dimnames(imp_posterior)[[4]] <- c("low25","low16","low05","median","high75","high84","high95")
  #------------------------------ start computing irfs  ---------------------------------------------------#
  start.comp <- Sys.time()
  if(multithread){
    numCores <- detectCores()
    registerDoParallel(cores=numCores)
    cat(paste("Start impulse response analysis on ", numCores, " cores", " (",thinsaves," stable draws in total).",sep=""),"\n")
    imp.obj <-foreach(irep=1:thinsaves) %dopar% {
      Ginv <- Ginv_large[irep,,]
      Fmat <- adrop(F_large[irep,,,,drop=FALSE],drop=1)
      Smat <- S_large[irep,,]
      imp.obj    <- irf(x=x,plag=plag,nhor=nhor,Ginv=Ginv,Fmat=Fmat,Smat=Smat,shock=shock,sign.constr=sign.constr,Global=Global,
                        MaxTries=MaxTries,shock.nr=shock.nr)
      if(!is.null(sign.constr)){
        if(!any(is.null(imp.obj$rot))){
          cat("\n",as.character(Sys.time()), "MCMC draw", irep, ": rotation found after ",imp.obj$icounter," tries", "\n")
        }else{
          cat("\n",as.character(Sys.time()), "MCMC draw", irep, ": no rotation found", "\n")
        }
      }
      return(list(impl=imp.obj$impl,rot=imp.obj$rot))
    }
    for(irep in 1:thinsaves){
      R_store[irep,,]    <- ifelse(is.null(imp.obj[[irep]]$rot),NA,imp.obj[[irep]]$rot)
      IRF_store[irep,,,] <- imp.obj[[irep]]$impl
    }
  }else{
    cat(paste("Start impulse response analysis on single core", " (",thinsaves," stable draws in total).",sep=""),"\n")
    for(irep in 1:thinsaves){
      Ginv <- Ginv_large[irep,,]
      Fmat <- adrop(F_large[irep,,,,drop=FALSE],drop=1)
      Smat <- S_large[irep,,]
      imp.obj    <- irf(x=x,plag=plag,nhor=nhor,Ginv=Ginv,Fmat=Fmat,Smat=Smat,shock=shock,sign.constr=sign.constr,Global=Global,
                        MaxTries=MaxTries,shock.nr=shock.nr)
      if(!is.null(sign.constr)){
        if(!any(is.na(imp.obj$rot))){
          cat("\n",as.character(Sys.time()), "MCMC draw", irep, ": rotation found after ",imp.obj$icounter," tries", "\n")
        }else{
          cat("\n",as.character(Sys.time()), "MCMC draw", irep, ": no rotation found", "\n")
        }
      }
      R_store[irep,,]    <- ifelse(is.null(imp.obj$rot),NA,imp.obj$rot)
      IRF_store[irep,,,] <- imp.obj$impl
    }
  }
  end.comp <- Sys.time()
  diff.comp <- difftime(end.comp,start.comp,units="mins")
  mins <- round(diff.comp,0); secs <- round((diff.comp-floor(diff.comp))*60,0)
  cat(paste("\nImpulse response analysis took ",mins," ",ifelse(mins==1,"min","mins")," ",secs, " ",ifelse(secs==1,"second.","seconds.\n"),sep=""))
  #------------------------------ post processing  ---------------------------------------------------#
  # re-set IRF object in case we have found only a few rotation matrices
  if(ident=="sign"){
    idx<-which(!is.na(apply(IRF_store,1,sum)))
    rot.nr<-paste("For ", length(idx), " draws out of ", thinsaves, " draws, a rotation matrix has been found.")
    if(length(idx)==0){
      stop("No rotation matrix found with imposed sign restrictions. Please respecify.")
    }
    print(rot.nr)
    # subset posterior draws
    IRF_store <- IRF_store[idx,,,,drop=FALSE]
    Ginv_large<-Ginv_large[idx,,,drop=FALSE]
    A_large   <- A_large[idx,,,drop=FALSE]
    S_large   <- S_large[idx,,,drop=FALSE]
    R_store   <- R_store[idx,,,drop=FALSE]
    thinsaves <- length(idx)
  }
  # Normalization
  if(thinsaves>0){
    if(!is.null(scal)){
      # do this for each identified shock (if sign restrictions, there might be more than one shock identified)
      if(Global){
        for(z in 1:shock.nr){
          Mean<-colMeans(IRF_store[,ss[[z]],1,z]) #in case we have multiple country shocks, scale as average on impact
          for(irep in 1:thinsaves){
            IRF_store[irep,,,z]<-(IRF_store[irep,,,z]/Mean[irep])*scal[z]
          }
        }
      }else{ #if no global shock has been carried out, but a local shock
        for(z in 1:shock.nr){
          Mean<-IRF_store[,ss[[z]],1,z]
          for(irep in 1:thinsaves){
            IRF_store[irep,,,z]<-(IRF_store[irep,,,z]/Mean[irep])*scal[z]
          }
        }
      }
    }
    
    for(i in 1:shock.nr){
      imp_posterior[,,i,"low25"]  <- apply(IRF_store[,,,i],c(2,3),quantile,0.25,na.rm=TRUE)
      imp_posterior[,,i,"low16"]  <- apply(IRF_store[,,,i],c(2,3),quantile,0.16,na.rm=TRUE)
      imp_posterior[,,i,"low05"]  <- apply(IRF_store[,,,i],c(2,3),quantile,0.05,na.rm=TRUE)
      imp_posterior[,,i,"median"] <- apply(IRF_store[,,,i],c(2,3),median,na.rm=TRUE)
      imp_posterior[,,i,"high75"] <- apply(IRF_store[,,,i],c(2,3),quantile,0.75,na.rm=TRUE)
      imp_posterior[,,i,"high84"] <- apply(IRF_store[,,,i],c(2,3),quantile,0.84,na.rm=TRUE)
      imp_posterior[,,i,"high95"] <- apply(IRF_store[,,,i],c(2,3),quantile,0.95,na.rm=TRUE)
    }
  }
  # calculate objects needed for HD and struc shock functions later---------------------------------------------
  # median quantitities
  A       <- apply(A_large,c(2,3),median)
  Fmat    <- apply(F_large,c(2,3,4),median)
  Ginv    <- apply(Ginv_large,c(2,3),median)
  Smat    <- apply(S_large,c(2,3),median)
  Sigma_u <- Ginv%*%Smat%*%t(Ginv)
  if(ident=="sign"){
    imp.obj    <- try(irf(x=x,plag=plag,nhor=nhor,Ginv=Ginv,Fmat=Fmat,Smat=Smat,shock=shock,
                          sign.constr=sign.constr,Global=Global,MaxTries=MaxTries,shock.nr=shock.nr),silent=TRUE)
    if(!is(imp.obj,"try-error")){
      Rmed<-imp.obj$rot
    }else{
      Rmed<-NULL
    }
  }
  struc.obj <- list(A=A,Fmat=Fmat,Ginv=Ginv,Smat=Smat,xglobal=xglobal,plag=plag,Rmed=Rmed)
  model.obj <- list(xglobal=xglobal,plag=plag)
  #--------------------------------- prepare output----------------------------------------------------------------------#
  out <- structure(list("posterior"   = imp_posterior,
                        "rot.nr"      = rot.nr,
                        "shock"       = shock,
                        "sign.constr" = sign.constr,
                        "ident"       = ident,
                        "struc.obj"   = struc.obj,
                        "model.obj"   = model.obj), 
                   class="bgvar.irf")
  if(save.store){
    out$IRF_store = IRF_store
  }
  cat(paste("\nSize of IRF object: ", format(object.size(out),unit="MB")))
  end.irf <- Sys.time()
  diff.irf <- difftime(end.irf,start.irf,units="mins")
  mins.irf <- round(diff.irf,0); secs.irf <- round((diff.irf-floor(diff.irf))*60,0)
  cat(paste("\nNeeded time for impulse response analysis: ",mins.irf," ",ifelse(mins.irf==1,"min","mins")," ",secs.irf, " ",ifelse(secs.irf==1,"second.","seconds.\n"),sep=""))
  return(out)
}

#' @name IRF.cf
#' @title Counterfactual Analysis
#' @description Function to perform counterfactual analysis. It enables to neutralize the response of a specific variable to a given shock.
#' @usage IRF.cf(obj,shockvar,resp,nhor=24,save.store=FALSE)
#' @param obj an object of class \code{bgvar}.
#' @param shockvar structural shock of interest.
#' @param resp response variable to neutralize.
#' @param nhor forecasting horizon.
#' @param save.store If set to \code{TRUE} the full posterior is returned. Default is set to \code{FALSE} in order to save storage.
#' @author Maximilian Boeck, Martin Feldkircher
#' @examples 
#' \donttest{
#' library(BGVAR)
#' data(eerData)
#' model.ssvs.eer<-bgvar(Data=eerData,W=W.trade0012,saves=100,burns=100,plag=1,prior="SSVS",
#'                       eigen=TRUE)
#' }
#' # very time-consuming
#' \dontrun{
#' IRF.cf(model.ssvs.eer,shockvar="US.stir",resp="US.rer",nhor=24)
#' }
#' @importFrom stats quantile
#' @export
IRF.cf <- function(obj,shockvar,resp,nhor=24,save.store=FALSE){
  start.irf <- Sys.time()
  if(!inherits(obj, "bgvar")) {stop("Please provide a `bgvar.irf` object.")}
  cat("\nStart counterfactual analysis of Bayesian Global Vector Autoregression.\n\n")
  #----------------get stuff-------------------------------------------------------#
  plag        <- obj$args$plag
  xglobal     <- obj$xglobal
  bigK        <- ncol(xglobal)
  A_large     <- obj$stacked.results$A_large
  F_large     <- obj$stacked.results$F_large
  S_large     <- obj$stacked.results$S_large
  Ginv_large  <- obj$stacked.results$Ginv_large
  F.eigen     <- obj$stacked.results$F.eigen
  thinsaves   <- length(F.eigen)
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
  cat(paste("Shock of interest: ",shockvar,".\n",sep=""))
  cat(paste("Response to neutralize: ",resp,".\n",sep=""))
  #--------------compute-----------------------------------------------------------#
  cat("Start computing...\n")
  IRF_store     <- array(NA, dim=c(thinsaves,bigK,bigK,nhor))
  dimnames(IRF_store)[[2]] <- dimnames(IRF_store)[[3]] <- colnames(xglobal)
  pb <- txtProgressBar(min = 0, max = thinsaves, style = 3)
  for(irep in 1:thinsaves){
    Sigma_u <- Ginv_large[irep,,]%*%S_large[irep,,]%*%t(Ginv_large[irep,,])
    irf<-Phi2<- .impulsdtrf(B=adrop(F_large[irep,,,,drop=FALSE],drop=1),
                           smat=t(chol(Sigma_u)),nstep=nhor)
    for(h in 1:nhor){
      aux<-NULL
      e0<-Phi2[neutR,,h]/irf[neutR,neutS,1] # shocks are vectorized, no loop here
      for(i in 1:bigK){ # loop over varibles / responses
        idx<-c(1:(nhor-h+1))
        aux<-rbind(aux,matrix(irf[i,neutS,idx],nrow=bigK,ncol=length(idx),byrow=TRUE)*e0)
      }
      dim(aux)<-c(bigK,bigK,length(idx));aux<-aperm(aux,c(2,1,3))
      Phi2[,,(h:nhor)]<-Phi2[,,(h:nhor),drop=FALSE]-aux
    }
    IRF_store[irep,,,] <- Phi2
    setTxtProgressBar(pb, irep)
  }
  #---------------------compute posterior----------------------------------------#
  imp_posterior <- array(NA, dim=c(bigK,bigK,nhor,7))
  dimnames(imp_posterior)[[1]] <- colnames(xglobal)
  dimnames(imp_posterior)[[2]] <- colnames(xglobal)
  dimnames(imp_posterior)[[3]] <- 1:nhor
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
  struc.obj <- list(A=A,Fmat=Fmat,Ginv=Ginv,Smat=Smat,xglobal=xglobal,plag=plag)
  model.obj <- list(xglobal=xglobal,plag=plag)
  #--------------------------------- prepare output----------------------------------------------------------------------#
  out <- structure(list("posterior"   = imp_posterior,
                        "struc.obj"   = struc.obj,
                        "model.obj"   = model.obj), 
                   class="bgvar.irf")
  if(save.store){
    out$IRF_store = IRF_store
  }
  cat(paste("\nSize of IRF object: ", format(object.size(out),unit="MB")))
  end.irf <- Sys.time()
  diff.irf <- difftime(end.irf,start.irf,units="mins")
  mins.irf <- round(diff.irf,0); secs.irf <- round((diff.irf-floor(diff.irf))*60,0)
  cat(paste("\nNeeded time for impulse response analysis: ",mins.irf," ",ifelse(mins.irf==1,"min","mins")," ",secs.irf, " ",ifelse(secs.irf==1,"second.","seconds.\n"),sep=""))
  return(out)
}

#' @name plot.bgvar.irf
#' @title Plot predictions of bgvar
#' @description  Plots the predictions of an object of class \code{bgvar.predict}.
#' @param x an object of class \code{bgvar.irf}.
#' @param ... additional arguments.
#' @param resp specify a variable to plot predictions.
#' @param shock.nr specify shock to be plotted.
#' @param cumulative whether cumulative impulse response functions should be plotted. Default is set to \code{FALSE}.
#' @author Maximilian Boeck, Martin Feldkircher
#' @examples
#' \donttest{
#' library(BGVAR)
#' data(eerData)
#' model.ssvs.eer<-bgvar(Data=eerData,W=W.trade0012,saves=100,burns=100,plag=1,prior="SSVS")
#' # US monetary policy shock
#' shocks<-list();shocks$var="stir";shocks$cN<-"US";shocks$ident="chol";shocks$scal=-100
#' irf.chol.us.mp<-IRF(obj=model.ssvs.eer,shock=shocks,nhor=24)
#' # plots an impulse response function
#' plot(irf.chol.us.mp,resp="US.y")
#' }
#' @seealso \code{\link{IRF}}
#' @importFrom graphics abline matplot polygon segments
#' @importFrom stats rnorm
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
plot.bgvar.irf<-function(x, ...,resp,shock.nr=1,cumulative=FALSE){
  if(!inherits(x, "bgvar.irf")) {stop("Please provide a `bgvar.irf` object.")}
  if(length(shock.nr)!=1){stop("Please select only one shock.")}
  posterior <- x$posterior
  varNames  <- dimnames(posterior)[[1]]
  varAll    <- varNames
  cN        <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x) x[1]))
  vars      <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x) x[2]))
  
  if(!is.null(resp)){
    resp.p <- strsplit(resp,".",fixed=TRUE)
    resp.c <- sapply(resp.p,function(x) x[1])
    resp.v <- sapply(resp.p,function(x) x[2])
    if(!all(unique(resp.c)%in%cN)){
      stop("Please provide country names corresponding to the ones of the 'bgvar.irf' object.")
    }
    cN       <- cN[cN%in%resp.c]
    varNames <- lapply(cN,function(x)varNames[grepl(x,varNames)])
    if(all(!is.na(resp.v))){
      if(!all(unlist(lapply(resp,function(r)r%in%varAll)))){
        stop("Please provide correct variable names corresponding to the ones in the 'bgvar' object.")
      }
      varNames <- lapply(varNames,function(l)l[l%in%resp])
    }
    max.vars <- unlist(lapply(varNames,length))
  }else{
    varNames <- lapply(cN,function(cc) varAll[grepl(cc,varAll)])
  }
  for(cc in 1:length(cN)){
    rows <- max.vars[cc]/2
    if(rows<1) cols <- 1 else cols <- 2
    if(rows%%1!=0) rows <- ceiling(rows)
    if(rows%%1!=0) rows <- ceiling(rows)
    par(mfrow=c(rows,cols),mar=bgvar.env$mar)
    for(kk in 1:max.vars[cc]){
      idx  <- grep(cN[cc],varAll)
      idx <- idx[varAll[idx]%in%varNames[[cc]]][kk]
      x<-posterior[idx,,shock.nr,c("low25","median","high75"),drop=TRUE] # first dimension is nr. of variables to be plotted
      y<-posterior[idx,,shock.nr,c("low16","median","high84"),drop=TRUE] # first dimension is nr. of variables to be plotted
      if(cumulative){x<-apply(x,2,cumsum);y<-apply(y,2,cumsum)}
      b <- range(x,y)
      b1<-b[1];b2<-rev(b)[1]
      matplot(x,type="l",col=c("black",bgvar.env$plot$col.50,"black"),yaxt="n",xaxt="n",
              lwd=bgvar.env$plot$lwd.line,ylab="",main=varAll[idx],cex.main=bgvar.env$plot$cex.main,
              cex.axis=bgvar.env$plot$cex.axis,cex.lab=bgvar.env$plot$cex.lab,lty=c(0,1,0),ylim=c(b1,b2))
      polygon(c(1:nrow(y),rev(1:nrow(y))),c(y[,1],rev(y[,3])),col=bgvar.env$plot$col.75,border=NA)
      polygon(c(1:nrow(x),rev(1:nrow(x))),c(x[,1],rev(x[,3])),col=bgvar.env$plot$col.68,border=NA)
      segments(x0=1,y0=0,x1=nrow(x),y1=0,col=bgvar.env$plot$col.zero,lty=bgvar.env$plot$lty.zero,lwd=bgvar.env$plot$lwd.zero)
      lines(x[,2],col=bgvar.env$plot$col.50,lwd=4)
      axis(2, at=seq(b1,b2,length.out=5), labels=format(seq(b1,b2,length.out=5),digits=1,nsmall=1),cex.axis=1.2,las=1)
      axisindex<-seq(1,nrow(x),by=4)
      axis(side=1, las=1,at=axisindex, labels=c(0:nrow(x))[axisindex], cex.axis=1.6,tick=FALSE)
      abline(v=axisindex,col=bgvar.env$plot$col.tick,lty=bgvar.env$plot$lty.tick)
    }
    if(cc<length(cN)) readline(prompt="Press enter for next country...")
  }
}

 