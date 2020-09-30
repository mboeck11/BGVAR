#' @name bgvar
#' @export
#' @title Estimation of Bayesian GVAR
#' @description Estimates a Bayesian GVAR with either the Stochastic Search Variable Selection (SSVS), the Minnesota prior (MN), or the Normal-Gamma prior. All specifications can be estimated with stochastic volatility.
#' @usage 
#' bgvar(Data, W, plag=1, draws=5000, burnin=5000, prior="NG", SV=TRUE, h=0, thin=1, 
#'       hyperpara=NULL, eigen=FALSE, variable.list=NULL, OE.weights=NULL, Wex.restr=NULL,
#'       Ex=NULL, trend=FALSE, save.country.store=FALSE, applyfun=NULL, cores=NULL, 
#'       verbose=TRUE)
#' @param Data Either a \itemize{
#' \item{\code{list object}}{ of length \code{N} that contains the data. Each element of the list refers to a country/entity. The number of columns (i.e., variables) in each country model can be different. The \code{T} rows (i.e., number of time observations), however, need to be the same for each country. Country and variable names are not allowed to contain a dot \code{.} (i.e., a dot) since this is our naming convention.}
#' \item{\code{matrix object}}{ of dimension \code{T} times \code{K}, with \code{K} denoting the sum of all endogenous variables of the system. The column names should consist of two parts, separated by a \code{.} (i.e., a dot). The first part should denote the country / entity name and the second part the name of the variable. Country and variable names are not allowed to contain a \code{.} (i.e., a dot).}
#' }
#' @param W An N times N weight matrix with 0 elements on the diagonal and row sums that sum up to unity or a list of weight matrices. 
#' @param plag Number of lags used (the same for domestic, exogenous and weakly exogenous variables.). Default set to \code{plag=1}.
#' @param draws Number of draws saved. Default set to \code{draws=5000}.
#' @param burnin Number of burn-ins. Default set to \code{burnin=5000}.
#' @param prior Either \code{SSVS} for the Stochastic Search Variable Selection prior, \code{MN} for the Minnesota prior or \code{NG} for the Normal-Gamma prior. See Details below.
#' @param SV If set to \code{TRUE}, models are fitted with stochastic volatility using the \code{stochvol} package. Due to storage issues, not the whole history of the \code{T} variance covariance matrices are kept, only the median. Consequently, the \code{BGVAR} package shows only one set of impulse responses (with variance covariance matrix based on mean sample point volatilities) instead of \code{T} sets. Specify \code{SV=FALSE} to turn SV off.
#' @param h Defines the hold-out sample. Default without hold-out sample, thus set to zero.
#' @param thin Is a thinning interval of the MCMC chain. As a rule of thumb, workspaces get large if draws/thin>500. Default set to \code{thin=1}.
#' @param hyperpara Is a list object that defines the hyperparameters when the prior is set to either \code{MN}, \code{SSVS} or \code{NG}. \itemize{
#' \item{\code{a_1}}{ is the prior hyperparameter for the inverted gamma prior (shape) (set a_1 = b_1 to a small value for the standard uninformative prior). Default is set to \code{a_1=0.01}.}
#' \item{\code{b_1}}{ is the prior hyperparameter for the inverted gamma prior (rate). Default is set to \code{b_1=0.01}.}
#' \item{\code{prmean}}{ Prior mean on the first lag of the autoregressive coefficients, standard value is \code{prmean=1} for non-stationary data. Prior mean for the remaining autoregressive coefficients automatically set to 0.}
#' \item{\code{bmu}}{ If \code{SV=TRUE}, this is the prior hyperparameter for the mean of the the mean of the log-volatilities. Default is \code{bmu=0}.}
#' \item{\code{Bmu}}{ If \code{SV=TRUE}, this is the prior hyperparameter for the variance of the mean of the log-volatilities. Default is \code{Bmu=100}.}
#' \item{\code{a0}}{ If \code{SV=TRUE}, this is the hyperparameter of the shape1 parameter for the Beta prior on the persistence parameter of the log-volatilities. Default is \code{a0=25}.}
#' \item{\code{b0}}{ If \code{SV=TRUE}, this is the hyperparameter of the shape2 parameter for the Beta prior on the persistence parameter of the log-volatilities. Default is \code{b0=1.5}.}
#' \item{\code{Bsigma}}{ If \code{SV=TRUE}, this is the hyperparameter for the Gamma prior on the variance of the log-volatilities. Default is set to \code{Bsigma=1}.}
#' \item{"MN"}{\itemize{
#'       \item{\code{shrink1}}{ Starting value of \code{shrink1}. Default set to 0.1.}
#'       \item{\code{shrink2}}{ Starting value of \code{shrink2}. Default set to 0.2.}
#'       \item{\code{shrink3}}{ Hyperparameter of \code{shrink3}. Default set to 100.}
#'       \item{\code{shrink4}}{ Starting value of \code{shrink4}. Default set to 0.1.}
#'       }}
#' \item{"SSVS"}{\itemize{
#'       \item{\code{tau0}}{ is the prior variance associated with the normal prior on the regression coefficients if a variable is NOT included (spike, tau0 should be close to zero).}
#'       \item{\code{tau1}}{ is the prior variance associated with the normal prior on the regression coefficients if a variable is  included (slab, tau1 should be large).}
#'       \item{\code{kappa0}}{ is the prior variance associated with the normal prior on the covariances if a covariance equals zero (spike, kappa0 should be close to zero).}
#'       \item{\code{kappa1}}{  is the prior variance associated with the normal prior on the covariances if a covariance is  unequal to zero (slab, kappa1 should be large).}
#'       \item{\code{p_i}}{ is the prior inclusion probability for each regression coefficient whether it is included in the model (default set to \code{p_i=0.5}).}
#'       \item{\code{q_ij}}{ is the prior inclusion probability for each covariance whether it is included in the model (default set to \code{q_ij=0.5}).}
#'       }}
#' \item{"NG":}{\itemize{
#'       \item{\code{e_lambda}}{ Prior hyperparameter for the Gamma prior on the lag-specific shrinkage components, standard value is \code{e_lambda=1.5}.}
#'       \item{\code{d_lambda}}{ Prior hyperparameter for the Gamma prior on the lag-specific shrinkage components, standard value is \code{d_lambda=1}.}
#'       \item{\code{a_start}}{ Parameter of the Normal-Gamma prior that governs the heaviness of the tails of the prior distribution. A value of a_start=1 would lead to the Bayesian LASSO. Default value differs per entity and set to \code{a_start=1/log(M)}, where \code{M} is the number of endogenous variables per entity.}
#'       \item{\code{sample_A}}{ If set to \code{TRUE} \code{a_start} is sampled.}
#'       }}
#'  }
#' @param eigen Set to TRUE if you want to compute the largest eigenvalue of the companion matrix for each posterior draw. If the modulus of the eigenvalue is significantly larger than unity, the model is unstable. Unstable draws exceeding an eigenvalue of one are then excluded. If \code{eigen} is set to a numeric value, then this corresponds to the maximum eigenvalue. The default is set to 1.05 (which excludes all posterior draws for which the eigenvalue of the companion matrix was larger than 1.05 in modulus).
#' @param variable.list In case \code{W} is a list of weight matrices, specify here which set of variables should be weighted by which weight matrix. See the help file on \code{getweights} for details. Default is \code{NULL}.
#' @param OE.weights Default value is set to \code{NULL}. Can be used to provide information of how to handle additional country models (other entities). Additional country models can be used to endogenously determine variables that are (weakly) exogenous for the majority of the other country models. As examples, one could think of an additional oil price model (see also Mohaddes and Raissi 2019) or a model for the joint euro area monetary policy (see also Georgiadis 2015; Feldkircher, Gruber and Huber (2020)). The data for these additional country models has to be contained in \code{Data}. The number of additional country models is unlimited. Each list entry of \code{OE.weights} has to be named similar to the name of the additional country model contained in \code{Data}. Each slot of \code{OE.weight} has to contain the following information: \itemize{
#' \item{\code{weights}}{ a vector of weights with names relating to the countries for which data should be aggregated. Can also relate to a subset of countries contained in the data.}
#' \item{\code{variables}}{ a vector of variables names that should be included in the additional country model. Variables that are not contained in the data slot of the extra country model are assumed to be weakly exogenous for the additional country model (aggregated with \code{weight}).}
#' \item{\code{exo}}{ a vector of variable names that should be fed into the other countries as (weakly) exogenous variables.}
#' }
#' @param Wex.restr A character vector that contains variables that should only be specified as weakly exogenous if not contained as endogenous variable in a particular country. An example that has often been used in the literature is to place these restrictions on nominal exchange rates. Default is \code{NULL} in which case all weakly exogenous variables are treated symmetrically. See function \code{getweights} for more details.
#' @param Ex For including truly exogenous variables to the model. Either a \itemize{
#' \item{\code{list object}}{ of maximum length \code{N} that contains the data. Each element of the list refers to a country/entity and has to match the country/entity names in \code{Data}. If no truly exogenous variables are added to the respective country/entity model, omit the entry. The \code{T} rows (i.e., number of time observations), however, need to be the same for each country. Country and variable names are not allowed to contain a dot \code{.} (i.e., a dot) since this is our naming convention.}
#' \item{\code{matrix object}}{ of dimension \code{T} times number of truly exogenous variables. The column names should consist of two parts, separated by a \code{.} (i.e., a dot). The first part should denote the country / entity name and the second part the name of the variable. Country and variable names are not allowed to contain a \code{.} (i.e., a dot).}
#' }
#' @param trend If set to \code{TRUE} a deterministic trend is added to the country models.
#' @param save.country.store If set to \code{TRUE} then function also returns the container of all draws of the individual country models. Significantly raises object size of output and default is thus set to \code{FALSE}.
#' @param applyfun Allows for user-specific apply function, which has to have the same interface than \code{lapply}. If \code{cores=NULL} then \code{lapply} is used, if set to a numeric either \code{parallel::parLapply()} is used on Windows platforms and \code{parallel::mclapply()} on non-Windows platforms.
#' @param cores Specifies the number of cores which should be used. Default is set to \code{NULL} and \code{applyfun} is used.
#' @param verbose If set to \code{FALSE} it suppresses printing messages to the console.
#' @details We provide three priors, the Minnesota labeled \code{MN}, the Stochastic Search Variable Selection prior labeled \code{SSVS} and the Normal-Gamma prior labeled \code{NG}. The first one has been implemented for global VARs in Feldkircher and Huber (2016) and the second one in Crespo Cuaresma et al. (2016), while the last one has been introduced to VAR modeling in Huber and Feldkircher (2019).
#'  Please consult these references for more details on the specification. In the following we will briefly explain the difference between the three priors. The Minnesota prior pushes the variables in the country-specific VAR towards their unconditional stationary mean, or toward a situation where there is at least one unit root present. The SSVS prior is a form of a 'spike' and 'slab' prior. Variable selection is based on the probability of assigning the corresponding regression coefficient to the 'slab' component. If a regression coefficient is non informative, the 'spike' component pushes the associated posterior estimate more strongly towards zero. Otherwise, the slab component resembles a non-informative prior that has little impact on the posterior. Following George et. al. (2008) we set the prior variances for the normal distribution in a semi-automatic fashion. This implies scaling the mixture normal with the OLS standard errors of the coefficients for the full model. The NG prior is a form of global-local shrinkage prior. Hence, the local component shrinks each coefficient towards zero if there is no information for the associated dependent variable. Otherwise, the prior exerts a fat-tail structure such that deviations from zero are possible. The global component is present for each lag, thus capturing the idea that higher lags should be shrunk more aggressively towards zero.
#' @author Maximilian Boeck, Martin Feldkircher, Florian Huber
#' @return Returns a list of class \code{bgvar} with the following elements: \itemize{
#' \item{\code{args}}{ is a list object that contains the arguments submitted to function \code{bgvar}.}
#' \item{\code{xglobal}}{ is a matrix object of dimension T times N (T # of observations, K # of variables in the system).}
#' \item{\code{gW}}{ is the global weight matrix. It is a list, with \code{N} entries, each of which contains the weight matrix of each country.}
#' \item{\code{country.res}}{ is a matrix that contains the posterior mean of the  country models' residuals. The residuals have been obtained as a running mean and thus always relate to the full set of posterior draws. This implies that in case you have opted for trimming the draws the residuals do not correspond to the posterior draws of the "trimmed" coefficients. This is a storage problem, rather than a statistical problem. Experiments, however, show that residual properties (autocorrelation, cross-sectional correlation) of trimmed and reported residuals are close.}
#' \item{\code{stacked results}}{\itemize{
#'       \item{\code{S_large}}{ is a three-dimensional array (K times K times draws) of the (block-diagonal) posterior variance covariance matrix.}
#'       \item{\code{F_large}}{ is a four-dimensional array (K times K times lags times draws) of the coefficients.}
#'       \item{\code{Ginv_large}}{ is a three-dimensional array (K times K times draws) of the inverse of the G matrix.}
#'       \item{\code{A_large}}{ is a three-dimensional array (K times K+1 times draws) of the posterior estimates for the K coefficients plus a global constant.}
#'       \item{\code{F.eigen}}{ in case \code{eigen="TRUE"}, returns a vector that contains for each posterior draw the modulus of the largest eigenvalue of the companion matrix.}
#'       \item{\code{trim.info}}{ is a character vector. Contains information regarding the nr. of stable draws out of total (thinned) draws. Experience shows that a maximum eigenvalue of \code{1.05} seems a reasonable choice when working with data in levels to generate stable impulse responses.}
#' }}
#' \item{\code{cc.results}}{ each entry of this list contains an list object of length \code{N}. Each entry in the list corresponds to one country model and contains one of the following posterior medians.
#' \itemize{
#'       \item{\code{coeffs}}{ contains in each entry the matrix with the posterior median of the estimated coefficients. Columns of the matrix correspond to an equation in the country model (i.e., the dependent variable) and rows to coefficient estimates of the explanatory variables.}
#'       \item{\code{sig}}{ contains in each entry the variance-covariance matrix for each point in time. If \code{SV=FALSE} all entries along the time dimension are the same.}
#'       \item{\code{theta}}{ contains in each entry the estimated prior variances for the coefficients. Explains how much shrinkage is induced on each coefficient depending on the prior setup.}
#'       \item{\code{res}}{ contains in each entry a matrix of dimension (T-p  times K) with the posterior median of the residuals of the cross-country models.}
#'       \item{\code{shrink}}{ in case \code{prior="MN"} each entry contains the estimated shrinkage parameters.}
#'       \item{\code{PIP}}{ in case \code{prior="SSVS"} returns a list object. The first slot in the list \code{PIP.cc}, is a list of length \code{N} and contains the posterior inclusion probabilities of the country models. The second slot in the list, named \code{PIP.avg} yields simple averages (over the country models where a particular variable has been included) of the posterior inclusion probabilities.}
#'       \item{\code{lambda2}}{ in case \code{prior="NG"} each entry contains the estimated global shrinkage parameters. It is a matrix of dimension (p+1 times 3). Columns refer to the endogenous, weakly exogenous and shrinkage parameters for the covariances. Rows correspond to different degree of shrinkage per lag of the variables starting with the contemporaneous lag (only for weakly exogenous variables). In case of the covariances just one global shrinkage parameter is estimated.}
#'       \item{\code{tau}}{ in case \code{prior="NG"} each entry contains the estimated parameter that governs the heaviness of the tails of the marginal prior distribution of the coefficients associated to endogenous variables. Structure is the same as \code{lambda2}.}
#' }}}
#' @examples
#' \dontshow{
#' library(BGVAR)
#' data(eerData)
#' cN<-c("EA","US","UK")
#' eerData<-eerData[cN]
#' W.trade0012<-apply(W.trade0012[cN,cN],2,function(x)x/rowSums(W.trade0012[cN,cN]))
#' hyperpara <- list(tau0=0.1,tau1=3,kappa0=0.1,kappa1=7,a_1=0.01,b_1=0.01,p_i=0.5,q_ij=0.5)
#' model.ssvs <- bgvar(Data=eerData,W=W.trade0012,plag=1,draws=100,burnin=100,
#'                     prior="SSVS",SV=FALSE,hyperpara=hyperpara,thin=1)
#' 
#' W.list<-lapply(W.list,function(l){l<-apply(l[cN,cN],2,function(x)x/rowSums(l[cN,cN]))})
#' variable.list<-list();variable.list$real<-c("y","Dp","tb");variable.list$fin<-c("stir","ltir","rer")
#' model.mn <- bgvar(Data=eerData, W=W.list[c("tradeW.0012","finW0711")],plag=1,draws=100, 
#'                   burnin=100,prior="MN",SV=TRUE,thin=2,variable.list=variable.list)
#' print(model.mn)
#' }
#' \donttest{
#' library(BGVAR)
#' # replicate Feldkircher and Huber (2016) using trade based weights
#' data(eerData)
#' hyperpara <- list(tau0=0.1,tau1=3,kappa0=0.1,kappa1=7,a_1=0.01,b_1=0.01,p_i=0.5,q_ij=0.5)
#' model.ssvs <- bgvar(Data=eerData,W=W.trade0012,plag=1,draws=100,burnin=100,
#'                     prior="SSVS",SV=FALSE,hyperpara=hyperpara,thin=1)
#' print(model.ssvs)
#' 
#' data("eerData")
#' variable.list<-list();variable.list$real<-c("y","Dp","tb");variable.list$fin<-c("stir","ltir","rer")
#' model.mn <- bgvar(Data=eerData, W=W.list[c("tradeW.0012","finW0711")], plag=1, draws=200, 
#'                   burnin=100,prior="MN",SV=TRUE,thin=2,variable.list=variable.list)
#' print(model.mn)
#' 
#' data(monthlyData)
#' EA.weights$variables <- c("EAstir","total.assets","M3","ciss","y","p")
#' OC.weights$variables <- c("poil","qoil","y")
#' OE.weights <- list(EB=EA.weights,OC=OC.weights)
#' hyperpara<-list(c_tau = 0.01, d_tau = 0.01,e_lambda=1.5,d_lambda=1, 
#'                 prmean=0,a_i=0.01,b_i=0.01,a_start=.6,sample_A=FALSE)
#' model.ssvs <- bgvar(Data=monthlyData,W=W,plag=2,draws=100,burnin=100,prior="SSVS",
#'                     hyperpara=hyperpara,eigen=TRUE,SV=TRUE,OE.weights=OE.weights)
#' print(model.ssvs)
#' }
#' @references 
#' Crespo Cuaresma, J., Feldkircher, M. and F. Huber (2016) Forecasting with Global Vector Autoregressive Models: A Bayesian Approach. \emph{Journal of Applied Econometrics}, Vol. 31(7), pp. 1371-1391.

#' Doan, T. R., Litterman, B. R. and C. A. Sims (1984) Forecasting and Conditional Projection Using Realistic Prior Distributions. \emph{Econometric Reviews}, Vol. 3, pp. 1-100.
#' 
#' Dovern, J., Feldkircher, M. and F. Huber (2016) Does joint modelling of the world economy pay off? Evaluating multivariate forecasts from a Bayesian GVAR. \emph{Journal of Economic Dynamics and Control}, Vol. 70, pp. 86-100.
#' 
#' Feldkircher, M. and F. Huber (2016) The International Transmission of US Shocks - Evidence from Bayesian Global Vector Autoregressions. \emph{European Economic Review}, Vol. 81, pp. 167-188.
#' 
#' Feldkircher, M. Gruber, T. and F. Huber (2020) International effects of a compression of euro area yield curves. \emph{Journal of Banking & Finance}, Vol. 113, pp. 11-14.
#' 
#' George, E.I., Sun, D. and S. Ni (2008) Bayesian stochastic search for var model restrictions. \emph{Journal of Econometrics}, Vol. 142, pp. 553-580.
#' 
#' Georgiadis, G. (2015) Examining asymmetries in the transmission of monetary policy in the euro area: Evidence from a mixed cross-section global VAR model. \emph{European Economic Review}, Vol. 75, pp. 195-215.
#' 
#' Huber, F. and M. Feldkircher (2016) Adaptive Shrinkage in Bayesian Vector Autoregressive Models. \emph{Journal of Business and Economic Statistics}, Vol. 37(1), pp. 27-39.
#' 
#' Mohaddes, K. and M. Raissi (2018). Compilation, Revision and Updating of the Global VAR (GVAR) Database, 1979Q2-2016Q4. University of Cambridge: Faculty of Economics (mimeo).
#' 
#' Mohaddes, K. and M. Raissi (2019) The US oil supply revolution and the global economy. \emph{Empirical Economics}, Vol. 57, pp. 515-546.
#' 
#' Pesaran, M.H., Schuermann T. and S.M. Weiner (2004) Modeling Regional Interdependencies Using a Global Error-Correcting Macroeconometric Model. \emph{Journal of Business and Economic Statistics}, Vol. 22, pp. 129-162.
#' 
#' Sims, C. A. (1992) Bayesian Inference for Multivariate Time Series with Trend. \emph{Mimeo}, presented at the American statistical Association meeting.
#' 
#' Sims, C.A. and T. Zha (1998) Bayesian Methods for Dynamic Multivariate Models. \emph{International Economic Review}, Vol. 39, pp. 949-968.
#' @importFrom abind adrop
#' @importFrom GIGrvg rgig
#' @importFrom parallel parLapply mclapply
#' @importFrom Rcpp evalCpp
#' @importFrom stats is.ts median time ts
#' @importFrom xts is.xts
#' @importFrom zoo coredata
bgvar<-function(Data,W,plag=1,draws=5000,burnin=5000,prior="NG",SV=TRUE,h=0,thin=1,hyperpara=NULL,eigen=FALSE,variable.list=NULL,OE.weights=NULL,Wex.restr=NULL,Ex=NULL,trend=FALSE,save.country.store=FALSE,applyfun=NULL,cores=NULL,verbose=TRUE){
  start.bgvar <- Sys.time()
  #--------------------------------- checks  ------------------------------------------------------#
  if(!is.list(Data) & !is.matrix(Data)){
    stop("Please provide the argument 'Data' either as 'list' or as 'matrix' object.")
  }
  if(!is.list(W) & !is.matrix(W)){
    stop("Please provide the argument 'W' either as 'list' or as 'matrix' object.")
  }
  if(!is.null(Ex)){
    if(!is.list(Ex) & !is.matrix(Ex)){
      stop("Please provide the argument 'Ex' either as 'list' or as 'matrix' object.")
    }
  }
  if(!is.numeric(plag)){
    stop("Please specify number of lags as numeric.")
  }
  if(any(is.na(plag))){
    stop("Please specify number of lags.")
  }
  if(length(plag)>1 || plag<1){
    stop("Please specify number of lags accordingly. One lag length parameter for the whole model.")
  }
  if(!is.numeric(draws) | !is.numeric(burnin)){
    stop("Please specify number of draws and burnin as numeric.")
  }
  if(length(draws)>1 || draws<0 || length(burnin)>1 || burnin<0){
    stop("Please specify number of draws and burnin accordingly. One draws and burnin parameter for the whole model.")
  }
  #-------------------------- construct arglist ----------------------------------------------------#
  args <- .construct.arglist(bgvar)
  if(verbose){
    cat("\nStart estimation of Bayesian Global Vector Autoregression.\n\n")
    cat(paste("Prior: ",ifelse(prior=="MN","Minnesota prior",ifelse(prior=="SSVS","Stochastic Search Variable Selection prior","Normal-Gamma prior")),".\n",sep=""))
    cat(paste("Lag order: ",plag,"\n",sep=""))
    cat(paste("Stochastic volatility: ", ifelse(SV,"enabled","disabled"),".\n",sep=""))
  }
  #------------------------------ user checks  ---------------------------------------------------#
  # check Data
  if(is.matrix(Data)){
    if(any(is.na(Data))){
      stop("The data you have submitted contains NAs. Please check the data.")
    }
    if(!all(grepl("\\.",colnames(Data)))){
      stop("Please separate country- and variable names with a point.")
    }
    cN <- unique(unlist(lapply(strsplit(colnames(Data),".",fixed=TRUE),function(l) l[1])))
    N  <- length(cN)
    if(!all(nchar(cN)>1)){
      stop("Please provide entity names with minimal two characters.")
    }
    temp <- list()
    for(cc in 1:N){
      temp[[cN[cc]]] <- Data[,grepl(cN[cc],colnames(Data))]
      colnames(temp[[cN[cc]]]) <- unlist(lapply(strsplit(colnames(temp[[cN[cc]]]),".",fixed=TRUE),function(l)l[2]))
    }
    Data <- temp
  }else if(is.list(Data)){
    if(any(unlist(lapply(Data,is.na)))){
      stop("The data you have submitted contains NAs. Please check the data.")
    }
    N <- length(Data)
    # check names
    if(is.null(names(Data))){
      names(Data)<-paste(c,1:length(Data),sep="")
    }
    cN <- names(Data)
    if(!all(nchar(cN)>1)){
      stop("Please provide entity names with minimal two characters..")
    }
    isTS  <- unlist(lapply(Data,function(l)is.ts(l)))
    isXTS <- unlist(lapply(Data,function(l)is.xts(l)))
    if(!all(isTS) & any(isTS)){
      stop("Please provide all list elements as time-series objects.")
    }
    if(!all(isXTS) & any(isXTS)){
      stop("Please provide all list elements as xts objects.")
    }
    isTS <- all(isTS); isXTS <- all(isXTS)
    Traw <- unique(unlist(lapply(Data,function(l)nrow(l))))
    if(length(Traw)>1){
      stop("Please provide same sample size for all countries.")
    }
    for(cc in 1:N){
      if(isTS || isXTS){
        temp       <- as.character(time(Data[[cc]]))
        years      <- unique(regmatches(temp,regexpr("^[0-9]{4}",temp)))
        months     <- temp
        for(kk in 1:length(years)) months <- gsub(paste(years[kk],"(\\.)?",sep=""),"",months)
        freq       <- length(unique(months))
        months     <- strtrim(months,3)
        startmonth <- ifelse(months[1]=="","01",ifelse(months[1]=="083","02",ifelse(months[1]=="166","03",ifelse(months[1]=="25","04",
                      ifelse(months[1]=="333","05",ifelse(months[1]=="416","06",ifelse(months[1]=="5","07",ifelse(months[1]=="583","08",
                      ifelse(months[1]=="666","09",ifelse(months[1]=="75","10",ifelse(months[1]=="833","11","12")))))))))))
        timeindex  <- seq.Date(from=as.Date(paste(years[1],"-",startmonth,"-01",sep=""), format="%Y-%m-%d"), 
                               by=ifelse(freq==12,"months","quarter"), length.out = Traw)
        Data[[cc]] <- ts(coredata(Data[[cc]]), start=c(as.numeric(years[1]),as.numeric(startmonth)),
                        frequency=freq)
      }else{
        timeindex  <- seq.Date(from=as.Date("1830-08-01", format="%Y-%m-%d"), by="month", length.out = Traw)
        temp       <- coredata(Data[[cc]])
        Data[[cc]] <- ts(temp, start=c(1830,8), frequency=12)
      }
    }
    args$time <- timeindex
    args$Traw <- length(timeindex)
  }
  args$Data <- Data
  # check Weight matrix if matrix
  if(is.matrix(W)){
    W.aux<-list();W.aux$W<-W;W<-W.aux;rm(W.aux) # convert W into a list
  }
  if(any(unlist(lapply(W,is.na)))){
    stop("The weight matrix you have provided contains NAs. Please check the weight matrix.")
  }
  for(ww in 1:length(W)){
    if(is.null(OE.weights)){
      if(!nrow(W[[ww]])==N){
        stop("Data and W matrix not of the same dimension.")
      }
      if(!all(cN%in%rownames(W[[ww]]))){
        stop("Please provide the same country names for the Data and W objects.")
      }
      # make sure that W and Data are in the same order
      W[[ww]]<-W[[ww]][cN,cN]
    }else{
      if(!(nrow(W[[ww]])+length(OE.weights))==N){
        stop("Data and W matrix plus additional weights for other entities are not of the same dimension.")
      }
      if(!all(cN%in%c(rownames(W[[ww]]),names(OE.weights)))){
        stop("Please provide the same country names for the Data and W matrix plus additional weights for other entities.")
      }
      W[[ww]] <- W[[ww]][cN[!cN%in%names(OE.weights)],cN[!cN%in%names(OE.weights)]]
    }
  }
  args$W <- W
  # check truly exogenous variables
  if(!is.null(Ex)){
    if(is.matrix(Ex)){
      if(any(is.na(Ex))){
        stop("The data for exogenous variables you have submitted contains NAs. Please check the data.")
      }
      if(nrow(Ex)!=args$Traw){
        stop("Provided data and truly exogenous data not equally long. Please check.")
      }
      if(!all(grepl("\\.",colnames(Ex)))){
        stop("Please separate country- and variable names with a point.")
      }
      ExcN <- unique(unlist(lapply(strsplit(colnames(Ex),".",fixed=TRUE),function(l) l[1])))
      if(!all(ExcN%in%cN)){
        stop("Provided country names in data and truly exogenous data not equal. Please check.")
      }
      ExN  <- length(ExcN)
      if(!all(nchar(ExcN)>1)){
        stop("Please provide entity names with minimal two characters.")
      }
      temp <- list()
      for(cc in 1:ExN){
        temp[[cN[cc]]] <- Ex[,grepl(ExcN[cc],colnames(Ex))]
        colnames(temp[[ExcN[cc]]]) <- unlist(lapply(strsplit(colnames(temp[[ExcN[cc]]]),".",fixed=TRUE),function(l)l[2]))
      }
      Ex <- temp
    }else if(is.list(Ex)){
      # check for NAs
      if(any(unlist(lapply(Ex,is.na)))){
        stop("The data for exogenous variables you have submitted contains NAs. Please check the data.")
      }
      ExN <- length(Ex)
      # check names
      if(is.null(names(Ex))){
        names(Ex)<-paste(c,1:length(Ex),sep="")
      }
      ExcN <- names(Ex)
      if(!all(nchar(ExcN)>1)){
        stop("Please provide entity names with minimal two characters..")
      }
      if(!all(ExcN%in%cN)){
        stop("Provided country names in data and truly exogenous data not equal. Please check.")
      }
      isTS  <- unlist(lapply(Ex,function(l)is.ts(l)))
      isXTS <- unlist(lapply(Ex,function(l)is.xts(l)))
      if(!all(isTS) & any(isTS)){
        stop("Please provide all list elements as time-series objects.")
      }
      if(!all(isXTS) & any(isXTS)){
        stop("Please provide all list elements as xts objects.")
      }
      isTS <- all(isTS); isXTS <- all(isXTS)
      ExTraw <- unique(unlist(lapply(Ex,function(l)nrow(l))))
      if(length(ExTraw)>1){
        stop("Please provide same sample size for all countries.")
      }
      if(ExTraw!=args$Traw){
        stop("Provided data and truly exogenous data not equally long. Please check.")
      }
    }
  }
  args$Ex <- Ex
  # check prior
  if(!prior%in% c("MN","SSVS","NG")){
    stop("Please selecte one of the following prior options: MN, SSVS or NG")
  }
  # check thinning factor
  if(thin<1){
    thin_mess <- paste("Thinning factor of ",thin," not possible. Adjusted to ",round(1/thin,2),".\n",sep="")
    thin <- round(1/thin,2)
  }
  if(draws%%thin!=0){
    thin_mess <- paste("Thinning factor of ",thin," no divisor of ",draws," (number of draws to save for posterior analysis).\n",sep="")
    div <- .divisors(draws,thin)
    thin <- min(div[which(abs(div-thin)==min(abs(div-thin)))])
    thin_mess <- paste(thin_mess,"New thinning factor: ", thin,". This means every", ifelse(thin==1,"",ifelse(thin==2,paste(" ",thin,"nd ",sep=""), ifelse(thin==3,paste(" ",thin,"rd ",sep=""),paste(" ",thin,"th ",sep="")))), "draw is saved.\n",sep="")
  }else{
    thin_mess <- paste("Thinning factor: ", thin,". This means every ",ifelse(thin==1,"",ifelse(thin==2,paste(thin,"nd ",sep=""),ifelse(thin==3,paste(thin,"rd ",sep=""),paste(thin,"th ",sep="")))),"draw is saved.\n",sep="")
  }
  if(verbose) cat(thin_mess)
  args$thindraws <- draws/thin
  # set default
  if(verbose) cat("Hyperparameter setup: \n")
  default_hyperpara <- list(a_1=0.01,b_1=0.01, prmean=0,# Gamma hyperparameter SIGMA (homoskedastic case) and mean
                            Bsigma=1, a0=25, b0=1.5, bmu=0, Bmu=100^2, # SV hyper parameter
                            shrink1=0.1,shrink2=0.2,shrink3=10^2,shrink4=0.1, # MN
                            tau0=.1,tau1=3,kappa0=0.1,kappa1=7,p_i=0.5,q_ij=0.5,   # SSVS
                            e_lambda=0.01,d_lambda=0.01,a_start=0.7,sample_A=FALSE,a_log=TRUE) # NG
  paras     <- c("a_1","b_1","prmean","Bsigma_sv","a0","b0","bmu","Bmu","shrink1","shrink2","shrink3",
                 "shrink4","tau0","tau1","kappa0","kappa1","p_i","q_ij","e_lambda","d_lambda","a_start","sample_A")
  if(is.null(hyperpara)){
    if(verbose) cat("\t No hyperparameters are chosen, default setting applied.\n")
  }
  if(!is.null(hyperpara)){
    for(para in names(hyperpara)){
      if(!para%in%paras){
        warning(paste0(para," no valid hyperparameter. Please check.\n"))
        next
      }
      default_hyperpara[para] <- hyperpara[para]
      if(para=="a_start") default_hyperpara["a_log"] <- FALSE
    }
    if(verbose) cat("Default values for chosen hyperparamters overwritten.\n")
  }
  #------------------------------ get weights -----------------------------------------------------------------#
  xglobal <- .getweights(Data=Data,W=W,OE.weights=OE.weights,Wex.restr=Wex.restr,variable.list=variable.list)
  
  exo.countries<-xglobal$exo.countries
  exo     <- xglobal$exo
  endo    <- xglobal$endo
  gW      <- xglobal$gW
  xglobal <- xglobal$bigx
  #---------------------------------hold out sample------------------------------------------------------------#
  args$yfull <- xglobal
  xglobal    <- xglobal[1:(nrow(xglobal)-h),,drop=FALSE]
  args$time  <- args$time[1:(length(args$time)-h)]
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
  #------------------------------ estimate BVAR ---------------------------------------------------------------#
  if(verbose) cat("\nEstimation of country models starts... ")
  start.estim <- Sys.time()
  globalpost <- applyfun(1:N, function(cc){
    .BVAR_linear_wrapper(cc=cc,cN=cN,xglobal=xglobal,gW=gW,prior=prior,plag=plag,draws=draws,burnin=burnin,trend=trend,SV=SV,thin=thin,default_hyperpara=default_hyperpara,Ex=Ex)
  })
  names(globalpost) <- cN
  end.estim <- Sys.time()
  diff.estim <- difftime(end.estim,start.estim,units="mins")
  mins <- round(diff.estim,0); secs <- round((diff.estim-floor(diff.estim))*60,0)
  if(verbose) cat(paste(" took ",mins," ",ifelse(mins==1,"min","mins")," ",secs, " ",ifelse(secs==1,"second.","seconds.\n"),sep=""))
  #--------------------------- stacking part for global model -----------------------------------------------------#
  if(is.logical(eigen)){
    if(eigen){trim<-1.05}else{trim<-NULL}
  }else{
    trim<-eigen;eigen<-TRUE
  }
  if(verbose) cat("Start stacking: \n")
  # insert stacking function here
  stacked.results <- .gvar.stacking.wrapper(xglobal=xglobal,plag=plag,globalpost=globalpost,draws=draws,thin=thin,trend=trend,eigen=eigen,trim=trim,verbose=verbose)
  if(!is.null(trim)) {args$thindraws <- length(stacked.results$F.eigen)}
  if(verbose) cat("\nStacking finished.\n")
  #--------------------------- prepare country models -------------------------------------------------------------#
  # country model residuals
  country.coeffs <- lapply(globalpost,function(l) l$post$A_post)
  country.sig    <- lapply(globalpost,function(l) l$post$SIGMA_post)
  country.theta  <- lapply(globalpost,function(l) l$post$theta_post)
  country.res    <- lapply(globalpost,function(l) l$post$res_post)
  varNames       <- lapply(gW,function(x) dimnames(x)[[1]])
  for(cc in 1:N){
    varx <- varNames[[cc]]
    endo <- grep(cN[cc],varx)
    exx  <- which(varx%in%names(exo))
    wex  <- seq(1,length(varx))[-c(endo,exx)]
    if(length(exx)>0){
      wex0 <- c(paste(varx[wex],"*",sep=""),paste(varx[exx],"**",sep=""))
    }else{
      wex0 <- c(paste(varx[wex],"*",sep=""))
    }
    wexL <- endoL <- c()
    for(pp in 1:plag){
      if(length(exx)>0){
        wexL <- c(wexL, paste(varx[wex],"*_lag",pp,sep=""), paste(varx[exx],"**_lag",pp,sep=""))
      }else{
        wexL   <- c(wexL, paste(varx[wex],"*_lag",pp,sep=""))
      }
      endoL  <- c(endoL,paste(varx[endo],"_lag",pp,sep=""))
    }
    if(cN[cc]%in%names(Ex)){
      tex <- colnames(Ex[[cN[cc]]])
    }else{tex<-NULL}
    names <- c(endoL,wex0,wexL,tex,"cons")
    if(trend) names <- c(names,"trend")
    
    rownames(country.coeffs[[cc]]) <- names
    dimnames(country.sig[[cc]])[[2]]<-dimnames(country.sig[[cc]])[[3]]<-varx[endo]
  }
  cc.results <- list(coeffs=country.coeffs,sig=country.sig,theta=country.theta,res=country.res)
  if(prior=="MN"){
    cc.results$shrink <- lapply(globalpost,function(l) l$post$shrink)
  }else if(prior=="SSVS"){
    country.shrink <- lapply(globalpost,function(l) l$post$PIP)
    for(cc in 1:N) rownames(country.shrink[[cc]]) <- rownames(country.coeffs[[cc]])
    cc.results$PIP <- .avg.shrink(country.shrink,prior="SSVS")
  }else if(prior=="NG"){
    cc.results$lambda2 <- lapply(globalpost,function(l) l$post$lambda2_post)
    cc.results$tau     <- lapply(globalpost,function(l) l$post$tau_post)
  }
  if(save.country.store){
    cc.results$store <- lapply(globalpost,function(l) l$store)
  }
  #---------------------- return output ---------------------------------------------------------------------------#
  out  <- structure(list("args"=args,
                         "xglobal"=xglobal,
                         "gW"=gW,
                         "stacked.results"=stacked.results,
                         "cc.results"=cc.results), class = "bgvar")
  end.bgvar <- Sys.time()
  diff.bgvar <- difftime(end.bgvar,start.bgvar,units="mins")
  mins.bgvar <- round(diff.bgvar,0); secs.bgvar <- round((diff.bgvar-floor(diff.bgvar))*60,0)
  if(verbose) cat(paste("\n Needed time for estimation of bgvar: ",mins.bgvar," ",ifelse(mins.bgvar==1,"min","mins")," ",secs.bgvar, " ",ifelse(secs.bgvar==1,"second.","seconds.\n"),sep=""))
  return(out)
}

#' @method print bgvar
#' @export
#' @importFrom utils object.size
print.bgvar<-function(x, ...){
  cat("---------------------------------------------------------------------------------------")
  cat("\n")
  cat("Model Info:")
  cat("\n")
  cat(paste("Prior: ",x$args$prior,sep=""))
  cat("\n")
  cat(paste("Nr. of lags: ",x$args$plag,sep=""))
  cat("\n")
  cat(paste("Nr. of posterior draws: ",x$args$draws,"/",x$args$thin,"=",floor(x$args$draws/x$args$thin),sep=""))
  cat("\n")
  cat(paste("Size of GVAR object: ",format(object.size(x),units="MB"),sep=""))
  cat("\n")
  cat(x$stacked.results$trim.info)
  cat("\n")
  cat("---------------------------------------------------------------------------------------")
  cat("\n")
  cat("Model specification:")
  cat("\n")
  
  endo <- lapply(x$cc.results$coeffs,colnames)
  exo  <- lapply(x$cc.results$coeffs,rownames)
  cN   <- names(endo)
  
  vars <- list()
  for(i in 1:length(endo)){
    vars[[i]] <- c(gsub(paste(cN[i],".",sep=""),"",endo[[i]]), exo[[i]][-grep("_lag",exo[[i]])])
    vars[[i]] <- vars[[i]][-charmatch("cons",vars[[i]])]
  }
  
  varNames <- lapply(vars,function(l) paste(l,collapse=", "))
  names(varNames) <- cN
  
  for(i in 1:length(varNames)){
    cat("\n")
    cat(paste0(names(varNames[i]),": ",varNames[[i]]))
  }
  invisible(x)
}

#' @name summary
#' @title Summary of Bayesian GVAR
#' @description Output gives model information as well as some descriptive statistics on convergence properties, likelihood, serial autocorrelation in the errors and the average pairwise autocorrelation of cross-country residuals.
#' @aliases summary summary.bgvar
#' @param object an object of class \code{bgvar}.
#' @param ... other arguments.
#' @return No return value.
#' @seealso
#' \code{\link{bgvar}} to estimate a \code{bgvar} object.
#' \code{\link{avg.pair.cc}} to compute average pairwise cross-country correlation of cross-country residuals separately.
#' \code{\link{resid.corr.test}} to compute F-test on first-order autocorrelation of cross-country residuals separately.
#' @author Maximilian Boeck
#' @examples
#' \dontshow{
#' library(BGVAR)
#' data(eerData)
#' cN<-c("EA","US")
#' eerData<-eerData[cN]
#' W.trade0012<-apply(W.trade0012[cN,cN],2,function(x)x/rowSums(W.trade0012[cN,cN]))
#' model.ssvs <- bgvar(Data=eerData,W=W.trade0012,plag=1,draws=50,burnin=50,
#'                     prior="SSVS",SV=TRUE,trend=TRUE)
#' summary(model.ssvs)
#' }
#' \donttest{
#' set.seed(571)
#' library(BGVAR)
#' data(eerData)
#' model.ssvs <- bgvar(Data=eerData,W=W.trade0012,plag=1,draws=100,burnin=100,
#'                     prior="SSVS",thin=1,SV=TRUE,trend=TRUE)
#' summary(model.ssvs)
#' }
#' @export
summary.bgvar <- function(object, ...){
  CD  <- conv.diag(object)
  res <- residual.corr.test(object,lag.cor=1,alpha=0.95)
  cross.corr <- avg.pair.cc(object)
  LL <- logLik(object)
  out  <- structure(list("object"=object,
                         "CD"=CD,
                         "res"=res,
                         "cross.corr"=cross.corr,
                         "logLik"=LL), class = "bgvar.summary")
  return(out)
}

#' @method print bgvar.summary
#' @export
print.bgvar.summary <- function(x, ...){
  cat("-------------------------------------------------------------------------------")
  cat("\n")
  cat("Model Info:")
  cat("\n")
  cat(paste("Prior: ",x$object$args$prior,sep=""))
  cat("\n")
  cat(paste("Nr. of lags: ",x$object$args$plag,sep=""))
  cat("\n")
  cat(paste("Nr. of posterior draws: ",x$object$args$draws,"/",x$object$args$thin,"=",x$object$args$draws/x$object$args$thin,sep=""))
  cat("\n")
  if(x$object$args$eigen){
    cat("Number of stable posterior draws: ",length(x$object$stacked.results$F.eigen))
    cat("\n")
  }
  cat(paste("Number of countries: ",length(x$object$gW),sep=""))
  cat("\n")
  cat("-------------------------------------------------------------------------------")
  cat("\n")
  cat("Convergence diagnostics")
  cat("\n")
  cat(paste("Geweke statistic:\n",x$CD$perc,sep=""))
  cat("\n")
  cat("-------------------------------------------------------------------------------")
  cat("\n")
  cat(paste("Global Likelihood: ",round(x$logLik,2),sep=""))
  cat("\n")
  cat("F-test, first order serial autocorrelation of cross-country residuals")
  cat("\n")
  cat("Summary statistics:")
  cat("\n")
  temp <- kable(x$res$p.res, row.names=TRUE, "rst")
  for(ii in 1:length(temp)){
    cat(paste0(temp[ii],"\n"))
  }
  cat("-------------------------------------------------------------------------------")
  cat("\n")
  cat("Average pairwise cross-country correlation of country model residuals")
  cat("\n")
  cat("Summary statistics:")
  cat("\n")
  temp <- kable(x$cross.corr$res.res, row.names=TRUE, "rst")
  for(ii in 1:length(temp)){
    cat(paste0(temp[ii],"\n"))
  }
  cat("-------------------------------------------------------------------------------")
  cat("\n")
  
  invisible(x)
}

#' @name residuals
#' @export
#' @title Extract Residuals of Bayesian GVAR
#' @description Calculate residuals of the global model and the country models.
#' @aliases residuals residuals.bgvar
#' @param object a fitted \code{bgvar} object.
#' @param ... other arguments.
#' @details This function calculates residuals of the global and the country models based on a \code{bgvar} object. Country models' residuals are equivalent to output generated by the \code{print.bgvar} function in case no trimming has been used. If trimming was invoked to discard unstable draws output of both functions might differ since \code{print.bgvar} calculates residuals as a running mean to save storage which is based on the \emph{whole} set of posterior draws (including discarded draws). In this case it is recommended to recalculate the residuals with \code{residuals.bgvar} and re-do the serial autocorrelation or average pairwise cross-correlation analysis using functions \code{resid.corr.test} and \code{avg.pair.cc}.
#' 
#' @return returns a list with the following arguments \itemize{
#' \item{\code{global}}{ is a (T-p) times K times draws/thin array containing the residuals of the global model.}
#' \item{\code{country}}{ is a (T-p) times K times draws/thin array containing the residuals of the country models.}
#' \item{\code{Data}}{ is a (T-p) times K matrix containing the data of the model.}
#' }
#' @author Maximilian Boeck, Martin Feldkircher
#' @seealso \code{\link{bgvar}} for estimation of a \code{bgvar} object.
#' @importFrom stats resid
#' @examples
#' \donttest{
#' library(BGVAR)
#' data(eerData)
#' model.ssvs <- bgvar(Data=eerData,W=W.trade0012,plag=1,draws=100,burnin=100,
#'                     prior="SSVS")
#' res <- residuals(model.ssvs)
#' }
residuals.bgvar <- function(object, ...){
  if(!inherits(object, "bgvar")) {stop("Please provide a `bgvar` object.")}
  G.mat   <- object$stacked.results$Ginv_large
  A.mat   <- object$stacked.results$A_large
  plag    <- object$args$plag
  draws   <- object$args$thindraws
  time    <- object$args$time
  trend   <- object$args$trend
  xglobal <- object$xglobal
  YY      <- xglobal[(plag+1):nrow(xglobal),]
  XX      <- cbind(.mlag(xglobal,plag),1)
  XX      <- XX[(plag+1):nrow(XX),]
  if(trend) XX <- cbind(XX,seq(1,nrow(XX)))
  
  rownames(YY) <- as.character(time[-c(1:plag)])
  res.array.country<-res.array.global<-array(0,dim=c(draws,dim(YY)))
  for(irep in 1:draws){
    res.array.global[irep,,]  <- (YY-XX%*%t(A.mat[irep,,]))
    res.array.country[irep,,] <- (res.array.global[irep,,]%*%t(solve(G.mat[irep,,])))
  }
  out <- structure(list(global=res.array.global,country=res.array.country,Data=YY),
                   class = "bgvar.resid")
  return(out)
}

#' @rdname residuals
#' @examples 
#' \donttest{
#' resid(model.ssvs)
#' }
#' @export
resid.bgvar <- residuals.bgvar

#' @name coef
#' @title Extract Model Coefficients of Bayesian GVAR
#' @description Extracts the global model coefficients for \code{bgvar} for certain quantiles of the posterior distribution. \code{coefficients} is an \emph{alias} for it.
#' @param object an object of class \code{bgvar}.
#' @param ... additional arguments.
#' @param quantile reported quantiles. Default is set to the median.
#' @return Returns an \code{q} times \code{K} times \code{K} times \code{p} array of the global coefficients, where \code{q} is the number of specified quantiles (this dimension is dropped if \code{q=1}), \code{K} the number of endogenous variables and \code{p} number of lags.
#' @export
#' @importFrom stats quantile
#' @examples
#' \donttest{
#' library(BGVAR)
#' data(eerData)
#' model.ng <- bgvar(Data=eerData,W=W.trade0012,plag=1,draws=100,burnin=100)
#' coef(model.ng)
#' }
coef.bgvar<-function(object, ..., quantile=.50){
  out <- apply(object$stacked.results$F_large,c(2,3,4),quantile,quantile,na.rm=TRUE)
  dimnames(out)[[2]] <- colnames(object$xglobal)
  return(out)
}

#' @rdname coef
#' @examples
#' \donttest{
#' coefficients(model.ng)
#' }
#' @importFrom stats coefficients
#' @export
coefficients.bgvar <- coef.bgvar

#' @name vcov
#' @title Extract Variance-covariance Matrix of Bayesian GVAR
#' @description Extracts the global variance-covariance matrix for \code{bgvar} for certain quantiles of the posterior distribution. 
#' @param object an object of class \code{bgvar}.
#' @param ... additional arguments.
#' @param quantile reported quantiles. Default is set to median.
#' @return Returns an \code{q} times \code{K} times \code{K} array of the global variance-covariance matrix, where \code{q} is the number of specified quantiles (this dimension is dropped if \code{q=1}) and  \code{K} the number of endogenous variables.
#' @importFrom stats vcov
#' @examples
#' \donttest{
#' library(BGVAR)
#' data(eerData)
#' model.ng <- bgvar(Data=eerData,W=W.trade0012,plag=1,draws=100,burnin=100)
#' vcov(model.ng)
#' }
#' @export
vcov.bgvar<-function(object, ..., quantile=.50){
  S_qu <- apply(object$stacked.results$S_large,c(2,3),quantile,quantile,na.rm=TRUE)
  Ginv_qu <- apply(object$stacked.results$Ginv_large,c(2,3),quantile,quantile,na.rm=TRUE)
  if(length(quantile)==1){
    out <- Ginv_qu%*%S_qu%*%t(Ginv_qu)
  }else{
    out <- sapply(1:length(quantile),function(qq)Ginv_qu[qq,,]%*%S_qu[qq,,]%*%t(Ginv_qu[qq,,]),simplify="array")
    out <- aperm(out,c(3,1,2))
  }
  return(out)
}

#' @name fitted
#' @title Extract Fitted Values of Bayesian GVAR
#' @description Extracts the fitted values for \code{bgvar}.
#' @param object an object of class \code{bgvar}.
#' @param ... additional arguments.
#' @param global if \code{TRUE} global fitted values are returned otherwise country fitted values.
#' @return Returns an \code{T} times \code{K} matrix, where \code{T} is the number of observations and \code{K} number of endogenous variables.
#' @importFrom stats fitted
#' @examples 
#' \donttest{
#' library(BGVAR)
#' data(eerData)
#' model.ng <- bgvar(Data=eerData,W=W.trade0012,plag=1,draws=100,burnin=100)
#' fitted(model.ng)
#' }
#' @export
fitted.bgvar<-function(object, ..., global=TRUE){
  plag     <- object$args$plag
  xglobal  <- object$xglobal
  trend    <- object$args$trend
  XX       <- .mlag(xglobal,plag)
  YY       <- xglobal[-c(1:plag),,drop=FALSE]
  XX       <- cbind(XX[-c(1:plag),,drop=FALSE],1)
  bigT     <- nrow(YY)
  if(trend) XX <- cbind(XX,seq(1,bigT))
  if(global){
    A_post <- apply(object$stacked.results$A_large,c(2,3),median)
    fit    <- XX%*%t(A_post)
  }else{
    fit <- YY-do.call("cbind",object$cc.results$res)
  }
  return(fit)
}

#' @name logLik
#' @title Extract Log-likelihood of Bayesian GVAR
#' @description Extract Log-Likelihood for \code{bgvar}.
#' @param object an object of class \code{bgvar}.
#' @param ... additional arguments.
#' @param quantile reported quantiles. Default is set to median.
#' @return Returns an vector of dimension \code{q} (number of specified quantiles) of global log-likelihoods.
#' @importFrom stats logLik
#' @examples 
#' \donttest{
#' library(BGVAR)
#' data(eerData)
#' model.ng <- bgvar(Data=eerData,W=W.trade0012,plag=1,draws=100,burnin=100)
#' logLik(model.ng)
#' }
#' @export
logLik.bgvar<-function(object, ..., quantile=.50){
  if(length(quantile)!=1){
    stop("Please provide only one quantile.")
  }
  temp <- object$args$logLik
  if(is.null(temp)){
    xglobal   <- object$xglobal
    plag      <- object$args$plag
    trend     <- object$args$trend
    bigT      <- nrow(xglobal)
    bigK      <- ncol(xglobal)
    thindraws <- object$args$thindraws
    X_large   <- cbind(.mlag(xglobal,plag),1)
    if(trend) X_large <- cbind(X_large,seq(1:bigT))
    Y_large   <- xglobal[(plag+1):bigT,,drop=FALSE]
    X_large   <- X_large[(plag+1):bigT,,drop=FALSE]
    A_large   <- object$stacked.results$A_large
    S_large   <- object$stacked.results$S_large
    Ginv_large<- object$stacked.results$Ginv_large
    globalLik <- c(globalLik(Y_in=Y_large,X_in=X_large,A_in=A_large,S_in=S_large,Ginv_in=Ginv_large,thindraws=thindraws)$globalLik)
    
    out <- quantile(globalLik,quantile,na.rm=TRUE)
    eval.parent(substitute(object$args$logLik<-out))
  }
  attributes(out) <- list(nall=bigT, nobs=bigT, df=bigK)
  class(out) <- "logLik"
  return(out)
}
