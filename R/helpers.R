#' @name avg.pair.cc
#' @export
#' @title Average Pairwise Cross-Sectional Correlations
#' @description Computes average pairwise cross-sectional correlations of the data and the country models' residuals.
#' @details If used for analyzing the country models' residuals, \code{avg.pair.cc} computes for each country and a given variable, the average cross-sectional correlation (either for the data or for the residuals). In theory, including foreign variables should soak up cross-sectional residual dependence and  correlation of the residuals should be small. Otherwise dynamic analysis, especially using GIRFs, might lead to invalid results. See Dees et al. (2007) for more details.
#' @usage avg.pair.cc(object, digits=3)
#' @param object Either an object of class \code{bgvar} or residuals of class \code{bgvar.res}.
#' @param digits Number of digits that should be used to print output to the console.
#' @return Returns a list with the following elements
#' \item{\code{data.cor}}{ is a matrix containing in the rows the cross-sections and in the columns the cross-sectional pairwise correlations of the data per variable.}
#' \item{\code{resid.cor}}{ is a matrix containing in the rows the cross-sections and in the columns the cross-sectional pairwise correlations of the country models' residuals per variable.}
#' \item{\code{resid.corG}}{ is a matrix containing in the rows the cross-sections and in the columns the cross-sectional pairwise correlations of the global models' residuals per variable. Only available when \code{avg.pair.cc} has been applied to a \code{bgvar.res} object from \code{residuals}.}
#' \item{\code{data.res}}{ is a summary object showing the number and percentage of correlations <0.1, between 0.1-0.2, 0.2-0.5 and <0.5 per variable of the data.}
#' \item{\code{res.res}}{ is a summary object showing the number and percentage of correlations <0.1, between 0.1-0.2, 0.2-0.5 and <0.5 per variable of the country models' residuals. This is also what is used by \code{print.bgvar}.}
#' \item{\code{res.resG}}{ is a summary object showing the number and percentage of correlations <0.1, between 0.1-0.2, 0.2-0.5 and <0.5 per variable of the global models' residuals. Only available when \code{avg.pair.cc} has been applied to a \code{bgvar. res} object from \code{residuals}.}
#' @author Martin Feldkircher
#' @references 
#' Dees, S., Di Mauro F., Pesaran, M. H. and Smith, L. V. (2007) \emph{Exploring the international linkages of the euro area: A global VAR analysis.} Journal of Applied Econometrics, Vol. 22, pp. 1-38.
#' @seealso 
#' \code{\link{bgvar}} for estimation of a \code{bgvar} object.
#' \code{\link{residuals}} for calculating the residuals from a \code{bgvar} object and creating a \code{bgvar.res} object.
#' @examples
#' \donttest{
#' library(BGVAR)
#' data(testdata)
#' model.mn <- bgvar(Data=testdata,W=W.test,plag=1,SV=TRUE,
#'                   draws=100,burnin=100,prior="MN")
#' avg.pair.cc(model.mn)
#' 
#' res <- residuals(model.mn)
#' avg.pair.cc(res)
#' }
#' @importFrom stats cor
avg.pair.cc=function(object, digits=3){
  if(inherits(object, "bgvar")){
    lags <- object$args$lags
    pmax <- max(lags)
    dat  <- object$xglobal[-c(1:pmax),]
    res  <- do.call("cbind",object$cc.results$res)
    res  <- res[,colnames(dat)] 
  }
  if(inherits(object, "bgvar.resid")){
    dat    <- object$Data
    res    <- apply(object$country,c(2,3),mean)  
    res.g  <- apply(object$global,c(2,3),mean) 
  }
  bigT     <- nrow(res)
  varNames <- colnames(dat)
  cN   <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x) x[1]))
  vars <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x) x[2]))
  
  idx<-lapply(as.list(1:length(vars)),function(x) grep(paste(".",vars[x],sep=""),colnames(dat)))
  names(idx)<-vars
  # kick out exo variables
  exo<-which(sapply(idx,length)==1)
  if(length(exo)>0){
    idx<-idx[-exo]
  }
  
  datL<-resL<-resL.g<-matrix("-",nrow=length(cN),ncol=length(idx))
  rownames(datL)<-rownames(resL)<-rownames(resL.g)<-cN
  colnames(datL)<-colnames(resL)<-colnames(resL.g)<-names(idx)
  
  if(inherits(object,"bgvar")){
    for(i in 1:length(idx)){
      aux.dat <- cor(dat[,idx[[i]]])
      aux.res <- cor(res[,idx[[i]]])
      diag(aux.dat)<-diag(aux.res)<-NA
      ii <- sapply(strsplit(rownames(aux.dat),".",fixed=TRUE),function(x) x[1]) 
      aux.dat <- round(rowMeans(aux.dat,na.rm=TRUE),digits=digits)
      aux.res <-round(rowMeans(aux.res,na.rm=TRUE),digits=digits)
      datL[ii,i]<-aux.dat
      resL[ii,i]<-aux.res
    }
  }
  if(inherits(object,"bgvar.resid")){  # include analysis based on residuals of the global model as well
    for(i in 1:length(idx)){
      aux.dat <- cor(dat[,idx[[i]]])
      aux.res <- cor(res[,idx[[i]]])
      aux.rg  <- cor(res.g[,idx[[i]]])
      diag(aux.dat)<-diag(aux.res)<-diag(aux.rg)<-NA
      ii <- sapply(strsplit(rownames(aux.dat),".",fixed=TRUE),function(x) x[1]) #should be the same for u
      aux.dat <- round(rowMeans(aux.dat,na.rm=TRUE),digits=digits)
      aux.res <- round(rowMeans(aux.res,na.rm=TRUE),digits=digits)
      aux.rg  <- round(rowMeans(aux.rg,na.rm=TRUE),digits=digits)
      datL[ii,i]   <- aux.dat
      resL[ii,i]   <- aux.res
      resL.g[ii,i] <- aux.rg
    }
  }
  
  # Generate summary-table
  pp<-suppressWarnings(apply(datL,2,as.numeric))
  rr<-suppressWarnings(apply(resL,2,as.numeric))
  rg<-suppressWarnings(apply(resL.g,2,as.numeric))
  dat.res<-res.res<-res.resG<-matrix(0,nrow=4,ncol=ncol(datL))
  for(i in 1:ncol(datL)){
    aux<-pp[,i];aux<-abs(aux[which(!is.na(aux))]);K<-length(aux)
    aux2<-rr[,i];aux2<-abs(aux2[which(!is.na(aux2))]);K2<-length(aux2)
    aux3<-rg[,i];aux3<-abs(aux3[which(!is.na(aux3))]);K3<-length(aux3)
    
    dat.res[1,i]<-paste(length(which(aux<=0.1))," (",round((length(which(aux<=0.1))/K)*100,2),"%)",sep="")
    res.res[1,i]<-paste(length(which(aux2<=0.1))," (",round((length(which(aux2<=0.1))/K2)*100,2),"%)",sep="")
    temp<-round((length(which(aux3<=0.1))/K3)*100,2)
    res.resG[1,i]<-paste(length(which(aux3<=0.1))," (",ifelse(is.nan(temp),0,1),"%)",sep="")
    
    dat.res[2,i]<-paste(length(which(aux>0.1&aux<=0.2))," (",round((length(which(aux>0.1&aux<=0.2))/K)*100,2),"%)",sep="")
    res.res[2,i]<-paste(length(which(aux2>0.1&aux2<=0.2))," (",round((length(which(aux2>0.1&aux2<=0.2))/K2)*100,2),"%)",sep="")
    temp<-round((length(which(aux3>0.1&aux3<=0.2))/K3)*100,2)
    res.resG[2,i]<-paste(length(which(aux3>0.1&aux2<=0.2))," (",ifelse(is.nan(temp),0,1),"%)",sep="")
    
    dat.res[3,i]<-paste(length(which(aux>0.2&aux<=0.5))," (",round((length(which(aux>0.2&aux<=0.5))/K)*100,2),"%)",sep="")
    res.res[3,i]<-paste(length(which(aux2>0.2&aux2<=0.5))," (",round((length(which(aux2>0.2&aux2<=0.5))/K2)*100,2),"%)",sep="")
    temp<-round((length(which(aux3>0.2&aux3<=0.5))/K3)*100,2)
    res.resG[3,i]<-paste(length(which(aux3>0.2&aux3<=0.5))," (",ifelse(is.nan(temp),0,temp),"%)",sep="")
    
    
    dat.res[4,i]<-paste(length(which(aux>0.5&aux<=1))," (",round((length(which(aux>0.5&aux<=1))/K)*100,2),"%)",sep="")
    res.res[4,i]<-paste(length(which(aux2>0.5&aux2<=1))," (",round((length(which(aux2>0.5&aux2<=1))/K2)*100,2),"%)",sep="")
    temp<-round((length(which(aux3>0.5&aux3<=1))/K3)*100,2)
    res.resG[4,i]<-paste(length(which(aux3>0.5&aux3<=1))," (",ifelse(is.nan(temp),0,temp),"%)",sep="")
    
  }
  colnames(dat.res) <- colnames(res.res) <- colnames(res.resG) <- colnames(datL)
  rownames(dat.res) <- rownames(res.res) <- rownames(res.resG) <- c("<0.1","0.1-0.2","0.2-0.5",">0.5")
  #dat.res  <- rbind(c("",colnames(datL)),cbind(c("<0.1","0.1-0.2","0.2-0.5",">0.5"),dat.res))
  #res.res  <- rbind(c("",colnames(datL)),cbind(c("<0.1","0.1-0.2","0.2-0.5",">0.5"),res.res))
  #res.resG <- rbind(c("",colnames(datL)),cbind(c("<0.1","0.1-0.2","0.2-0.5",">0.5"),res.resG))
  
  avg.cc<-list(data.cor=datL,resid.cor=resL,resid.corG=resL.g,dat.res=dat.res,res.res=res.res,res.resG=res.resG)
  return(avg.cc)
}

#' @name conv.diag
#' @export
#' @title MCMC Convergence Diagnostics
#' @description This function computes Geweke's Convergence diagnostic making use of the \code{coda} package.
#' @usage conv.diag(object, crit.val=1.96)
#' @param object A fitted \code{bgvar} object.
#' @param crit.val Critical value used for test statistic.
#' @details Geweke (1992) proposed a convergence diagnostic for Markov chains based on a test for equality of the means of the first and last part of a Markov chain (by default we use the first 10\% and the last 50\%). If the samples are drawn from the stationary distribution of the chain, the two means are equal and Geweke's statistic has an asymptotically standard normal distribution. The test statistic is a standard Z-score: the difference between the two sample means divided by its estimated standard error. The standard error is estimated from the spectral density at zero and so takes into account any autocorrelation.
#' @return Returns an object of class \code{bgvar.CD}. This is a list with \describe{
#' \item{\code{geweke.z}}{ Z-scores for a test of equality of means between the first and last parts of the chain. A separate statistic is calculated for each variable in each chain.}
#' \item{\code{perc}}{ is the percentage of Z-scores exceeding \code{crit.val} (in absolute terms).}
#' }
#' @seealso 
#' \code{\link[coda]{geweke.diag}} in the \code{coda} package.
#' \code{\link{bgvar}} for estimation of a \code{bgvar} object.
#' @author Martin Feldkircher
#' @examples
#' \donttest{
#' library(BGVAR)
#' data(testdata)
#' model.mn <- bgvar(Data=testdata,W=W.test,plag=1,draws=200,burnin=200,prior="MN")
#' geweke <- conv.diag(model.mn)
#' }
#' @references 
#' Geweke, J. (1992) Evaluating the accuracy of sampling-based approaches to calculating posterior moments. \emph{Bayesian Statistics} 4 (edited by JM Bernado, JO Berger, AP Dawid and AFM Smith). Clarendon Press, Oxford, UK.
#' @importFrom coda mcmc geweke.diag
conv.diag<-function(object, crit.val=1.96){
  if(!inherits(object, "bgvar")) {stop("Please provide a `bgvar` object.")}
  
  ALPHA <- object$stacked.results$A_large
  draws <- dim(ALPHA)[3]
  d1    <- dim(ALPHA)[1]
  d2    <- dim(ALPHA)[2]
  K     <- d1*d2
  
  geweke.z<-NULL
  for(i in 1:d1){
    for(j in 1:d2){
      mcmc.obj<-mcmc(ALPHA[i,j,])
      geweke<-try(geweke.diag(mcmc.obj),silent=TRUE)
      if(!is(geweke,"try-error")){
        geweke.z<-c(geweke.z,as.numeric(geweke$z))
      }else{
        K<-K-1
      }
    }
  }
  idx<-which(abs(geweke.z)>crit.val) # if z is smaller or greater than 1.96 there is evidence that the means of both distributions are different
  xx<-paste(length(idx), " out of ",K, " variables' z-values exceed the 1.96 threshold", " (", round(length(idx)/K*100,2),"%).",sep="")
  return <- structure(list(geweke.z=geweke.z,perc=xx), class="bgvar.CD")
  return(return)
}

#' @method print bgvar.CD
#' @export
print.bgvar.CD <- function(x, ...){
  cat(x$perc)
}

#' @name resid.corr.test
#' @export resid.corr.test
#' @title Residual Autocorrelation Test
#' @description An F-test for serial autocorrelation in the residuals.
#' @usage resid.corr.test(obj, lag.cor=1, alpha=0.95, dig1=5, dig2=3)
#' @param obj An object of class \code{bgvar}.
#' @param lag.cor The order of serial correlation to be tested for. Default is set to \code{lag.cor=1}.
#' @param alpha Significance level of test. Default is set to \code{alpha=0.95}.
#' @param dig1 Number of digits to display F-statistics and its critical values.
#' @param dig2 Number of digits to display p-values.
#' @details It is the F-test of the familiar Lagrange Multiplier (LM) statistic (see Godfrey 1978a, 1978b), also known as the 'modified LM' statistic. The null hypothesis is that \eqn{rho}, the autoregressive parameter on the residuals, equals 0 indicating absence of serial autocorrelation. For higher order serial correlation, the null is that all \eqn{rho}'s jointly are 0. The test is implemented as in Vanessa Smith's and Alessandra Galesi's ''GVAR toolbox 2.0 User Guide'', page 129.
#' @return Returns a list with the following objects \describe{
#' \item{\code{Fstat}}{ contains a list of length \code{N} with the associated F-statistic for each variable in each country.}
#' \item{\code{resTest}}{ contains a matrix of size 2N times K+3, with the F-statistics for each country and each variable.}
#' \item{\code{p.res}}{ contains a table which summarizes the output.}
#' \item{\code{pL}}{ contains a list of length \code{N} with the associated p-values for each variable in each country.}
#' }
#' @author Martin Feldkircher
#' @seealso \code{\link{bgvar}} for estimation of a \code{bgvar} object.
#' @references 
#' Godfrey, L.G. (1978a) \emph{Testing Against General Autoregressive and Moving Average Error Models When the Regressors Include Lagged Dependent Variables.} Econometrica, 46, pp. 1293-1302.
#' Godfrey, L.G. (1978b) \emph{Testing for Higher Order Serial Correlation in Regression Equations When the Regressors Include Lagged Dependent Variables.} Econometrica, 46, pp. 1303-1310.
#' Smith, L. V. and A. Galesi (2014) \emph{GVAR Toolbox 2.0 User Guide}, available at \url{https://sites.google.com/site/gvarmodelling/gvar-toolbox}.
#' @examples
#' \donttest{
#' library(BGVAR)
#' data(testdata)
#' model.ng <- bgvar(Data=testdata,W=W.test,draws=100,burnin=100)
#' resid.corr.test(model.ng)
#' }
#' @importFrom stats pf qf
resid.corr.test=function(obj, lag.cor=1, alpha=0.95, dig1=5, dig2=3){
  # Residual correlation test
  # Tests the residuals of the country VECM models for serial autocorrelation
  # Input arguments:
  #                  x......GVAR object
  #                  lag....the degree of serial autocorrelation that should be tested for
  #                  alpha..significance level
  #                  dig1 / dig2...nr. of digits for the test statistic / p-value
  if(!inherits(obj, "bgvar")) {stop("Please provide a `bgvar` object.")}
  # get data and arguments - note each second column of V has sign switched -> does not impact on results of F test so keep it as it is
  xglobal  <- obj$xglobal    
  res      <- obj$cc.results$res
  lags     <- obj$args$lags
  pmax     <- max(lags)
  bigT     <- nrow(xglobal)-pmax
  pidx     <- 1:lag.cor
  varNames <- colnames(xglobal)
  cN       <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x) x[1]))
  vars     <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x) x[2]))
  
  # helper function to construct W projector matrix
  w.t<-function(x,lag=1){
    x.n<-c(rep(0,lag),rev(rev(x)[-c(1:lag)]))
    return(x.n)
  }
  
  # Calculate F-Statistic in a loop
  Fstat<-critL<-pL<-dofL<-list() # list objects since not for every country some nr. of regressors
  for(cc in 1:length(cN)){
    idx   <- grep(paste("^",cN[cc],".",sep=""),varNames) 
    X.dat <- xglobal[-c(1:pmax),idx,drop=FALSE]
    r.dat <- res[[cN[cc]]]
    ki    <- ncol(X.dat)
    dof   <- (bigT-ki-lag.cor)
    M     <- diag(bigT)-tcrossprod(X.dat%*%chol2inv(chol(crossprod(X.dat))),X.dat)
    # construct W matrix
    w.array <- array(0,dim=c(bigT,ki,lag.cor))
    faux<-critV<-pV<-NULL
    for(j in 1:ki){ # for all variables
      w<-NULL
      for(p in 1:lag.cor){
        w<-cbind(w,w.t(r.dat[,j],lag=p))
      }
      aux  <- bigT*((crossprod(r.dat[,j],w)%*%solve(t(w)%*%M%*%w)%*%(t(w)%*%r.dat[,j]))/crossprod(r.dat[,j]))
      faux <- c(faux,(dof/lag.cor)*(aux/(bigT-aux)))
      pV   <- c(pV,c(1-pf(aux,lag.cor,dof)))
    }
    names(faux) <- sapply(strsplit(colnames(X.dat),".",fixed=TRUE),function(x) x[[2]])
    Fstat[[cc]] <- faux
    critL[[cc]] <- critV <-c(critV,qf(alpha, lag.cor,dof))
    pL[[cc]]    <- pV
    dofL[[cc]]  <- dof
  }
  names(Fstat) <- names(pL) <- cN
  
  # Generate Output
  resTest<-array("-",dim=c(length(cN)*2,(length(vars)+3)))
  colnames(resTest)<-c("Country","DoF",paste("F-crit."," (",alpha,")",sep=""),vars)
  arrayIdx<-(1:(length(cN)*2))[c(TRUE, FALSE)]
  for(i in 1:length(arrayIdx)){
    resTest[arrayIdx[i],1:3]<-c(cN[i],paste("F(",lag.cor," ,",dofL[[i]],")",sep=""),format(round(critL[[i]],dig1)))
    resTest[arrayIdx[i],names(Fstat[[i]])]<-c(format(round(Fstat[[i]],dig1)))
    resTest[arrayIdx[i]+1,names(Fstat[[i]])]<-paste("(",format(round(pL[[i]],dig2)),")",sep="")
  }
  # Generate p-table
  pp<-unlist(pL);K<-length(pp)
  p.res<-matrix(0,nrow=4,ncol=2)
  p.res[1,]<-c(length(which(pp>0.10)),paste(round((length(which(pp>0.10))/K)*100,2),"%",sep=""))
  p.res[2,]<-c(length(which(pp>0.05& pp<=0.1)),paste(round((length(which(pp>0.05& pp<=0.1))/K)*100,2),"%",sep=""))
  p.res[3,]<-c(length(which(pp>0.01& pp<=0.05)),paste(round((length(which(pp>0.01& pp<=0.05))/K)*100,2),"%",sep=""))
  p.res[4,]<-c(length(which(pp<=0.01)),paste(round((length(which(pp<=0.01))/K)*100,2),"%",sep=""))
  rownames(p.res)<-c(">0.1","0.05-0.1","0.01-0.05","<0.01")
  colnames(p.res)<-c("# p-values", "in %")
  
  res<-list(Fstat=Fstat,resTest=resTest,p.res=p.res,pL=pL)
  return(res)
}

#' @name matrix_to_list
#' @export
#' @title Convert Input Matrix to List
#' @description Converts a big input matrix to an appropriate list for use of \code{bgvar}. 
#' @usage matrix_to_list(datamat)
#' @details Note the naming convention. Columns should indicate entity and variable name, separated by a dot, e.g. \code{US.y}.
#' @param datamat A matrix of size T times K, where T are time periods and K total amount of variables.
#' @return returns a list of length \code{N} (number of entities).
#' @author Maximilian Boeck
#' @seealso \code{\link{bgvar}} for estimation of a \code{bgvar} object.
#' @importFrom stats time
matrix_to_list <- function(datamat){
  if(any(is.na(datamat))){
    stop("The data you have submitted contains NAs. Please check the data.")
  }
  if(!all(grepl("\\.",colnames(datamat)))){
    stop("Please seperate country- and variable names with a point.")
  }
  cN <- unique(unlist(lapply(strsplit(colnames(datamat),".",fixed=TRUE),function(l) l[1])))
  N  <- length(cN)
  if(!all(nchar(cN)>1)){
    stop("Please provide entity names with minimal two characters.")
  }
  datalist <- list()
  for(cc in 1:N){
    datalist[[cN[cc]]] <- datamat[,grepl(cN[cc],colnames(datamat)),drop=FALSE]
    colnames(datalist[[cN[cc]]]) <- unlist(lapply(strsplit(colnames(datalist[[cN[cc]]]),".",fixed=TRUE),function(l)l[2]))
  }
  return(datalist)
}

#' @name list_to_matrix
#' @export
#' @title Convert Input List to Matrix
#' @description Converts a list to an appropriate input matrix for use of \code{bgvar}. 
#' @usage list_to_matrix(datalist)
#' @details Note the naming convention. Columns should indicate entity and variable name, separated by a dot, e.g. \code{US.y}.
#' @param datalist A list of length \code{N} which contains each a matrix of size T times k, where T are time periods and k variables per entity.
#' @return Returns a matrix of size T times K (number of time periods times number of total variables).
#' @author Maximilian Boeck
#' @seealso \code{\link{bgvar}} for estimation of a \code{bgvar} object.
#' @importFrom stats time
list_to_matrix <- function(datalist){
  if(any(unlist(lapply(datalist,function(l)any(is.na(l)))))){
    stop("The data you have submitted contains NAs. Please check the data.")
  }
  if(!all(nchar(names(datalist))>1)){
    stop("Please provide entity names with minimal two characters.")
  }
  cN <- names(datalist)
  N  <- length(cN)
  datamat<-NULL
  cc<-NULL
  for(i in 1:N){
    datamat<-cbind(datamat,datalist[[i]])
    cc<-c(cc,paste0(cN[i],".",colnames(datalist[[i]]),sep=""))
  }
  colnames(datamat)<-cc
  return(datamat)
}

#' @name excel_to_list
#' @export
#' @title Read Data from Excel
#' @description Reads a spreadsheet from excel and converts it to a list for use of \code{bgvar}.
#' @usage excel_to_list(file, first_column_as_time=TRUE, skipsheet=NULL, ...)
#' @details Note that each sheet has to be named for a respective country. Column names are used as variable names. Reader uses the \code{readxl} R package, hence additional arguments can be passed to the function. Furthermore, if \code{first_column_as_time=TRUE} then the column name has also to be time.
#' @param file A path to the file.
#' @param first_column_as_time Logical indicating whether the first column indicates the time.
#' @param skipsheet If one or more sheets should be skipped for reading, this can be provided with this argument. Either a vector of numeric indices or a vector of strings.
#' @param ... Additional arguments.
#' @return Returns a list of length \code{N} which contains each a matrix of size T times k, where T are time periods and k variables per entity.
#' @author Maximilian Boeck
#' @seealso \code{\link{bgvar}} for estimation of a \code{bgvar} object.
#' @importFrom readxl excel_sheets excel_format read_xls read_xlsx
#' @importFrom xts xts
excel_to_list <- function(file, first_column_as_time=TRUE, skipsheet=NULL, ...){
  if(!file.exists(file))
    stop("The provided file does not exist.")
  if(!grepl("(xls|xlsx)$",file))
    stop("Please provide a path to an excel filesheet (ending with xls/xlsx).")
  
  skiptype <- typeof(skipsheet)
  if(!(skiptype %in% c("numeric","character","NULL")))
    stop("Please provide skipsheet argument in right format.")
  cN <- excel_sheets(file)
  if(skiptype == "character")
    cN <- cN[!cN%in%skipsheet]
  if(skiptype == "numeric")
    cN <- cN[-skipsheet]
  
  format <- excel_format(file)
  datalist <- list()
  for(cc in cN){
    if(format == "xls"){
      temp <- read_xls(path = file, sheet = cc, col_names = TRUE, ...)
    }else if(format == "xlsx"){
      temp <- read_xlsx(path = file, sheet = cc, col_names = TRUE, ...)
    }
    if(first_column_as_time){
      if(typeof(as.matrix(temp[,1])) != "character")
        stop(paste0("Please provide as first column in sheet ",cc," as time in character format."))
      time <- as.Date(c(as.matrix(temp[,1])))
      temp <- as.matrix(temp[,2:ncol(temp)])
      temp <- xts(temp, order.by=time)
    }else{
      temp <- as.matrix(temp)
    }
    
    datalist[[cc]] <- temp
  }
  
  return(datalist)
}

