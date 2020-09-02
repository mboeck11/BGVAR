#' @name plot.bgvar
#' @title Plotting function for fitted values
#' @description Plots the fitted values in red of either the country VARs or the GVAR (default) along with the original data.
#' @param x an object of class \code{bgvar}.
#' @param ... additional arguments.
#' @param global if \code{TRUE} global fitted values are plotted, otherwise country fitted values.
#' @param resp if only a subset of variables or countries should be plotted. If set to default value \code{NULL} all countries/variables are plotted.
#' @return No return value.
#' @export
#' @examples
#' \donttest{
#' library(BGVAR)
#' data(eerData)
#' model.ssvs <- bgvar(Data=eerData,W=W.trade0012,plag=1,draws=100,burnin=100,
#'                     prior="SSVS")
#' summary(model.ssvs)
#' plot(model.ssvs, resp="EA")
#' }
#' @importFrom graphics axis lines par plot abline
#' @importFrom stats median plot.ts
plot.bgvar <- function(x, ..., global=TRUE, resp=NULL){
  # reset user par settings on exit
  oldpar   <- par(no.readonly=TRUE)
  on.exit(par(oldpar))
  plag     <- x$args$plag
  xglobal  <- x$xglobal
  trend    <- x$args$trend
  XX       <- .mlag(xglobal,plag)
  YY       <- xglobal[-c(1:plag),,drop=FALSE]
  XX       <- cbind(XX[-c(1:plag),,drop=FALSE],1)
  bigT     <- nrow(YY)
  if(trend) XX <- cbind(XX,seq(1,bigT))
  time     <- .timelabel(x$args$time)
  varNames <- dimnames(xglobal)[[2]]
  varAll   <- varNames
  cN       <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x) x[1]))
  vars     <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x) x[2]))
  max.vars <- unlist(lapply(cN,function(x)length(grep(x,varNames))))
  if(global){
    A_post <- apply(x$stacked.results$A_large,c(2,3),median)
    fit    <- XX%*%t(A_post)
  }else{
    fit <- YY-do.call("cbind",x$cc.results$res)
  }
  if(!is.null(resp)){
    resp.p <- strsplit(resp,".",fixed=TRUE)
    resp.c <- sapply(resp.p,function(x) x[1])
    resp.v <- sapply(resp.p,function(x) x[2])
    if(!all(unique(resp.c)%in%cN)){
      stop("Please provide country names corresponding to the ones of the 'bgvar' object.")
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
    varNames <- lapply(cN,function(cc)varAll[grepl(cc,varAll)])
  }
  for(cc in 1:length(cN)){
    rows <- max.vars[cc]/2
    if(rows<1) cols <- 1 else cols <- 2
    if(rows%%1!=0) rows <- ceiling(rows)
    if(rows%%1!=0) rows <- ceiling(rows)
    # update par settings
    par(mar=bgvar.env$mar,mfrow=c(rows,cols))
    for(kk in 1:max.vars[cc]){
      idx  <- grep(cN[cc],varAll)
      idx <- idx[varAll[idx]%in%varNames[[cc]]][kk]
      lims <- c(min(fit[,idx],YY[,idx]),max(fit[,idx],YY[,idx]))
      plot.ts(fit[,idx], type="l", xlab="", ylab="", main = varAll[idx], ylim=lims,
              xaxt="n",yaxt="n", cex.main=bgvar.env$plot$cex.main, cex.lab=bgvar.env$plot$cex.lab, 
              lwd=3)
      lines(YY[,idx],col="grey40", lwd=3, lty=2)
      axisindex <- round(seq(1,bigT,length.out=8))
      axis(1, at=axisindex, labels=time[axisindex], las=2, cex.axis=0.6, cex.lab=2.5)
      axis(2, cex.axis=0.6, cex.lab=2.5)
      abline(v=axisindex,col=bgvar.env$plot$col.tick,lty=bgvar.env$plot$lty.tick)
    }
    if(cc<length(cN)) readline(prompt="Press enter for next country...")
  }
}

#' @name plot
#' @title Plotting function for residuals
#' @description Either plots country-residuals or the global-residuals. 
#' @param x an object of class \code{bgvar.res}.
#' @param ... additional arguments.
#' @param global if \code{TRUE} global residuals are plotted, otherwise country residuals.
#' @param resp default to \code{NULL}. Either specify a single country or a group of variables to be plotted.
#' @return No return value.
#' @export
#' @examples
#' \donttest{
#' library(BGVAR)
#' data(eerData)
#' model.ssvs <- bgvar(Data=eerData,W=W.trade0012,plag=1,draws=100,burnin=100,
#'                     prior="SSVS")
#' summary(model.ssvs)
#' res <- residuals(model.ssvs)
#' plot(res, resp="EA")
#' }
#' @importFrom graphics abline axis lines par plot
#' @importFrom stats quantile plot.ts
plot.bgvar.resid <- function(x, ..., global=TRUE, resp=NULL){
  # reset user par settings on exit
  oldpar   <- par(no.readonly=TRUE)
  on.exit(par(oldpar))
  bigT     <- nrow(x$Data)
  time     <- .timelabel(rownames(x$Data))
  varNames <- dimnames(x$Data)[[2]]
  varAll   <- varNames
  cN       <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x) x[1]))
  vars     <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x) x[2]))
  max.vars <- unlist(lapply(cN,function(x)length(grep(x,varNames))))
  if(global){
    res <- apply(x$global,c(2,3),median)
  }else{
    res <- apply(x$country,c(2,3),median)
  }
  if(!is.null(resp)){
    resp.p <- strsplit(resp,".",fixed=TRUE)
    resp.c <- sapply(resp.p,function(x) x[1])
    resp.v <- sapply(resp.p,function(x) x[2])
    if(!all(unique(resp.c)%in%cN)){
      stop("Please provide country names corresponding to the ones of the 'bgvar' object.")
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
    varNames <- lapply(cN,function(cc)varAll[grepl(cc,varAll)])
  }
  for(cc in 1:length(cN)){
    rows <- max.vars[cc]/2
    if(rows<1) cols <- 1 else cols <- 2
    if(rows%%1!=0) rows <- ceiling(rows)
    if(rows%%1!=0) rows <- ceiling(rows)
    # update par settings
    par(mar=bgvar.env$mar,mfrow=c(rows,cols))
    for(kk in 1:max.vars[cc]){
      idx  <- grep(cN[cc],varAll)
      idx <- idx[varAll[idx]%in%varNames[[cc]]][kk]
      lims <- c(min(res[,idx]),max(res[,idx]))
      plot.ts(res[,idx], type="l", xlab="", ylab="",main = varAll[idx], ylim=lims,
              cex.main=bgvar.env$plot$cex.main, cex.lab=bgvar.env$plot$cex.lab, lwd=3, xaxt="n", yaxt="n")
      axisindex <- seq(1,bigT,length.out=8)
      axis(1, at=axisindex, labels=time[axisindex], las=2, cex.axis=0.6, cex.lab=2.5)
      axis(2, cex.axis=0.6, cex.lab=2.5)
      abline(v=axisindex,col=bgvar.env$plot$col.tick,lty=bgvar.env$plot$lty.tick)
    }
    if(cc<length(cN)) readline(prompt="Press enter for next country...")
  }
}

#' @name plot
#' @title Plot predictions of \code{bgvar}
#' @description  Plots the predictions of an object of class \code{bgvar.predict}.
#' @param x an object of class \code{bgvar.predict}.
#' @param ... additional arguments.
#' @param resp specify a variable to plot predictions.
#' @param Cut length of series to be plotted before prediction begins.
#' @return No return value.
#' @author Maximilian Boeck, Martin Feldkircher
#' @examples
#' \dontshow{
#' library(BGVAR)
#' data(eerData)
#' cN<-c("EA","US","UK")
#' eerData<-eerData[cN]
#' W.trade0012<-apply(W.trade0012[cN,cN],2,function(x)x/rowSums(W.trade0012[cN,cN]))
#' model.ssvs.eer<-bgvar(Data=eerData,W=W.trade0012,draws=100,burnin=100,plag=1,prior="SSVS",
#'                       eigen=TRUE)
#' fcast <- predict(model.ssvs.eer,n.ahead=8,save.store=TRUE)
#' plot(fcast, resp="US.Dp", Cut=20)
#' }
#' \donttest{
#' library(BGVAR)
#' data(eerData)
#' model.ssvs.eer<-bgvar(Data=eerData,W=W.trade0012,draws=100,burnin=100,plag=1,prior="SSVS",
#'                       eigen=TRUE)
#' fcast <- predict(model.ssvs.eer,n.ahead=8,save.store=TRUE)
#' plot(fcast, resp="US.Dp", Cut=20)
#' }
#' @importFrom graphics abline matplot polygon
#' @importFrom stats rnorm
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
plot.bgvar.pred<-function(x, ..., resp=NULL,Cut=40){
  # reset user par settings on exit
  oldpar<- par(no.readonly=TRUE)
  on.exit(par(oldpar))
  fcast <- x$fcast
  Xdata <- x$xglobal
  hstep <- x$n.ahead
  thin<-nrow(Xdata)-hstep
  if(thin>Cut){
    Xdata<-Xdata[(nrow(Xdata)-Cut+1):nrow(Xdata),]
  }
  varNames  <- colnames(Xdata)
  varAll    <- varNames
  cN        <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x) x[1]))
  vars      <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x) x[2]))
  if(!is.null(resp)){
    resp.p <- strsplit(resp,".",fixed=TRUE)
    resp.c <- sapply(resp.p,function(x) x[1])
    resp.v <- sapply(resp.p,function(x) x[2])
    if(!all(unique(resp.c)%in%cN)){
      stop("Please provide country names corresponding to the ones of the 'bgvar.predict' object.")
    }
    cN       <- cN[cN%in%resp.c]
    varNames <- lapply(cN,function(x)varNames[grepl(x,varNames)])
    if(all(!is.na(resp.v))){
      if(!all(unlist(lapply(resp,function(r)r%in%varAll)))){
        stop("Please provide correct variable names corresponding to the ones in the 'bgvar.predict' object.")
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
    # update par settings
    par(mar=bgvar.env$mar,mfrow=c(rows,cols))
    for(kk in 1:max.vars[cc]){
      idx  <- grep(cN[cc],varAll)
      idx <- idx[varAll[idx]%in%varNames[[cc]]][kk]
      x <- rbind(cbind(NA,Xdata[,idx],NA),fcast[idx,,c("low25","median","high75")])
      y <- rbind(cbind(NA,Xdata[,idx],NA),fcast[idx,,c("low16","median","high84")])
      b <- range(x,y, na.rm=TRUE)
      b1<-b[1];b2<-rev(b)[1]
      matplot(x,type="l",col=c("black","black","black"),xaxt="n",lwd=4,ylab="",main=varAll[idx],yaxt="n",
              cex.main=bgvar.env$plot$cex.main,cex.axis=bgvar.env$plot$cex.axis,
              cex.lab=bgvar.env$plot$cex.lab,lty=c(0,1,0),ylim=c(b1,b2))
      polygon(c(1:nrow(y),rev(1:nrow(y))),c(y[,1],rev(y[,3])),col=bgvar.env$plot$col.75,border=NA)
      polygon(c(1:nrow(x),rev(1:nrow(x))),c(x[,1],rev(x[,3])),col=bgvar.env$plot$col.68,border=NA)
      lines(c(rep(NA,Cut),x[seq(Cut+1,Cut+hstep),2]),col=bgvar.env$plot$col.50,lwd=4)
      
      axisnames <- c(rownames(Xdata),paste("t+",1:hstep,sep=""))
      axisindex <- c(round(seq(1,Cut,length.out=8)),seq(Cut+1,Cut+hstep))
      axis(side=1, at=axisindex, labels=axisnames[axisindex], cex.axis=0.6,tick=FALSE,las=2)
      axis(side=2, cex.axis=0.6)
      abline(v=axisindex,col=bgvar.env$plot$col.tick,lty=bgvar.env$plot$lty.tick)
    }
    if(cc<length(cN)) readline(prompt="Press enter for next country...")
  }
}

#' @name plot
#' @title Plot predictions of bgvar
#' @description  Plots the predictions of an object of class \code{bgvar.predict}.
#' @param x an object of class \code{bgvar.irf}.
#' @param ... additional arguments.
#' @param resp specify a variable to plot predictions.
#' @param shock.nr specify shock to be plotted.
#' @param cumulative whether cumulative impulse response functions should be plotted. Default is set to \code{FALSE}.
#' @return No return value.
#' @author Maximilian Boeck, Martin Feldkircher
#' @examples
#' \donttest{
#' library(BGVAR)
#' data(eerData)
#' model.ssvs.eer<-bgvar(Data=eerData,W=W.trade0012,draws=100,burnin=100,plag=1,prior="SSVS",
#'                       eigen=TRUE)
#' # US monetary policy shock
#' shocks<-list();shocks$var="stir";shocks$cN<-"US";shocks$ident="chol";shocks$scal=-100
#' irf.chol.us.mp<-irf(model.ssvs.eer,shock=shocks,nhor=24)
#' # plots an impulse response function
#' plot(irf.chol.us.mp,resp="US.y")
#' }
#' @seealso \code{\link{irf}}
#' @importFrom graphics abline matplot polygon segments
#' @importFrom stats rnorm
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
plot.bgvar.irf<-function(x, ...,resp=NULL,shock.nr=1,cumulative=FALSE){
  # restore user par settings on exit
  oldpar <- par(no.readonly=TRUE)
  on.exit(par(oldpar))
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
    # update par settings
    par(mar=bgvar.env$mar,mfrow=c(rows,cols))
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

#' @name plot
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
#' model.ssvs.eer<-bgvar(Data=eerData,W=W.trade0012,draws=100,burnin=100,plag=1,
#'                       prior="SSVS",thin=1,eigen=TRUE)
#'                       
#' # US monetary policy shock
#' shocks<-list();shocks$var="stir";shocks$cN<-"US";shocks$ident="chol";shocks$scal=-100
#' irf.chol.us.mp<-irf(obj=model.ssvs.eer,shock=shocks,nhor=48)
#' 
#' # calculates FEVD for variables US.Dp and EA.y
#' fevd.us.mp=fevd(obj=irf.chol.us.mp,var.slct=c("US.Dp","EA.y"))
#' 
#' plot(fevd.us.mp, ts="US.Dp", k.max=10)
#' }
#' @seealso \code{\link{fevd}} \code{\link{gfevd}}
#' @importFrom graphics abline matplot polygon
#' @importFrom MASS Null
#' @importFrom stats rnorm
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
plot.bgvar.fevd<-function(x, ..., ts, k.max=10){
  # restore user par settings on exit
  oldpar <- par(no.readonly=TRUE)
  on.exit(par(oldpar))
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
    par(mar=bgvar.env$mar,mfrow=c(rows,cols))
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