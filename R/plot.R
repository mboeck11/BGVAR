#' @name plot
#' @title Graphical Summary of Output Created with \code{bgvar}
#' @description Plotting function for fitted values, residuals, predictions, impulse responses and forecast error variance decompositions created with the \code{BGVAR} package.
#' @param x Either an object of class \code{bgvar}, \code{bgvar.res}, \code{bgvar.irf}, \code{bgvar.predict} or \code{bgvar.fevd}.
#' @param ... Additional arguments; set graphical parameters.
#' @param resp If only a subset of variables or countries should be plotted. If set to default value \code{NULL} all countries/variables are plotted.
#' @param global If \code{TRUE} global fitted values are plotted, otherwise country fitted values.
#' @param quantiles Numeric vector with posterior quantiles. Default is set to plot median along with 68\%/80\% confidence intervals.
#' @return No return value.
#' @author Maximilian Boeck, Martin Feldkircher
#' @export
#' @examples
#' \donttest{
#' library(BGVAR)
#' data(testdata)
#' model.ssvs <- bgvar(Data=testdata,W=W.test,plag=1,draws=100,burnin=100,
#'                     prior="SSVS",eigen=TRUE)
#'
#' # example for class 'bgvar'
#' plot(model.ssvs, resp=c("EA.y","US.Dp"))
#' }
#' @importFrom graphics axis lines par plot abline matplot polygon segments
#' @importFrom stats median quantile plot.ts
plot.bgvar <- function(x, ..., resp=NULL, global=TRUE){
  # reset user par settings on exit
  oldpar   <- par(no.readonly=TRUE)
  on.exit(par(oldpar))
  plag     <- x$args$plag
  xglobal  <- x$xglobal
  trend    <- x$args$trend
  XX       <- .mlag(xglobal,plag[1])
  YY       <- xglobal[-c(1:plag[1]),,drop=FALSE]
  XX       <- cbind(XX[-c(1:plag[1]),,drop=FALSE],1)
  bigT     <- nrow(YY)
  if(trend) XX <- cbind(XX,seq(1,bigT))
  time     <- .timelabel(x$args$time)
  varNames <- dimnames(xglobal)[[2]]
  cN   <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x) x[1]))
  vars <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x) x[2]))
  bigK <- length(vars)
  Ki   <- unlist(lapply(cN,function(x)length(grep(x,varNames))))
  if(global){
    A_post <- apply(x$stacked.results$A_large,c(1,2),median)
    fit    <- XX%*%t(A_post)
  }else{
    fit <- YY-do.call("cbind",x$cc.results$res)
  }
  # adapt styles
  bgvar.env$plot$cex.axis = 1.1 # adjust for this particular plot
  args <- list(...)
  args.env <- names(bgvar.env$plot)
  if(length(args)>0){
    for(aa in args.env){
      if(aa%in%names(args)) bgvar.env$plot[[aa]] = args[[aa]]
    }
  }
  if(is.null(resp)){
    nrc  <- lapply(Ki,function(k).get_nrc(k))
    for(cc in 1:length(nrc)){
      par(mar=bgvar.env$mar,mfrow=c(nrc[[cc]][1],nrc[[cc]][2]))
      for(kk in 1:bigK){
        idx <- which(paste0(cN[cc],".",vars[kk])==varNames)
        if(length(idx) == 0) next
        lims <- c(min(fit[,idx],YY[,idx]),max(fit[,idx],YY[,idx]))
        plot.ts(fit[,idx], type="l", xlab="", ylab="", main = varNames[idx], ylim=lims,
                xaxt="n",yaxt="n", cex.main=bgvar.env$plot$cex.main, cex.lab=bgvar.env$plot$cex.lab, 
                lwd=3)
        lines(YY[,idx],col="grey40", lwd=3, lty=2)
        axisindex <- round(seq(1,bigT,length.out=8))
        axis(1, at=axisindex, labels=time[axisindex], las=2, cex.axis=bgvar.env$plot$cex.axis, cex.lab=bgvar.env$plot$cex.lab)
        axis(2, cex.axis=bgvar.env$plot$cex.axis, cex.lab=bgvar.env$plot$cex.lab)
        abline(v=axisindex,col=bgvar.env$plot$col.tick,lty=bgvar.env$plot$lty.tick)
      }
    }
  }else if(all(resp%in%cN)){
    cidx <- which(cN%in%resp)
    cN   <- cN[cidx]; Ki <- Ki[cidx]
    nrc  <- lapply(Ki,function(k).get_nrc(k))
    for(cc in 1:length(nrc)){
      par(mar=bgvar.env$mar,mfrow=c(nrc[[cc]][1],nrc[[cc]][2]))
      for(kk in 1:bigK){
        idx <- which(paste0(cN[cc],".",vars[kk])==varNames)
        if(length(idx) == 0) next
        lims <- c(min(fit[,idx],YY[,idx]),max(fit[,idx],YY[,idx]))
        plot.ts(fit[,idx], type="l", xlab="", ylab="", main = varNames[idx], ylim=lims,
                xaxt="n",yaxt="n", cex.main=bgvar.env$plot$cex.main, cex.lab=bgvar.env$plot$cex.lab, 
                lwd=3)
        lines(YY[,idx],col="grey40", lwd=3, lty=2)
        axisindex <- round(seq(1,bigT,length.out=8))
        axis(1, at=axisindex, labels=time[axisindex], las=2, cex.axis=bgvar.env$plot$cex.axis, cex.lab=bgvar.env$plot$cex.lab)
        axis(2, cex.axis=bgvar.env$plot$cex.axis, cex.lab=bgvar.env$plot$cex.lab)
        abline(v=axisindex,col=bgvar.env$plot$col.tick,lty=bgvar.env$plot$lty.tick)
      }
    }
  }else if(all(resp%in%vars)){
    vidx <- which(vars%in%resp)
    vars <- vars[vidx]; Ki <- rep(length(cN),length(vidx))
    nrc  <- lapply(Ki,function(k).get_nrc(k))
    for(vv in 1:length(nrc)){
      par(mar=bgvar.env$mar,mfrow=c(nrc[[vv]][1],nrc[[vv]][2]))
      for(kk in 1:Ki[vv]){
        idx <- which(paste0(cN[kk],".",vars[vv])==varNames)
        if(length(idx)==0) next
        lims <- c(min(fit[,idx],YY[,idx]),max(fit[,idx],YY[,idx]))
        plot.ts(fit[,idx], type="l", xlab="", ylab="", main = varNames[idx], ylim=lims,
                xaxt="n",yaxt="n", cex.main=bgvar.env$plot$cex.main, cex.lab=bgvar.env$plot$cex.lab, 
                lwd=3)
        lines(YY[,idx],col="grey40", lwd=3, lty=2)
        axisindex <- round(seq(1,bigT,length.out=8))
        axis(1, at=axisindex, labels=time[axisindex], las=2, cex.axis=bgvar.env$plot$cex.axis, cex.lab=bgvar.env$plot$cex.lab)
        axis(2, cex.axis=bgvar.env$plot$cex.axis, cex.lab=bgvar.env$plot$cex.lab)
        abline(v=axisindex,col=bgvar.env$plot$col.tick,lty=bgvar.env$plot$lty.tick)
      }
    }
  }else if(all(resp%in%varNames)){
    ridx <- which(varNames%in%resp)
    Ki <- length(ridx)
    nrc <- .get_nrc(Ki)
    par(mar=bgvar.env$mar,mfrow=c(nrc[1],nrc[2]))
    for(kk in 1:Ki){
      idx <- ridx[kk]
      lims <- c(min(fit[,idx],YY[,idx]),max(fit[,idx],YY[,idx]))
      plot.ts(fit[,idx], type="l", xlab="", ylab="", main = varNames[idx], ylim=lims,
              xaxt="n",yaxt="n", cex.main=bgvar.env$plot$cex.main, cex.lab=bgvar.env$plot$cex.lab, 
              lwd=3)
      lines(YY[,idx],col="grey40", lwd=3, lty=2)
      axisindex <- round(seq(1,bigT,length.out=8))
      axis(1, at=axisindex, labels=time[axisindex], las=2, cex.axis=bgvar.env$plot$cex.axis, cex.lab=bgvar.env$plot$cex.axis)
      axis(2, cex.axis=bgvar.env$plot$cex.axis, cex.lab=bgvar.env$plot$cex.lab)
      abline(v=axisindex,col=bgvar.env$plot$col.tick,lty=bgvar.env$plot$lty.tick)
    }
  }else{
    stop("Please specify 'resp' either as one or more specific variable names in the dataset, as general variable name or as unit name, but not as a combination therof. Respecify.")
  }
  return(invisible(x))
}

#' @name plot
#' @param global If \code{global=TRUE} global residuals are plotted, otherwise country residuals.
#' @param resp Default to \code{NULL}. Either specify a single country or a group of variables to be plotted.
#' @export
#' @examples
#' \donttest{
#' # example for class 'bgvar.resid'
#' res <- residuals(model.ssvs)
#' plot(res, resp="EA.y")
#' }
plot.bgvar.resid <- function(x, ..., resp=NULL, global=TRUE){
  # reset user par settings on exit
  oldpar   <- par(no.readonly=TRUE)
  on.exit(par(oldpar))
  bigT     <- nrow(x$Data)
  time     <- .timelabel(rownames(x$Data))
  varNames <- dimnames(x$Data)[[2]]
  cN       <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x) x[1]))
  vars     <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x) x[2]))
  bigK     <- length(vars)
  Ki       <- unlist(lapply(cN,function(x)length(grep(x,varNames))))
  if(global){
    res <- apply(x$global,c(2,3),median)
  }else{
    res <- apply(x$country,c(2,3),median)
  }
  # adapt styles
  args <- list(...)
  args.env <- names(bgvar.env$plot)
  if(length(args)>0){
    for(aa in args.env){
      if(aa%in%names(args)) bgvar.env$plot[[aa]] = args[[aa]]
    }
  }
  if(is.null(resp)){
    nrc  <- lapply(Ki,function(k).get_nrc(k))
    for(cc in 1:length(nrc)){
      par(mar=bgvar.env$mar,mfrow=c(nrc[[cc]][1],nrc[[cc]][2]))
      for(kk in 1:bigK){
        idx <- which(paste0(cN[cc],".",vars[kk])==varNames)
        if(length(idx) == 0) next
        lims <- c(min(res[,idx]),max(res[,idx]))
        plot.ts(res[,idx], type="l", xlab="", ylab="", main = varNames[idx], ylim=lims,
                xaxt="n",yaxt="n", cex.main=bgvar.env$plot$cex.main, cex.lab=bgvar.env$plot$cex.lab, 
                lwd=3)
        axisindex <- round(seq(1,bigT,length.out=8))
        axis(1, at=axisindex, labels=time[axisindex], las=2, cex.axis=bgvar.env$plot$cex.axis, cex.lab=bgvar.env$plot$cex.lab)
        axis(2, cex.axis=bgvar.env$plot$cex.axis, cex.lab=bgvar.env$plot$cex.lab)
        abline(v=axisindex,col=bgvar.env$plot$col.tick,lty=bgvar.env$plot$lty.tick)
      }
    }
  }else if(all(resp%in%cN)){
    cidx <- which(cN%in%resp)
    cN   <- cN[cidx]; Ki <- Ki[cidx]
    nrc  <- lapply(Ki,function(k).get_nrc(k))
    for(cc in 1:length(nrc)){
      par(mar=bgvar.env$mar,mfrow=c(nrc[[cc]][1],nrc[[cc]][2]))
      for(kk in 1:bigK){
        idx <- which(paste0(cN[cc],".",vars[kk])==varNames)
        if(length(idx) == 0) next
        lims <- c(min(res[,idx]),max(res[,idx]))
        plot.ts(res[,idx], type="l", xlab="", ylab="", main = varNames[idx], ylim=lims,
                xaxt="n",yaxt="n", cex.main=bgvar.env$plot$cex.main, cex.lab=bgvar.env$plot$cex.lab, 
                lwd=3)
        axisindex <- round(seq(1,bigT,length.out=8))
        axis(1, at=axisindex, labels=time[axisindex], las=2, cex.axis=bgvar.env$plot$cex.axis, cex.lab=bgvar.env$plot$cex.lab)
        axis(2, cex.axis=bgvar.env$plot$cex.axis, cex.lab=bgvar.env$plot$cex.lab)
        abline(v=axisindex,col=bgvar.env$plot$col.tick,lty=bgvar.env$plot$lty.tick)
      }
    }
  }else if(all(resp%in%vars)){
    vidx <- which(vars%in%resp)
    vars <- vars[vidx]; Ki <- rep(length(cN),length(vidx))
    nrc  <- lapply(Ki,function(k).get_nrc(k))
    for(vv in 1:length(nrc)){
      par(mar=bgvar.env$mar,mfrow=c(nrc[[vv]][1],nrc[[vv]][2]))
      for(kk in 1:Ki[vv]){
        idx <- which(paste0(cN[kk],".",vars[vv])==varNames)
        if(length(idx)==0) next
        lims <- c(min(res[,idx]),max(res[,idx]))
        plot.ts(res[,idx], type="l", xlab="", ylab="", main = varNames[idx], ylim=lims,
                xaxt="n",yaxt="n", cex.main=bgvar.env$plot$cex.main, cex.lab=bgvar.env$plot$cex.lab, 
                lwd=3)
        axisindex <- round(seq(1,bigT,length.out=8))
        axis(1, at=axisindex, labels=time[axisindex], las=2, cex.axis=bgvar.env$plot$cex.axis, cex.lab=bgvar.env$plot$cex.lab)
        axis(2, cex.axis=bgvar.env$plot$cex.axis, cex.lab=bgvar.env$plot$cex.lab)
        abline(v=axisindex,col=bgvar.env$plot$col.tick,lty=bgvar.env$plot$lty.tick)
      }
    }
  }else if(all(resp%in%varNames)){
    ridx <- which(varNames%in%resp)
    Ki <- length(ridx)
    nrc <- .get_nrc(Ki)
    par(mar=bgvar.env$mar,mfrow=c(nrc[1],nrc[2]))
    for(kk in 1:Ki){
      idx <- ridx[kk]
      lims <- c(min(res[,idx]),max(res[,idx]))
      plot.ts(res[,idx], type="l", xlab="", ylab="", main = varNames[idx], ylim=lims,
              xaxt="n",yaxt="n", cex.main=bgvar.env$plot$cex.main, cex.lab=bgvar.env$plot$cex.lab, 
              lwd=3)
      axisindex <- round(seq(1,bigT,length.out=8))
      axis(1, at=axisindex, labels=time[axisindex], las=2, cex.axis=bgvar.env$plot$cex.axis, cex.lab=bgvar.env$plot$cex.lab)
      axis(2, cex.axis=bgvar.env$plot$cex.axis, cex.lab=bgvar.env$plot$cex.lab)
      abline(v=axisindex,col=bgvar.env$plot$col.tick,lty=bgvar.env$plot$lty.tick)
    }
  }else{
    stop("Please specify 'resp' either as one or more specific variable names in the dataset, as general variable name or as unit name, but not as a combination therof. Respecify.")
  }
  return(invisible(x))
}

#' @name plot
#' @param resp Specify a variable to plot predictions.
#' @param cut Length of series to be plotted before prediction begins.
#' @examples
#' \donttest{
#' # example for class 'bgvar.pred'
#' fcast <- predict(model.ssvs,n.ahead=8)
#' plot(fcast, resp="y", cut=20)
#' }
#' @export
plot.bgvar.pred<-function(x, ..., resp=NULL, cut=40, quantiles=c(.10,.16,.50,.84,.90)){
  # reset user par settings on exit
  oldpar<- par(no.readonly=TRUE)
  on.exit(par(oldpar))
  fcast <- x$fcast
  Xdata <- x$xglobal
  hstep <- x$n.ahead
  if(!all(paste0("Q",quantiles*100)%in%dimnames(fcast)[[3]])){
    stop("Please provide available quantiles.")
  }
  thin<-nrow(Xdata)-hstep
  if(thin>cut){
    Xdata<-Xdata[(nrow(Xdata)-cut+1):nrow(Xdata),]
  }
  varNames  <- colnames(Xdata)
  cN        <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x) x[1]))
  vars      <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x) x[2]))
  K         <- length(vars)
  Ki        <- unlist(lapply(cN,function(x)length(grep(x,varNames))))
  Q         <- length(quantiles)
  if((Q %% 2) == 0){
    stop("Please provide odd numbers of quantiles: median along with intervals.")
  }
  # adapt styles
  args <- list(...)
  args.env <- names(bgvar.env$plot)
  if(length(args)>0){
    for(aa in args.env){
      if(aa%in%names(args)) bgvar.env$plot[[aa]] = args[[aa]]
    }
  }
  if(is.null(resp)){
    nrc  <- lapply(Ki,function(k).get_nrc(k))
    for(cc in 1:length(nrc)){
      par(mar=bgvar.env$mar,mfrow=c(nrc[[cc]][1],nrc[[cc]][2]))
      for(kk in 1:K){
        idx <- which(paste0(cN[cc],".",vars[kk])==varNames)
        if(length(idx) == 0) next
        x <- rbind(cbind(matrix(NA,nrow(Xdata),floor(Q/2)),Xdata[,idx],matrix(NA,nrow(Xdata),floor(Q/2))),fcast[idx,,paste0("Q",quantiles*100)])
        b <- range(x,na.rm=TRUE); b1<-b[1];b2<-rev(b)[1]
        plot.ts(x[,median(seq(Q))], col=bgvar.env$plot$col.50, lty=1, yaxt="n", xaxt="n",
                lwd=bgvar.env$plot$lwd.line,ylab="",main=varNames[idx],cex.main=bgvar.env$plot$cex.main,
                cex.axis=bgvar.env$plot$cex.axis,cex.lab=bgvar.env$plot$cex.lab,ylim=c(b1,b2))
        for(qq in 1:floor(Q/2)){
          polygon(c(1:nrow(x),rev(1:nrow(x))),c(x[,qq],rev(x[,Q-qq+1])),col=bgvar.env$plot$col.unc[qq],border=NA)
        }
        lines(c(rep(NA,cut),x[seq(cut+1,cut+hstep),median(seq(Q))]),col=bgvar.env$plot$col.50,lwd=4)
        axisnames <- c(rownames(Xdata),paste("t+",1:hstep,sep=""))
        axisindex <- round(seq(1,length(axisnames),length.out=8))
        axis(side=1, at=axisindex, labels=axisnames[axisindex], cex.axis=bgvar.env$plot$cex.axis,tick=FALSE,las=2)
        axis(side=2, cex.axis=bgvar.env$plot$cex.axis)
        abline(v=axisindex,col=bgvar.env$plot$col.tick,lty=bgvar.env$plot$lty.tick)
      }
    }
  }else if(all(resp%in%cN)){
    cidx <- which(cN%in%resp)
    cN   <- cN[cidx]; Ki <- Ki[cidx]
    nrc  <- lapply(Ki,function(k).get_nrc(k))
    for(cc in 1:length(nrc)){
      par(mar=bgvar.env$mar,mfrow=c(nrc[[cc]][1],nrc[[cc]][2]))
      for(kk in 1:K){
        idx <- which(paste0(cN[cc],".",vars[kk])==varNames)
        if(length(idx) == 0) next
        x <- rbind(cbind(matrix(NA,nrow(Xdata),floor(Q/2)),Xdata[,idx],matrix(NA,nrow(Xdata),floor(Q/2))),fcast[idx,,paste0("Q",quantiles*100)])
        b <- range(x,na.rm=TRUE); b1<-b[1];b2<-rev(b)[1]
        plot.ts(x[,median(seq(Q))], col=bgvar.env$plot$col.50, lty=1, yaxt="n", xaxt="n",
                lwd=bgvar.env$plot$lwd.line,ylab="",main=varNames[idx],cex.main=bgvar.env$plot$cex.main,
                cex.axis=bgvar.env$plot$cex.axis,cex.lab=bgvar.env$plot$cex.lab,ylim=c(b1,b2))
        for(qq in 1:floor(Q/2)){
          polygon(c(1:nrow(x),rev(1:nrow(x))),c(x[,qq],rev(x[,Q-qq+1])),col=bgvar.env$plot$col.unc[qq],border=NA)
        }
        lines(c(rep(NA,cut),x[seq(cut+1,cut+hstep),median(seq(Q))]),col=bgvar.env$plot$col.50,lwd=4)
        axisnames <- c(rownames(Xdata),paste("t+",1:hstep,sep=""))
        axisindex <- round(seq(1,length(axisnames),length.out=8))
        axis(side=1, at=axisindex, labels=axisnames[axisindex], cex.axis=bgvar.env$plot$cex.axis, tick=FALSE, las=2)
        axis(side=2, cex.axis=bgvar.env$plot$cex.axis)
        abline(v=axisindex,col=bgvar.env$plot$col.tick,lty=bgvar.env$plot$lty.tick)
      }
    }
  }else if(all(resp%in%vars)){
    vidx <- which(vars%in%resp)
    vars <- vars[vidx]; Ki <- rep(length(cN),length(vidx))
    nrc  <- lapply(Ki,function(k).get_nrc(k))
    for(vv in 1:length(nrc)){
      par(mar=bgvar.env$mar,mfrow=c(nrc[[vv]][1],nrc[[vv]][2]))
      for(kk in 1:Ki[vv]){
        idx <- which(paste0(cN[kk],".",vars[vv])==varNames)
        if(length(idx)==0) next
        x <- rbind(cbind(matrix(NA,nrow(Xdata),floor(Q/2)),Xdata[,idx],matrix(NA,nrow(Xdata),floor(Q/2))),fcast[idx,,paste0("Q",quantiles*100)])
        b <- range(x,na.rm=TRUE); b1<-b[1];b2<-rev(b)[1]
        plot.ts(x[,median(seq(Q))], col=bgvar.env$plot$col.50, lty=1, yaxt="n", xaxt="n",
                lwd=bgvar.env$plot$lwd.line,ylab="",main=varNames[idx],cex.main=bgvar.env$plot$cex.main,
                cex.axis=bgvar.env$plot$cex.axis,cex.lab=bgvar.env$plot$cex.lab,ylim=c(b1,b2))
        for(qq in 1:floor(Q/2)){
          polygon(c(1:nrow(x),rev(1:nrow(x))),c(x[,qq],rev(x[,Q-qq+1])),col=bgvar.env$plot$col.unc[qq],border=NA)
        }
        lines(c(rep(NA,cut),x[seq(cut+1,cut+hstep),median(seq(Q))]),col=bgvar.env$plot$col.50,lwd=4)
        axisnames <- c(rownames(Xdata),paste("t+",1:hstep,sep=""))
        axisindex <- round(seq(1,length(axisnames),length.out=8))
        axis(side=1, at=axisindex, labels=axisnames[axisindex], cex.axis=bgvar.env$plot$cex.axis,tick=FALSE,las=2)
        axis(side=2, cex.axis=bgvar.env$plot$cex.axis)
        abline(v=axisindex,col=bgvar.env$plot$col.tick,lty=bgvar.env$plot$lty.tick)
      }
    }
  }else if(all(resp%in%varNames)){
    ridx <- which(varNames%in%resp)
    Ki <- length(ridx)
    nrc <- .get_nrc(Ki)
    par(mar=bgvar.env$mar,mfrow=c(nrc[1],nrc[2]))
    for(kk in 1:Ki){
      idx <- ridx[kk]
      x <- rbind(cbind(matrix(NA,nrow(Xdata),floor(Q/2)),Xdata[,idx],matrix(NA,nrow(Xdata),floor(Q/2))),fcast[idx,,paste0("Q",quantiles*100)])
      b <- range(x,na.rm=TRUE); b1<-b[1];b2<-rev(b)[1]
      plot.ts(x[,median(seq(Q))], col=bgvar.env$plot$col.50, lty=1, yaxt="n", xaxt="n",
              lwd=bgvar.env$plot$lwd.line,ylab="",main=varNames[idx],cex.main=bgvar.env$plot$cex.main,
              cex.axis=bgvar.env$plot$cex.axis,cex.lab=bgvar.env$plot$cex.lab,ylim=c(b1,b2))
      for(qq in 1:floor(Q/2)){
        polygon(c(1:nrow(x),rev(1:nrow(x))),c(x[,qq],rev(x[,Q-qq+1])),col=bgvar.env$plot$col.unc[qq],border=NA)
      }
      lines(c(rep(NA,cut),x[seq(cut+1,cut+hstep),median(seq(Q))]),col=bgvar.env$plot$col.50,lwd=4)
      axisnames <- c(rownames(Xdata),paste("t+",1:hstep,sep=""))
      axisindex <- round(seq(1,length(axisnames),length.out=8))
      axis(side=1, at=axisindex, labels=axisnames[axisindex], cex.axis=bgvar.env$plot$cex.axis,tick=FALSE,las=2)
      axis(side=2, cex.axis=bgvar.env$plot$cex.axis)
      abline(v=axisindex,col=bgvar.env$plot$col.tick,lty=bgvar.env$plot$lty.tick)
    }
  }else{
    stop("Please specify 'resp' either as one or more specific variable names in the dataset, as general variable name or as unit name, but not as a combination therof. Respecify.")
  }
  return(invisible(x))
}

#' @name plot
#' @param resp Specify either a specific variable, a specific country or a specific variable in a specific country which should be plotted. If set to \code{NULL} all countries is plotted.
#' @param shock Specify the shock which should be plotted.
#' @param cumulative Default is set to \code{FALSE}. If \code{cumulative=TRUE} cumulative impulse response functions are plotted.
#' @examples
#' \donttest{
#' # example for class 'bgvar.irf'
#' shockinfo <- get_shockinfo("chol")
#' shockinfo$shock <- "US.stir"; shockinfo$scale <- +1
#' irf.chol<-irf(model.ssvs, n.ahead=24, shockinfo=shockinfo)
#' plot(irf.chol, resp="US")
#' }
#' @export
plot.bgvar.irf<-function(x, ...,resp=NULL, shock=1, quantiles=c(.10,.16,.50,.84,.90), cumulative=FALSE){
  # restore user par settings on exit
  oldpar <- par(no.readonly=TRUE)
  on.exit(par(oldpar))
  if(length(shock)!=1){
    stop("Please select only one shock.")
  }
  posterior <- x$posterior
  if(!all(paste0("Q",quantiles*100)%in%dimnames(posterior)[[4]])){
    stop("Please provide available quantiles.")
  }
  varNames  <- dimnames(posterior)[[1]]
  cN        <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x) x[1]))
  vars      <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x) x[2]))
  Ki        <- unlist(lapply(cN,function(x)length(grep(x,varNames))))
  K         <- length(vars)
  Q         <- length(quantiles)
  if((Q %% 2) == 0){
    stop("Please provide odd numbers of quantiles: median along with intervals.")
  }
  # adapt styles
  args <- list(...)
  args.env <- names(bgvar.env$plot)
  if(length(args)>0){
    for(aa in args.env){
      if(aa%in%names(args)) bgvar.env$plot[[aa]] = args[[aa]]
    }
  }
  irf_list = list(); count<-0
  if(is.null(resp)){
    nrc  <- lapply(Ki,function(k).get_nrc(k))
    for(cc in 1:length(nrc)){
      par(mar=bgvar.env$mar,mfrow=c(nrc[[cc]][1],nrc[[cc]][2]))
      for(kk in 1:K){
        
        plot_varname = paste0(cN[cc],".",vars[kk])
        idx          = which(plot_varname==varNames)
        if(length(idx) == 0) next
        
        # get plot data
        x<-posterior[idx,,shock,paste0("Q",quantiles*100),drop=TRUE] 
        if(cumulative){x<-apply(x,2,cumsum)}
        # save plot data
        irf_list[[paste0("IRF.",plot_varname)]] = x
        count = count+1
        
        # do plot
        b <- range(x);b1<-b[1];b2<-rev(b)[1]
        plot.ts(x[,median(seq(Q))], col=bgvar.env$plot$col.50, lty=1, yaxt="n", xaxt="n",
                lwd=bgvar.env$plot$lwd.line,ylab="",xlab="",main=varNames[idx],cex.main=bgvar.env$plot$cex.main,
                cex.axis=bgvar.env$plot$cex.axis,cex.lab=bgvar.env$plot$cex.lab,ylim=c(b1,b2))
        for(qq in 1:floor(Q/2)){
          polygon(c(1:nrow(x),rev(1:nrow(x))),c(x[,qq],rev(x[,Q-qq+1])),col=bgvar.env$plot$col.unc[qq],border=NA)
        }
        lines(x[,median(seq(Q))],col=bgvar.env$plot$col.50,lwd=4)
        segments(x0=1,y0=0,x1=nrow(x),y1=0,col=bgvar.env$plot$col.zero,lty=bgvar.env$plot$lty.zero,lwd=bgvar.env$plot$lwd.zero)
        axis(2, at=seq(b1,b2,length.out=5), labels=format(seq(b1,b2,length.out=5),digits=1,nsmall=1),cex.axis=bgvar.env$plot$cex.axis,las=1)
        axisindex<-seq(1,nrow(x),by=4)
        axis(side=1, las=1,at=axisindex, labels=axisindex-1, cex.axis=bgvar.env$plot$cex.axis,tick=FALSE)
        abline(v=axisindex,col=bgvar.env$plot$col.tick,lty=bgvar.env$plot$lty.tick)
      }
    }
  }else if(all(resp%in%cN)){
    cidx <- which(cN%in%resp)
    cN   <- cN[cidx]; Ki <- Ki[cidx]
    nrc  <- lapply(Ki,function(k).get_nrc(k))
    for(cc in 1:length(nrc)){
      par(mar=bgvar.env$mar,mfrow=c(nrc[[cc]][1],nrc[[cc]][2]))
      vars.cc <- sapply(strsplit(varNames[grep(cN[cc],varNames)],".",fixed=TRUE),function(x) x[2])
      for(kk in 1:K){
        
        plot_varname = paste0(cN[cc],".",vars[kk])
        idx <- which(plot_varname==varNames)
        if(length(idx) == 0) next
        
        # get plot data
        x<-posterior[idx,,shock,paste0("Q",quantiles*100),drop=TRUE] 
        if(cumulative){x<-apply(x,2,cumsum)}
        # save plot data
        irf_list[[paste0("IRF.",plot_varname)]] = x
        count = count+1
        
        # do plot
        b <- range(x);b1<-b[1];b2<-rev(b)[1]
        plot.ts(x[,median(seq(Q))], col=bgvar.env$plot$col.50, lty=1, yaxt="n", xaxt="n",
                lwd=bgvar.env$plot$lwd.line,ylab="",xlab="",main=varNames[idx],cex.main=bgvar.env$plot$cex.main,
                cex.axis=bgvar.env$plot$cex.axis,cex.lab=bgvar.env$plot$cex.lab,ylim=c(b1,b2))
        for(qq in 1:floor(Q/2)){
          polygon(c(1:nrow(x),rev(1:nrow(x))),c(x[,qq],rev(x[,Q-qq+1])),col=bgvar.env$plot$col.unc[qq],border=NA)
        }
        lines(x[,median(seq(Q))],col=bgvar.env$plot$col.50,lwd=4)
        segments(x0=1,y0=0,x1=nrow(x),y1=0,col=bgvar.env$plot$col.zero,lty=bgvar.env$plot$lty.zero,lwd=bgvar.env$plot$lwd.zero)
        axis(2, at=seq(b1,b2,length.out=5), labels=format(seq(b1,b2,length.out=5),digits=1,nsmall=1),cex.axis=bgvar.env$plot$cex.axis,las=1)
        axisindex<-seq(1,nrow(x),by=4)
        axis(side=1, las=1,at=axisindex, labels=axisindex-1, cex.axis=bgvar.env$plot$cex.axis,tick=FALSE)
        abline(v=axisindex,col=bgvar.env$plot$col.tick,lty=bgvar.env$plot$lty.tick)
      }
    }
  }else if(all(resp%in%vars)){
    vidx <- which(vars%in%resp)
    vars <- vars[vidx]; Ki <- rep(length(cN),length(vidx))
    nrc  <- lapply(Ki,function(k).get_nrc(k))
    for(vv in 1:length(nrc)){
      par(mar=bgvar.env$mar,mfrow=c(nrc[[vv]][1],nrc[[vv]][2]))
      for(kk in 1:Ki[vv]){
        plot_varname = paste0(cN[kk],".",vars[vv])
        idx <- which(plot_varname==varNames)
        if(length(idx)==0) next
        
        # get plot data
        x<-posterior[idx,,shock,paste0("Q",quantiles*100),drop=TRUE] 
        if(cumulative){x<-apply(x,2,cumsum)}
        # save plot data
        irf_list[[paste0("IRF.",plot_varname)]] = x
        count = count+1
        
        # do plot
        b <- range(x);b1<-b[1];b2<-rev(b)[1]
        plot.ts(x[,median(seq(Q))], col=bgvar.env$plot$col.50, lty=1, yaxt="n", xaxt="n",
                lwd=bgvar.env$plot$lwd.line,ylab="",xlab="",main=varNames[idx],cex.main=bgvar.env$plot$cex.main,
                cex.axis=bgvar.env$plot$cex.axis,cex.lab=bgvar.env$plot$cex.lab,ylim=c(b1,b2))
        for(qq in 1:floor(Q/2)){
          polygon(c(1:nrow(x),rev(1:nrow(x))),c(x[,qq],rev(x[,Q-qq+1])),col=bgvar.env$plot$col.unc[qq],border=NA)
        }
        lines(x[,median(seq(Q))],col=bgvar.env$plot$col.50,lwd=4)
        segments(x0=1,y0=0,x1=nrow(x),y1=0,col=bgvar.env$plot$col.zero,lty=bgvar.env$plot$lty.zero,lwd=bgvar.env$plot$lwd.zero)
        axis(2, at=seq(b1,b2,length.out=5), labels=format(seq(b1,b2,length.out=5),digits=1,nsmall=1),cex.axis=bgvar.env$plot$cex.axis,las=1)
        axisindex<-seq(1,nrow(x),by=4)
        axis(side=1, las=1,at=axisindex, labels=axisindex-1, cex.axis=bgvar.env$plot$cex.axis,tick=FALSE)
        abline(v=axisindex,col=bgvar.env$plot$col.tick,lty=bgvar.env$plot$lty.tick)
      }
    }
  }else if(all(resp%in%varNames)){
    ridx <- which(varNames%in%resp)
    Ki <- length(ridx)
    nrc <- .get_nrc(Ki)
    par(mar=bgvar.env$mar,mfrow=c(nrc[1],nrc[2]))
    for(kk in 1:Ki){
      idx <- ridx[kk]
      
      # get plot data
      x<-posterior[idx,,shock,paste0("Q",quantiles*100),drop=TRUE] 
      if(cumulative){x<-apply(x,2,cumsum)}
      # save plot data
      irf_list[[paste0("IRF.",resp[kk])]] = x
      count = count+1
      
      # do plot
      b <- range(x);b1<-b[1];b2<-rev(b)[1]
      plot.ts(x[,median(seq(Q))], col=bgvar.env$plot$col.50, lty=1, yaxt="n", xaxt="n",
              lwd=bgvar.env$plot$lwd.line,ylab="",xlab="",main=varNames[idx],cex.main=bgvar.env$plot$cex.main,
              cex.axis=bgvar.env$plot$cex.axis,cex.lab=bgvar.env$plot$cex.lab,ylim=c(b1,b2))
      for(qq in 1:floor(Q/2)){
        polygon(c(1:nrow(x),rev(1:nrow(x))),c(x[,qq],rev(x[,Q-qq+1])),col=bgvar.env$plot$col.unc[qq],border=NA)
      }
      lines(x[,median(seq(Q))],col=bgvar.env$plot$col.50,lwd=4)
      segments(x0=1,y0=0,x1=nrow(x),y1=0,col=bgvar.env$plot$col.zero,lty=bgvar.env$plot$lty.zero,lwd=bgvar.env$plot$lwd.zero)
      axis(2, at=seq(b1,b2,length.out=5), labels=format(seq(b1,b2,length.out=5),digits=1,nsmall=1),cex.axis=bgvar.env$plot$cex.axis,las=1)
      axisindex<-seq(1,nrow(x),by=4)
      axis(side=1, las=1,at=axisindex, labels=axisindex-1, cex.axis=bgvar.env$plot$cex.axis,tick=FALSE)
      abline(v=axisindex,col=bgvar.env$plot$col.tick,lty=bgvar.env$plot$lty.tick)
    }
  }else{
    stop("Please specify 'resp' either as one or more specific variable names in the dataset, as general variable name or as unit name, but not as a combination therof. Respecify.")
  }
  return(invisible(irf_list))
}

#' @name plot
#' @param k.max plots the k series with the highest for the decomposition of \code{resp}.
#' @examples
#' \donttest{
#' # example for class 'bgvar.fevd'
#' fevd.us=fevd(irf.chol,var.slct=c("US.stir"))
#' plot(fevd.us, resp="US.stir", k.max=10)
#' }
#' @export
plot.bgvar.fevd<-function(x, ..., resp, k.max=10){
  # restore user par settings on exit
  oldpar <- par(no.readonly=TRUE)
  on.exit(par(oldpar))
  fevd      <- x[[1]]
  xglobal   <- x$xglobal
  varNames  <- colnames(xglobal)
  varAll    <- varNames
  cN        <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x) x[1]))
  vars      <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x) x[2]))
  
  resp0     <- paste("Decomp. of ",resp,sep="")
  if(length(resp0)>1){
    stop("Please provide just one time series in 'resp'.")
  }
  if(!(resp0%in%dimnames(fevd)[[2]])){
    stop("Please provide time series present in dataset.")
  }
  if(is.numeric(k.max)){
    mean <- apply(fevd,c(1,2),mean)
    resp <- names(sort(mean[,resp0], decreasing=TRUE))[1:k.max]
  }
  varNames <- list()
  for(kk in 1:ceiling(k.max/10)){
    if(kk*10>k.max) kk.max <- k.max else kk.max <- kk*10
    varNames[[kk]] <- resp[((kk-1)*10+1):kk.max]
  }
  # adapt styles
  args <- list(...)
  args.env <- names(bgvar.env$plot)
  if(length(args)>0){
    for(aa in args.env){
      if(aa%in%names(args)) bgvar.env$plot[[aa]] = args[[aa]]
    }
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
      x<-fevd[idx,resp0,]
      b<-range(x); b1<-b[1]; b2<-b[2]
      plot.ts(x,col=bgvar.env$plot$col.50,xaxt="n",yaxt="n",lwd=bgvar.env$plot.lwd.line,ylab="",xlab="",
              main=varAll[idx],cex.main=bgvar.env$plot.cex.main,cex.axis=bgvar.env$plot$cex.axis,
              cex.lab=bgvar.env$plot$cex.lab,lty=1,ylim=c(b1,b2))
      
      axis(2, at=seq(b1,b2,length.out=5), labels=format(seq(b1,b2,length.out=5),digits=2,nsmall=1),cex.axis=bgvar.env$plot$cex.axis,las=1)
      axisindex<-seq(1,length(x),by=4)
      axis(side=1, las=1,at=axisindex, labels=c(0:length(x))[axisindex], cex.axis=bgvar.env$plot$cex.axis,tick=FALSE)
      abline(v=axisindex,col=bgvar.env$plot$col.tick,lty=bgvar.env$plot$lty.tick)
    }
  }
  return(invisible(x))
}
