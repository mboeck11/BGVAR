% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bgvar.R
\name{summary}
\alias{summary}
\alias{summary.bgvar}
\title{Summary of Bayesian GVAR}
\usage{
\method{summary}{bgvar}(object, ...)
}
\arguments{
\item{object}{an object of class \code{bgvar}.}

\item{...}{other arguments.}
}
\value{
No return value.
}
\description{
Output gives model information as well as some descriptive statistics on convergence properties, likelihood, serial autocorrelation in the errors and the average pairwise autocorrelation of cross-country residuals.
}
\examples{
\dontshow{
library(BGVAR)
data(eerData)
cN<-c("EA","US")
eerData<-eerData[cN]
W.trade0012<-apply(W.trade0012[cN,cN],2,function(x)x/rowSums(W.trade0012[cN,cN]))
model.ssvs <- bgvar(Data=eerData,W=W.trade0012,plag=1,draws=50,burnin=50,
                    prior="SSVS",SV=TRUE,trend=TRUE)
summary(model.ssvs)
}
\donttest{
set.seed(571)
library(BGVAR)
data(eerData)
model.ssvs <- bgvar(Data=eerData,W=W.trade0012,plag=1,draws=100,burnin=100,
                    prior="SSVS",thin=1,SV=TRUE,trend=TRUE)
summary(model.ssvs)
}
}
\seealso{
\code{\link{bgvar}} to estimate a \code{bgvar} object.
\code{\link{avg.pair.cc}} to compute average pairwise cross-country correlation of cross-country residuals separately.
\code{\link{resid.corr.test}} to compute F-test on first-order autocorrelation of cross-country residuals separately.
}
\author{
Maximilian Boeck
}