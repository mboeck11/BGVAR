#' @name BGVAR-package
#' @title BGVAR: Bayesian Global Vector Autoregressions
#' @description The Bayesian Global Vector Autoregression (BGVAR) package allows to estimate Global Vector Autoregressions and consists of various tools for predicting and doing structural analysis.
#' @details It provides a fully Bayesian implementation of Global Vector Autoregressions. It utilizes Markov chain Monte Carlo (MCMC) samplers to conduct inference by obtaining draws from the posterior distribution of parameters. One of the main advantages is the implementation of different shrinkage prior setups for estimating the model. The packages consists thus of various post-processing functions to carry out predictions or structural analysis. It is possible to perform structural identification via short-run or sign/zero restrictions. The available structural tools comprise impulse response functions, historical decompositions and forecast error variance decompositions. For all the aforementioned tools plotting functions are implemented. Furthermore, various functions of the package are intended to inspect the convergence properties of the MCMC chain and to do model evaluation. The main focus of this paper is to show the functionality of \code{BGVAR}. In addition, it provides a brief mathematical description of the model, an overview of the implemented sampling scheme, and several illustrative examples using global macroeconomic data.
#' @importFrom utils data
#' @docType package
#' @seealso 
#' \code{\link{bgvar}} for estimating a Bayesian GVAR.
#' \code{\link{predict}} for doing predictions with a Bayesian GVAR.
#' \code{\link{irf}} for doing impulse response analysis with a Bayesian GVAR.
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @useDynLib BGVAR, .registration=TRUE
NULL

#' @title Example data set to replicate Feldkircher and Huber (2016)
#' @description This data set contains 76 quarterly observations by country, spanning the period from 1995Q1 to 2013Q4. The country coverage is 43 countries and the Euro area (EA) as a regional aggregate.
#' @format The data loads two objects \code{eerData}, which is a list object of length \code{N} (i.e, the number of countries) and \code{W.trade0012}, which is an \code{N} times \code{N} weight matrix with rowsums summing up to unity and zero elements on its diagonal. The global variable, oil prices, is included in the US country model as e.g., in Dees et al. (2007). The countries are abbreviated using ISO-2 codes. The weight matrix corresponds to average annual bilateral trade flows (including services) over the period from 2000 to 2012.\code{eerData} contains the country data, for more details, see below:
#' \describe{
#'   \item{\code{W.trade0012}}{\code{N} times \code{N} weight matrix based on trade flows, rowsums equal unity.}
#'   \item{\code{W.list}}{A list of 10 weight matrices, described in Feldkircher and Huber (2016).}
#'   \item{\code{eerData}}{ is a list object of length \code{N} containing \itemize{
#'   \item{\code{y}}{ Real GDP, average of 2005=100. Seasonally adjusted, in logarithms.}
#'   \item{\code{Dp}}{ Consumer prices (period-on-period). CPI seasonally adjusted, in logarithm.}
#'   \item{\code{stir}}{ Short-term interest rate, typically 3-months money market rate.}
#'   \item{\code{ltir}}{ Long-term interest rates, typically 10-year government bond yields.}
#'   \item{\code{reer}}{ Real effective exchange rate, deflated by consumer prices.}
#'   \item{\code{tb}}{ Trade balance (ratio of real exports to real imports).}
#'   \item{\code{poil}}{ Price of oil, seasonally adjusted, in logarithms.}
#'   }}
#' }
#' @aliases W.list W.trade0012
#' @docType data
"eerData"

#' @title Example data set to show functionality of the package
#' @description This data set is a subset of \code{eerData} containing just three countries with 76 quarterly observations, spanning the period from 1995Q1 to 2013Q4. The country coverage are the United States, the United Kingdom and the Euro area (EA) as a regional aggregate.
#' @format The data loads two objects \code{eerDatasmall}, which is a list object of length \code{N} (i.e, the number of countries) and \code{W.trade0012}, which is an \code{N} times \code{N} weight matrix with rowsums summing up to unity and zero elements on its diagonal. The global variable, oil prices, is included in the US country model as e.g., in Dees et al. (2007). The countries are abbreviated using ISO-2 codes. The weight matrix corresponds to average annual bilateral trade flows (including services) over the period from 2000 to 2012.\code{eerDatasmall} contains the country data, for more details, see below:
#' \describe{
#'   \item{\code{W.trade0012.small}}{\code{N} times \code{N} weight matrix based on trade flows, rowsums equal unity.}
#'   \item{\code{eerDatasmall}}{ is a list object of length \code{N} containing \itemize{
#'         \item{\code{y}}{ Real GDP, average of 2005=100. Seasonally adjusted, in logarithms.}
#'         \item{\code{Dp}}{ Consumer prices (period-on-period). CPI seasonally adjusted, in logarithm.}
#'         \item{\code{stir}}{ Short-term interest rate, typically 3-months money market rate.}
#'         \item{\code{ltir}}{ Long-term interest rates, typically 10-year government bond yields.}
#'         \item{\code{reer}}{ Real effective exchange rate, deflated by consumer prices.}
#'         \item{\code{tb}}{ Trade balance (ratio of real exports to real imports).}
#'         \item{\code{poil}}{ Price of oil, seasonally adjusted, in logarithms.}
#'   }}
#' }
#' @aliases W.trade0012.small
#' @docType data
"eerDatasmall"

#' @title eerData extended with expectations data
#' @description This data set contains 76 quarterly observations by country, spanning the period from 1995Q1 to 2013Q4. The country coverage is 43 countries + the euro area (EA) as a regional aggregate. Additionally, the US country dataset is extended with four quarter ahead expectation data on output, prices and short-term interest rates from the Survey of Professional Forecasters.
#' @format The data loads two objects \code{eerData}, which is a list object of length \code{N} (i.e, the number of countries) and \code{W.trade0012}, which is an \code{N} times \code{N} weight matrix with rowsums summing up to unity and zero elements on its diagonal. The global variable, oil prices, is included in the US country model as e.g., in Dees et al. (2007). The countries are abbreviated using ISO-2 codes. The weight matrix corresponds to average annual bilateral trade flows (including services) over the period from 2000 to 2012.\code{eerData} contains the country data, for more details, see below:
#' \describe{
#'   \item{\code{W.trade0012spf}}{\code{N} times \code{N} weight matrix based on trade flows, rowsums equal unity.}
#'   \item{\code{eerDataspf}}{ is a list object of length \code{N} containing \itemize{
#'   \item{\code{y_t+4}}{ four quarter ahead expectation of Real GDP growth.}
#'   \item{\code{Dp_t+4}}{ four quarter ahead expectation of consumer price inflation.}
#'   \item{\code{stir_t+4}}{ four quarter ahead expectation of short-term interest rates.}
#'   \item{\code{y}}{ Real GDP growth.}
#'   \item{\code{Dp}}{ Consumer price inflation (period-on-period).}
#'   \item{\code{stir}}{ Short-term interest rate, typically 3-months money market rate.}
#'   \item{\code{ltir}}{ Long-term interest rates, typically 10-year government bond yields.}
#'   \item{\code{reer}}{ Real effective exchange rate, deflated by consumer prices.}
#'   \item{\code{tb}}{ Trade balance (ratio of real exports to real imports).}
#'   \item{\code{poil}}{ Price of oil, seasonally adjusted, in logarithms.}
#'   }}
#' }
#' 
#' @aliases W.trade0012.spf
#' @docType data
"eerDataspf"

#' @title Monthly EU / G8 countries macroeconomic dataset
#' @description This data set contains monthly observations on industrial production, consumer price indices, short- and long-term interest rates, real effective exchange rates and equity prices. The time period covered is from January 2000 to December 2015 and the country coverage amounts to  28 countries -- roughly corresponding to EU member states + G-8 countries and a country model to model common monetary policy in the euro area.
#' @format The data loads three objects \code{monthlyData}, which is a list object of length \code{N+1} (i.e, the number of countries + the ECB country model), \code{W}, which is an \code{N} times \code{N} weight matrix with rowsums summing up to unity and zero elements on its diagonal. The countries are abbreviated using ISO-2 codes. The weight matrix corresponds to average annual input output flows for the \code{N} countries over the period from 2000 to 2014. The data are from the world input output table database (\url{http://www.wiod.org/home}) and are fully described in Timmerman et al. (2015). \code{monthlyData} contains the country data. Per default, variables that should affect all countries (global variables) are treated as endogenous variables in the US country model (\code{poil}, \code{pcom}, \code{vix}). Akin to Georgiadis (2015), interest setting in the euro area is modeled by a Taylor rule that includes ppp-weighted output and prices of euro area countries. The euro area interest rate enters other country models as an additional exogenous variable. For more details, see below:
#' \itemize{
#' \item{W.} {\code{N} times \code{N} weight matrix, rowsums equal unity.}
#' \item{monthlyData} { is a list object of length \code{N} containing \itemize{
#' \item{\code{y}}{ Industrial production index, in real terms, logarithmic transform and seasonally adjusted.}
#' \item{\code{p}}{ Harmonized Consumer Price Index (HCPI) for EU member states, for other countries Consumer Price Index. Data in logarithmic transform and seasonally adjusted.}
#' \item{\code{stir}}{ Short-term interest rate, typically 3 months money market rate.}
#' \item{\code{EAstir}}{ Short-term interest rate, typically 3 months money market rate (3 months euribor).}
#' \item{\code{ltir}}{ Long term interest rates, typically 10-year government bond yields.}
#' \item{\code{er}}{ Real effective exchange rate index, deflated by consumer prices.}
#' \item{\code{eq}}{ Equity price index, in logarithmic transform.}
#' \item{\code{poil}}{ Price of oil, seasonally adjusted, in logarithms.}
#' \item{\code{pcom}}{ Commodity price index, seasonally adjusted, in logarithms.}
#' \item{\code{vix}}{ Volatility index, in logarithms.}
#' }}}
#' @aliases EA.weights OC.weights W
#' @docType data
"monthlyData"

#' @title pesaranData
#' @description This data set contains 151 quarterly observations by country, spanning the period from 1979Q2 to 2016Q4. It can be downloaded from \url{https://sites.google.com/site/gvarmodelling/gvar-toolbox}. The country coverage is 33 countries.
#' @format The data loads three objects \code{pesaranData}, which is a list object of length \code{N} (i.e, the number of countries) and \code{W.1316}, which is an \code{N} times \code{N} weight matrix with rowsums summing up to unity and zero elements on its diagonal. The global variable, oil prices, is included in the US country model as e.g., in Dees et al. (2007). The countries are abbreviated using ISO-2 codes. The weight matrix corresponds to average annual bilateral trade flows over the period from 2013 to 2016.\code{peseranData} contains the country data in logarithms, for more details, see below:
#' \describe{
#'   \item{\code{W.1316}}{\code{N} times \code{N} weight matrix based on trade flows, rowsums equal unity.}
#'   \item{\code{tA}}{\code{N} times \code{N} times \code{T}  array that contains the yearly, bilateral trade flows, which were used to construct \code{W.1316}.}
#'   \item{\code{peseranData}}{ is a list object of length \code{N} containing \itemize{
#'   \item{\code{y}}{ Real GDP.}
#'   \item{\code{Dp}}{ Consumer price inflation.}
#'   \item{\code{eq}}{ Equity prices.}
#'   \item{\code{stir}}{ Short-term interest rate, typically 3-months money market rate.}
#'   \item{\code{ltir}}{ Long-term interest rates, typically 10-year government bond yields.}
#'   \item{\code{poil}}{ Price of oil.}
#'   \item{\code{pmetal}}{ Price of metals.}
#'   \item{\code{pmat}}{ Price of agricultural products.}
#'   }}
#' }
#' @aliases W.1316 tA
#' @docType data
"pesaranData"

bgvar.env <- new.env()
bgvar.env$plot <- list(
  cex.main = 2,
  cex.axis = 2.5,
  cex.lab  = 2.5,
  col.unc  = c("grey60","grey40","grey20"),
  col.50   = "black",
  col.tick = "lightgrey",
  col.zero = "red",
  lty.zero = 2,
  lty.tick = 3,
  lwd.line = 4,
  lwd.zero = 3
)
bgvar.env$mar <- c(4.3,4.3,2.3,2.3)

.onAttach <- function(lib, pkg) {
  if(interactive() || getOption("verbose")){
    packageStartupMessage(sprintf("Package %s %s attached. To cite, see citation(\"%s\").", pkg, utils::packageDescription(pkg)$Version, pkg))
  }
}

.onUnload <- function (libpath) {
  library.dynam.unload("BGVAR", libpath)
}
