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
#'   \item{\code{W.trade0012}}{ Weight matrix based on trade flows, rowsums equal unity.}
#'   \item{\code{W.list}}{ List of ten weight matrices, described in Feldkircher and Huber (2016).}
#'   \item{\code{eerData}}{ is a list object of length \code{N} containing \itemize{
#'   \item{\code{y}}{ Real GDP, average of 2005=100. Seasonally adjusted, in logarithms.}
#'   \item{\code{Dp}}{ Consumer prices (period-on-period). CPI seasonally adjusted, in logarithm.}
#'   \item{\code{stir}}{ Short-term interest rate, typically 3-months money market rate.}
#'   \item{\code{ltir}}{ Long-term interest rates, typically 10-year government bond yields.}
#'   \item{\code{reer}}{ Real effective exchange rate, deflated by consumer prices.}
#'   \item{\code{tb}}{ Trade balance (ratio of real exports to real imports).}
#'   \item{\code{poil}}{ Price of oil, seasonally adjusted, in logarithms.}
#'   }}
#'   \item{\code{USexpectations}}{ is a time series object containing US expectations data: \itemize{
#'   \item{\code{y_t+4}}{ Four-quarter ahead expectation of Real GDP growth.}
#'   \item{\code{Dp_t+4}}{ Four-quarter ahead expectation of consumer price inflation.}
#'   \item{\code{stir_t+4}}{ Four-quarter ahead expectation of short-term interest rates.}
#'   }}
#' }
#' @aliases W.list W.trade0012 USexpectations
#' @docType data
"eerData"

#' @title Example data set to show functionality of the package
#' @description This data set is a subset of \code{eerData} containing just three countries with 76 quarterly observations, spanning the period from 1995Q1 to 2013Q4. The country coverage are the United States, the United Kingdom and the Euro area (EA) as a regional aggregate.
#' @format The data loads two objects \code{eerDatasmall}, which is a list object of length \code{N} (i.e, the number of countries) and \code{W.trade0012}, which is an \code{N} times \code{N} weight matrix with rowsums summing up to unity and zero elements on its diagonal. The global variable, oil prices, is included in the US country model as e.g., in Dees et al. (2007). The countries are abbreviated using ISO-2 codes. The weight matrix corresponds to average annual bilateral trade flows (including services) over the period from 2000 to 2012.\code{eerDatasmall} contains the country data, for more details, see below:
#' \describe{
#'   \item{\code{W.test}}{ Weight matrix based on trade flows, rowsums equal unity.}
#'   \item{\code{testdata}}{ List object of length \code{N} containing \itemize{
#'         \item{\code{y}}{ Real GDP, average of 2005=100. Seasonally adjusted, in logarithms.}
#'         \item{\code{Dp}}{ Consumer prices (period-on-period). CPI seasonally adjusted, in logarithm.}
#'         \item{\code{stir}}{ Short-term interest rate, typically 3-months money market rate.}
#'         \item{\code{ltir}}{ Long-term interest rates, typically 10-year government bond yields.}
#'         \item{\code{reer}}{ Real effective exchange rate, deflated by consumer prices.}
#'         \item{\code{tb}}{ Trade balance (ratio of real exports to real imports).}
#'         \item{\code{poil}}{ Price of oil, seasonally adjusted, in logarithms.}
#'   }}
#' }
#' @aliases W.test
#' @docType data
"testdata"

#' @title Monthly EU / G8 countries macroeconomic dataset
#' @description This data set contains monthly observations on industrial production, consumer price indices, short- and long-term interest rates, real effective exchange rates and equity prices. The time period covered is from January 2000 to December 2015 and the country coverage amounts to  28 countries -- roughly corresponding to EU member states + G-8 countries and a country model to model common monetary policy in the euro area.
#' @format The data loads three objects \code{monthlyData}, which is a list object of length \code{N+1} (i.e, the number of countries + the ECB country model), \code{W}, which is an \code{N} times \code{N} weight matrix with rowsums summing up to unity and zero elements on its diagonal. The countries are abbreviated using ISO-2 codes. The weight matrix corresponds to average annual input output flows for the \code{N} countries over the period from 2000 to 2014. The data are from the world input output table database (\url{http://www.wiod.org/home}) and are fully described in Timmerman et al. (2015). \code{monthlyData} contains the country data. Per default, variables that should affect all countries (global variables) are treated as endogenous variables in the US country model (\code{poil}, \code{pcom}, \code{vix}). Akin to Georgiadis (2015), interest setting in the euro area is modeled by a Taylor rule that includes ppp-weighted output and prices of euro area countries. The euro area interest rate enters other country models as an additional exogenous variable. For more details, see below:
#' \itemize{
#' \item{W} {Weight matrix, rowsums equal unity.}
#' \item{monthlyData} { List object containing \itemize{
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
#' @description This data set contains quarterly observations by country, spanning the period from 1979Q2 to 2019Q4. It can be downloaded from \url{https://www.mohaddes.org/gvar}. The country coverage is 33 countries.
#' @format The data loads \code{pesaranData}, which is a list object of length \code{N} (i.e, the number of countries) and contains the country-level data as described in Mohaddes and Raissi (2020). \code{EAData} contains the same data but euro area countries (AT, BE, DE, ES, FI, FR, IT, NL) are aggregated into a euro area block using purchasing power parities (\code{ppp}). The countries are abbreviated using ISO-2 codes. Furthermore, we also provide two datasets with first differences of some variables in \code{pesarnDiff} and \code{EAdiff}. \code{dominant} contains data that is considered global. \code{tA} is a three-dimensional array that contains \code{N} times \code{N} annual trade flow matrices over the period from 1980 to 2016. This array can be used to construct weight matrices. For more details, see below:
#' \describe{
#'   \item{\code{W.8016}}{ Weight matrix for the \code{pesaran.level} and \code{pesaran.diff} data sets, based on averaged trade flows covering the period 1980 to 2016 (based on \code{tA}).}
#'   \item{\code{tA}}{ Three-dimensional array that contains the yearly, bilateral trade flows, which were used to construct \code{W.8016}.}
#'   \item{\code{peseranData}}{ List object of length \code{N} containing \itemize{
#'     \item{\code{y}}{ Real GDP.}
#'     \item{\code{Dp}}{ Consumer price inflation.}
#'     \item{\code{r}}{ Short-term interest rate, typically 3-months money market rate.}
#'     \item{\code{lr}}{ Long-term interest rate.}
#'     \item{\code{eq}}{ Equity prices.}
#'     \item{\code{ep}}{ Exchange rate vis a vis the US dollar, deflated by the domestic CPI.}}}
#'  \item{\code{pesaranDiff}}{ List object of length \code{N} containing \itemize{
#'     \item{\code{y}}{ Growth rate of real GDP.}
#'     \item{\code{Dp}}{ First differences of consumer price inflation.}
#'     \item{\code{r}}{ First differences of short-term interest rate, typically 3-months money market rate.}
#'     \item{\code{lr}}{ Long-term interest rate.}
#'     \item{\code{eq}}{ Equity prices.}
#'     \item{\code{ep}}{ Exchange rate vis a vis the US dollar, deflated by the domestic CPI.}}}
#'  \item{\code{dominant}}{ Data set containing global variables: \itemize{
#'     \item{\code{poil}}{ Oil prices.}
#'     \item{\code{pmetal}}{ Metal price index.}
#'     \item{\code{pmat}}{ Agricultural price index.}}}
#' }
#' @aliases pesaranDiff EAData EADiff W.8016 W.EA.9916 tA tA.EA ppp dominant
#' @references 
#' Mohaddes, K. and M. Raissi (2018). Compilation, Revision and Updating of the Global VAR (GVAR) Database, 1979Q2-2016Q4. University of Cambridge: Faculty of Economics (mimeo).
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
