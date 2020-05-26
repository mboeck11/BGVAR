# BGVAR: Bayesian Global Vector Autoregressions

<!-- badges: start -->
<!-- badges: end -->

Estimation of Bayesian Global Vector Autoregressions with different prior setups and the possibility to introduce stochastic volatility. Built-in priors include the SIMS, SSVS and NG prior. Post-processing functions allow for doing predictions, structurally identify the model with short-run or sign-restrictions and compute impulse response function, historical decompositions and forecast error variance decompositions. Plotting functions are also available.

## Installation

You can install the released version of BGVAR from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("BGVAR")
```

## Usage

The core function of the package is `bgvar()` to estimate Bayesian Global Vector Autoregressions with different shrinkage prior setups. Calls can be heavily customized with respect to the specification details of the model, the MCMC chain, hyperparameter setup and various extra features. The output of the estimation can then be used for a variety of tools implemented for the **BGVAR** package.

Predictions are invoked with `predict()`, impulse responses are computed with `IRF()`, forecast error variance decompositions can be called with `fevd.decomp` and historical decompositions with `hd.decomp()`. Furthermore, counterfactual impulse responses are computed with `IRF.cf()` and conditional forecasts with `cond.pred()`. 

The package comes with standard methods to ease the analysis. The estimation output can be inspected with `print()`, `summary()`, `fitted()`, `coef()`, `vcov()` and `residuals()`. Default `plot()` is available for most outputs. Various other helper functions to ease analysis are also available.

## References

Crespo Cuaresma, J., Feldkircher, M. and F. Huber (2016) Forecasting with Global Vector Autoregressive Models: A Bayesian Approach. *Journal of Applied Econometrics*, Vol. 31(7), pp. 1371-1391.

Doan, T. R., Litterman, B. R. and C. A. Sims (1984) Forecasting and Conditional Projection Using Realistic Prior Distributions. *Econometric Reviews*, Vol. 3, pp. 1-100.

George, E.I., Sun, D. and S. Ni (2008) Bayesian stochastic search for var model restrictions. *Journal of Econometrics*, Vol. 142, pp. 553-580.

Huber, F. and M. Feldkircher (2016) Adaptive Shrinkage in Bayesian Vector Autoregressive Models. *Journal of Business and Economic Statistics*, Vol. 37(1), pp. 27-39.

Pesaran, M.H., Schuermann T. and S.M. Weiner (2004) Modeling Regional Interdependencies Using a Global Error-Correcting Macroeconometric Model. *Journal of Business and Economic Statistics*, Vol. 22, pp. 129-162.

