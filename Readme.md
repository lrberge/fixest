
# fixest: Fast and user-friendly fixed-effects estimation

<a href="https://cran.r-project.org/web/checks/check_results_fixest.html"><img src="https://cranchecks.info/badges/flavor/release/fixest" alt="CRAN status"></a>
<a href="https://CRAN.R-project.org/package=fixest"><img src="http://www.r-pkg.org/badges/version/fixest" alt="Version"> </a>
<a href="https://ipub.com/dev-corner/apps/r-package-downloads/"> <img src="https://cranlogs.r-pkg.org/badges/fixest" alt = "Downloads"> </a>

The `fixest` package offers a family of functions to perform estimations with multiple fixed-effects in both an OLS and a GLM context. Please refer to the [introduction](https://CRAN.R-project.org/package=fixest/vignettes/fixest_walkthrough.html) for a walk-through.

At the time of writing of this page (February 2020), `fixest` is the fastest existing method to perform fixed-effects estimations, often by orders of magnitude. See below for a benchmarking with the fastest alternative software. 

## Benchmarking

Here is a comparison of the performance of `fixest` functions to other state of the art methods to perform estimations with multiple fixed-effects. The results are reported in the five figures below. Package `fixest` (black lines) is consistently faster in all situations.

![](https://github.com/lrberge/fixest/blob/master/vignettes/images/benchmark_gaussian.png?raw=true)

![](https://github.com/lrberge/fixest/blob/master/vignettes/images/benchmark_difficult.png?raw=true)

![](https://github.com/lrberge/fixest/blob/master/vignettes/images/benchmark_poisson.png?raw=true)

![](https://github.com/lrberge/fixest/blob/master/vignettes/images/benchmark_negbin.png?raw=true)

![](https://github.com/lrberge/fixest/blob/master/vignettes/images/benchmark_logit.png?raw=true)

### Setup

The benchmarking was performed as follows: In the OLS context, we estimate the following equation:

<!-- $$y_{ijk} = \alpha_i + \beta_j + \gamma_k + \delta x_{ijk} + \epsilon_{ijk}$$ -->
![](https://github.com/lrberge/fixest/blob/master/vignettes/images/equation.PNG?raw=true)
 
The same functional form (one variable, three fixed-effects) is estimated for the Poisson, the Negative Binomial and the Logit cases (with ad hoc modifications to fit each model). See Berge (2018) for more details on the setup.

For the "difficult" benchmark (OLS only), the data is generated in a way that makes the convergence of the fixed-effects slow. The phenomenon of slow convergence is frequent for real micro-level data sets involving employee and firm fixed-effects for instance.

Each estimation is replicated 10 times and the average computing time is reported in the figures. 

The alternative methods used for comparison are:

* OLS: felm (R: package [lfe](https://github.com/sgaure/lfe)), [reghdfe](https://github.com/sergiocorreia/reghdfe) (Stata) and [FixedEffectModels](https://github.com/FixedEffects/FixedEffectModels.jl) (Julia)
* Poisson: glmmboot (R: package [glmmML](https://cran.r-project.org/package=glmmML)), feglm (R: package [alpaca](https://github.com/amrei-stammann/alpaca)) and [ppmlhdfe](https://github.com/sergiocorreia/ppmlhdfe) (Stata)
* Negative Binomial: glm.nb (R: package MASS) and nbreg (Stata)
* Logit: glmmboot (R: package [glmmML](https://cran.r-project.org/package=glmmML)), feglm (R: package [alpaca](https://github.com/amrei-stammann/alpaca)) and logit (Stata)

All the aforementioned packages were updated at the benchmarking date: February 2020. 

The code and data for the benchmarking can be found [in this folder](https://drive.google.com/drive/folders/1-1M_vLGduByk5P3qHl7AMBBl85RzHpv7?usp=sharing).

## Acknowledgements

Of course the development of `fixest` has been inspired and pushed forward by (almost all) these (great) packages used in the benchmarking and I am deeply indebted to their authors. Although `fixest` contains many features, some are still uncovered and you should definitely have a look at these packages. 



