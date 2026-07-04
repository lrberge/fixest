
# fixest: fast and user-friendly econometric estimations

<a href="https://cran.r-project.org/web/checks/check_results_fixest.html"><img src="https://badges.cranchecks.info/worst/fixest.svg" alt="CRAN status"></a>
<a href="https://fastverse.r-universe.dev"><img src="https://fastverse.r-universe.dev/badges/fixest" alt="Version_ROpensci"></a>
<a href="https://CRAN.R-project.org/package=fixest"><img src="https://www.r-pkg.org/badges/version/fixest" alt="Version"> </a>
<a href="https://ipub.com/dev-corner/apps/r-package-downloads/"> <img src="https://cranlogs.r-pkg.org/badges/fixest" alt = "Downloads"> </a>


`fixest` is an R package offering fast and flexible econometric estimations. It provides a unified
framework for applied research, with comprehensive support for a diverse class of models:
ordinary least squares (OLS), instrumental variables (IV), generalized linear models (GLM), maximum
likelihood (ML), and difference-in-differences (DiD).

## Features

A few features:
- super fast
- large array of models (OLS, IV, GLM, ML) with a consistent syntax throughout
- highly optimised handling of fixed-effects
- 6 built-in family of VCOVs for inference
- run multiple estimations with minimal typing, at dazzling speed
- integrated panel features (lead/lag/diff)
- top notch error handling to facilitate the user's life
- programmatic formula manipulation with built-in interpolation and macro variables
- compatible with > 20 `base`/`stats` methods augmented with numerous options (`predict`, etc)
- methods to extract parts of the, possibly multiple, results (`coeftable`, `se`, `pvalue`, `tstat`)
- easy reporting of multiple results on the console, or in publication-ready Latex tables
- coefficient plots, including event-study graphs

If you are new to `fixest`, you may be interested in the `fixest` article (https://arxiv.org/abs/2601.21749) or the [introductory walk-through](https://CRAN.R-project.org/package=fixest/vignettes/fixest_walkthrough.html). 

If you are coming from Stata, the stata2R website is highly relevant: https://stata2r.github.io/fixest/.

For more details, there are four dedicated guides: to the [standard-errors](https://CRAN.R-project.org/package=fixest/vignettes/standard_errors.html), to [collinearity](https://CRAN.R-project.org/package=fixest/vignettes/collinearity.html), to [multiple estimations](https://CRAN.R-project.org/package=fixest/vignettes/multiple_estimations.html), to the [exportation of tables](https://CRAN.R-project.org/package=fixest/vignettes/exporting_tables.html). Quickly find how to specify the formula/VCOV with this [cheat sheet](https://lrberge.github.io/fixest/articles/cheat-sheet.html).

## Installation

```R
# To install from CRAN:
install.packages("fixest")

# To install the latest stable development release:
install.packages("fixest", repos = ropensci = 'https://fastverse.r-universe.dev')
```

## Quickstart

In this example, we use the built-in `airquality` data set to estimate the effect of solar radiations and temperature on ozone concentration.

```R
library(fixest)
data(airquality)

# OLS model: feols
feols(Ozone ~ Solar.R + Temp, airquality)
#> NOTE: 42 observations removed because of NA values (LHS: 37, RHS: 7).
#> OLS estimation, Dep. Var.: Ozone
#> Observations: 111
#> Standard-errors: IID 
#>                Estimate Std. Error  t value   Pr(>|t|)    
#> (Intercept) -145.703155  18.446718 -7.89860 2.5293e-12 ***
#> Solar.R        0.057110   0.025719  2.22053 2.8471e-02 *  
#> Temp           2.278467   0.245996  9.26222 2.2156e-15 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> RMSE: 23.2   Adj. R2: 0.501249
```

We add as controls the month and the day. We treat these variables as fixed-effects, that is categorical variables for which we do not report the estimated coefficients.

```R
feols(Ozone ~ Solar.R + Temp | Month + Day, airquality)
#> NOTES: 42 observations removed because of NA values (LHS: 37, RHS: 7).
#>        0/2 fixed-effect singletons were removed (2 observations).
#> OLS estimation, Dep. Var.: Ozone
#> Observations: 109
#> Fixed-effects: Month: 5,  Day: 29
#> Standard-errors: IID 
#>         Estimate Std. Error t value   Pr(>|t|)    
#> Solar.R 0.040986   0.027355 1.49828 1.3831e-01    
#> Temp    2.782978   0.381101 7.30246 2.6723e-10 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> RMSE: 16.5     Adj. R2: 0.643603
#>              Within R2: 0.48899 
```

We now use the multiple estimations syntax to run the two previous models at once and display the results with `fixest`'s `etable`:

```R
# sw0(): 'sw' means "stepwise", '0' means starting with the empty element
est_multi = feols(Ozone ~ Solar.R + Temp | sw0(Month + Day), airquality)
etable(est_multi)
#>                               x.1               x.2
#> Dependent Var.:             Ozone             Ozone
#>                                                    
#> Constant        -145.7*** (18.45)                  
#> Solar.R          0.0571* (0.0257)   0.0410 (0.0274)
#> Temp            2.278*** (0.2460) 2.783*** (0.3811)
#> Fixed-Effects:  ----------------- -----------------
#> Month                          No               Yes
#> Day                            No               Yes
#> _______________ _________________ _________________
#> S.E. type                     IID               IID
#> Observations                  111               109
#> R2                        0.51032           0.75580
#> Within R2                      --           0.48899
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Let's report the same results but with clustered standard-errors:

```R
# we use the `vocv` argument to change the standard-errors
etable(est_multi, vcov = vcov_cluster("Month"))
#>                      est_multi.1      est_multi.2
#> Dependent Var.:            Ozone            Ozone
#>                                                  
#> Constant         -145.7* (32.70)                 
#> Solar.R          0.0571 (0.0289)  0.0410 (0.0403)
#> Temp            2.278** (0.3428) 2.783** (0.3395)
#> Fixed-Effects:  ---------------- ----------------
#> Month                         No              Yes
#> Day                           No              Yes
#> _______________ ________________ ________________
#> S.E.: Clustered        by: Month        by: Month
#> Observations                 111              109
#> R2                       0.51032          0.75580
#> Within R2                     --          0.48899
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Let's split the sample by month and represent the coefficients for `Solar.R` and `Temp`.

```R
est_month = feols(Ozone ~ Solar.R + Temp, airquality, split = ~Month)
# we use the argument `drop` to drop the Constant
coefplot(est_month, drop = "Constant")
# we add a legend representing the months
legend("topleft", col = 1:5, lty = 1, lwd = 2, legend = paste0("Month = ", 5:9))
## Note: You could obtain the legend programmatically with
# stringmagic::sma("{sample.var} = {sample}", .data = models(est_month))
```

![](https://github.com/lrberge/fixest/blob/master/vignettes/images/readme/coefplot.png?raw=true)

## Benchmarks

Benchmarks should never be taken at face value given that the timings depend on many parameters: the hardware, the software versions installed, and the specific numerical setups. However, with these caveats in mind, they can still be informative. Below is a set of benchmarks on simulated data with simple/difficult fixed-effects convergence properties, run on 12 January 2026 with up to date software at that time:

![](https://github.com/lrberge/fixest/blob/master/vignettes/images/readme/bench_ols_simple.png?raw=true)
![](https://github.com/lrberge/fixest/blob/master/vignettes/images/readme/bench_ols_difficult.png?raw=true)
![](https://github.com/lrberge/fixest/blob/master/vignettes/images/readme/bench_poisson_simple.png?raw=true)
![](https://github.com/lrberge/fixest/blob/master/vignettes/images/readme/bench_poisson_difficult.png?raw=true)

Please refer to "Section 8 Benchmarks" of the [fixest paper](https://arxiv.org/abs/2601.21749) for details on the setup. The code to run the benchmarks is [on zenodo](https://zenodo.org/records/20704145). See "Section 8.2 NYC taxi data" for a benchmark highlighting `fixest`-specific features with a large real-world data set.

Across the board, it is safe to say that `fixest` is one of the fastest tools out there to perform econometric estimations with fixed-effects. Julia's [`FixedEffectModels`](https://github.com/FixedEffects/FixedEffectModels.jl) is very close and might be faster on specific contexts. Even if python's [`pyfixest`](https://github.com/py-econometrics/pyfixest) is slower on these benchmarks, this is likely to change in the close future: Recently the `pyfixest` team, led by Alex Fisher, has been very active in the development of new demeaning algorithms and have come up with new, and very creative, solutions. In particular for 3+ fixed-effects with difficult convergence properties, you should [check it out](https://pyfixest.org/how-to/demeaner-backends.html) if interested.

## Acknowledgements

The development of `fixest` has been inspired and pushed forward by many existing software. In particular, the formula syntax with the pipe to separate the fixed-effects from the independent variables is not from me but from Simen Gaure's [`lfe`](https://cran.r-project.org/package=lfe). I learned a lot on R-language manipulation by perusing the source code of Achim Zeileis' [`Formula`](https://cran.r-project.org/package=Formula) package. I was spurred to improve `fixest`'s algorithms by a sane (yet intense!) competition: first by `lfe`, then by Julia's `FixedEffectModels`.

I am truly indebted to the community who, with countless -- and often of very high quality -- bug reports, has tremendously improved the robustness of the software.

Big thanks to the direct contributors (esp. Grant, Kyle, Sebastian) for taking the time to contribute to this project with their expertise.


