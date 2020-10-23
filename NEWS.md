
# News for the R Package `fixest`

## Changes in version 0.7.0

#### Bugs

 - Major bug when fixed-effects were combined with `^` and they contained NAs (thanks to @poliquin [#35](https://github.com/lrberge/fixest/issues/35)).
 
 - Bug when using lead/lags in estimations. The bug was due to a bug in a dependency ([dreamerr](https://cran.r-project.org/package=dreamerr)) and was fixed. Now fixest requires **dreamerr** version >= 1.2.1. Bug spotted by @seunghoon001 ([#44](https://github.com/lrberge/fixest/issues/44)).
 
 - Major bug when n_obs x n_vars > 2B or n_obs x n_fixed-effects > 2B. In such cases estimations could just not be done, even leading R to crash when using nthreads > 1. The algorithm was fixed to allow datasets with up to 2B observations to be estimated in all circumstances. Bug reported, and many help for checking provided, by Howard Zihao Zhang.
 
 - `coefplot`: Problem regarding interactions when observations, and hence coefficients, were removed from the estimation. Now the coefficients are removed from the plot. Bug reported by @phisherblack [#45](https://github.com/lrberge/fixest/issues/45).
 
 - `coefplot`: Corrected various bugs when asked for the plotting of several estimations. 
 
 - Fix the stack imbalance warning (report by @shoonlee, [#46](https://github.com/lrberge/fixest/issues/46)).
 
#### Internal improvements

 - Brand new internal algorithm which now uses closed form solutions when dealing with variables with varying slopes. This means that when variables with varying slopes are present, the algorithm is incomparably faster and more accurate.

 - Two deep copies of some data are now avoided in the demeaning function. This improves the performance in terms of memory footprint, and also makes the algorithm faster. 
 
#### Standard-errors, important changes
  
 - New default values for standard-errors (only concerns multiway clustering). They become similar to `reghdfe` to increase cross-software comparability. Computing the standard-errors the old way is still possible using the argument `dof`. See the dedicated vignette: [On standard errors](https://cran.r-project.org/package=fixest/vignettes/standard_errors.html).
 
 - Name change in `summary`/`vcov`/`etable`: To get heteroskedasticity-robust standard-errors, `se = "hetero"` now replaces `se = "white"` to enhance clarity. Note that `se = "white"` still works.
 
#### New function: `fitstat` 

  - New function `fitsat` that computes various fit statistics. It is integrated with `etable` and can be invoked with the argument `fitstat`. So far only two fit statistics are included, but more will come.
 
#### New features in `interact()`

  - You can now use `i(var)` to treat the variable `var` as a factor. You can select which values to drop/keep with the respective arguments. 
  
  - Using `i(var)` leads to a special treatment of these variables in the functions `coefplot` and `etable`.
 
#### New features in `etable`

  - New argument `placement` to define the position of the float in Latex (suggestion by Caleb Kwon).
  
  - New argument `drop.section`, with which you can drop a) the fixed-effects, b) the variables with varying slopes, or c) the statistics, sections (suggestion by Caleb Kwon).
  
  - Fix glitch in help pages regarding the use of the '%' (percentage) character in regular expressions.
  
  - Two new arguments `.vcov` and `.vcov_args` to compute the standard-errors with custom functions.
  
  - The number of observations (`n`) is now treated as a regular statistic and can be placed where one wants.
  
  - The statistics can now have custom aliases using the argument `dict`.
  
  - The overdispersion becomes a regular fit statistic that can be included (or not) using `fitstat`.
  
  - The dictionnary now applies to the factors of interactions, and the values of factors.

#### User visible changes

 - Argument `nthreads`:
 
  * The new default of argument `nthreads` is 50% of all available threads. 
  * Accepts new values: a) 0 means all available threads, b) a number strictly between 0 and 1 will represent the fraction of all threads to use.
  
  - When setting formula macros:
  
    * The functions `xpd` and `setFixest_fml` now accept character vectors and numeric scalars on top of formulas.
  
  - `demean`:
  
    * speed improvement.
  
  - `coefplot`:
  
    * The argument `group` now accepts a special character `"^^"`, when used, it cleans the beginning of the coefficient name. Very useful for, e.g., factors although factors created with `i()` need not that.
    
     * When `horiz = TRUE`, the order of the coefficients is not reversed any more.

 - Improved display of numbers in `print` method.
 
 - Added variables names to `X_demeaned` from `feols`.
 
 - Lagging functions:
 
  * Now `time.step = NULL` by default, which means that the choice of how to lag is automatically set. This means that the default behavior for time variables equal to Dates or character values should be appropriate.
  
  * New operator `d` which is the difference operator.
 
 - In all estimations: 
 
    * new argument `mem.clean`: internally, intermediary objects are removed as much as possible and `gc()` is called before each memory intensive C++ section. Only useful when you're at the edge of reaching the memory limit.
    * new output: `collin.min_norm`, this value informs on the possible presence of collinearity in the system of variables.
    * new arguments `only.env` and `env`:
      * The first, `only.env`, allows to recover only the environment used to perform the estimation (i.e. all the preprocessing done before the estimation).
      * The second, `env`, accepts a fixest environment created by `only.env`, and performs the estimation using this environment--all other arguments are ignored. 
      * These changes are a prerequisite to the efficient implementation of bootstraping (since, by applying modifications directly in `env`, we cut all preprocessing).
    

 - In non-linear estimations: 
 
    * non-numeric variables can now be used. 
    * argument `NL.start` now accepts numeric scalars, initializing all coefficients to the same value (avoids the use of the other argument `NL.start.init`).
    
  - `summary.fixest`:
  
    * argument `.vcov` now accepts functions that compute the vcov. This ensures convenient compatibility with the `sandwich` package (compatibility is still not full though: bootstraped SEs don't work yet).
    
  - `update.fixest`:
  
    * new argument `evaluate` to ensure consistency with the `update` method from stats.
    
  - `feols` & `feglm`:
  
    * the Cholesky decomposition now checks for user interrupts (matters for models with MANY variables to estimate).
  
#### Deprecation

 - Argument `na_inf.rm` has been removed. It was present for historical reasons, and is removed to increase code clarity.

## Changes in version 0.6.0 (13-07-2020)

#### Bugs

 - In `vcov`, the degree-of-freedom in the small sample correction correction was fixed to "nested" and couldn't be modified, now corrected. Further, "nested" was not properly accounted for, now corrected.
 
 - In `etable`, `fitsat = FALSE` or `fitsat = NA` led to a bug.
 
 - `r2`: bug when the estimation contained only fixed effects (thanks to Luis Fonseca [#27](https://github.com/lrberge/fixest/issues/27)).
 
 - Now the `BIC` of `feglm` is similar to the one of `glm`.
 
 - Bug in the log-likelihood in the presence of weights, now corrected.
 
 - Bug in `coefplot` when some interacted variables were removed because of collinearity. Now corrected.
 
#### New vignettes

 - On standard-errors: how are the SEs computed in fixest and how to replicate the SEs from other software.
 
 - Exporting estimation tables: how to use efficiently `etable`, in particular how to customize the tables.
 
#### Major changes: etable

  - New arguments: `group`, `extraline`, `notes`, `tablefoot`. 
  
    - `group` allows to eliminate variables (like `drop`) and adds an extra line with TRUE/FALSE if the model contained those variables.
    
    - `extraline` allows to add extra lines with any content.
    
    - `notes` allows to add notes after the table (suggestion by @bgchamps [#25](https://github.com/lrberge/fixest/issues/25)).
    
    - `tablefoot` controls whether the table footer, containing the type of standard-errors and the significance codes, should be displayed.
    
  - Renaming: `yesNoFixef` => `yesNo`.
  
  - Most default values can be set globally with the new function `setFixest_etable`.
 
#### Major changes: dof

 - Function `dof`, used to adjust the small sample corrections, is now much more complete and allows to replicate a large set of estimation results from alternative software.

#### User visible changes

 - You can now provide custom VCOVs to summary by using the argument `.vcov`.
 
 - A warning is now prompted when the maximum number of iterations of the algorithm is reached (suggestion by @clukewatson  [#24](https://github.com/lrberge/fixest/issues/24)]). 
 
 - The types of standard-errors can now be set globally with the function `setFixest_se` (suggestion by @dlindzee [#28](https://github.com/lrberge/fixest/issues/28))
 
 - New `feols` argument `demeaned`. If `TRUE`, then the centered variables are returned (`y_demeaned` and `X_demeaned`). (Suggestion by Linus Holtermann.)
 
 - `interact` gains two new arguments: `drop` and `keep` (suggestion by @SuperMayo [#23](https://github.com/lrberge/fixest/issues/23)).
 
#### New methods

 - `hatvalues` has been implemented for feols and feglm estimations.
 
 - the `estfun` from `sandwich` has been implemented.

## Changes in version 0.5.1 (18-06-2020)

#### Hotfix
 
 - Fixed bug introduced in the previous update (memory access error). Does not affect any of the results but could lead R to crash unexpectedly (odds were low though since access was adjacent).
 
#### Bugs

 - Fix image link to the equation in the README.md.
 - Fix bug R2 and logLik when observations were removed because of NA values. Due to the update in `residuals.fixest`.
 
#### User visible change

 - Rewriting of the internal algorithm computing the VCOV. 1) About 30% performance gain for estimations with many variables. 2) The code is much less memory hungry. 

#### Major update of `etable`

 - New argument `style` which allows to set many elements of the output table. 
 - (minor) `signifCode` can be equal to `"letters"` to display letters instead of stars.
 
#### Other

 - `setFixest_nthreads` now respects the `OMP_THREAD_LIMIT` environment variable.
 - Rd links are now made to the proper htlm files.


## Changes in version 0.5.0 (10-06-2020)

#### Bug fixes
        
 - Bug with estimations with varying slopes if the fixed-effect relative to the slope is not in its decreasing order (thanks to Davide Proserpio).
 - Bug when interacting two variables with the `var::fe` syntax with `confirm = TRUE` and no reference.
 - Bug in `etable` when the standard-errors where `NA`.
 - Fixed very minor bug when computing the SEs (1e-6 difference).
 - Standard-errors in `feglm` for non-poisson, non-binomial families, are now correct (minor differences).
 - `fixef` did not work when the slope was an integer, now corrected (thanks to @clerousset [#20](https://github.com/lrberge/fixest/issues/20)).

#### New functionality: formula macros
        
 - You can use macros in formulas.
 - To set a macro variable, use e.g., `setFixest_fml(..ctrl = ~ var1 + var2)`. Here the macro variable `..ctrl` has been set to the value `"var1 + var2"`.
 - Now you can use this macro variable in any `fixest` estimation: e.g. `data(airquality) ; setFixest_fml(..ctrl = ~ Temp + Day) ; feols(Ozone ~ Wind + ..ctrl, airquality)`.
 - You can use macros in non-fixest estimations with `xpd`, which expands formulas. E.g. `lm(xpd(Ozone ~ Wind + ..ctrl), airquality)`.
 
#### New functions

 - `to_integer`: user-level version of the internal algorithm transforming any kind of vector (or combination of vectors) into an integer ranging from 1 to the number of unique elements of the vector. Very fast.
 - `demean`: user-level version of the demeaning algorithm used in `feols`.

#### Major user-visible changes
        
 - New internal algorithm to estimate OLS (applies to both `feols` and `feglm`):
 
    1. It is numerically more stable.
    
    2. Incomparably faster when factors are to be estimated (and not explicitly used as fixed-effects).
    
    3. Collinear variables are removed on the fly.

#### User-visible changes
        
 - Interactions in `var::fe(ref)` now accept multiple references (i.e. `ref` can be a vector).
 - In `etable`, the variable names of non-Latex output can now be changed.
 - You can use the argument `n` when applying summary to choose the number of coefficients to display.
 - Argument `confirm` has been removed from the function `interact`.
 - `r2` allows more flexibility in the keywords it accepts.
 - Function `dof` gains a new argument `adj` which allows to make different types of common small sample corrections. Its other arguments have been renamed for clarity (`fixef` => `fixef.K`, `exact` => `fixef.exact`, `cluster` => `cluster.adj`).
 - Now t-statistics are used for `feols` and non-poisson, non-binomial models in `feglm`. For all other models, z-statistics are used. This complies with the default's R-stats behavior.

#### New Methods
        
 - The `residuals` method has been substantially improved, now allowing different types.
 - New stats methods: sigma, deviance, weights.

#### Vignette and Readme
        
 - Typos corrected.
 - Images in the Readme set to 1200px.
 
#### Issue found: convergence problems with multiples variables with varying slopes
        
 - Convergence problems may arise in the presence of **multiple** variables with varying slopes. Theoretical work helped find a solution to this problem, but the implementation in R is proving not instantaneous.
 - In the meantime, now a warning is prompted when the algorithm suspects a convergence problem leading to poor precision of the estimated coefficients.

#### Error-handling
        
 - Improved error-handling with [dreamerr](https://cran.r-project.org/package=dreamerr)'s functions.
 
#### Other

 -  Dependency to MASS has been removed.

## Changes in version 0.4.1 (2020-04-13)

#### Bug fixes
  
 - Major bug leading R to crash when using non-linear-in-parameters right-hand-sides in feNmlm. Only occured when some observations were removed from the data set (due to NAness or to perfect fit). [Thanks to @marissachilds, GH issue [#17]](https://github.com/lrberge/fixest/issues/17).]
 - In the `collinearity` help pages: an example could lead to an error (due to random data generation). It has been removed.
 - In `collinearity`, corrected the problem of display of the intercept in some situations.
 - Defaults for the arguments `cex` and `lwd` in `coefplot` have been changed to 1 and 1 (instead of par("cex") and par("lwd")). Otherwise this led to the creation of `Rplots.pdf` in the working directory (thanks to Kurt Hornik).
 - Corrected a typo in the article's title in the vignette.

#### Help
  
 - Rewriting of sections, correction of small mistakes (wrong argument names), dropping completely the 'cluster' terminology (meant for fixed-effects), addition of what is contained in fixest objects.

#### Other
  
 - Adding a README.md.
 - Small corrections in the vignette.

## Changes in version 0.4.0 (2020-03-27)

#### User visible changes: Latex export
        
 - Better Latex special character escaping (errors reported by @dlindzee [#15](https://github.com/lrberge/fixest/issues/15)).
 - New argument `fixef_sizes.simplify`, which provides the sizes of the fixed-effects in parentheses when there is no ambiguity.
 - You can suppress the line with the significance codes with `signifCode = NA`.
 - New argument `float` which decides whether to embed the table into a table environment. By default it is set to `TRUE` if a `title` or `label` is present.
 - New argument `keep` to select the variables to keep in the table.
 - New way to keep/drop/order variables with the special argument "\%". If you use "\%var", then it makes reference to the original variable name, not the aliased one (which is the default).
 - New argument `coefstat` defining what should be shown below the coefficients (standard-errors, t-stats or confidence intervals). Suggestion by @d712 [#16](https://github.com/lrberge/fixest/issues/16).
 - Better rendering of significant digits.

#### User visible changes: coefplot
        
 - Argument `horiz`. The coefficients can now be displayed horizontally instead of vertically.
 - The coefficient labels, when in the x-axis, can now be displayed in three different ways thanks to the new argument `lab.fit`: "simple", the classic axis, "multi", the labels appear across multiple lines to avoid collision, and "tilted" for tilted labels.
 - The margins now automatically fit.
 - Argument `style` allows you to set styles with the function `setFixest_coefplot`, you can then summon the style in `coefplot` with this argument.
 - Use the ampersand to set dictionnary variables specific to `coefplot` (e.g. setFixest_dict(c(var1 = "a $\\times$ b", "&var1" = "&a\%*\%b")), the variable `var1` will be rendered differently in `etable` and in `coefplot`).
 - Better display of groups (with the arguments `group` and `group.par`).

#### New methods
        
 - `terms.fixest` giving the terms of the estimation.

#### Other
        
 - All `donttest` sections were removed from help pages.


## Changes in version 0.3.1 (2020-02-09)

#### Major bug fix
        
 - [panel] Fixed faulty memory access when taking the lead of a variable.

#### Other bug fixes
        
 - [esttable/esttex] These two functions were replaced by the function `etable`. In the process, some of their arguments were "lost", this is now corrected.
 - [etable] Better escaping of special characters.
 - [estimations] Bug when particular non-numeric vectors were used in explanatory variables.

#### New features
        
 - [coefplot] The function `coefplot` now accepts lists of estimations.


## Changes in version 0.3.0 (2020-02-01)

#### New feature: Lagging
         
 - You can now add lags and leads in any `fixest` estimations. You only need to provide the panel identifiers with the new argument `panel.id`, then you're free to use the new functions `l()` for lags and `f()` for leads.
 
 - You can also set up a panel data set using the function `panel` which allows you to use the lagging functions without having to provide the argument `panel.id`, and which dispose of more options for setting the panel.

#### New feature: Interactions
         
 - You can now add interactions in formulas with a new syntax: `var::fe(ref)`
 
 - The command `var::fe(ref)` interacts the variable `var` with each value of `fe` and sets `ref` as a reference. Note that if you don't use the argument `ref`, the command `var::fe` is identical to `var:factor(fe)`.
 
 - Using `var::fe(ref)` to write interactions opens up a special treatment of such variables in the exporting function `etable` and in the coefficient plotting function `coefplot`.

#### New feature: `coefplot`
         
 - You can plot coefficients and their associated confidence intervals with the function `coefplot`.
 
 - `coefplot` dispose of many options, whose default values can be set with the function `setFixest_coefplot`.
 
 - As for the function `etable`, you can easily rename/drop/order the coefficients.
 
 - `coefplot` detects when interactions have been used and offers a special display for it.

#### New functions
         
 - [etable] Estimations table: new function to export the results of multiple estimations. Replaces the two functions `esttex` and `esttable` (the two functions still exist but they will be deprecated in the future).
 
 - [Lagging] New functions related to lagging: `l`, `f`, `panel`, `unpanel` and `[.fixest_panel`.
 
 - [Utilities] A set of small utility functions has been added. They allow to extract part a coefficient table or parts of it (like the t-statistics of the standard-error) from an estimation. These functions are `coeftable`, `ctable` (an alias to `coeftable`), `se`, `tstat` and `pvalue`.
 
 - [coefplot] The functions `coefplot` and `setFixest_coefplot`.
 
 - [dof] New function to set the type of degree of freedom adjustment when computing the variance-covariance matrix. You can permanently set the type of DoF adjustment with the new function setFixest_dof().

#### User visible changes
         
 - [all estimations] A key pre-processing step has been paralellized => algorithm faster in general and much faster for multi-FEs.
 - [predict & fitted] Predict and fitted now returns vectors of the length equal to the one of original data.
 - [standard-errors] New ways to compute the standard-errors have been implemented. In particular, now it account for the "nestedness" of the fixed-effects in the clusters by default. You can freely change how to compute the degrees of freedom correction with the function dof().
 - [r2] Computation of the within-R2 for feglm models is now self-contained.
 - [all estimations] New, more accurate, stopping criterion for 2+ fixed-effects.
 - [feols] Estimations are slightly faster.
 - [etable/esttex] When there are interactions, R may change the order of the interactions, making two interactions in two different estimations look different while they are in fact the same (e.g. x3:x2 and x2:x3). Now esstable automatically reorders the interactions when needed for comparison across estimations.
 - [etable/esttable] The type of standard errors is now always shown.
 - [etable/esttex] The aliases provided by 'dict' are also applied within interactions. For example: `dict=c(x1="Wind", x2="Rain")`, with an estimation with the following variables 'x1', 'x2', 'x1:x2' will lead to the following aliases in Latex 'Wind', 'Rain' and 'Wind $\\times$ Rain'.
 - [etable/esttex] Interactions of similar values but of different order (e.g. x1:x2 and x2:x1) are reorderd to appear in the same lines.
 - [etable/esttex] The i) type of standard errors and ii) the significance codes, are now displayed in two separate lines (otherwise the line would be too wide).
 - [etable/esttex] Argument `yesNoFixef` can be of length one, defaulting the second element to the empty string.
 - [etable/esttex] Escaping of Latex special characters is now much more robust.

#### Bug correction
         
 - [all estimations] Fixed: bug when functions in the formula returned matrices.
 - [update] Fixed: error message when the data is missing.
 - [feglm] Fixed: bug double estimation when family not equal to poisson or logit
 - [feglm] Fixed: severe bug occurring for families not equal to poisson or logit
 - [predict] Fixed: bug when the estimation contained combined FEs.
 - [summary] Regarding small sample only: now Student t distribution is used instead of the Normal to compute the pvalue.
 - [esttex] Different variables with the same aliases (given by the argument 'dict') now appear in the same row.
 - [esttex] Arguments 'drop' and 'order' are now applied post aliasing (alias given by the argument 'dict').
 - [esttex] But when exporting multi-way standard errors.
 - [r2] Small bug regarding objects obtained with `did_estimate_yearly_effects`.
 - [estimations] bug when using weights in `feglm`.

## Changes in version 0.2.1 (2019-11-22)

#### Major bug correction
        
 - lag.formula: Bug introduced from previous update which could lead to wrong results. Now fixed.

#### Major user visible changes
         
 - [All estimation methods] Significant speed improvement when the fixed-effects variables (i.e. the identifiers) are string vectors.


## Changes in version 0.2.0 (2019-11-19)

#### New function
         
 -[did_means] New function `did_means` to conveniently compare means of groups of observations (both treat/control and pre/post). Contains tools to easily export in Latex.



#### Major user visible changes
        
 - [All estimation methods] Significant speed improvement when the fixed-effects variables (i.e. the identifiers) are of type integer or double.
 - [esttex, esttable] New argument 'fitstat' to select which fit statistic to display. The default adapts to the models. Old arguments (loglik, bic, aic, sq.cor) are dropped.
 - [esttable] Significantly better rendering of SE types.
 - [r2] Now NA is returned for R2s that have no theoretical justification (e.g. within R2 when no FEs, or 'regular' R2 for ML models).



#### Minor user visible changes
        
 - [did_plot_yearly_effects] Now the name of the dependent variable appears on the y-axis.
 - [esttex] Usage of the `sym` macro in Latex is dropped.



#### Bug correction
        
 - [fixef.fixest] bug could appear when using varying slopes coefficients in some specific circumstances (when the slope FEs were different from the regular FEs).
 - [fixef.fixest] bug when many regular FEs jointly with varying slopes.
 - [fixef.fixest] regarding slope coefficients: now the algorithm also evaluates functions of variables.
 - [esttable] Width of the "separating lines" now appropriately set for long dependent variable names.
 - [esttex] Spelling mistake corrected.
 - [estimations] Bug could occur for extremely small data sets (< 10 observations).



#### Error handling
        
 - [esttex, esttable] More informative error messages in functions esttex and esttable.


## Changes in version 0.1.2 (2019-10-04)
#### Major bug correction
        
 - lag.formula: When the data was not in a particular format, the results could be wrong. Now corrected.


## Changes in version 0.1.1 (2019-09-20)
#### Major bug correction
	    
 - feglm: bug when a) the deviance at initialization was higher than the deviance of the first iteration of the IRWLS and b) the step-halving was unable to find a lower deviance. This led the estimation to fail with an error although it should have been performed properly.
 - did_estimate_yearly_effects: bug when the estimation involved periods with negative values.

#### Minor bug correction
        
 - esttex: bug regarding the number of digits of negative coefficients to be displayed
 - esttex: now properly escaping the percentage and the underscore for exports in Latex
 - esttex: bug when changing the names of the dependent variables using a dictionnary
 - vcov: some warning messages were misleading
 - update: bug update when using the argument nframes
 - update: bug when updating the function fepois

#### Error handling
        
 - Better error messages for: did_estimate_yearly_effects, main estimation functions, setFixest_dict, fepois and fenegbin.

## Version 0.1.0 (2019-09-03)

#### First version
	    
 - This package is an effort to create a family of fast and user-friendly functions to perform estimations with multiple fixed-effects (F.E.).

 - Estimations with fixed-effects (or call it factor variables) is a staple in social science. Hence having a package gathering many methods with fast execution time is of prime importance. At the time of this version, this is the fastest existing method to perform F.E. estimations (often by orders of magnitude, compared to the most efficient alternative methods [both in R and Stata]). The underlying method to obtain the F.E. is based on Berge 2018, and the workhorse of the code is in c++ parallelized via OpenMP (btw thanks Rcpp for simplifying coders' life!).

 - This package is the follow up of the (now deprecated) package `FENmlm` which performed fixed-effects estimations but for only four likelihood families. Package `fixest` completely supersedes `FENmlm` by extending the method to regular OLS and all GLM families, and adding new utility functions. Further, the design of the functions has been completely overhauled and extended towards much more user-friendliness. Massive effort has been put into providing a set of informative error messages to the user for quick debugging of her workflow (e.g. one of the functions contains over 100 different errors messages).


