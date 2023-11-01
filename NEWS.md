
# fixest 0.11.2

## Bugs

- fix bug "t value" displaying in lieu of "z value". Thanks to issue on [easystats/parameters #892](https://github.com/easystats/parameters/issues/892). Very ancient bug!

- fix bug in coefplot: now the confidence intervals use the Student t when appropriate (before that, only the Normal law was used). Thanks to Grant McDermott, [#409](https://github.com/lrberge/fixest/pull/408)

# fixest 0.11.1

## Documentation bug

- fix bug in the help of `bin`

## C++ code

- replace stdint with cstdint

# fixest 0.11.0

## Bug fixes

 - fix bug in function `coef()` leading to methods to throw errors in R devel. Thanks to @vincentarelbundock for reporting (#291).

 - fix bug in the `predict` method when applied to objects estimated with `feNmlm`. Thanks again to @vincentarelbundock for reporting (#292)! 
 
 - fix missing variable names in the VCOV matrix of `feNmlm` models. Thanks (yet again!) to @vincentarelbundock for reporting (#293). Comme on dit : jamais deux sans trois !
 
 - fix display bug in cluster names in `etable`. 
 
 - fix bug in IV estimations with no exogenous variable and no fixed-effect (thanks to Kyle Butts, #296).
 
 - fix bug panel vs panel.id behaving differently in terms of default type of VCOV when the estimation did not contained lags. 
 
 - fix bug in `confint.fixest` when only one variable was estimated (thanks to @joachim-gassen, #296).
 
 - fix several bugs in predict when using `i()`, in particular when used in combination with a factor or `poly()` (reported by @rfbressan, #301).
 
 - fix bug in `etable` relating to ampersands not being correctly escaped.
 
 - fix bug in `sunab` when the time variable is exactly named `t` (reported by Florian Hollenbach, #330).
 
 - fix bug in `feols.fit` when `vcov` was supplied and the estimation did not contain fixed-effects (reported by @grlju, #341).
 
 - fix bug in `sample_df` when the name of the variable was too long.
 
 - fix rare bug regarding an error message when a missing variable did exist as a function in the environment. 
 
 - fix bug preventing the use of binning with formulas reported by @tlcaputi, #359).
 
 - fix various errors in the documentation (thanks to Ed Rubin and others!).
 
 - fix bug in warning message in peculiar case of divergence in GLM (reported by @pachadotdev, #315).
 
 - fix bug preventing the use of the global data set in wrapper functions (`fepois`, `fenegbin`, etc). Reported by @turbanisch, #343.
 
 - fix bug preventing the use of `split` in non-GLM, non-OLS estimations (reported by @bberger94, #333).
 
## Multiple estimations

 - new internal algorithm leading to an object very much like a plain list, much easier to interact with.
 
 - new function `models` to extract the matrix of reporting which model has been estimated.
```R
base = setNames(iris, c("y", "x1", "x2", "x3", "species"))
mult_est = feols(y ~ csw(x.[,1:3]), base)
models(mult_est)
#>   id          rhs
#> 1  1           x1
#> 2  2      x1 + x2
#> 3  3 x1 + x2 + x3
```
 - in multiple estimations: all warnings are turned to notes and all notes are delayed and stacked.
 
 - `coef.fixest_multi`: Now reports the model of each estimation in the first columns. Also gains the arguments `collin`, `long` (to display the results in a long format) and `na.rm`.
```R
coef(mult_est)
#>   id          rhs (Intercept)         x1       x2         x3
#> 1  1           x1    6.526223 -0.2233611       NA         NA
#> 2  2      x1 + x2    2.249140  0.5955247 0.471920         NA
#> 3  3 x1 + x2 + x3    1.855997  0.6508372 0.709132 -0.5564827

# Now in long format
coef(mult_est, long = TRUE)
#>    id          rhs coefficient   estimate
#> 1   1           x1 (Intercept)  6.5262226
#> 2   1           x1          x1 -0.2233611
#> 5   2      x1 + x2 (Intercept)  2.2491402
#> 6   2      x1 + x2          x1  0.5955247
#> 7   2      x1 + x2          x2  0.4719200
#> 9   3 x1 + x2 + x3 (Intercept)  1.8559975
#> 10  3 x1 + x2 + x3          x1  0.6508372
#> 11  3 x1 + x2 + x3          x2  0.7091320
#> 12  3 x1 + x2 + x3          x3 -0.5564827
```
 
 - new methods: `coeftable.fixest_multi`, `se.fixest_multi`, `tstat.fixest_multi`, `pvalue.fixest_multi` to easily extract the results from multiple estimations.
 
```R
coeftable(mult_est)
#>   id          rhs coefficient   Estimate Std. Error   t value     Pr(>|t|)
#> 1  1           x1 (Intercept)  6.5262226 0.47889634 13.627631 6.469702e-28
#> 2  1           x1          x1 -0.2233611 0.15508093 -1.440287 1.518983e-01
#> 3  2      x1 + x2 (Intercept)  2.2491402 0.24796963  9.070224 7.038510e-16
#> 4  2      x1 + x2          x1  0.5955247 0.06932816  8.589940 1.163254e-14
#> 5  2      x1 + x2          x2  0.4719200 0.01711768 27.569160 5.847914e-60
#> 6  3 x1 + x2 + x3 (Intercept)  1.8559975 0.25077711  7.400984 9.853855e-12
#> 7  3 x1 + x2 + x3          x1  0.6508372 0.06664739  9.765380 1.199846e-17
#> 8  3 x1 + x2 + x3          x2  0.7091320 0.05671929 12.502483 7.656980e-25
#> 9  3 x1 + x2 + x3          x3 -0.5564827 0.12754795 -4.362929 2.412876e-05
```

- new method `confint.fixest_multi` to extract the confidence intervals of multiple estimations.
 
## xpd

 - new argument `add` to facilitate adding elements to the formula.
 
 - new argument `frame` to tell where to fetch the values of the variables expanded with the dot square bracket operator.
 
 - empty strings or empty elements expanded with `.[]` are now set to be equal to `1` (the neutral element in formulas):
```R
x = ""
xpd(y ~ .[x] + .[NULL])
#> y ~ 1 + 1

```
 
 - regex values can be negated: just start with a `!`:
```R
xpd(am ~ ..("!^am"), data = mtcars)
#> am ~ mpg + cyl + disp + hp + drat + wt + qsec + vs + gear + carb
```
 
 - auto-completion of variables names is now enabled with the '..' suffix.
```R
base = setNames(iris, c("y", "x1", "x2", "x3", "species"))
xpd(y ~ x.., data = base)
#> y ~ x1 + x2 + x3
feols(y ~ x.., base)
#> OLS estimation, Dep. Var.: y
#> Observations: 150 
#> Standard-errors: IID 
#>              Estimate Std. Error  t value   Pr(>|t|)    
#> (Intercept)  1.855997   0.250777  7.40098 9.8539e-12 ***
#> x1           0.650837   0.066647  9.76538  < 2.2e-16 ***
#> x2           0.709132   0.056719 12.50248  < 2.2e-16 ***
#> x3          -0.556483   0.127548 -4.36293 2.4129e-05 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> RMSE: 0.310327   Adj. R2: 0.855706
```

 - when using `xpd` in non-fixest functions, the algorithm tries to guess the `data` so that calls to `..("regex")` or auto-completion can be used seamlessly.
```R
lm(xpd(y ~ x..), base)
#> Call:
#> lm(formula = xpd(y ~ x..), data = base)
#> 
#> Coefficients:
#> (Intercept)           x1           x2           x3  
#>      1.8560       0.6508       0.7091      -0.5565 
```

 - the dot-square-bracket operator in `xpd` also expands one-sided formulas:
```R
x_all = ~sepal + petal
xpd(color ~ .[x_all])
#> color ~ sepal + petal
```

## etable

 - in the argument `fitstat`, the formula is now automatically expanded with `xpd`. This means that you can set fit statistics macro which can be summoned from `etable`. Useful to set default fit statistics for: IVs, GLMs, etc.
 ```R
base = setNames(iris, c("y", "x1", "x2", "x3", "species"))
est = feols(y ~ csw(x.[,1:3]), base)

# setting the macro
setFixest_fml(..fit_ols = ~ n + ar2 + my)

# summoning it
etable(est, fitstat = ~..fit_ols)
#>                             est.1              est.2               est.3
#> Dependent Var.:                 y                  y                   y
#>                                                                         
#> Constant        6.526*** (0.4789)  2.249*** (0.2480)   1.856*** (0.2508)
#> x1               -0.2234 (0.1551) 0.5955*** (0.0693)  0.6508*** (0.0667)
#> x2                                0.4719*** (0.0171)  0.7091*** (0.0567)
#> x3                                                   -0.5565*** (0.1275)
#> _______________ _________________ __________________ ___________________
#> S.E. type                     IID                IID                 IID
#> Observations                  150                150                 150
#> Adj. R2                   0.00716            0.83800             0.85571
#> Dep. Var. mean             5.8433             5.8433              5.8433
 ```
 - now there is support for models with no coefficient (only fixed-effects).
 
 - the application of markdown markup is now more robust and can also be escaped with a backslash. The escaping has been ported to c++.
 
## coeftable

 - it gains the argument `list`. If `TRUE`, then the result is returned in a list form. Useful in Rmarkdown documents for quick reference to specific values.
```R
est = feols(mpg ~ cyl + drat + wt, mtcars)
ct = coeftable(est, list = TRUE)
ct$constant$coef
#> Estimate 
#> 39.76766 
ct$wt$se
#> Std. Error 
#>  0.8293065
```
 
## All estimations

 - arguments `split` and `fsplit` gain the `%keep%` and `%drop%` operators which allow to split the sample only on a subset of elements. All estimations also gain the arguments `split.keep` and `split.drop` which do the same thing as the previous operators.
```R
base = setNames(iris, c("y", "x1", "x2", "x3", "species"))
est = feols(y ~ x.[1:3], base, fsplit = ~species %keep% c("set", "vers"))
etable(est)
#>                              model 1            model 2            model 3
#> Sample (species)         Full sample             setosa         versicolor
#> Dependent Var.:                    y                  y                  y
#>                                                                           
#> (Intercept)        1.856*** (0.2508)  2.352*** (0.3929)  1.896*** (0.5071)
#> x1                0.6508*** (0.0667) 0.6548*** (0.0925)   0.3869. (0.2045)
#> x2                0.7091*** (0.0567)    0.2376 (0.2080) 0.9083*** (0.1654)
#> x3               -0.5565*** (0.1275)    0.2521 (0.3469)   -0.6792 (0.4354)
#> ________________ ___________________ __________________ __________________
#> S.E. type                        IID                IID                IID
#> Observations                     150                 50                 50
#> R2                           0.85861            0.57514            0.60503
#> Adj. R2                      0.85571            0.54743            0.57927
```
 
## Dictionary

 - new way to create dictionaries with `as.dict`:
```R
x = "
# Main vars
mpg: Miles per gallon
hp: Horsepower

# Categorical variables
cyl: Number of cylinders; vs: Engine"

as.dict(x)
#>                   mpg                    hp                   cyl                    vs 
#>    "Miles per gallon"          "Horsepower" "Number of cylinders"              "Engine" 

# setFixest_dict works directly with x
setFixest_dict(x)
```
- `setFixest_dict`: i) now the dictionary only grows, ii) you can define variables directly in the arguments of `setFixest_dict`, iii) `as.dict` is applied to the dictionary if relevant, iv) there's a new argument `reset`.
  
## New functions

 - new function `degrees_freedom_iid` which is a more user-friendly version of `degrees_freedom`.
 
 - new function `fdim` to print the dimension of a data set in an user-readable way.

## Other

 - remove warnings when a binomial family is used with weights in `feglm`.
 
 - add the arguments `y`, `X`, `weights`, `endo`, `inst` to the function `est_env` to make it more user-friendly.
 
 - fix documentation typos (thanks to Caleb Kwon).
 
 - `etable` now returns a `data.frame` whose first column is the variables names (before this was contained in the row names).
 
 - fix environment problems when `lean = TRUE`, leading to large objects when saved on disk.
 
 - `print.fixest` now displays the information on the sample/subset/offset/weights.


# fixest 0.10.4

## Hot fix

 - fix major bug related to the extraction of fixed-effects (function `fixef`) when there are 3+ fixed-effects. This bug led to, in some specific circumstances, wrong values for the fixed-effects coefficients. Thanks a lot to @pachadotdev (#286) for finding this out!

## Other bug fixes

 - fix bug in `confint` when `sunab` was used (thanks to Sarah Hofmann).
 
 - fix an important "documentation bug" on the Sun and Abraham method (thanks to Kyle Butts, #287).
 
 - fix bugs regarding `view`/`markdown` features of `etable`.
 
## Other

 - added compatibility with `car::deltaMethod` following Grant McDermott's suggestion.
 
 - new function `lag_fml` which is an alias to `lag.formula`. The latter being easily stomped by other function names from other packages.

# fixest 0.10.3

 - fix bug linked to the proper identification of estimations with only fixed-effects.
 
 - remove the use of `anyNA.data.frame` leading to a dependency to R 3.6.3 (reported by @MichaelChirico, #261).

# fixest 0.10.2

## Bug fixes

  - fix bug in stepwise estimations when two (stepwised) explanatory variables have exactly the same NAs values.

  - fix display bug regarding factors in `etable` when `dict` was present.
  
  - fix bug `tablefoot.value` not working any more (reported by @resulumit, #224).
  
  - fix possible environment problem when estimating non linear functions outside of the global environment.
  
  - fix bug in the stepwise functions `sw` and `csw` when they contained only one variable.
  
  - fix bug in `etable` preventing automatic headers to be displayed.
  
  - fix bug in `n_unik` preventing the auto completion of variable names.
  
  - fix bug in `fitstat` for the KPR statistic (reported by @etiennebacher, #161).
  
  - fix bug in `i` when two factor variables were interacted and one specific value of one variable was to be set as a reference.
  
  - fix bug in `model.matrix` when no variable was used in the estimation (reported by @kylebutts, #229).
  
  - `model.matrix` now returns the variables in the same order as in the estimation -- a discrepancy could happen in stepwise estimations with interactions in which the interactions were put *before* fixed covariates (related to @sergiu-burlacu, #231).
  
  - fix bugs in `feglm.fit` prevented the VCOV to be computed (reported by @etiennebacher and @edrubin, #237)
  
  - fix bug in `predict` with variables created with `i()` leading to a prediction even for values not included in the original estimation (reported by @vincentarelbundock, #235).
  
  - fix bug in multiple estimations when the data contains weights and there are missing values in the y's or X's (reported by @sahilchinoy, #263).
  
  - fix bug multiverse stepwise when the estimation contains fixed-effects or IVs (reported by @resulumit, #260).
  
  - fix bug in the startup message trigger.
  
  - increase the robustness of the code leading to the startup message (reported by @flycattt, #262).
  
  - improve the robustness of the algorithm parsing the fixed-effects (linked to issue, #253).
  
  - fix minor bug in the Cragg-Donald statistic.
  
  - fix peculiar problem on load when directories names end with ".R" (thanks to @kyleam, #271).
  
  - remove remaining large items from GLM estimations with `lean = TRUE`.
  
  - fix bug in removing the singletons from several fixed-effects (reported by @johannesbubeck, #244).
  
  - in rep.fixest: replace argument cluster with argument vcov to enable the use of any VCOV (related to, #258 by @ShunsukeMatsuno).
  
  - fix bug in predict, which automatically discarded NA values (reported by @ColinTB, #273).
  
## etable

#### New arguments

 - new argument `view` to display the latex table in the viewer pane (suggestion by Or Avishay-Rizi, #227). You need to a) have a working distribution of pdflatex, imagemagick and ghostscript, or b) have the R packages pdftools and tinytex installed, for this feature to work.
 
 - new argument `view.cache` in `setFixest_etable`: whether to cache the PNGs generated.
 
 - new argument `export` to export the Latex table in PNG to a file.
 
 - new (experimental) argument `markdown`: Latex tables can be automatically integrated in the non-Latex markdown document in PNG format.
 
 - new argument `div.class`. Linked to the `markdown` argument. In Rmarkdown documents, the table in PNG format is embedded in a `<div>` container. The class of the div is `div.class`, which is by default `"etable"`.
 
 - new argument `tpt` to nest the table in a `threeparttable` environment. Notes are then nested into the `tablenotes` environment.
 
 - in `style.tex`: new argument `notes.tpt.intro` to insert code right after the `tablenotes` environment and before any note (useful to set the font size of notes globally for instance).
 
 - new argument `arraystretch` to set the height of the table rows.
 
 - new argument `fontsize` which applies Latex font sizes to the table.
 
 - new argument `adjustbox`: `adjustbox = TRUE` nests the tabular into an `adjustbox` environment with `width = \\textwidth, center` as default option. Use `adjustbox = x` with `x` a number giving the text-width. Use `adjustbox = "x th"` with `x` a number giving the text-height. Finally you can use a character string, as in `adjustbox = "my options"`, that will be passed verbatim as a an `adjustbox` option.
 
 - new argument `highlight` to highlight the coefficients with a frame or by changing the row/cell color.
 
 - new argument `coef.style` to apply an arbitrary style to one or several coefficients.
 
 - in `style.tex`: new argument `rules_width` to easily set the width of the `booktabs` rules.
 
 - in `style.tex`: new argument `caption.after` to insert code right after the caption.
 
 - in `style.tex`: new argument `no_border` to remove the borders on the sides of the table.

#### New features

 - the quality of the tex output has been substantially improved.
 
 - `signif.code` now replaces `signifCode` (retro compatibility ensured).
 
 - `signifCode` is removed from `setFixest_etable`, and `signif.code` is added to both `style.tex` and `style.df` so that each style can have its own significance code defined globally.
 
 - the object returned by `etable` are now of class `etable_tex` (when `tex = TRUE`) or `etable_df`, both types having their own printing method.
  
 - the significance codes are now displayed under the table when the output is a `data.frame`.
 
 - in `headers`/`extralines`: `cmidrule` does not show up for empty column names any more.
 
 - new markup: markdown-style markup (e.g. `**text**`) can be used to put text in italic/bold in almost anything in the table.
 
 - `notes` can be set in the dictionary: useful for notes (like source for example) that gets repeated across tables.
 
 - `line.top` and `line.bottom` now admit the values `simple` and `double`. The argument `line.bottom` now affects the "effective" end of table, irrespective of the value of `tablefoot`. This is more in line with intuition.
 
 - improve the use of `tabularx`.
 
 - automatic support `makecell`: any new lines found in names within the table will be translated with `makecell`. For example: `"The \n long \n varname"` is automatically translated into `\makecell{The \\ long \\ varname}`.
 
## dsb

 - completely new function `dsb()` to manipulate strings. Applies many low level string operations very easily. The syntax may be a bit disturbing at first, but, unlike French grammar, there's some logic behind!
 
 - there are over 30 basic string operations available! Do complex string manipulations in a single call!

```R
# At first sight, it's impossible to understand what's going on.
# But I assure you, it's pretty logical! 
# Type dsb("--help") to get some help.

dollar = 6
reason = "glory"
dsb("Why do you develop packages? For .[`dollar`*c!$]?",
    "For money? No... for .[U,''s, c?reason]!", sep = "\n")
#> Why do you develop packages? For $$$$$$?
#> For money? No... for G L O R Y!
```

 - the dot square bracket operator in formulas now calls `dsb` when the calls are nested:
```R
xpd(~ sw(.[, "disp:.[/mpg, cyl]"]))
#> ~sw(disp:mpg, disp:cyl)
```

## New argument in all estimations

 - new argument `only.coef` in all estimation. If `TRUE`, then only the estimated coefficients are returned, which can be useful for MC experiments.
 
## New functions

- new function `est_env` to estimate a model from a `fixest` environment. Mostly useful to cut overheads in simulations.
```R
# First we get the environment (the estimation is not performed!)
env = feols(mpg ~ disp + drat, mtcars, only.env = TRUE)

# Then we estimate: we get the reult from feols(mpg ~ disp + drat, mtcars)

est_env(env)
#> OLS estimation, Dep. Var.: mpg
#> Observations: 32 
#> Standard-errors: IID 
#>              Estimate Std. Error  t value   Pr(>|t|)    
#> (Intercept) 21.844880   6.747971  3.23725 3.0167e-03 ** 
#> disp        -0.035694   0.006653 -5.36535 9.1914e-06 ***
#> drat         1.802027   1.542091  1.16856 2.5210e-01    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> RMSE: 3.07661   Adj. R2: 0.712458

# Why doing that? You can modify the env w/t incurring overheads

assign("weights.value", mtcars$wt, env)
# New estimation with weights
est_env(env)
#> OLS estimation, Dep. Var.: mpg
#> Observations: 32 
#> Standard-errors: IID 
#>              Estimate Std. Error  t value   Pr(>|t|)    
#> (Intercept) 21.967576   6.320006  3.47588 1.6241e-03 ** 
#> disp        -0.032922   0.005884 -5.59478 4.8664e-06 ***
#> drat         1.505517   1.470671  1.02369 3.1444e-01    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> RMSE: 5.08781   Adj. R2: 0.709392
```

- new function `ref` which allows to re-factor variables on-the-fly. This function always returns a factor and relocates the values given in the argument as the first factor levels. It also allows to bin values, similarly to the function `bin`:
```R
# We want to place 5 in the first place
ref(1:5, 5)
#> [1] 1 2 3 4 5
#> Levels: 5 1 2 3 4

# You can also bin at the same time
ref(1:5, .("4:5" = 4:5))
#> [1] 1   2   3   4:5 4:5
#> Levels: 4:5 1 2 3
```

 
## Other

 - `bin`: `cut::` now ignores white spaces, so that `cut:: q1 ] q3 [` works appropriately.
 
 - speed of stepwise estimations (using `sw` [not `csw`]) has been improved.
 
 - recursive formula macro definitions are allowed (feature request by @turbanisch, #234).
 
 - the startup message does not pop in Rmarkdown documents any more.
 
 - function `sample_df` gains the argument `previous` which recovers the previous draw.

# fixest 0.10.1

## Bug fixes

 - remove new R native piping test `|>` which led to errors in R < 4.1.0 despite conditional testing.
 
 - fix bug in `etable` `headers` when one wants to include several lines and the first line contains only one element repeated across columns.
 
 - fix bugs in predict: a) when variables are created with functions of the data, and b) when the new data contains single level factors (relates to issues #200 and #180 by @steffengreup and @IsadoraBM).
 
 - fix bug in `etable` non-clustered standard errors not displaying properly in footers.
 
 - fix bug in `etable` regarding the escaping of `fixef_sizes` (reported by Apoorva Lal, #201).
 
 - fix bug introduced in 0.10.0 preventing the estimation of IV models with interacted fixed-effects (reported by @etiennebacher, #203).
 
 - fix bug in IV estimations when: a) no exogenous variables were present AND the IV part contained at lags; and b) the endogenous variables contained at least two lags. Reported by Robbie Minton.
 
 - fix bug in the `.fit` methods when the argument `vcov` wasn't `NULL`.
 
 - fix bug in `summary.fixest_multi`: when the variance was NA and internal bug could pop in some circumstances.
 
 - fix bug `plot.fixef` not working for `fepois` (reported by @statzhero, #213).
 
 - fix error message when the (wrong) argument `X` is used in `feols`.
 
## Dot square bracket operator

 - add a comma first, like in `.[,stuff]`, to separate variables with commas (instead of separating them with additions):
```R
lhs_vars = c("var1", "var2")
xpd(c(.[,lhs_vars]) ~ csw(x.[,1:3]))
#> c(var1, var2) ~ csw(x1, x2, x3)
```

 - new function `dsb`: applies the dot square bracket operator to character strings.
 
 - in the function `dsb`, you can add a string literal in first or last position in `.[]` to "collapse" the character string in question. The way the collapse is performed depends on the position:
```R
name = c("Juliet", "Romeo")

# default behavior => vector
dsb("hello .[name], what's up?")
#> [1] "hello Juliet, what's up?" "hello Romeo, what's up?" 

# string literal in first position
dsb("hello .[' and ', name], what's up?")
#> [1] "hello Juliet and Romeo, what's up?"

# string literal in last position
dsb("hello .[name, ' and '], what's up?")
#> [1] "hello Juliet and hello Romeo, what's up?"

```

## bin

 - `bin`: numeric vectors can be 'cut' with the new special value `'cut::q3]p90]'`, check it out!
```R
data(iris)
plen = iris$Petal.Length

# 3 parts of (roughly) equal size
table(bin(plen, "cut::3"))
#> 
#> [1.0; 1.9] [3.0; 4.9] [5.0; 6.9] 
#>         50         54         46 

# Three custom bins
table(bin(plen, "cut::2]5]"))
#> 
#> [1.0; 1.9] [3.0; 5.0] [5.1; 6.9] 
#>         50         58         42 

# .. same, excluding 5 in the 2nd bin
table(bin(plen, "cut::2]5["))
#> 
#> [1.0; 1.9] [3.0; 4.9] [5.0; 6.9] 
#>         50         54         46 

# Using quartiles
table(bin(plen, "cut::q1]q2]q3]"))
#> 
#> [1.0; 1.6] [1.7; 4.3] [4.4; 5.1] [5.2; 6.9] 
#>         44         31         41         34 

# Using percentiles
table(bin(plen, "cut::p20]p50]p70]p90]"))
#> 
#> [1.0; 1.5] [1.6; 4.3] [4.4; 5.0] [5.1; 5.8] [5.9; 6.9] 
#>         37         38         33         29         13 

# Mixing all
table(bin(plen, "cut::2[q2]p90]"))
#> 
#> [1.0; 1.9] [3.0; 4.3] [4.4; 5.8] [5.9; 6.9] 
#>         50         25         62         13

# Adding custom names
table(bin(plen, c("cut::2[q2]p90]", "<2", "]2; Q2]", NA, ">90%")))
#>         <2    ]2; Q2] [4.4; 5.8]       >90% 
#>         50         25         62         13 
```

 - `bin` also accepts formulas, e.g. `bin = list("<2" = ~ x < 2)` (`x` must be the only variable).
 
 - `bin` accepts the use of `.()` for `list()`.
 
 - you can add the location of the element using `@d` in the name. Useful to rearrange factors:
```R
base = setNames(iris, c("y", "x1", "x2", "x3", "species"))
table(base$species)
#>     setosa versicolor  virginica 
#>         50         50         50

table(bin(base$species, .("@3" = "seto", "@1 VIRGIN" = "virg")))
#>     VIRGIN versicolor     setosa 
#>         50         50         50 
```
 
 
## etable

 - the tex output is now "nicely" formatted.
 
 - argument `extralines` replaces the argument `extraline` to increase coherence. Hence function `extraline_register` becomes `extralines_register` (the change is done without deprecation since I guess this function must be only rarely used).
 
 - arguments `extralines` and `headers` accept `.()` for `list()`.
```R
base = setNames(iris, c("y", "x1", "x2", "x3", "species"))

```
 
## New function
 
 - `check_conv_feols`: checks the convergence of the fixed-effects in `feols` models by looking at the first-order conditions.
 
## New functions, unrelated but possibly useful

Although a bit unrelated to the purpose of this package, these functions are so extensively used in the author's research that he decided to leverage his author privileges to include them in `fixest` to make them easier to share with co-authors.

 - `osize`: simple function returning a formatted object size.
 
 - `n_unik`: simple but flexible function returning the number of unique elements from variables in one or several data sets. Useful for checking keys.
 
 - `sample_df`: simple function to extract random lines from a `data.frame`.
 
## Other new features

 - when computing Newey-West standard-errors for time series, the bandwidth is now selected thanks to the [bwNeweyWest](https://sandwich.r-forge.r-project.org/reference/NeweyWest.html) function from the [sandwich](https://sandwich.r-forge.r-project.org/index.html) package. This function implements the method described in Newey and West 1994.
 
 - add `type = "se_long"` to `summary.fixest_multi` which yields all coefficients and SEs for all estimations in a "long" format.
 
 - only in `fixest` estimations, using a "naked" dot square bracket variable in the left-hand-side includes them as multiple left hand sides. Regular expressions can also be used in the LHS.
```R
base = setNames(iris, c("y", "x1", "x2", "x3", "species"))
y = c("y", "x1")
feols(.[y] ~ x2, base)
#> Standard-errors: IID 
#> Dep. var.: y
#>             Estimate Std. Error t value  Pr(>|t|)    
#> (Intercept) 4.306603   0.078389 54.9389 < 2.2e-16 ***
#> x2          0.408922   0.018891 21.6460 < 2.2e-16 ***
#> ---
#> Dep. var.: x1
#>              Estimate Std. Error  t value   Pr(>|t|)    
#> (Intercept)  3.454874   0.076095 45.40188  < 2.2e-16 ***
#> x2          -0.105785   0.018339 -5.76845 4.5133e-08 ***


etable(feols(..("x") ~ y + i(species), base))
#>                                  model 1            model 2            model 3
#> Dependent Var.:                       x1                 x2                 x3
#>                                                                               
#> (Intercept)            1.677*** (0.2354) -1.702*** (0.2301) -0.4794** (0.1557)
#> y                     0.3499*** (0.0463) 0.6321*** (0.0453) 0.1449*** (0.0306)
#> species = versicolor -0.9834*** (0.0721)  2.210*** (0.0705) 0.9452*** (0.0477)
#> species = virginica   -1.008*** (0.0933)  3.090*** (0.0912)  1.551*** (0.0617)
#> ____________________ ___________________ __________________ __________________
#> S.E. type                            IID                IID                IID
#> Observations                         150                150                150
#> R2                               0.56925            0.97489            0.93833
#> Adj. R2                          0.56040            0.97438            0.93706
```
 
## Other

 - improve error messages when `subset` does not select any element.
 
 - in `xpd` and `fixest` estimations, variables can be "grepped" from the data set with `regex("regex")`.
 
 - add inheritance of the default style in `iplot` when the style is set globally with `setFixest_coefplot`.
 
 - improve error messages in general by prompting additional error calls (when appropriate).
 
 - the dictionaries now ignore white spaces in coefficient names (thanks to Caleb Kwon).
 
 - the package startup messages have been improved (they should pop up less often).
 
 - to comply with CRAN policies, the startup message doesn't write on the .Renviron file any more.

# fixest 0.10.0

## Bugs fixes

 - Fix bug occurring in IV models with multiple instruments and with multithreading on. That bug could lead to the wrong imputation of the IV residuals, hence affecting the standard-errors (although the order of magnitude of the variation should be minor). Thanks to @whitfillp, #182.
 
 - Fix minor, rare, bug occurring in `feglm` when the model was badly specified and VAR(Y) >>>> VAR(X) and there were only one variable.
 
 - `model.matrix` did not work with `type = "fixef"` (thanks to @kylebutts, #172).
 
 - In nonlinear estimations:`fixef.rm = "none"` or `fixef.rm = "singleton"` did not work as expected (thanks to @kre32, #171).
 
 - Fix bug that could occur when observations had to be removed on several fixed-effects dimensions (had no impact on the estimates though).
 
 - Fix bug in `etable` when `file` is provided and `tex = FALSE` (thanks to @roussanoff, #169).
 
 - Fix bug when: i) a `fixest_panel` is used as a data set in an estimation, ii) NA values are to be removed and iii) fixed-effects are used. Thanks to Nicola Cortinovis for the report!
 
 - Fix bug in `to_integer` when converting multiple vectors and sorting is required, without items.
 
 - Fix bug in `feols.fit` when the matrix of regressors was only partially named (reported by @leucothea, #176).
 
 - Fix bug in the value of the fixed-effects coefficients in IV estimations (thanks to @tappek, #190).
 
 - Fix bug in `coefplot` when `lean = TRUE` in the estimation (reported by @adamaltmejd, #195).
 
 - Fix bug in `iplot` when IVs contained interactions.
 
 - Fix bug in `iplot` preventing some variables to be removed (reported by @roussanoff, #164).
 
 - Fix 0 right-padding of numbers displayed in estimation results that could be confusing (reported by @hjuerges, #197).
 
## Major changes

 - New argument `vcov`:
 
    * greatly simplifies and extends how to specify the SEs
 
    * completely replaces the arguments `se` and `cluster` (they still work)
    
    * accepts functions
    
    * see the documentation in the [vignette](https://lrberge.github.io/fixest/articles/fixest_walkthrough.html#the-vcov-argument-1)
    
 - New built-in VCOVs:

    * Newey-West (1987) for serially correlated errors
    
    * Driscoll-Kraay (1998) for cross-sectionally and serially correlated errors
    
    * Conley (1999) for spatially correlated errors
    
 - you can summon variables from the environment directly in the formula using the new dot square bracket (DSB) operator. The DSB operator can be used to create many variables at once, and can also be using within regular expressions. One example:
```R
base = setNames(iris, c("y", "x1", "x2", "x3", "species"))
i = 2:3
z = "i(species)"
feols(y ~ x.[i] + .[z], base)
#> OLS estimation, Dep. Var.: y
#> Observations: 150 
#> Standard-errors: IID 
#>                      Estimate Std. Error   t value   Pr(>|t|)    
#> (Intercept)          3.682982   0.107403 34.291343  < 2.2e-16 ***
#> x2                   0.905946   0.074311 12.191282  < 2.2e-16 ***
#> x3                  -0.005995   0.156260 -0.038368 9.6945e-01    
#> species::versicolor -1.598362   0.205706 -7.770113 1.3154e-12 ***
#> species::virginica  -2.112647   0.304024 -6.948940 1.1550e-10 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> RMSE: 0.333482   Adj. R2: 0.832221
```
    
 - in `i` and `sunab` you can now bin the variables on the fly with the new argument `bin`. The new function [`bin`](https://lrberge.github.io/fixest/reference/bin.html) is also available at the user-level.
 
 - Function `dof` has been renamed into `ssc` (stands for small sample correction) to improve clarity. Retro compatibility is partially ensured but the function `dof` will be removed at some point.
    
## Breaking changes

 - Functions `setFixest_dof` and `setFixest_se` have been renamed into `setFixest_ssc` and `setFixest_vcov`. **No retro compatibility ensured.**
 
 - Removal of the `var::factor` operator to interact a continuous variable with a variable treated as a factor.
 
## New features

 - `feglm` now accepts partially matched character shortcuts for families: "poisson", "logit", "probit" are now valid `family` arguments.
 
 - `predict.fixest` accepts the new argument `fixef` which, if `TRUE`, returns a data.frame of the fixed-effects coefficients for each observation, with the same number of columns as the number of fixed-effects (feature requests #144 and #175 by @pp2382 and @cseveren).
 
 - offsets present in the formula are now accepted.
 
 - VCOV aliases (both Grant McDermott's suggestions: thanks Grant!): 
 
    * the default standard-errors is now `"iid"` (former keywords still work)
    
    * the keyword `hc1` can be used to summon heteroskedasticity-robust SEs
    
  - *Argument sliding*: the argument `vcov` can be called implicitly when `data` is set up globally:
  ```R
  base = setNames(iris = c("y", "x1", "x2", "x3", "species"))
  
  # Setting up the data
  setFixest_estimation(data = base)
  
  # Now vcov can be used without using vcov = stuff:
  feols(y ~ x1 + x2, ~species)
  
  # => same as feols(y ~ x1 + x2, vcov = ~species)
  ```
  
  - piping the data now works:
  ```R
  mtcars |> feols(cyl ~ mpg)
  # => same as feols(cyl ~ mpg, mtcars)
  ```
  
  - the user can now specify custom degrees of freedom to compute the t-tests in `ssc()` (feature request by Kyle F. Butts).
  
  - the predict method gains the new argument `se.fit` and `interval` which computes the SEs/CI of the predicted variable. This only works for OLS models without fixed-effects. Feature request by Gábor Békés, #193.
  
  - `iplot` gains the argument `i.select` to navigate through the different variables created with `i()` (provided there are more than one of course).
  
## etable

 - new argument `interaction.order` to control the order in which interacted variables are displayed (feature request by @inkrement, #120).
 
 - new argument `i.equal` to control how the values taken by factor variables created with `i()` are displayed.
 
 - new `meta.XX` family of arguments when exporting to Latex. They include various type of information as comments before the table (suggestion of adding the time by Apoorva Lal, #184). So far the new arguments are: `meta.time`, `meta.author`, `meta.sys`, `meta.call`, `meta.comment`. The argument `meta` is a shortcut to all these.
 
 - default values can be saved at the project level using the argument `save = TRUE` in the function `setFixest_etable`. This means that default values will be automatically set without having to call `setFixest_etable` at the startup of the R session. For example, if you want to permanently add the creation time in your Latex exports, just use `setFixest_etable(meta.time = TRUE, save = TRUE)`, and you won't need to bother about it anymore in any future session in your current project.
 
 - some changes in `extraline`, now: *i)* it accepts raw vectors, *ii)* it accepts lines without title, *iii)* the elements are recycled across models, *iv)* it accepts elements of the form `list("item1" = #item1, "item2" = #item2, etc)`, and *v)* the elements are Latex-escaped.
 
 - in `style.tex`: all Latex-escaping is removed.
 
 - the argument `subtitle` has been renamed into `headers` (retro-compatibility is ensured).
 
 - on the new argument `headers`:
   * it accepts named lists of numbers, where the names represent the values in the cell and the numbers represent the span. For example `headers = list("Gender" = list("M" = 3, "F" = 4))` will create a line with 3 times "M" and 4 times "F". 
   * adding the special tag `":_:"` in the row name will add a rule for each column group (in the previous example `":_:Gender"` would do). Suggestion by @nhirschey, #173.
   * you can control the placement of the header line by using as first character the following special tags: "^" (top), "-" (mid, default), "_" (bottom). Ex: `headers = list("_Gender" = list("M" = 3, "F" = 4))` will place the header line at the very bottom of the headers.
   * by default all values are Latex-escaped. You can disable escaping by adding `":tex:"` in the row title.
   
 - argument `sdBelow` has been renamed into `se.below` (retro-compatibility is ensured).
 
 - new argument `se.row` to control whether the row displaying the standard-errors should be displayed (clarification requested by @waynelapierre, #127).
 
 - `dict` now directly modifies the entries in the global dictionary instead of creating a brand new one. For example if you have `setFixest_dict(c(cyl="Cylinder"))` and then use `etable` with `dict=c(mpg="miles per gallon")`, you end up with both the names `cyl` and `mpg` to be modified. To disable this behavior, you can add `"reset"` as the first element, like in `dict=c("reset", mpg="miles per gallon")`.
 
## Other

 - in multiple estimations, it is now made explicit that the information regarding NA values only concern the variables in common across all models (formerly, it was implicit and hence confusing).


# fixest 0.9.0

## Bugs

 - Major bug, leading R to crash, occurring when the same variable was used with several different slopes (thanks to @Oravishayrizi, #119). 
 
 - Major bug, leading R to crash, occurring when 3+ fixed-effects are to be combined.
 
 - Major bug, leading R to crash, occurring when multiple LHS are estimated with the option `fixef.rm = "singleton"` (thanks to Ole Rogeberg).
 
 - Major bug, leading R to crash, occurring when *many* fixed-effects have to be removed because of only 0/1 outcomes (thanks to @mangelett #146 and @ChristianDueben #157).
 
 - Fix bug occurring for undefined covariances with only one regressor (thanks to @joseph-richard-martinez, #118).
 
 - Fix bug in IV estimations regarding the Wald statistic of the first stage when `lean = TRUE` and the VCOV computation is done post estimation.  
 
 - Fix bug in the Wald test in IV estimations when variables are removed because of collinearity (thanks to @pei-huang, #117).
 
 - Fix bug regarding multiple estimations when the multiple fixed-effects contained variables with varying slopes.
 
 - Fix various display bugs in `fitstat`.
 
 - Fix bug in `etable`: using split sample estimations prevented the argument `title` to render correctly. 
 
 - Fix incorrect information message when observations are removed because of infinite values (in some circumstances the removal was wrongly attributed to NAness).
 
 - Fix bug in `etable` when checking the argument `coefstat` (thanks to @waynelapierre, #121).
 
 - Fix bug in `feols` when IV estimations contained fixed-effects and `lean = TRUE` (thanks to @adamaltmejd, #123).
 
 - Fix bug in IV estimations when an endogenous regressor was removed because of collinearity.
 
 - Fix bug estimation without intercept not working when lags are present in the formula (thanks to @nreigl, #126). 
 
 - Fix various bugs when using `subset` in the estimations (reported by @noahmbuckley and @Oravishayrizi, #129 and #131).

 - Fix error message when data cannot be fetched (reported by @Oravishayrizi, #134).
 
 - Fix bug getting the "G" statistic in `fitstat`.
 
 - Fix bug in `predict` when a `poly()` term was used and the formula was long (reported by @XiangLiu-github, #135). 
 
 - fix bug for extracting sub statistics of `"ivwald"` and `"ivf"` in `fitstat`.
 
 - fix bug when `i()` was used without intercept.
 
 - Fix display bug in `etable` when Tex output is requested and interactions composed of identical variables with different interacted orders are present (reported by @Oravishayrizi, #148).
 
 - Fix bug in `etabe` when `fixef.group` is used and fixed-effects are renamed (reported by @jamesfeigenbaum).
 
 - Fix bug when `fplit` is used with subset.
 
 - Fix bug when using `cluster` with `subset` and NA values are removed (reported by @adamaltmejd, #154).
 
 - Fix bug argument `lean` not working on `summary` when applied to an existing `summary` and only the argument `lean` was present (reported by @adamaltmejd).
 
 - Fix bug when using multiple LHS with lags in the formula (reported by @Nicolas Reigl, #158).
 
 - Fix bug regarding the intercept-only likelihood when weights are provided (only with Poisson and logit models), reported by @fostermeijer, #155.
 
## Breaking changes: new i() function

 - the function `i()`, used to create factors or interactions has been tidied up, leading to breaking changes.
 
 - the first two arguments have been swapped! such that now the first argument will always be treated as a factor. 
 
 - the new syntax is `i(factor_var, var, ref, keep, ref2, keep2)` where `var` can be either continuous or factor-like (the argument `f2`, for interaction with factors, has been removed).
 
 - Fix rare bug when the number of parameters is greater than the number of observations and the GLM family has a dispersion parameter.
 
 
## Breaking changes: new default family for feglm

 - to be in line with R stats's `glm`, the new default family for `feglm` is `gaussian` (previously it was Poisson, if you were using it, please now use the function `fepois` instead).
 
## Breaking changes: coefplot is now split in two

 - the function `coefplot` has been split in two: 
 
   - `coefplot`: always plots *all* the coefficients. 
   
   - `iplot`: plots only interactions or factors created with the function `i()`.
  
 - the function `iplot` hence replaces `coefplot`'s former argument `only.inter` which controlled whether or not to focus on interactions.

## etable

  - `group` and `extraline`: enhanced and simplified control of the placement of the new rows. Now only two special characters at the beginning of the row name decide its location.
  
  - new argument `fixest.group`. If `TRUE`, then fixed-effects appearing always jointly across models will be grouped in a single row. The user can alternatively specify a list to declare which fixed-effect to group and customize the row name. 
  
  - `sdBelow` now works when `tex = FALSE` (request by Sasha Indarte).
  
  - `extraline` can now be equal to a formula containing `extraline` macros or valid `fitstat` types.
  
    - When a list, `extraline` can contain functions (returning a scalar) that will be applied to each model.  
  
    - When a list, `extraline` can contain formulas containing `extraline` macros or valid `fitstat` types. 
  
  - You can register `extraline` macros with the new function `extraline_register`.
  
  - When `tex = TRUE`, n-way clustering now always leads to the name of clustered SEs (n-way is not shown any more).
  
  - Add the argument `coef.just` that controls the justification of the coefficients and standard-errors. Only works when `tex = FALSE` (i.e. a `data.frame` is requested).
  
  
## Sun and Abraham staggered DiD method

 - new function `sunab` that simplifies the implementation of the SA method. 
 
 - just type `sunab(cohort, period)` in a `fixest` estimation and it works!
  
## fixest_multi methods

Common methods have been extended to `fixest_multi` objects.

 - `coef.fixest_multi`: re-arranges the coefficients of multiple estimations into a matrix.
 
 - `resid.fixest_multi`: re-arranges the residuals of multiple estimations into a matrix.
 
 
## fitstat: New fit statistics

  - `kpr`: Kleibergen-Paap rank test for IV estimations.
  
  - `cd`: Cragg-Donald F statistic for IV estimations.
  
  - `my`: gives the mean of the dependent variable.
  
  
## New functions

  - `degrees_freedom`: to access the DoFs of the models (sometimes that can be intricate).
  
  - `feols.fit`: fit method for feols. 
  
  - `obs`: to obtain the observations used in the estimation.
  
## New features

  - All `fixest` estimation now accept scalars from the global environment (variables are still not allowed!).
  
  - Better handling of the DoFs in `fitstat` (in particular when the VCOV is clustered).
  
  - `model.matrix`: 
  
    * the endogenous and exogenous regressors, and the instruments from IV estimations can now be easily extracted.
    
    * new arguments `as.matrix` and `as.df` to coerce the result to a particular format.
  
  - `.fit` methods (`feols.fit` and `feglm.fit`) now handle multiple dependent variables.
  
  - `to_integer` now sorts appropriately any kind of vectors (not just numeric/character/factors).
  
  - substantial speed improvement when combining several vectors with *many* cases (> millions).
  
  - The number of threads to use can now be set permanently at the project level with the new argument `save` in the function `setFixest_nthreads`.
  
  
## Minor breaking changes

 - in `model.matrix` when `type = "lhs"` the return object is now a vector (previously it was a data.frame).
  
## Other changes

  - Improve error messages. 
  
  - `hatvalues.fixest`: now returns an error instead of a message when fixed-effects are present (it makes the interplay with `sadnwich` 'nicer').

# fixest 0.8.4

## Bugs

 - Fix bug `depvar = FALSE` not working when tex output was requested (thanks to @apoorvalal and @pbaylis, #104).
 
 - Fix bug in naming when `i()` led to only one variable being retained (thanks to @
colejharvey, #106). 

 - Fix bug display when only degrees of freedom are selected in `fitstat`.

 - Fix bug when `lean = TRUE` in IV estimations with fixed-effects (a large object was still present, thanks to @zozotintin).
 
 - Fix bug display of `etable` in Rmarkdown (thanks to @kdzhang, #93, and @nikolassch, #112)
 
 - Improve error messages in `fitstat` when selecting statistics components.
 
 - Fix bug in `predict` when `poly()` was used in the estimation (thanks to @tholdaway, #109).

 - Fix bug in `predict`: an error message would not pop when combined fixed-effects are used with `combine.quick = TRUE` (thanks to @benzipperer, #115).
 
 - Fix bug to properly account for the nestedness of combined fixed-effects when clustered standard-errors are requested (thanks to @Oravishayrizi , #116).
 
 - Fig major bug in `model.matrix` that could make it very slow. Led the function `aggregate` to be very slow (thanks to Benny Goldman).
 
 - Fix bug that prevented `aggregate` to effectively use weights (thanks to Benny Goldman).
 
## New features

 - `model.matrix` gains the new argument `subset` which allows the creation of the design matrix for a subset of variables only.
 
 - `drop.section` now works for `etable` when `tex = FALSE`.
 
 - The argument `panel.id` used in all estimations can be set globally with the function `setFixest_estimation`.

## Other changes

  - `i`: Factor variables with only the values of 0 and 1 are treated as numeric. 
  
  - `fitstat`: The statistic `G` is now equal to the degrees of freedom used in the t-test of coefficients testing.
  
  - `esttable` and `esttex` are not deprecated any more: they are now pure aliases of `etable`.
  
  - `aggregate`: for weighted regressions, `use_weights` controls whether to use the weights to perform the aggregation.

# fixest 0.8.3 (2021-03-01)

## Bugs

 - Remove test that leads to a (uber odd) bug on fedora devel.
 
 - Fix bug in IV estimation when using factors as instrumented variables (thanks to @adamaltmejd, #99).
 
 - Fix bug when using at least two fixed-effects and varying slopes with singletons (thanks to @adamtheising, #89). 
 
## New features

 - In `xpd`, macros are parsed even after creating the formula with the `lhs` and `rhs` arguments.  

# fixest 0.8.2 (2021-02-11)

## Bugs

  - Fix bug in IV estimations when `lean = TRUE` (thanks to @reifjulian, #88).
  
  - Fix various bugs related to the use of `summary` when `lean = TRUE` in the estimation.
  
  - Fix bug preventing `se = "cluster"` to be used in `etable` (thanks to Caleb Kwon).
  
  - Fix bug `etable` not escaping variable names properly when `sdBelow = FALSE` (thanks to Jeppe Viero).
  
  - Fix bug in IV estimation with `lean = TRUE`.
  
  - Fix bug preventing the return of demeaned variables in IV estimations (thanks to @amarbler, #94).

## Other

 - `i()` now automatically converts its first argument to numeric if it was of type logical. The user can still pass logicals to the argument `f2` if the expected behavior is really to treat it as a logical.
 
 - Improve `fitstat` help and error messages.

# fixest 0.8.1 (2021-01-13)

## Bugs

 - Bug in `etable` when the default value of `fitstat` was set with `setFixest_etable`.
 
 - Bug in `model.matrix` when the model contained fixed-effects and the RHS was requested: the intercept was wrongfully added.
 
 - Fix rare bug when `i()` was called within a very specific set of functions.
 
 - Fix bug in R old release due to `anyNA.data.frame`.
 
 - Fix bug regarding `panel` data sets when variables were created in a `data.table` within functions (thanks to @tcovert, #76).
 
 - Add extra elements to be removed when `lean = TRUE` to keep the object as small as possible (reported by @zozotintin, #81).
 
 - Fix bug in fixed-effects estimations with multiple LHS and different number of observations per estimation that prevented to get the default behavior for standard-errors to work.
 
 - Fix occasional bug when using `split` with fixed-effects.
 
 - `xpd` now appropriately returns a two sided formula when a one sided formula is fed in and the argument `lhs` is provided.
 
 - Fix bug in `coefplot` preventing the proper scaling of the x-axis for interactions when multiple models are displayed.
 
 - Fix occasional bug in the ordering of sub-selections of multiple estimations. 
 
## Sun and Abraham method for staggered DiD

 - For staggered difference-in-difference analyzes: the method of Sun and Abraham (forthcoming, Journal of Econometrics) has been implemented. 
 
 - After having used `i()` to interact cohort dummies with time to treatment dummies, use the function `aggregate` to recover the yearly treatment effects.
 
 - So far the way to do it, although easy, is a bit arcane but the next versions of the software will include a user-friendly way. 
 
 - For details, check out the help page of the function `aggregate` or the staggered DiD section in the vignette [fixest walkthrough](https://CRAN.R-project.org/package=fixest/vignettes/fixest_walkthrough.html).
 
## New features

 - Function `i()` now has the new arguments `f2`, `drop2` and `keep2` which allows the interaction of two factors (useful for staggered DiD estimations).
 
 - Argument `dof`, used to compute the standard-errors, can now be used at estimation time.
 
 - In `etable`, the argument `digits` can now accepts a character value specifying the way the decimals should be displayed. For example if `digits = "r2"` this means that all numbers will be rounded at two decimals and these two decimals will always be displayed. The default behavior is to display significant digits. Follows feature request #82 by @lyifa.
 
 - `etable` also gains the argument `digits.stats` which monitors how the fit statistics decimals should be displayed.
 
 - Argument `split` now accepts variable names.
 

## Other

  - More coherence regarding the use of `summary` applied to models for which the SEs were computed at estimation time. Now there is a memory of how the SEs were computed, so that, for example, if only the argument `dof` is passed to `summary`, then the SEs will be clustered in the same way as estimation time and only `dof` will change.
  
  - Now an error is raised when `i()` is used in the fixed-effects part of the formula. The appropriate way is indicated (related to #77 by 
@rrichmond).
  
  - Improved default setting of standard-errors.
  
  - Improved error messages.
  
  - In multiple estimations, models returning full NA coefficients are not returned (instead of raising an error).

# fixest 0.8.0 (2020-12-14)

## Bugs

 - Major bug when predict was used in the presence of fixed-effects (thanks to @jurojas5, #54). Introduced in version 0.7.

  - When using variable names to cluster the standard-errors inside functions, summary may not fetch the data in the right frame (thanks to @chenwang, #52). Now a completely new internal mechanic is in place.
  
  - When using variables with varying slopes and the number of iterations is greater than 300, a bug occurred in the function checking the convergence was right (thanks to @kendonB, #53). 
  
  - Fix bug in the demeaning algorithm when two variables with varying slopes were identical.
  
  - Fix bug in femlm/feNmlm when factor variables are removed due to the removal of some observations.
  
  - In `summary`, fix bug when the argument `cluster` was equal to a formula with expressions and not a variable name (thanks to @edrubin, #55).
  
  - Fix bug when integers are present in the RHS (thanks to @zozotintin, #56).
  
  - Fix bug when nb_FE >= 2 and the data was large (thanks to @zozotintin, #56). 
  
  - Fix bug display of how the standard-errors were clustered in `etable`.
  
  - Fix bug occurring when lags were used in combination with combined fixed-effects (i.e. fe1 ^ fe2) (thanks to @SuperMayo, #59).
  
  - Fix bug `coefplot` when representing multiple estimations and coefficient names are numbers.
 
## IV

 - IV estimations are now supported. It is summoned by adding a formula defining the endogenous regressors and the instruments after a pipe.
```R
base = iris
names(base) = c("y", "x1", "x_endo", "x_inst", "species")
base$endo_bis = 0.5 * base$y + 0.3 * base$x_inst + rnorm(150)
base$inst_bis = 0.2 * base$x_endo + 0.3 * base$endo_bis + rnorm(150)

# The endo/instrument is defined in a formula past a pipe
res_iv1 = feols(y ~ x1 | x_endo ~ x_inst, base)

# Same with the species fixed-effect
res_iv2 = feols(y ~ x1 | species | x_endo ~ x_inst, base)

# To add multiple endogenous regressors: embed them in c()
res_iv3 = feols(y ~ x1 | c(x_endo, x_endo_bis) ~ x_inst + x_inst_bis, base)

```

## fit statistics

  - The `fitstat` function has been significantly enhanced. 
  
  - Now the following types are supported:
  
    * Likelihood ratios
    
    * F-tests
    
    * Wald tests
    
    * IV related tests (F/Wald/Sargan)
    
    * common stats like the R2s, the RMSE, Log-likelihood, etc
    
  - You can register your own fit statistics. These can then be seamlessly summoned in `etable` via the argument `fitstat`.
  
  - The `print.fixest` function now supports the `fitstat` argument. This means that you can display your own desired fit statistics when printing `fixest` objects. This is especially useful in combination with the `setFixest_print` function that allows to define the default fit statistics to display once and for all. See the example in the "Instrumental variables" section of the Walkthrough vignette.
  
  - The new function `wald` computes basic Wald tests.
 
## Multiple estimations

  - New arguments `split` and `fsplit`: you can now perform split sample estimations (`fsplit` adds the full sample).
    
  - Estimations for multiple left-hand-sides can be done at once by wrapping the variables in `c()`.
    
  - In the right-hand-side and the fixed-effects parts of the formula, stepwise estimations can be performed with the new stepwise functions (`sw`, `sw0`, `csw` and `csw0`).
    
  - The object returned is of class `fixest_multi`. You can easily navigate through the results with its subset methods.

```R
aq = airquality[airquality$Month %in% 5:6, ]
est_split = feols(c(Ozone, Solar.R) ~ sw(poly(Wind, 2), poly(Temp, 2)),
                 aq, split = ~ Month)
                 
# By default: sample is the root
etable(est_split)

# Let's reorder, by considering lhs the root
etable(est_split[lhs = TRUE])

# Selecting only one LHS and RHS
etable(est_split[lhs = "Ozone", rhs = 1])

# Taking the first root (here sample = 5)
etable(est_split[I = 1])

# The first and last estimations
etable(est_split[i = c(1, .N)])
```

## Formula macros

  - The algorithm now accepts regular expressions with the syntax `..("regex")`:
  ```R
  data(longley)
  # All variables containing "GNP" or "ployed" in their names are fetched
  feols(Armed.Forces ~ Population + ..("GNP|ployed"), longley)
  ```
  
    
## New features in etable

  - New `style.tex` and `style.df` arguments that define the look of either Latex tables or the output data.frames. 
  
    * it can be set with the new functions `style.tex` and `style.df` that contain their own documentation. 
    
    * some `etable` arguments have been ported to the `style` functions (`yesNo`, `tablefoot`).
    
  - New `postprocess.tex` and `postprocess.df` arguments which allow the automatic postprocessing of the outputs. See the dedicated vignette on exporting tables for an illustration.
    
  - new `tabular` arguments which allows to create `tabular*` tables (suggestion by @fostermeijer, #51).

  - polynomials and powers are automatically renamed to facilitate comparison across models. You can set their style with the argument `poly_dict`.
  
  - the labeling of models is enhanced when `rep.fixest` is used with different standard-errors (the model names are now "model INDEX.SUB-INDEX").
  
  - the argument `subtitles` has been improved, and now automatically displays the samples when split sample estimations are performed.
 
## Other new features

 - In all estimations:
    
    * `subset`: regular subset (long overdue).
    
    * `split`, `fsplit`: to perform split sample estimations.
    
    * `se`, `cluster`: to cluster the standard-errors during the call.
    
    * `lean`: if `TRUE`, then summary is applied and any large object is removed from the result. To save memory => but many methods won't work afterwards.
    
    * `fixef.rm`: argument that accepts `none`, `perfect`, `singleton`, `both`. Controls the removal of fixed-effects from the observation.
    
    * auto parsing of powers. Now you don't need to use `I()` to have powers of variables in the RHS, it is automatically done for you (i.e. `x^3` becomes `I(x^3)`):
    ```R
    base = iris
    names(base) = c("y", "x1", "x2", "x3", "species")

    # The multiple estimation below works just fine
    feols(y ~ csw(x, x^2, x^3), base)
    ```
  - Estimation options can be set globally with `setFixest_estimation()`.
  
  - The `demean` function has been enhanced (with the contribution of Sebastian Krantz).
  
## Improvements of the internal algorithm

 - Internal demeaning algorithm: some copies of the data are avoided when using `feglm`.
 
 - Internal algorithm of `to_integer` (used in all estimations): one copy of the input data is now avoided.
 
 - All estimations: smarter handling of the intercept, thus avoiding the reconstruction of the design matrix.

# fixest 0.7.1 (2020-10-27)

## Hotfixes

 - Fix bug int overflow in estimations with only one variable.
 
 - Fix bug in tests occurring in R old release.
 
 - Fix bug in examples occurring in R old release.

## Improvements

 - Function `i()` now behaves as `factor()`, setting automatically a reference when appropriate.
 
 - Internal algorithm of `i()` is much faster.
 
## New features

 - In `etable`, the user can now provide a type of clustering for each model.
 
 - New method `rep.fixest` to replicate fixest objects, mostly useful in `etable` when several SEs for the same models are to be reported.
 
 - Automatic fix when the variance is not positive definite.


# fixest 0.7.0 (2020-10-24)

## Bugs

 - Major bug when fixed-effects were combined with `^` and they contained NAs (thanks to @poliquin, #35).
 
 - Bug when using lead/lags in estimations. The bug was due to a bug in a dependency ([dreamerr](https://cran.r-project.org/package=dreamerr)) and was fixed. Now fixest requires **dreamerr** version >= 1.2.1. Bug spotted by @seunghoon001 (#44).
 
 - Major bug when n_obs x n_vars > 2B or n_obs x n_fixed-effects > 2B. In such cases estimations could just not be done, even leading R to crash when using nthreads > 1. The algorithm was fixed to allow datasets with up to 2B observations to be estimated in all circumstances. Bug reported, and many help for checking provided, by Howard Zihao Zhang.
 
 - `coefplot`: Problem regarding interactions when observations, and hence coefficients, were removed from the estimation. Now the coefficients are removed from the plot. Bug reported by @phisherblack, #45.
 
 - `coefplot`: Corrected various bugs when asked for the plotting of several estimations. 
 
 - Fix the stack imbalance warning (report by @shoonlee, #46).
 
## Internal improvements

 - Brand new internal algorithm which now uses closed form solutions when dealing with variables with varying slopes. This means that when variables with varying slopes are present, the algorithm is incomparably faster and more accurate.

 - Two deep copies of some data are now avoided in the demeaning function. This improves the performance in terms of memory footprint, and also makes the algorithm faster. 
 
## Standard-errors, important changes
  
 - New default values for standard-errors (only concerns multiway clustering). They become similar to `reghdfe` to increase cross-software comparability. Computing the standard-errors the old way is still possible using the argument `dof`. See the dedicated vignette: [On standard errors](https://cran.r-project.org/package=fixest/vignettes/standard_errors.html).
 
 - Name change in `summary`/`vcov`/`etable`: To get heteroskedasticity-robust standard-errors, `se = "hetero"` now replaces `se = "white"` to enhance clarity. Note that `se = "white"` still works.
 
## New function: `fitstat` 

  - New function `fitsat` that computes various fit statistics. It is integrated with `etable` and can be invoked with the argument `fitstat`. So far only two fit statistics are included, but more will come.
 
## New features in `interact()`

  - You can now use `i(var)` to treat the variable `var` as a factor. You can select which values to drop/keep with the respective arguments. 
  
  - Using `i(var)` leads to a special treatment of these variables in the functions `coefplot` and `etable`.
 
## New features in `etable`

  - New argument `placement` to define the position of the float in Latex (suggestion by Caleb Kwon).
  
  - New argument `drop.section`, with which you can drop a) the fixed-effects, b) the variables with varying slopes, or c) the statistics, sections (suggestion by Caleb Kwon).
  
  - Fix glitch in help pages regarding the use of the '%' (percentage) character in regular expressions.
  
  - Two new arguments `.vcov` and `.vcov_args` to compute the standard-errors with custom functions.
  
  - The number of observations (`n`) is now treated as a regular statistic and can be placed where one wants.
  
  - The statistics can now have custom aliases using the argument `dict`.
  
  - The overdispersion becomes a regular fit statistic that can be included (or not) using `fitstat`.
  
  - The dictionnary now applies to the factors of interactions, and the values of factors.

## User visible changes

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
  
## Deprecation

 - Argument `na_inf.rm` has been removed. It was present for historical reasons, and is removed to increase code clarity.

# fixest 0.6.0 (2020-07-13)

## Bugs

 - In `vcov`, the degree-of-freedom in the small sample correction correction was fixed to "nested" and couldn't be modified, now corrected. Further, "nested" was not properly accounted for, now corrected.
 
 - In `etable`, `fitsat = FALSE` or `fitsat = NA` led to a bug.
 
 - `r2`: bug when the estimation contained only fixed effects (thanks to Luis Fonseca, #27).
 
 - Now the `BIC` of `feglm` is similar to the one of `glm`.
 
 - Bug in the log-likelihood in the presence of weights, now corrected.
 
 - Bug in `coefplot` when some interacted variables were removed because of collinearity. Now corrected.
 
## New vignettes

 - On standard-errors: how are the SEs computed in fixest and how to replicate the SEs from other software.
 
 - Exporting estimation tables: how to use efficiently `etable`, in particular how to customize the tables.
 
## Major changes: etable

  - New arguments: `group`, `extraline`, `notes`, `tablefoot`. 
  
    - `group` allows to eliminate variables (like `drop`) and adds an extra line with TRUE/FALSE if the model contained those variables.
    
    - `extraline` allows to add extra lines with any content.
    
    - `notes` allows to add notes after the table (suggestion by @bgchamps, #25).
    
    - `tablefoot` controls whether the table footer, containing the type of standard-errors and the significance codes, should be displayed.
    
  - Renaming: `yesNoFixef` => `yesNo`.
  
  - Most default values can be set globally with the new function `setFixest_etable`.
 
## Major changes: dof

 - Function `dof`, used to adjust the small sample corrections, is now much more complete and allows to replicate a large set of estimation results from alternative software.

## User visible changes

 - You can now provide custom VCOVs to summary by using the argument `.vcov`.
 
 - A warning is now prompted when the maximum number of iterations of the algorithm is reached (suggestion by @clukewatson , #24]). 
 
 - The types of standard-errors can now be set globally with the function `setFixest_se` (suggestion by @dlindzee, #28)
 
 - New `feols` argument `demeaned`. If `TRUE`, then the centered variables are returned (`y_demeaned` and `X_demeaned`). (Suggestion by Linus Holtermann.)
 
 - `interact` gains two new arguments: `drop` and `keep` (suggestion by @SuperMayo, #23).
 
## New methods

 - `hatvalues` has been implemented for feols and feglm estimations.
 
 - the `estfun` from `sandwich` has been implemented.

# fixest 0.5.1 (2020-06-18)

## Hotfix
 
 - Fixed bug introduced in the previous update (memory access error). Does not affect any of the results but could lead R to crash unexpectedly (odds were low though since access was adjacent).
 
## Bugs

 - Fix image link to the equation in the README.md.
 - Fix bug R2 and logLik when observations were removed because of NA values. Due to the update in `residuals.fixest`.
 
## User visible change

 - Rewriting of the internal algorithm computing the VCOV. 1) About 30% performance gain for estimations with many variables. 2) The code is much less memory hungry. 

## Major update of `etable`

 - New argument `style` which allows to set many elements of the output table. 
 - (minor) `signifCode` can be equal to `"letters"` to display letters instead of stars.
 
## Other

 - `setFixest_nthreads` now respects the `OMP_THREAD_LIMIT` environment variable.
 - Rd links are now made to the proper htlm files.


# fixest 0.5.0 (2020-06-10)

## Bug fixes
        
 - Bug with estimations with varying slopes if the fixed-effect relative to the slope is not in its decreasing order (thanks to Davide Proserpio).
 - Bug when interacting two variables with the `var::fe` syntax with `confirm = TRUE` and no reference.
 - Bug in `etable` when the standard-errors where `NA`.
 - Fixed very minor bug when computing the SEs (1e-6 difference).
 - Standard-errors in `feglm` for non-poisson, non-binomial families, are now correct (minor differences).
 - `fixef` did not work when the slope was an integer, now corrected (thanks to @clerousset, #20).

## New functionality: formula macros
        
 - You can use macros in formulas.
 - To set a macro variable, use e.g., `setFixest_fml(..ctrl = ~ var1 + var2)`. Here the macro variable `..ctrl` has been set to the value `"var1 + var2"`.
 - Now you can use this macro variable in any `fixest` estimation: e.g. `data(airquality) ; setFixest_fml(..ctrl = ~ Temp + Day) ; feols(Ozone ~ Wind + ..ctrl, airquality)`.
 - You can use macros in non-fixest estimations with `xpd`, which expands formulas. E.g. `lm(xpd(Ozone ~ Wind + ..ctrl), airquality)`.
 
## New functions

 - `to_integer`: user-level version of the internal algorithm transforming any kind of vector (or combination of vectors) into an integer ranging from 1 to the number of unique elements of the vector. Very fast.
 - `demean`: user-level version of the demeaning algorithm used in `feols`.

## Major user-visible changes
        
 - New internal algorithm to estimate OLS (applies to both `feols` and `feglm`):
 
    1. It is numerically more stable.
    
    2. Incomparably faster when factors are to be estimated (and not explicitly used as fixed-effects).
    
    3. Collinear variables are removed on the fly.

## User-visible changes
        
 - Interactions in `var::fe(ref)` now accept multiple references (i.e. `ref` can be a vector).
 - In `etable`, the variable names of non-Latex output can now be changed.
 - You can use the argument `n` when applying summary to choose the number of coefficients to display.
 - Argument `confirm` has been removed from the function `interact`.
 - `r2` allows more flexibility in the keywords it accepts.
 - Function `dof` gains a new argument `adj` which allows to make different types of common small sample corrections. Its other arguments have been renamed for clarity (`fixef` => `fixef.K`, `exact` => `fixef.exact`, `cluster` => `cluster.adj`).
 - Now t-statistics are used for `feols` and non-poisson, non-binomial models in `feglm`. For all other models, z-statistics are used. This complies with the default's R-stats behavior.

## New Methods
        
 - The `residuals` method has been substantially improved, now allowing different types.
 - New stats methods: sigma, deviance, weights.

## Vignette and Readme
        
 - Typos corrected.
 - Images in the Readme set to 1200px.
 
## Issue found: convergence problems with multiples variables with varying slopes
        
 - Convergence problems may arise in the presence of **multiple** variables with varying slopes. Theoretical work helped find a solution to this problem, but the implementation in R is proving not instantaneous.
 - In the meantime, now a warning is prompted when the algorithm suspects a convergence problem leading to poor precision of the estimated coefficients.

## Error-handling
        
 - Improved error-handling with [dreamerr](https://cran.r-project.org/package=dreamerr)'s functions.
 
## Other

 -  Dependency to MASS has been removed.

# fixest 0.4.1 (2020-04-13)

## Bug fixes
  
 - Major bug leading R to crash when using non-linear-in-parameters right-hand-sides in feNmlm. Only occured when some observations were removed from the data set (due to NAness or to perfect fit). [Thanks to @marissachilds, GH issue #17.]
 - In the `collinearity` help pages: an example could lead to an error (due to random data generation). It has been removed.
 - In `collinearity`, corrected the problem of display of the intercept in some situations.
 - Defaults for the arguments `cex` and `lwd` in `coefplot` have been changed to 1 and 1 (instead of par("cex") and par("lwd")). Otherwise this led to the creation of `Rplots.pdf` in the working directory (thanks to Kurt Hornik).
 - Corrected a typo in the article's title in the vignette.

## Help
  
 - Rewriting of sections, correction of small mistakes (wrong argument names), dropping completely the 'cluster' terminology (meant for fixed-effects), addition of what is contained in fixest objects.

## Other
  
 - Adding a README.md.
 - Small corrections in the vignette.

# fixest 0.4.0 (2020-03-27)

## User visible changes: Latex export
        
 - Better Latex special character escaping (errors reported by @dlindzee, #15).
 - New argument `fixef_sizes.simplify`, which provides the sizes of the fixed-effects in parentheses when there is no ambiguity.
 - You can suppress the line with the significance codes with `signifCode = NA`.
 - New argument `float` which decides whether to embed the table into a table environment. By default it is set to `TRUE` if a `title` or `label` is present.
 - New argument `keep` to select the variables to keep in the table.
 - New way to keep/drop/order variables with the special argument "\%". If you use "\%var", then it makes reference to the original variable name, not the aliased one (which is the default).
 - New argument `coefstat` defining what should be shown below the coefficients (standard-errors, t-stats or confidence intervals). Suggestion by @d712, #16.
 - Better rendering of significant digits.

## User visible changes: coefplot
        
 - Argument `horiz`. The coefficients can now be displayed horizontally instead of vertically.
 - The coefficient labels, when in the x-axis, can now be displayed in three different ways thanks to the new argument `lab.fit`: "simple", the classic axis, "multi", the labels appear across multiple lines to avoid collision, and "tilted" for tilted labels.
 - The margins now automatically fit.
 - Argument `style` allows you to set styles with the function `setFixest_coefplot`, you can then summon the style in `coefplot` with this argument.
 - Use the ampersand to set dictionary variables specific to `coefplot`.
 - Better display of groups (with the arguments `group` and `group.par`).

## New methods
        
 - `terms.fixest` giving the terms of the estimation.

## Other
        
 - All `donttest` sections were removed from help pages.


# fixest 0.3.1 (2020-02-09)

## Major bug fix
        
 - [panel] Fixed faulty memory access when taking the lead of a variable.

## Other bug fixes
        
 - [esttable/esttex] These two functions were replaced by the function `etable`. In the process, some of their arguments were "lost", this is now corrected.
 - [etable] Better escaping of special characters.
 - [estimations] Bug when particular non-numeric vectors were used in explanatory variables.

## New features
        
 - [coefplot] The function `coefplot` now accepts lists of estimations.


# fixest 0.3.0 (2020-02-01)

## New feature: Lagging
         
 - You can now add lags and leads in any `fixest` estimations. You only need to provide the panel identifiers with the new argument `panel.id`, then you're free to use the new functions `l()` for lags and `f()` for leads.
 
 - You can also set up a panel data set using the function `panel` which allows you to use the lagging functions without having to provide the argument `panel.id`, and which dispose of more options for setting the panel.

## New feature: Interactions
         
 - You can now add interactions in formulas with a new syntax: `var::fe(ref)`
 
 - The command `var::fe(ref)` interacts the variable `var` with each value of `fe` and sets `ref` as a reference. Note that if you don't use the argument `ref`, the command `var::fe` is identical to `var:factor(fe)`.
 
 - Using `var::fe(ref)` to write interactions opens up a special treatment of such variables in the exporting function `etable` and in the coefficient plotting function `coefplot`.

## New feature: `coefplot`
         
 - You can plot coefficients and their associated confidence intervals with the function `coefplot`.
 
 - `coefplot` dispose of many options, whose default values can be set with the function `setFixest_coefplot`.
 
 - As for the function `etable`, you can easily rename/drop/order the coefficients.
 
 - `coefplot` detects when interactions have been used and offers a special display for it.

## New functions
         
 - [etable] Estimations table: new function to export the results of multiple estimations. Replaces the two functions `esttex` and `esttable` (the two functions still exist but they will be deprecated in the future).
 
 - [Lagging] New functions related to lagging: `l`, `f`, `panel`, `unpanel` and `[.fixest_panel`.
 
 - [Utilities] A set of small utility functions has been added. They allow to extract part a coefficient table or parts of it (like the t-statistics of the standard-error) from an estimation. These functions are `coeftable`, `ctable` (an alias to `coeftable`), `se`, `tstat` and `pvalue`.
 
 - [coefplot] The functions `coefplot` and `setFixest_coefplot`.
 
 - [dof] New function to set the type of degree of freedom adjustment when computing the variance-covariance matrix. You can permanently set the type of DoF adjustment with the new function setFixest_dof().

## User visible changes
         
 - [all estimations] A key pre-processing step has been paralellized => algorithm faster in general and much faster for multi-FEs.
 - [predict & fitted] Predict and fitted now returns vectors of the length equal to the one of original data.
 - [standard-errors] New ways to compute the standard-errors have been implemented. In particular, now it account for the "nestedness" of the fixed-effects in the clusters by default. You can freely change how to compute the degrees of freedom correction with the function dof().
 - [r2] Computation of the within-R2 for feglm models is now self-contained.
 - [all estimations] New, more accurate, stopping criterion for 2+ fixed-effects.
 - [feols] Estimations are slightly faster.
 - [etable/esttex] When there are interactions, R may change the order of the interactions, making two interactions in two different estimations look different while they are in fact the same (e.g. x3:x2 and x2:x3). Now esstable automatically reorders the interactions when needed for comparison across estimations.
 - [etable/esttable] The type of standard errors is now always shown.
 - [etable/esttex] The aliases provided by 'dict' are also applied within interactions. For example: `dict=c(x1="Wind", x2="Rain")`, with an estimation with the following variables 'x1', 'x2', 'x1:x2' will lead to the following aliases in Latex 'Wind', 'Rain' and 'Wind $times$ Rain'.
 - [etable/esttex] Interactions of similar values but of different order (e.g. x1:x2 and x2:x1) are reorderd to appear in the same lines.
 - [etable/esttex] The i) type of standard errors and ii) the significance codes, are now displayed in two separate lines (otherwise the line would be too wide).
 - [etable/esttex] Argument `yesNoFixef` can be of length one, defaulting the second element to the empty string.
 - [etable/esttex] Escaping of Latex special characters is now much more robust.

## Bug correction
         
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

# fixest 0.2.1 (2019-11-22)

## Major bug correction
        
 - lag.formula: Bug introduced from previous update which could lead to wrong results. Now fixed.

## Major user visible changes
         
 - [All estimation methods] Significant speed improvement when the fixed-effects variables (i.e. the identifiers) are string vectors.


# fixest 0.2.0 (2019-11-19)

## New function
         
 -[did_means] New function `did_means` to conveniently compare means of groups of observations (both treat/control and pre/post). Contains tools to easily export in Latex.



## Major user visible changes
        
 - [All estimation methods] Significant speed improvement when the fixed-effects variables (i.e. the identifiers) are of type integer or double.
 - [esttex, esttable] New argument 'fitstat' to select which fit statistic to display. The default adapts to the models. Old arguments (loglik, bic, aic, sq.cor) are dropped.
 - [esttable] Significantly better rendering of SE types.
 - [r2] Now NA is returned for R2s that have no theoretical justification (e.g. within R2 when no FEs, or 'regular' R2 for ML models).



## Minor user visible changes
        
 - [did_plot_yearly_effects] Now the name of the dependent variable appears on the y-axis.
 - [esttex] Usage of the `sym` macro in Latex is dropped.



## Bug correction
        
 - [fixef.fixest] bug could appear when using varying slopes coefficients in some specific circumstances (when the slope FEs were different from the regular FEs).
 - [fixef.fixest] bug when many regular FEs jointly with varying slopes.
 - [fixef.fixest] regarding slope coefficients: now the algorithm also evaluates functions of variables.
 - [esttable] Width of the "separating lines" now appropriately set for long dependent variable names.
 - [esttex] Spelling mistake corrected.
 - [estimations] Bug could occur for extremely small data sets (< 10 observations).



## Error handling
        
 - [esttex, esttable] More informative error messages in functions esttex and esttable.


# fixest 0.1.2 (2019-10-04)
## Major bug correction
        
 - lag.formula: When the data was not in a particular format, the results could be wrong. Now corrected.


# fixest 0.1.1 (2019-09-20)
## Major bug correction
	    
 - feglm: bug when a) the deviance at initialization was higher than the deviance of the first iteration of the IRWLS and b) the step-halving was unable to find a lower deviance. This led the estimation to fail with an error although it should have been performed properly.
 - did_estimate_yearly_effects: bug when the estimation involved periods with negative values.

## Minor bug correction
        
 - esttex: bug regarding the number of digits of negative coefficients to be displayed
 - esttex: now properly escaping the percentage and the underscore for exports in Latex
 - esttex: bug when changing the names of the dependent variables using a dictionnary
 - vcov: some warning messages were misleading
 - update: bug update when using the argument nframes
 - update: bug when updating the function fepois

## Error handling
        
 - Better error messages for: did_estimate_yearly_effects, main estimation functions, setFixest_dict, fepois and fenegbin.

# fixest 0.1.0 (2019-09-03)

## First version
	    
 - This package is an effort to create a family of fast and user-friendly functions to perform estimations with multiple fixed-effects (F.E.).

 - Estimations with fixed-effects (or call it factor variables) is a staple in social science. Hence having a package gathering many methods with fast execution time is of prime importance. At the time of this version, this is the fastest existing method to perform F.E. estimations (often by orders of magnitude, compared to the most efficient alternative methods [both in R and Stata]). The underlying method to obtain the F.E. is based on Berge 2018, and the workhorse of the code is in c++ parallelized via OpenMP (btw thanks Rcpp for simplifying coders' life!).

 - This package is the follow up of the (now deprecated) package `FENmlm` which performed fixed-effects estimations but for only four likelihood families. Package `fixest` completely supersedes `FENmlm` by extending the method to regular OLS and all GLM families, and adding new utility functions. Further, the design of the functions has been completely overhauled and extended towards much more user-friendliness. Massive effort has been put into providing a set of informative error messages to the user for quick debugging of her workflow (e.g. one of the functions contains over 100 different errors messages).


