## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, comment = "#>")

if(requireNamespace("pander", quietly = TRUE)) library(pander)

require_pander_ON = function(){
  if(!requireNamespace("pander", quietly = TRUE)){
    knitr::opts_chunk$set(eval = FALSE)
    cat("Evaluation of the next chunks requires 'pander', which is not present.")
  }
}

require_pander_OFF = function(){
  knitr::opts_chunk$set(eval = TRUE)
}

library(fixest)
setFixest_notes(FALSE)
setFixest_etable(digits = 3)

## ---- eval = TRUE, results = "hide"-------------------------------------------
library(fixest)
data(airquality)

# Setting a dictionary 
setFixest_dict(c(Ozone = "Ozone (ppb)", Solar.R = "Solar Radiation (Langleys)",
                 Wind = "Wind Speed (mph)", Temp = "Temperature"))

## -----------------------------------------------------------------------------
# On multiple estimations: see the dedicated vignette
est = feols(Ozone ~ Solar.R + sw0(Wind + Temp) | csw(Month, Day), 
            airquality, cluster = ~Day)

## -----------------------------------------------------------------------------
etable(est)

## -----------------------------------------------------------------------------
etable(est, style.df = style.df(depvar.title = "", fixef.title = "", 
                                 fixef.suffix = " fixed effect", yesNo = "yes"))

## ---- echo = FALSE------------------------------------------------------------
require_pander_ON()

## -----------------------------------------------------------------------------
library(pander)

etable(est, postprocess.df = pandoc.table.return, style = "rmarkdown")

## -----------------------------------------------------------------------------
my_style = style.df(depvar.title = "", fixef.title = "", 
                    fixef.suffix = " fixed effect", yesNo = "yes")
setFixest_etable(style.df = my_style, postprocess.df = pandoc.table.return)

## -----------------------------------------------------------------------------
etable(est[rhs = 2], style = "rmarkdown", caption = "New default values")

## ---- echo = FALSE------------------------------------------------------------
require_pander_OFF()

## -----------------------------------------------------------------------------
est_slopes = feols(Ozone ~ Solar.R + Wind | Day + Month[Temp], airquality)

## ---- results = 'hide'--------------------------------------------------------
etable(est, est_slopes, tex = TRUE)

## ---- include = FALSE---------------------------------------------------------
# etable(est, est_slopes, file = "../_VIGNETTES/vignette_etable.tex", replace = TRUE)
# etable(est, est_slopes, file = "../_VIGNETTES/vignette_etable.tex", style.tex = style.tex("aer"), fitstat = ~ r2 + n, signifCode = NA)

## ---- results = 'hide'--------------------------------------------------------
etable(est, est_slopes, style.tex = style.tex("aer"), 
       signifCode = NA, fitstat = ~ r2 + n, tex = TRUE)

## -----------------------------------------------------------------------------
set_rules = function(x, heavy, light){
  # x: the character vector returned by etable
  
  tex2add = ""
  if(!missing(heavy)){
    tex2add = paste0("\\setlength\\heavyrulewidth{", heavy, "}\n")
  }
  if(!missing(light)){
    tex2add = paste0(tex2add, "\\setlength\\lightrulewidth{", light, "}\n")
  }
  
  if(nchar(tex2add) > 0){
    x[x == "%start:tab\n"] = tex2add
  }
  
  x
}

## ---- results = 'hide'--------------------------------------------------------
etable(est, est_slopes, postprocess.tex = set_rules, heavy = "0.14em", tex = TRUE)

## -----------------------------------------------------------------------------
setFixest_etable(style.tex = style.tex("aer"), postprocess.tex = set_rules, 
                 fitstat = ~ r2 + n, signifCode = NA)

## -----------------------------------------------------------------------------
etable(est, est_slopes, heavy = "0.14em", tex = TRUE)

## -----------------------------------------------------------------------------
fitstat_register(type = "p_s", alias = "pvalue (standard)",
                 fun = function(x) pvalue(x, se = "s")["Solar.R"])

fitstat_register(type = "p_h", alias = "pvalue (Heterosk.)",
                 fun = function(x) pvalue(x, se = "h")["Solar.R"])

fitstat_register(type = "p_day", alias = "pvalue (Day)",
                 fun = function(x) pvalue(x, cluster = "Day")["Solar.R"])

fitstat_register(type = "p_month", alias = "pvalue (Month)",
                 fun = function(x) pvalue(x, cluster = "Month")["Solar.R"])

# We first reset the default values set in the previous sections
setFixest_etable(reset = TRUE)
# Now we display the results with the new fit statistics
etable(est, fitstat = ~ . + p_s + p_h + p_day + p_month)

## ---- eval = FALSE------------------------------------------------------------
#  summary(.l(est, est_slopes), cluster = ~ Month)

