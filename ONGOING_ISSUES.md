
## Ongoing issues

#### Compatibility with the sandwich package

So far `fixest` objects are compatible with variances computed with the package `sandwich` with the following exceptions:

 - `vcovHC` with `type = "HC4"`: The compatibility is only partial, i.e. it works for models without fixed-effects but breaks with fixed-effects estimations. This is because the hatvalues, needed for this kind of VCOV correction, are not defined when fixed-effects are used (actually I could make it compatible but the model would need to be reestimated with dummy variables, and this hardly squares with the idea of efficient fixed-effects estimation).
 
 - `vcovBS`, i.e. bootstraped standard-errors: Absolutely not compatible so far, but it will become compatible later (hopefully).




