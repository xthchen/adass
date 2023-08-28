# adass
This R package is used for selecting a subset of best populations from a set of normally distributed populations, with “best” being defined as the one with the largest true location parameter. The main function subset_selection() takes as input the population samples, and returns the R-values of each population as output. Two adaptive methods, Adaptive Selection (AS) and Randomised Adaptive Selection (RAS) are available to reduce the conservativeness of the result. For the tuning parameter lambda, an automatic tuning method lamda_tuning() is also available in this package

In order to run this package, the stats package, as well as [chi](https://github.com/dkahle/chi) and [cubature](https://github.com/bnaras/cubature) R packages are required.

## Example of subset_selection usage

_r_value = subset_selection(pop = samples, var_type = “known”, var = 1, method = “RAS”, lambda = “MSE”)_

Here samples is a list, with each item being a vector of samples from a population. The variance is assumed to be known with value 1. The Randomised Adaptive Selection method is used, as well as the automatic tuning method.

The unknown but equal variance and unknown and unequal variance cases are available by changing the input of var_type to “eq_var” or “uneq_var” respectively, and ignoring the var parameter.

Basic subset selection method without the adaptive approach is available via setting method = “base”, and the AS method is available via method = “AS”.

If the user prefer to manually set the tuning parameter lambda, simply give a numeric input to the parameter lambda. We suggest an input of smaller than 0.5 to avoid excess anticonservativeness when the data is close to the least favourable condition.

For further information and examples refer to the manual typing ?subset_selection in R.


