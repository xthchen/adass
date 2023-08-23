#' Automatic tuning of lambda through MSE and bootstrap
#'
#' This function performs the AS method to estimate the number of effective populations K
#' @param R A numeric vector of R-values
#' @param method Factor, method used for obtaining the number of effective populations. Can be "AS" or "RAS".
#' @param S Numeric, number of bootstrap replicates, defaults to 100.
#' @param k Numeric, total number of populations
#' @return the number of effective populations K
#' @export
#' 

lambda_tuning = function(R, method, k, S = 100){
  lambda_choice = seq(0.05,0.95,0.05)
  #candidates of lambda
  
    lambda_K = k
    lambda_dash = 0
    for (i in lambda_choice){
      if(method(R, i) < lambda_K){
        lambda_K = AS(R, i)
        lambda_dash = i
      }
    }
    MSE = 10^6
    #some large dummy number
    for (j in lambda_choice){
      MSE_hold = 0
      for (s in 1:S){
        MSE_hold = MSE_hold + (method(sample(R, length(R), replace = TRUE),j) - lambda_K)^2
      }
      if ((MSE_hold/S) < MSE){
        lambda = j
        MSE = MSE_hold/S
      }
    }
  return(lambda)
}
