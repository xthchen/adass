#' Randomised Adaptive selection method
#'
#' This function performs the AS method to estimate the number of effective populations K
#' @param R A numeric vector of R-values
#' @param lambda Numeric, the tuning parameter
#' @return The number of effective populations K
#' @export
#' 
RAS = function(R, lambda){
  c = 0
  c_old = 0
  #computation of c through maximisation
  for (l in seq(0,1,0.01)){
    c_new = 0
    for (j in 1:length(R)){
      if (R[j] >= l){
        c_new = c_new+lambda
      }else if(R[j] <= (lambda*l)){
        c_new = c_new+1
      }
    }
    if (c_new > c_old){
      c = l
      c_old = c_new
    }
  }
  
  R_rand = c()
  #randomising R-values
  for (j in 1:length(R)){
    if (R[j] >= c){
      R_rand = c(R_rand, runif(1,0,1))
    }else{
      R_rand = c(R_rand, R[j]/c)
    }
  }
  return(max(sum(R_rand > lambda)/(1-lambda),2))
}


