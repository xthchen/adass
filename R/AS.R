#' Adaptive selection method
#'
#' This function performs the AS method to estimate the number of effective populations K
#' @param R A numeric vector of R-values
#' @param lambda Numeric, the tuning parameter
#' @return The number of effective populations K
#' @export
#' 
AS = function(R, lambda){
  return(max(sum(R > lambda)/(1-lambda),2))
}
