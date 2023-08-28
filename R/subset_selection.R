#' subset selection method
#'
#' This function performs the subset selection method, given normally distributed data, to produce R-values for all populations
#' @param pop A list. Each list item contains the sample vector of each population. 
#' @param var_type Factor, can be either "known", for the known variance case, if so specify the variance at parameter "var". "eq_var", for the unknown but equal variance case, and "uneq_var", for the unknown and unequal variance case.
#' @param var Numeric, value of variance when variance is known. If unknown variance, ignore this parameter.
#' @param method Factor, can be either "base" for the basic subset selection, AS for the Adaptive method, or RAS for the randomised Adaptive method.
#' @param lambda Factor, can be either "MSE" which uses the automatic lambda tuning method, or a predetermined number between 0 and 1.
#' @return A numeric vector of R-values, one for each population in the order of the input data.
#' @export
#' 
subset_selection = function(pop, var_type, var, method, lambda){
  k = length(pop)
  if (var_type == "known"){
    #known variance case
    samsize = lengths(pop)
      #equal sample size case
      pop_mean = c()
      for (j in 1:k){
        pop_mean = c(pop_mean, mean(pop[[j]]))
      }
      
      integrand = function(x){
        pnorm(x+d)^(k-1)*dnorm(x)
      }
      if (length(unique(lengths(pop)))==1){
        statistic = (max(pop_mean) - pop_mean)*sqrt(lengths(pop)[1])/sqrt(var)
      }else{
        statistic = (max(pop_mean) - pop_mean)/sqrt(sqrt(var)*0.5*((1/samsize[which(pop_mean == max(pop_mean))])+1/samsize))
      }
      R = c()
      for (d in statistic){
        R = c(R, 1 - integrate(integrand, -Inf, Inf)$value)
      }
      if (is.character(method)){
        if (method == "base"){
          #basic subset selection
          return (R)
        }
      }else{
        #adaptive method
        if (lambda == "MSE"){
          #lambda tuning
          lambda = lambda_tuning(R, method, k = k)
        }
        k_h = method(R, lambda)
        integrand_h = function(x){
          pnorm(x+d)^(k_h-1)*dnorm(x)
        }
        R_h = c()
        for (d in statistic){
          R_h = c(R_h, 1 - integrate(integrand_h, -Inf, Inf)$value)
        }
        return(R_h)
      }
  }else if (var_type == "eq_var"){
    #unknown but equal variance case
    samsize = lengths(pop)
    #equal sample size case
    pop_mean = c()
    var = 0
    for (j in 1:k){
      pop_mean = c(pop_mean, mean(pop[[j]]))
      var = var + sum((pop[[j]] - pop_mean[j])^2)/(samsize[j]-1)
    }
    df = sum(samsize) - k
    var = var/df
    
    integrand = function(arg){
      u = arg[1]
      x = arg[2]
      pnorm(d*u+x)^(k-1)*dnorm(x)*dchi(u*sqrt(df),df)*sqrt(df)
    }
    
    if (length(unique(lengths(pop)))==1){
      statistic = (max(pop_mean) - pop_mean)*sqrt(lengths(pop)[1])/sqrt(var)
    }else{
      statistic = (max(pop_mean) - pop_mean)/sqrt(sqrt(var)*0.5*((1/samsize[which(pop_mean == max(pop_mean))])+1/samsize))
    }
    
    R = c()
    for (d in statistic){
      R = c(R, 1-pcubature(integrand, lowerLimit = c(0,-Inf), upperLimit = c(Inf, Inf))$integral)
    }
    if (is.character(method)){
      if (method == "base"){
        #basic subset selection
        return (R)
      }
    }else{
      #adaptive method
      if (lambda == "MSE"){
        #lambda tuning
        lambda = lambda_tuning(R, method, k = k)
      }
      k_h = method(R, lambda)
      
      integrand_h = function(arg){
        u = arg[1]
        x = arg[2]
        pnorm(d*u+x)^(k_h-1)*dnorm(x)*dchi(u*sqrt(df),df)*sqrt(df)
      }
      
      R_h = c()
      for (d in statistic){
        R_h = c(R_h, 1 - pcubature(integrand_h, lowerLimit = c(0,-Inf), upperLimit = c(Inf, Inf))$integral)
      }
      return(R_h)
    }
  }else if (var_type == "uneq_var"){
    #unknown and unequal variance case
    samsize = lengths(pop)
    #equal sample size case
    pop_mean = c()
    var = c()
    for (j in 1:k){
      pop_mean = c(pop_mean, mean(pop[[j]]))
      var = c(var, sum((pop[[j]] - pop_mean[j])^2)/(samsize[j]-1))
    }
    
    integrand = function(arg){
      u = arg[1]
      x = arg[2]
      pnorm(d*u+x)^(k-1)*dnorm(x)*dchi(u*sqrt(df),df)*sqrt(df)
    }
    
    max_mean_samsize = samsize[which(pop_mean == max(pop_mean))]
    max_mean_var = var[which(pop_mean == max(pop_mean))]
    
    statistic = (max(pop_mean) - pop_mean)/(sqrt(0.5*((max_mean_var/max_mean_samsize)+var/samsize)))
    
    R = c()
    for (y in 1:k){
      d = statistic[y]
      df = (max_mean_var/max_mean_samsize + var[y]/samsize[y])^2/
        (max_mean_var^2/((max_mean_samsize^2)*(max_mean_samsize-1)) + var[y]^2/((samsize[y]^2)*samsize[y]-1))
      
      R = c(R, 1-pcubature(integrand, lowerLimit = c(0,-Inf), upperLimit = c(Inf, Inf))$integral)
      
    }
    
    if (is.character(method)){
      if (method == "base"){
        #basic subset selection
        return (R)
      }
    }else{
      #adaptive method
      if (lambda == "MSE"){
        #lambda tuning
        lambda = lambda_tuning(R, method, k = k)
      }
      k_h = method(R, lambda)
      
      integrand_h = function(arg){
        u = arg[1]
        x = arg[2]
        pnorm(d*u+x)^(k_h-1)*dnorm(x)*dchi(u*sqrt(df),df)*sqrt(df)
      }
      R_h = c()
      for (y in 1:k){
        d = statistic[y]
        df = (max_mean_var/max_mean_samsize + var[y]/samsize[y])^2/
          (max_mean_var^2/((max_mean_samsize^2)*(max_mean_samsize-1)) + var[y]^2/((samsize[y]^2)*samsize[y]-1))
        #df = sum(samsize)-k
        R_h = c(R_h, 1 - pcubature(integrand_h, lowerLimit = c(0,-Inf), upperLimit = c(Inf, Inf))$integral)
      }
      return(R_h)
    }
    }
  
  }
