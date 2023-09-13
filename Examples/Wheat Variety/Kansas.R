Kansas_2022 = read.csv("./Wheat Variety/Kansas Winter Wheat productivity - 2022.csv")
Kansas_2022 = Kansas_2022[rowSums(!is.na(Kansas_2022)) > 6,]
#The variety needs to be present at more than 6 farms to be included

farm_mean = c()
for (i in 2:ncol(Kansas_2022)){
  farm_mean = c(farm_mean, mean(as.numeric(na.omit(Kansas_2022[-1,i]))))
}
################################################################################
#Option 1, unknown equal variance unequal sample size

num_pop = nrow(Kansas_2022) - 1
df = sum(!is.na(as.matrix(Kansas_2022[-1,-1]))) - num_pop
samsize_all = c()
mean_all = c()
var = 0
for (i in 2:nrow(Kansas_2022)){
  variety = as.numeric(Kansas_2022[i,-1])
  variety = variety - farm_mean
  variety = variety[!is.na(variety)]
  mean_all = c(mean_all, mean(variety))
  samsize_all = c(samsize_all, length(variety))
  var = var + sum((mean(variety) - variety)^2)
}
var = var/df
#Obtaining all required data
integrand = function(arg){
  u = arg[1]
  x = arg[2]
  pnorm(d*u+x)^(num_pop-1)*dnorm(x)*dchi(u*sqrt(df),df)*sqrt(df)
}
max_loc = which(mean_all == max(mean_all))
statistic = (max(mean_all) - mean_all)/(sqrt(var)*0.5*sqrt(1/samsize_all[max_loc] + 1/samsize_all))
#Test statistic
R = c()
for (d in statistic){
  R = c(R, 1-pcubature(integrand, lowerLimit = c(0,-Inf), upperLimit = c(Inf, Inf))$integral)
}
length(R[R>0.05])
#11 varieties are best.
Kansas_2022[order(R, decreasing  = TRUE)[1:11],1]
#The best are: "KIVARI_AX", "WB4401", "AP18_AX",  "GUARDIAN", "CP7266AX",     
#"STRAD_CL_PLUS", "AP_EVERROCK", "CRESCENT_AX", "LCS_RUNNER", "BOB_DOLE"     
#"LCS_VALIANT" 
#adaptive using effective method
h = 0.5
k_h = max(sum(R > 1-h)/h,2)
#8 effective populations
df_h = sum(!is.na(as.matrix(Kansas_2022[order(R, decreasing  = TRUE)[1:k_h],-1]))) - k_h
#degrees of freedom for the 8 effective populations
integrand_h = function(arg){
  u = arg[1]
  x = arg[2]
  pnorm(d*u+x)^(k_h-1)*dnorm(x)*dchi(u*sqrt(df_h),df_h)*sqrt(df_h)
}

R_h = c()
for (d in statistic){
  R_h = c(R_h, 1-pcubature(integrand_h, lowerLimit = c(0,-Inf), upperLimit = c(Inf, Inf))$integral)
}
length(R_h[R_h>0.05])
#7 are best
Kansas_2022[order(R_h, decreasing  = TRUE)[1:7],1]
#name of best
#"KIVARI_AX", "WB4401", "AP18_AX",  "GUARDIAN", "CP7266AX", "STRAD_CL_PLUS"
#"AP_EVERROCK"
################################################################################
num_pop = nrow(Kansas_2022) - 1
df = sum(!is.na(as.matrix(Kansas_2022[-1,-1]))) - num_pop

samsize_all = c()
mean_all = c()
var_all = c()
for (i in 2:nrow(Kansas_2022)){
  variety = as.numeric(Kansas_2022[i,-1])
  variety = variety - farm_mean
  variety = variety[!is.na(variety)]
  mean_all = c(mean_all, mean(variety))
  samsize_all = c(samsize_all, length(variety))
  var_all = c(var_all, (sum((mean(variety) - variety)^2)/(length(variety)-1)))
}
#Obtaining all required data
integrand = function(arg){
  u = arg[1]
  x = arg[2]
  pnorm(d*u+x)^(num_pop-1)*dnorm(x)*dchi(u*sqrt(df),df)*sqrt(df)
}
max_loc = which(mean_all == max(mean_all))
statistic = (max(mean_all) - mean_all)/(0.5*sqrt(var_all[max_loc]/samsize_all[max_loc] + var_all/samsize_all))

#Test statistic
R = c()
for (d in statistic){
  R = c(R, 1-pcubature(integrand, lowerLimit = c(0,-Inf), upperLimit = c(Inf, Inf))$integral)
}
length(R[R>0.05])
#12 varieties are best
Kansas_2022[order(R, decreasing  = TRUE)[1:12],1]
#name of best
#"KIVARI_AX", "WB4401", "AP18_AX",  "GUARDIAN", "STRAD_CL_PLUS"
#"CP7266AX",  "LCS_RUNNER", "BOB_DOLE", "CRESCENT_AX",  "AP_EVERROCK"  
#"WB4595",  "LCS_VALIANT"  
h = 0.5
k_h = max(sum(R > 1-h)/h,2)
#8 effective populations
df_h = sum(!is.na(as.matrix(Kansas_2022[order(R, decreasing  = TRUE)[1:k_h],-1]))) - k_h
#degrees of freedom for the 8 effective populations
integrand_h = function(arg){
  u = arg[1]
  x = arg[2]
  pnorm(d*u+x)^(k_h-1)*dnorm(x)*dchi(u*sqrt(df_h),df_h)*sqrt(df_h)
}

R_h = c()
for (d in statistic){
  R_h = c(R_h, 1-pcubature(integrand_h, lowerLimit = c(0,-Inf), upperLimit = c(Inf, Inf))$integral)
}
length(R_h[R_h>0.05])
#10 best
Kansas_2022[order(R_h, decreasing  = TRUE)[1:10],1]
#name of best
#"KIVARI_AX", "WB4401", "AP18_AX",  "GUARDIAN", "STRAD_CL_PLUS",
#"CP7266AX",  "LCS_RUNNER", "BOB_DOLE", "CRESCENT_AX",  "AP_EVERROCK" 