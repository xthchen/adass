#If "devtools" not already installed
install.packages("devtools")
#install "adass" package
install_github("xthchen/adass")
#load package
#The "adass" package should also automatically load the packages "chi" and "cubature".
library("adass")
library("stats")
################################################################################
#different subset selection methods
methods = c("base", AS, RAS)
#colours for the plot
colourblind = c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#Global simulation function
sim = function(data, var_type, var, method, lambda, loc){
  R = c()
  R_AS = c()
  R_RAS = c()
  for (i in 1:10^5){
    #10^5 simulations
    pop = data(i)
    #data
    R = c(R, subset_selection(pop = pop, var_type = var_type, var = var, method = method[[1]])[loc])
    #no adaptive
    R_AS = c(R_AS, subset_selection(pop = pop, var_type = var_type, var = var, method = methods[[2]], lambda = lambda)[loc])
    #adaptive selection
    R_RAS = c(R_RAS, subset_selection(pop = pop, var_type = var_type, var = var, method = methods[[3]], lambda = lambda)[loc])
    #randomised adaptive selection
  }
  prob = c()
  prob_as = c()
  prob_ras = c()
  for (j in seq(0,0.5,0.001)){
    #corresponding 1 - PCS for the level of significance up to 0.5
    prob = c(prob, sum(R < j)/length(R))
    prob_as = c(prob_as, sum(R_AS < j)/length(R_AS))
    prob_ras = c(prob_ras, sum(R_RAS < j)/length(R_RAS))
  }
  return(cbind(prob,prob_as,prob_ras))
}
###############################Figure 1#########################################
#Figure 1_1, far from LFC case with equal sample size and known and equal variance.
sim_data = function(seed = seed){
  set.seed(seed)
  return(c(rnorm(10,mean = 0),rnorm(90,mean = -2)))
}
#1 - PCS for the 3 methods
PCS = sim(data = sim_data, var_type = "known", var = 1, method = methods, lambda = 0.5, loc = 1)
#plot
plot(1, type = "n", xlim = c(0,0.5), ylim = c(0,0.5),
     xlab = expression(level~of~significance(alpha)),
     ylab = "1 - PCS of best populations")
lines(seq(0,0.5,0.001), PCS[,1], type = "l", col = colourblind[2], lty = "longdash")
lines(seq(0,0.5,0.001), PCS[,2], type = "l", col = colourblind[4], lty = "dotted")
lines(seq(0,0.5,0.001), PCS[,3], type = "l", col = colourblind[8], lty = "dotdash")
lines(c(0,0.5), c(0,0.5), type = "l", col = colourblind[1])
legend(0,0.5,legend = c("No Adaptive", expression(AS~lambda==0.5),
                         expression(RAS~lambda==0.5)),
       lty = c("longdash","dotted","dotdash"), col = colourblind[c(2,4,8)])
#Figure 1_2, LFC case with equal sample size and known and equal variance.
sim_data = function(seed = seed){
  set.seed(seed)
  return(c(rnorm(100,mean = 0),rnorm(0,mean = -2)))
}
#1 - PCS for the 3 methods
PCS = sim(data = sim_data, var_type = "known", var = 1, method = methods, lambda = 0.5, loc = 1)
#plot
plot(1, type = "n", xlim = c(0,0.5), ylim = c(0,0.5),
     xlab = expression(level~of~significance(alpha)),
     ylab = "1 - PCS of best populations")
lines(seq(0,0.5,0.001), PCS[,1], type = "l", col = colourblind[2], lty = "longdash")
lines(seq(0,0.5,0.001), PCS[,2], type = "l", col = colourblind[4], lty = "dotted")
lines(seq(0,0.5,0.001), PCS[,3], type = "l", col = colourblind[8], lty = "dotdash")
lines(c(0,0.5), c(0,0.5), type = "l", col = colourblind[1])
legend(0,0.5,legend = c("No Adaptive", expression(AS~lambda==0.5),
                        expression(RAS~lambda==0.5)),
       lty = c("longdash","dotted","dotdash"), col = colourblind[c(2,4,8)])
###############################Figure 2#########################################
#different percentage of true best
PCS = c()
for (alpha in seq(5,100,5)){
  #data
  sim_data = function(seed = seed){
    set.seed(seed)
    return(c(rnorm(alpha,mean = 0),rnorm((100-alpha),mean = -2)))
  }
  PCS = rbind(PCS, sim(data = sim_data, var_type = "known", var = 1, method = methods, lambda = 0.5, loc = 1)[51,])
  #fixing the level of significance to 0.05 = seq(0,0.5,0.001)[51]
}
plot(1, type = "n", xlim = c(0,100), ylim = c(0,0.07),
     xlab = "true number of best",
     ylab = expression(1-PCS~of~best~populations~alpha==0.05))
lines(seq(5,100,5), PCS[,1], type = "l", col = colourblind[2], lty = "longdash")
lines(seq(5,100,5), PCS[,2], type = "l", col = colourblind[4], lty = "dotted")
lines(seq(5,100,5), PCS[,3], type = "l", col = colourblind[8], lty = "dotdash")
lines(c(5,100), c(0.05,0.05), type = "l", col = colourblind[1])
legend(0,0.07,legend = c("No Adaptive", expression(AS~lambda==0.5),
                        expression(RAS~lambda==0.5)),
       lty = c("longdash","dotted","dotdash"), col = colourblind[c(2,4,8)])
###############################Figure 3#########################################
#Figure 3_1, far from LFC case with equal sample size and known and equal variance.
sim_data = function(seed = seed){
  set.seed(seed)
  return(c(rnorm(10,mean = 0),rnorm(90,mean = -2)))
}
#1 - PCS for the 3 methods
PCS = sim(data = sim_data, var_type = "known", var = 1, method = methods, lambda = 0.5, loc = 100)
#plot
plot(1, type = "n", xlim = c(0,0.5), ylim = c(0,1),
     xlab = expression(level~of~significance(alpha)),
     ylab = "power")
lines(seq(0,0.5,0.001), PCS[,1], type = "l", col = colourblind[2], lty = "longdash")
lines(seq(0,0.5,0.001), PCS[,2], type = "l", col = colourblind[4], lty = "dotted")
lines(seq(0,0.5,0.001), PCS[,3], type = "l", col = colourblind[8], lty = "dotdash")
lines(c(0,0.5), c(0,0.5), type = "l", col = colourblind[1])
legend(0,1,legend = c("No Adaptive", expression(AS~lambda==0.5),
                        expression(RAS~lambda==0.5)),
       lty = c("longdash","dotted","dotdash"), col = colourblind[c(2,4,8)])
#Figure 3_2, varying percentage of best, equal sample size and known and equal variance scenario
PCS = c()
for (alpha in seq(0,95,5)){
  #data
  sim_data = function(seed = seed){
    set.seed(seed)
    return(c(rnorm(alpha,mean = 0),rnorm((100-alpha),mean = -2)))
  }
  PCS = rbind(PCS, sim(data = sim_data, var_type = "known", var = 1, method = methods, lambda = 0.5, loc = 100)[51,])
  #fixing the level of significance to 0.05 = seq(0,0.5,0.001)[51]
}
plot(1, type = "n", xlim = c(0,100), ylim = c(0,0.6),
     xlab = "true number of best",
     ylab = "power")
lines(seq(0,95,5), PCS[,1], type = "l", col = colourblind[2], lty = "longdash")
lines(seq(0,95,5), PCS[,2], type = "l", col = colourblind[4], lty = "dotted")
lines(seq(0,95,5), PCS[,3], type = "l", col = colourblind[8], lty = "dotdash")
legend(60,0.3,legend = c("No Adaptive", expression(AS~lambda==0.5),
                         expression(RAS~lambda==0.5)),
       lty = c("longdash","dotted","dotdash"), col = colourblind[c(2,4,8)])
###############################Figure 4#########################################
#Figure 4_1, far from LFC case with equal sample size and known and equal variance
#comparison with lambda automatically tuned when using the RAS approach.
sim_data = function(seed = seed){
  set.seed(seed)
  return(c(rnorm(10,mean = 0),rnorm(90,mean = -2)))
}
#1 - PCS for the fix lambda case
PCS = sim(data = sim_data, var_type = "known", var = 1, method = methods, lambda = 0.5, loc = 1)
#for the automatic tuning case
PCS_lambda = sim(data = sim_data, var_type = "known", var = 1, method = methods, lambda = "MSE", loc = 1)
#plot
plot(1, type = "n", xlim = c(0,0.5), ylim = c(0,0.5),
     xlab = expression(level~of~significance(alpha)),
     ylab = "1 - PCS of best populations")
lines(seq(0,0.5,0.001), PCS[,1], type = "l", col = colourblind[2], lty = "longdash")
lines(seq(0,0.5,0.001), PCS[,3], type = "l", col = colourblind[4], lty = "dotted")
lines(seq(0,0.5,0.001), PCS_lambda[,3], type = "l", col = colourblind[8], lty = "dotdash")
lines(c(0,0.5), c(0,0.5), type = "l", col = colourblind[1])
legend(0,0.5,legend = c("No Adaptive", expression(RAS~lambda==0.5),
                        expression(RAS~lambda:MSE~tuned)),
       lty = c("longdash","dotted","dotdash"), col = colourblind[c(2,4,8)])
#Figure 4_2, LFC case with equal sample size and known and equal variance.
#comparison with lambda automatically tuned when using the RAS approach.
sim_data = function(seed = seed){
  set.seed(seed)
  return(c(rnorm(100,mean = 0),rnorm(0,mean = -2)))
}
#1 - PCS for the fixed lambda case
PCS = sim(data = sim_data, var_type = "known", var = 1, method = methods, lambda = 0.5, loc = 1)
#for the automatic tuning case
PCS_lambda = sim(data = sim_data, var_type = "known", var = 1, method = methods, lambda = "MSE", loc = 1)
#plot
plot(1, type = "n", xlim = c(0,0.5), ylim = c(0,0.5),
     xlab = expression(level~of~significance(alpha)),
     ylab = "1 - PCS of best populations")
lines(seq(0,0.5,0.001), PCS[,1], type = "l", col = colourblind[2], lty = "longdash")
lines(seq(0,0.5,0.001), PCS[,3], type = "l", col = colourblind[4], lty = "dotted")
lines(seq(0,0.5,0.001), PCS_lambda[,3], type = "l", col = colourblind[8], lty = "dotdash")
lines(c(0,0.5), c(0,0.5), type = "l", col = colourblind[1])
legend(0,0.5,legend = c("No Adaptive", expression(RAS~lambda==0.5),
                        expression(RAS~lambda:MSE~tuned)),
       lty = c("longdash","dotted","dotdash"), col = colourblind[c(2,4,8)])
###############################Figure 5#########################################
#Figure 5_1, far from LFC case with unequal sample size randomly chosen between 5 to 10
#known and equal variance.
sim_data = function(seed = seed){
  set.seed(seed)
  samsize = sample(5:10,100, replace = TRUE)
  samvar = rep(1,100)
  #sample size of the 100 populations
  pop = list()
  for (j in 1:10){
    pop = append(pop, list(rnorm(samsize[j], mean = 0, sd = sqrt(samvar[j]))))
  }
  #best populations
  for (j in 11:100){
    pop = append(pop, list(rnorm(samsize[j], mean = -2, sd = sqrt(samvar[j]))))
  }
  #non-best
  return(pop)
}
#1 - PCS for the 3 methods
PCS = sim(data = sim_data, var_type = "known", var = 1, method = methods, lambda = 0.1, loc = 1)
#plot
plot(1, type = "n", xlim = c(0,0.5), ylim = c(0,0.5),
     xlab = expression(level~of~significance(alpha)),
     ylab = "1 - PCS of best populations")
lines(seq(0,0.5,0.001), PCS[,1], type = "l", col = colourblind[2], lty = "longdash")
lines(seq(0,0.5,0.001), PCS[,2], type = "l", col = colourblind[4], lty = "dotted")
lines(seq(0,0.5,0.001), PCS[,3], type = "l", col = colourblind[8], lty = "dotdash")
lines(c(0,0.5), c(0,0.5), type = "l", col = colourblind[1])
legend(0,0.5,legend = c("No Adaptive", expression(AS~lambda==0.1),
                        expression(RAS~lambda==0.1)),
       lty = c("longdash","dotted","dotdash"), col = colourblind[c(2,4,8)])
#Figure 5_2, LFC case with unequal sample size randomly chosen between 5 to 10.
#known and equal variance.
sim_data = function(seed = seed){
  set.seed(seed)
  samsize = sample(5:10,100, replace = TRUE)
  samvar = rep(1,100)
  #sample size of the 100 populations
  pop = list()
  for (j in 1:100){
    pop = append(pop, list(rnorm(samsize[j], mean = 0, sd = sqrt(samvar[j]))))
  }
  #best populations
  return(pop)
}
#1 - PCS for the 3 methods
PCS = sim(data = sim_data, var_type = "known", var = 1, method = methods, lambda = 0.1, loc = 1)
#plot
plot(1, type = "n", xlim = c(0,0.5), ylim = c(0,0.5),
     xlab = expression(level~of~significance(alpha)),
     ylab = "1 - PCS of best populations")
lines(seq(0,0.5,0.001), PCS[,1], type = "l", col = colourblind[2], lty = "longdash")
lines(seq(0,0.5,0.001), PCS[,2], type = "l", col = colourblind[4], lty = "dotted")
lines(seq(0,0.5,0.001), PCS[,3], type = "l", col = colourblind[8], lty = "dotdash")
lines(c(0,0.5), c(0,0.5), type = "l", col = colourblind[1])
legend(0,0.5,legend = c("No Adaptive", expression(AS~lambda==0.1),
                        expression(RAS~lambda==0.1)),
       lty = c("longdash","dotted","dotdash"), col = colourblind[c(2,4,8)])
###############################Figure 6#########################################
#Figure 6_1, far from LFC case with unequal sample size randomly chosen between 5 to 10.
#unknown and unequal variance chosen between 0.5 and 1.5.
sim_data = function(seed = seed){
  set.seed(seed)
  samsize = sample(5:10,100, replace = TRUE)
  samvar = sample(seq(0.5,1.5,0.1), 100, replace = TRUE)
  #sample size of the 100 populations
  pop = list()
  for (j in 1:10){
    pop = append(pop, list(rnorm(samsize[j], mean = 0, sd = sqrt(samvar[j]))))
  }
  #best populations
  for (j in 11:100){
    pop = append(pop, list(rnorm(samsize[j], mean = -2, sd = sqrt(samvar[j]))))
  }
  #non-best
  return(pop)
}
#1 - PCS for the 3 methods
PCS = sim(data = sim_data, var_type = "uneq_var", method = methods, lambda = 0.1, loc = 1)
#plot
plot(1, type = "n", xlim = c(0,0.5), ylim = c(0,0.5),
     xlab = expression(level~of~significance(alpha)),
     ylab = "1 - PCS of best populations")
lines(seq(0,0.5,0.001), PCS[,1], type = "l", col = colourblind[2], lty = "longdash")
lines(seq(0,0.5,0.001), PCS[,2], type = "l", col = colourblind[4], lty = "dotted")
lines(seq(0,0.5,0.001), PCS[,3], type = "l", col = colourblind[8], lty = "dotdash")
lines(c(0,0.5), c(0,0.5), type = "l", col = colourblind[1])
legend(0,0.5,legend = c("No Adaptive", expression(AS~lambda==0.1),
                        expression(RAS~lambda==0.1)),
       lty = c("longdash","dotted","dotdash"), col = colourblind[c(2,4,8)])
#Figure 6_2, LFC case with unequal sample size randomly chosen between 5 to 10.
#unknown and unequal variance chosen between 0.5 and 1.5.
sim_data = function(seed = seed){
  set.seed(seed)
  samsize = sample(5:10,100, replace = TRUE)
  samvar = sample(seq(0.5,1.5,0.1), 100, replace = TRUE)
  #sample size of the 100 populations
  pop = list()
  for (j in 1:100){
    pop = append(pop, list(rnorm(samsize[j], mean = 0, sd = sqrt(samvar[j]))))
  }
  #best populations
  return(pop)
}
#1 - PCS for the 3 methods
PCS = sim(data = sim_data, var_type = "uneq_var", method = methods, lambda = 0.1, loc = 1)
#plot
plot(1, type = "n", xlim = c(0,0.5), ylim = c(0,0.5),
     xlab = expression(level~of~significance(alpha)),
     ylab = "1 - PCS of best populations")
lines(seq(0,0.5,0.001), PCS[,1], type = "l", col = colourblind[2], lty = "longdash")
lines(seq(0,0.5,0.001), PCS[,2], type = "l", col = colourblind[4], lty = "dotted")
lines(seq(0,0.5,0.001), PCS[,3], type = "l", col = colourblind[8], lty = "dotdash")
lines(c(0,0.5), c(0,0.5), type = "l", col = colourblind[1])
legend(0,0.5,legend = c("No Adaptive", expression(AS~lambda==0.1),
                        expression(RAS~lambda==0.1)),
       lty = c("longdash","dotted","dotdash"), col = colourblind[c(2,4,8)])
