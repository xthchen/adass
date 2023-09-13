library("ACER")
library("cubature")
#If "devtools" not already installed
install.packages("devtools")
#install "haplotest" package
install_github("xthchen/haplotest")
#install "adass" package
install_github("xthchen/adass")
#install "ACER" package
install_github("MartaPelizzola/ACER")

#load package
library("haplotest")
library("adass")
library("stats")
library("ACER")
##########################with replicate scenario###############################
subset_sel = c()
for (i in 1:300){
  set.seed(i^2)
  print(i)
  subset_itr = c()
True_hap <- readRDS("./Evolve & Resequence/True_hap.rds")
True_hap = True_hap[,sample.int(180,10)]
while (!any(rowSums(True_hap)==2)){
  True_hap = True_hap[,sample(1:180,10)]
}

repli = 3

benef = benef_sim(True_hap, n_benef = 1, min = 0.05, max = 0.05, fix_sel = 2, repli = repli)
hap_freq = Frequency_sim(True_hap, t = 6, tdelta = 10, benef, rand_start = FALSE,
                         Ne = 500, repli = repli)
snp_freq = True_hap %*% hap_freq[[1]]
#allele frequency matrix
p_val = adapted.cmh.test(snp_freq,matrix(100000, nrow = nrow(snp_freq), ncol = ncol(snp_freq)),
                           Ne = rep(500,repli), gen = seq(0,60,10), IntGen = TRUE, repl = 1:repli)
#adapted chi square test
mp_val = p.adjust(p_val, method = "BH")
#multiple testing adjustment
mp_pos = which(mp_val < 0.05)
#position of potentially selected SNP
subset_itr = c(subset_itr, length(mp_pos))
subset_itr = c(subset_itr, (benef[[2]][[1]] %in% mp_pos))
snp_freq = snp_freq[mp_pos,]

snp_freq_diff = c()
for (j in 1:repli){
  snp_freq_diff = cbind(snp_freq_diff, (snp_freq[,j*7] - snp_freq[,(j-1)*7+1]))
}
#frequency difference
mean_all = rowSums(snp_freq_diff)/repli
#means
df = (repli-1)*nrow(snp_freq_diff)
var = sum((snp_freq_diff - mean_all)^2)/df
#assuming equal variance across SNPs (not true)
#var = rowSums((snp_freq_diff - mean_all)^2)/(repli-1)
#assuming unequal variance across SNPs

###########non-adaptive method##############
integrand = function(arg){
  u = arg[1]
  x = arg[2]
  pnorm(d*u+x)^(nrow(snp_freq_diff)-1)*dnorm(x)*dchi(u*sqrt(df),df)*sqrt(df)
}


statistic = (max(mean_all) - mean_all)*sqrt(nrow(snp_freq_diff))/sqrt(var)
R = c()
for (d in statistic){
  R = c(R, 1-pcubature(integrand, lowerLimit = c(0,-Inf), upperLimit = c(Inf, Inf))$integral)
}


subset_itr = c(subset_itr, length(R[R>0.05]))
subset_itr = c(subset_itr, (benef[[2]][[1]] %in% mp_pos[R>0.05]))
############adaptive method###############
k_h = AS(R, 0.5)
df_h = (repli-1)*k_h
integrand_h = function(arg){
  u = arg[1]
  x = arg[2]
  pnorm(d*u+x)^(k_h-1)*dnorm(x)*dchi(u*sqrt(df_h),df_h)*sqrt(df_h)
}
R_h = c()
for (d in statistic){
  R_h = c(R_h, 1-pcubature(integrand_h, lowerLimit = c(0,-Inf), upperLimit = c(Inf, Inf))$integral)
}
subset_itr = c(subset_itr, length(R_h[R_h>0.05]))
subset_itr = c(subset_itr, benef[[2]][[1]] %in% mp_pos[R_h>0.05])
############RAS#######################################
k_h = RAS(R, 0.5)
df_h = (repli-1)*k_h
integrand_h = function(arg){
  u = arg[1]
  x = arg[2]
  pnorm(d*u+x)^(k_h-1)*dnorm(x)*dchi(u*sqrt(df_h),df_h)*sqrt(df_h)
}
R_h = c()
for (d in statistic){
  R_h = c(R_h, 1-pcubature(integrand_h, lowerLimit = c(0,-Inf), upperLimit = c(Inf, Inf))$integral)
}
subset_itr = c(subset_itr, length(R_h[R_h>0.05]))
subset_itr = c(subset_itr, benef[[2]][[1]] %in% mp_pos[R_h>0.05])

subset_sel = rbind(subset_sel, subset_itr)
}
#Figure7_1, adapted CMH test case
hist(subset_sel[,1], xlab = "size of subset", main = "adapted CMH test", breaks = 20)
#Figure7_2, non-adaptive subset selection case
hist(subset_sel[,3], xlab = "size of subset", main = "subset selection", breaks = 20)
#Figure7_3, adaptive subset selection case, using AS method
hist(subset_sel[,5], xlab = "size of subset", main = "adaptive subset selection", breaks = 20)
#RAS method omitted due to its similar performance to the AS method.
