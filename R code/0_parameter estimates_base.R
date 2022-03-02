##################################################
##
## FUNCTIONS FOR PLOTS
##
##
library(adegenet)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

IRS_cleaner_f_base = function(data1,k0){
  
  n_t=data1$Ntotalfemalemosq_IRS
  d_t=data1$Ntotaldied_IRS
  fed_t=round(k0 * data1$Nbloodfed_IRS * (1 - data1$Ntotaldied_IRS/data1$Ntotalfemalemosq_IRS),0)
  deterrence_IRS = ifelse(c(data1$Ntotalfemalemosq_C-data1$Ntotalfemalemosq_IRS)<0,0,
                          c(data1$Ntotalfemalemosq_C-data1$Ntotalfemalemosq_IRS))
  deterrence_total = c(data1$Ntotalfemalemosq_IRS+data1$Ntotalfemalemosq_C)
  time=data1$Months_since_IRS*30
  
  return(list(N=nrow(data1),
              n_t=n_t,
              d_t=d_t,
              fed_t=fed_t,
              deterrence_IRS = deterrence_IRS,
              deterrence_total = deterrence_total,
              time=time) ) 
}
##Preparing data using
#Ngufor data.r
#Moore data clean v2.r
#Data cleaning v1.r
#Fitting base and random effect of study.R

## Load the raw data for each chemistry
######################
##
## ACTELLIC 300 CS
##
######################
k0 = 0.699

act_dat = read.csv("data/data_summaryActellic.csv",header=TRUE)

unique(act_dat$Study)
act_dat$rep_of_studies = ifelse(act_dat$Study == "Agossa2014", 1,
                                ifelse(act_dat$Study == "Tchicaya2014", 2,
                                       ifelse(act_dat$Study == "Rowland2013", 3,
                                              ifelse(act_dat$Study == "Unpublished_Muller_Sumi_trial",4,
                                                     ifelse(act_dat$Study == "Rowland2013mud",5,6)))))

act_dat$rep_of_studlab = act_dat$rep_of_studies
act_dat2 = subset(act_dat,act_dat$Study != "UNPUBLISHED_Sarah_Moore")
act_dat3 = subset(act_dat2,act_dat2$Ntotalfemalemosq_IRS > 5)
act_dattest_base=IRS_cleaner_f_base(act_dat3,k0)
act_dat3$timedays = act_dat3$time*30

######################
##
## PYRETHRIODS
##
######################

pyr_dat = read.csv("data/data_summaryPyrethroids.csv",header=TRUE)

pyr_dat$rep_of_studies = ifelse(pyr_dat$Study == "Agossa2014", 1,
                                ifelse(pyr_dat$Study == "UNPUBLISHED_Data_Corbel", 2,
                                       ifelse(pyr_dat$Study == "Rowland 2013", 3,
                                              ifelse(pyr_dat$Study == "Tchicaya2014", 4,
                                                     ifelse(pyr_dat$Study == "Agossa2015PLoSOne", 6,
                                                            ifelse(pyr_dat$Study == "Ngufor_2016", 5,7))))))

pyr_dat$rep_of_studlab = ifelse(pyr_dat$Study == "Agossa2014", 1,
                                ifelse(pyr_dat$Study == "UNPUBLISHED_Data_Corbel", 6,
                                       ifelse(pyr_dat$Study == "Rowland 2013", 3,
                                              ifelse(pyr_dat$Study == "Tchicaya2014", 2,
                                                     ifelse(pyr_dat$Study == "Agossa2015PLoSOne", 7,
                                                            ifelse(pyr_dat$Study == "Ngufor_2016", 8,9))))))

pyr_dat2 = subset(pyr_dat,pyr_dat$rep_of_studies < 7)
pyr_dat3 = subset(pyr_dat2,pyr_dat2$Ntotalfemalemosq_IRS > 5)

pyreth_dattest_base = IRS_cleaner_f_base(pyr_dat3,k0)
pyr_dat3$timedays = pyr_dat3$time*30

######################
##
## BENDIOCARB
##
######################
ben_dat = read.csv("data/data_summaryBendiocarb.csv",header=TRUE)

ben_dat$rep_of_studies = ifelse(ben_dat$Study == "Agossa2014", 2,
                                ifelse(ben_dat$Study == "Akogbeto2010", 4,
                                       ifelse(ben_dat$Study == "UNPUBLISHED_Data_Corbel", 1,
                                              ifelse(ben_dat$Study == "Rowland2013", 3,5))))

ben_dat$rep_of_studlab = ifelse(ben_dat$Study == "Agossa2014", 1,
                                ifelse(ben_dat$Study == "Akogbeto2010", 9,
                                       ifelse(ben_dat$Study == "UNPUBLISHED_Data_Corbel", 6,
                                              ifelse(ben_dat$Study == "Rowland2013", 3,10))))

##sorting here is based on decisions from authors
ben_dat2 = subset(ben_dat,ben_dat$Study == "UNPUBLISHED_Sarah_Moore")
ben_dat2 = subset(ben_dat,ben_dat$Study != "UNPUBLISHED_Sarah_Moore" & ben_dat$Study != "Akogbeto2010")
dim(ben_dat2)
ben_dat2 = subset(ben_dat,ben_dat$rep_of_studies < 3)
ben_dat3 = subset(ben_dat2,ben_dat2$Ntotalfemalemosq_IRS > 5)
ben_dat3$timedays = ben_dat3$Months_since_IRS*30

ben_dattest_base = IRS_cleaner_f_base(ben_dat2,k0)

######################
##
## SUMISHIELD
##
######################
sum_dat = read.csv("data/data_summarySumishield.csv",header=TRUE)
sum_dat$rep_of_studies = ifelse(sum_dat$Study == "UNPUBLISHED_Data_Corbel", 1,
                                ifelse(sum_dat$Study == "UNPUBLISHED_PieMullerCDIvoire", 2,3))

sum_dat$rep_of_studlab = ifelse(sum_dat$Study == "UNPUBLISHED_Data_Corbel", 6,
                                ifelse(sum_dat$Study == "UNPUBLISHED_PieMullerCDIvoire", 4,10))

sum_dat2 = subset(sum_dat,sum_dat$rep_of_studies == 3)
sum_dat2 = subset(sum_dat,sum_dat$rep_of_studies < 3 & sum_dat$concn > 0.3)
dim(sum_dat2)
sum_dattest_base = IRS_cleaner_f_base(sum_dat2,k0)
sum_dat2$timedays = sum_dat2$Months_since_IRS*30



##Base information
#working out quantiles
stan_base <- stan(file="R code\\stan models\\probability_estimates_lq and kq_random_effect_mean.stan", 
                  data=act_dattest_base, 
                  warmup=1000,
                  control = list(adapt_delta = 0.8,
                                 max_treedepth = 20),
                  iter=2000, chains=4)
base <- extract(stan_base)
#Actellic base only
mean(base$alpha1);mean(base$alpha2) ##  2.025098;-0.00889995
mean(base$beta1);mean(base$beta2)   ## -2.465118; 0.00695401
mean(base$omega1);mean(base$omega2) ## -1.232258; -0.002467599

ACTELLIC_ls_theta = sample(base$alpha1,50,replace=FALSE)
ACTELLIC_ls_gamma = sample(base$alpha2,50,replace=FALSE)
ACTELLIC_ks_theta = sample(base$beta1,50,replace=FALSE)
ACTELLIC_ks_gamma = sample(base$beta2,50,replace=FALSE)
ACTELLIC_ms_theta = sample(base$omega1,50,replace=FALSE)
ACTELLIC_ms_gamma = sample(base$omega2,50,replace=FALSE)


stan_base2 <- stan(file="R code\\stan models\\probability_estimates_lq and kq_random_effect_mean.stan", 
                   data=ben_dattest_base, 
                   warmup=1000,
                   control = list(adapt_delta = 0.8,
                                  max_treedepth = 20),
                   iter=2000, chains=4)
base2 <- extract(stan_base2)
#Actellic base only
mean(base2$alpha1);mean(base2$alpha2) ##  1.494069;-1.368968
mean(base2$beta1);mean(base2$beta2)   ## -1.89659; 1.503587
mean(base2$omega1);mean(base2$omega2) ## -1.701177; 1.464068

BEND_ls_theta = sample(base2$alpha1,50,replace=FALSE)
BEND_ls_gamma = sample(base2$alpha2,50,replace=FALSE)
BEND_ks_theta = sample(base2$beta1,50,replace=FALSE)
BEND_ks_gamma = sample(base2$beta2,50,replace=FALSE)
BEND_ms_theta = sample(base2$omega1,50,replace=FALSE)
BEND_ms_gamma = sample(base2$omega2*0,50,replace=FALSE)

stan_base3 <- stan(file="R code\\stan models\\probability_estimates_lq and kq_random_effect_mean.stan", 
                   data=sum_dattest_base, 
                   warmup=1000,
                   control = list(adapt_delta = 0.8,
                                  max_treedepth = 20),
                   iter=2000, chains=4)
base3 <- extract(stan_base3)
#Actellic base only
mean(base3$alpha1);mean(base3$alpha2) ##  1.494069;-1.368968
mean(base3$beta1);mean(base3$beta2)   ## -1.89659; 1.503587
mean(base3$omega1);mean(base3$omega2) ## -1.701177; 1.464068

SUM_ls_theta = sample(base3$alpha1,50,replace=FALSE)
SUM_ls_gamma = sample(base3$alpha2,50,replace=FALSE)
SUM_ks_theta = sample(base3$beta1,50,replace=FALSE)
SUM_ks_gamma = sample(base3$beta2,50,replace=FALSE)
SUM_ms_theta = sample(base3$omega1,50,replace=FALSE)
SUM_ms_gamma = sample(base3$omega2,50,replace=FALSE)


stan_base4 <- stan(file="R code\\stan models\\probability_estimates_lq and kq_random_effect_mean.stan", 
                   data=pyreth_dattest_base, 
                   warmup=1000,
                   control = list(adapt_delta = 0.8,
                                  max_treedepth = 20),
                   iter=2000, chains=4)
base4 <- extract(stan_base4)
#Actellic base only
mean(base4$alpha1);mean(base4$alpha2) ##  1.494069;-1.368968
mean(base4$beta1);mean(base4$beta2)   ## -1.89659; 1.503587
mean(base4$omega1);mean(base4$omega2) ## -1.701177; 1.464068

PYR_ls_theta = sample(base4$alpha1,50,replace=FALSE)
PYR_ls_gamma = sample(base4$alpha2,50,replace=FALSE)
PYR_ks_theta = sample(base4$beta1,50,replace=FALSE)
PYR_ks_gamma = sample(base4$beta2,50,replace=FALSE)
PYR_ms_theta = sample(base4$omega1,50,replace=FALSE)
PYR_ms_gamma = sample(base4$omega2,50,replace=FALSE)


Samples = data.frame(ACTELLIC_ls_theta,ACTELLIC_ls_gamma,ACTELLIC_ks_theta,ACTELLIC_ks_gamma,ACTELLIC_ms_theta,ACTELLIC_ms_gamma,
                     BEND_ls_theta,BEND_ls_gamma,BEND_ks_theta,BEND_ks_gamma,BEND_ms_theta,BEND_ms_gamma,
                     SUM_ls_theta,SUM_ls_gamma,SUM_ks_theta,SUM_ks_gamma,SUM_ms_theta,SUM_ms_gamma,
                     PYR_ls_theta,PYR_ls_gamma,PYR_ks_theta,PYR_ks_gamma,PYR_ms_theta,PYR_ms_gamma)

