library(adegenet)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

IRS_cleaner_f_base = function(data1,k0){
  
  n_t=data1$Ntotalfemalemosq_IRS
  d_t=data1$Ntotaldied_IRS
  fed_t=round(data1$Nbloodfed_IRS * (1 - data1$Ntotaldied_IRS/data1$Ntotalfemalemosq_IRS),0)
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
k0 = 0.699 ##the prob of feeding in absence of interventions

act_dat = read.csv("data/data_summaryActellic.csv",header=TRUE)

unique(act_dat$Study)
act_dat$rep_of_studies = ifelse(act_dat$Study == "Agossa2014", 1,
                                ifelse(act_dat$Study == "Tchicaya2014", 2,
                                       ifelse(act_dat$Study == "Rowland2013", 3,
                                              ifelse(act_dat$Study == "Unpublished_Muller_Sumi_trial",4,
                                                     ifelse(act_dat$Study == "Rowland2013mud",5,
                                                            ifelse(act_dat$Study == "Oxborough2014MalJ",6,7))))))
act_dat$rep_of_studlab = ifelse(act_dat$Study == "Agossa2014", 1,
                                ifelse(act_dat$Study == "Tchicaya2014", 5,
                                       ifelse(act_dat$Study == "Rowland2013", 4,
                                              ifelse(act_dat$Study == "Unpublished_Muller_Sumi_trial",12,
                                                     ifelse(act_dat$Study == "Rowland2013mud",13,
                                                            ifelse(act_dat$Study == "Oxborough2014MalJ",9,20))))))

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
unique(pyr_dat$Study)
pyr_dat$rep_of_studies = ifelse(pyr_dat$Study == "Agossa2014", 1,
                                ifelse(pyr_dat$Study == "UNPUBLISHED_Data_Corbel", 2,
                                       ifelse(pyr_dat$Study == "Rowland 2013", 3,
                                              ifelse(pyr_dat$Study == "Tchicaya2014", 4,
                                                     ifelse(pyr_dat$Study == "Agossa2015PLoSOne" & pyr_dat$GROUP == 3, 5,
                                                            ifelse(pyr_dat$Study == "Agossa2015PLoSOne" & pyr_dat$GROUP == 6, 6,
                                                                   ifelse(pyr_dat$Study == "Agossa2015PLoSOne" & pyr_dat$GROUP == 7, 7,
                                                                          ifelse(pyr_dat$Study == "Ngufor_2016", 8,
                                                                                 ifelse(pyr_dat$Study == "Ngufor_2017_PLoSOne" & pyr_dat$GROUP == 16, 9,
                                                                                        ifelse(pyr_dat$Study == "Ngufor_2017_PLoSOne" & pyr_dat$GROUP == 17, 10,
                                                                                               ifelse(pyr_dat$Study == "Ngufor_2017_PLoSOne" & pyr_dat$GROUP == 18, 11,
                                                                                                      ifelse(pyr_dat$Study == "Akogbeto2010" & pyr_dat$GROUP == 23,12,
                                                                                                             ifelse(pyr_dat$Study == "Akogbeto2010" & pyr_dat$GROUP == 24,13,
                                                                                                                    ifelse(pyr_dat$Study == "Agossa_2018ParaVec",14,20))))))))))))))

pyr_dat$rep_of_studlab = ifelse(pyr_dat$Study == "Agossa2014", 1,
                                ifelse(pyr_dat$Study == "UNPUBLISHED_Data_Corbel", 11,
                                       ifelse(pyr_dat$Study == "Rowland 2013", 4,
                                              ifelse(pyr_dat$Study == "Tchicaya2014", 5,
                                                     ifelse(pyr_dat$Study == "Agossa2015PLoSOne" & pyr_dat$GROUP == 3, 2,#alpha
                                                            ifelse(pyr_dat$Study == "Agossa2015PLoSOne" & pyr_dat$GROUP == 6, 14,#delta
                                                                   ifelse(pyr_dat$Study == "Agossa2015PLoSOne" & pyr_dat$GROUP == 7, 15,#lambda
                                                                          ifelse(pyr_dat$Study == "Ngufor_2016", 3,
                                                                                 ifelse(pyr_dat$Study == "Ngufor_2017_PLoSOne" & pyr_dat$GROUP == 16, 7, 
                                                                                        ifelse(pyr_dat$Study == "Ngufor_2017_PLoSOne" & pyr_dat$GROUP == 17, 17, 
                                                                                               ifelse(pyr_dat$Study == "Ngufor_2017_PLoSOne" & pyr_dat$GROUP == 18, 18, 
                                                                                                      ifelse(pyr_dat$Study == "Akogbeto2010" & pyr_dat$GROUP == 23,6,#delta
                                                                                                             ifelse(pyr_dat$Study == "Akogbeto2010" & pyr_dat$GROUP == 24,16,#alpha
                                                                                                                    ifelse(pyr_dat$Study == "Agossa_2018ParaVec",10,20))))))))))))))



pyr_dat2 = subset(pyr_dat,pyr_dat$rep_of_studies < 18)
pyr_dat3 = subset(pyr_dat2,pyr_dat2$Ntotalfemalemosq_IRS > 5)

pyreth_dattest_base = IRS_cleaner_f_base(pyr_dat3,k0)
pyr_dat3$timedays = pyr_dat3$time*30

######################
##
## BENDIOCARB
##
######################
ben_dat = read.csv("data/data_summaryBendiocarb.csv",header=TRUE)
unique(ben_dat$Study)
ben_dat$rep_of_studies = ifelse(ben_dat$Study == "Agossa2014", 1,
                                ifelse(ben_dat$Study == "Akogbeto2010", 2,
                                       ifelse(ben_dat$Study == "UNPUBLISHED_Data_Corbel", 3,
                                              ifelse(ben_dat$Study == "Djenontin2010",4,20))))

ben_dat$rep_of_studlab = ifelse(ben_dat$Study == "Agossa2014", 1,
                                ifelse(ben_dat$Study == "Akogbeto2010", 6,
                                       ifelse(ben_dat$Study == "UNPUBLISHED_Data_Corbel", 11,
                                              ifelse(ben_dat$Study == "Djenontin2010", 8,20))))


#ben_dat2 = subset(ben_dat,ben_dat$Study != "UNPUBLISHED_Sarah_Moore")
#ben_dat2 = subset(ben_dat,ben_dat$Study != "UNPUBLISHED_Sarah_Moore" & ben_dat$Study != "Akogbeto2010")
ben_dat2 = subset(ben_dat,ben_dat$rep_of_studies < 5)
ben_dat3 = subset(ben_dat2,ben_dat2$Ntotalfemalemosq_IRS > 5)
ben_dat3$timedays = ben_dat3$Months_since_IRS*30
dim(ben_dat2)
dim(ben_dat3)

#ben_dat2$timedays = ben_dat2$Months_since_IRS*30

ben_dattest_base = IRS_cleaner_f_base(ben_dat3,k0)

######################
##
## SUMISHIELD
##
######################
sum_dat = read.csv("data/data_summarySumishield.csv",header=TRUE)
dim(sum_dat)
unique(sum_dat$Study)
sum_dat$rep_of_studies = ifelse(sum_dat$Study == "UNPUBLISHED_Data_Corbel", 1,
                                ifelse(sum_dat$Study == "UNPUBLISHED_PieMullerCDIvoire", 2,
                                       ifelse(sum_dat$Study == "Ngufor_2017_PLoSOne",3,
                                              ifelse(sum_dat$Study == "Agossa_2018ParaVec",4,20))))

sum_dat$rep_of_studlab = ifelse(sum_dat$Study == "UNPUBLISHED_Data_Corbel", 11,
                                ifelse(sum_dat$Study == "UNPUBLISHED_PieMullerCDIvoire", 12,
                                       ifelse(sum_dat$Study == "Ngufor_2017_PLoSOne",7,
                                              ifelse(sum_dat$Study == "Agossa_2018ParaVec",10,20))))

summary(sum_dat)
sum_dat2 = subset(sum_dat,sum_dat$rep_of_studies == 3)
sum_dat2 = subset(sum_dat,sum_dat$rep_of_studies < 5 & sum_dat$concn > 0.1)
dim(sum_dat2)
sum_dattest_base = IRS_cleaner_f_base(sum_dat2,k0)
sum_dat2$timedays = sum_dat2$Months_since_IRS*30

act_dat3 = act_dat3[order(act_dat3$rep_of_studies),]
pyr_dat3 = pyr_dat3[order(pyr_dat3$rep_of_studies),]
ben_dat3 = ben_dat3[order(ben_dat3$rep_of_studies),]
sum_dat3 = sum_dat2[order(sum_dat2$rep_of_studies),]


IRS_cleaner_f2 = function(data1,spray,no_of_studies,rep_of_studies){
  
  
  n_t=data1$Ntotalfemalemosq_IRS
  d_t=data1$Ntotaldied_IRS
  fed_t=round(data1$Nbloodfed_IRS * (1 - data1$Ntotaldied_IRS/data1$Ntotalfemalemosq_IRS),0)
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
              time=time,
              N_IRS = no_of_studies,
              IRS = rep_of_studies) ) 
}

act_dattest2=IRS_cleaner_f2(act_dat3,"Actellic",6,act_dat3$rep_of_studies)
pyr_dattest2=IRS_cleaner_f2(pyr_dat3,"pyrethroid",14,pyr_dat3$rep_of_studies)
ben_dattest2=IRS_cleaner_f2(ben_dat3,"Bendiocarb",4,ben_dat3$rep_of_studies)
sum_dattest2=IRS_cleaner_f2(sum_dat3,"Sumishield",4,sum_dat2$rep_of_studies)

###########################
##
## 1_estimating the uncertainty 
## around IRS impacts

## Take the minimum estimate for each time
## point through the year

act_dattest2$mortality = act_dattest2$d_t/act_dattest2$n_t 

unique(act_dattest2$time)

minimum_mort = minimum_d_t = minimum_n_t = minimum_fed_t = times =
  minimum_deterrence_IRS = minimum_deterrence_total = minimum_time = 
  maximum_mort = maximum_d_t = maximum_n_t = maximum_fed_t = 
  maximum_deterrence_IRS = maximum_deterrence_total = maximum_time  = numeric(length(unique(act_dattest2$time))) 
for(i in 1:length(unique(act_dattest2$time))){
  minimum_mort[i] = min(act_dattest2$mortality[act_dattest2$time == unique(act_dattest2$time)[i]])
  maximum_mort[i] = max(act_dattest2$mortality[act_dattest2$time == unique(act_dattest2$time)[i]])

  
  }
for(i in 1:length(unique(act_dattest2$time))){
  act_dattest2$minimum_d_t             = ifelse(act_dattest2$mortality == minimum_mort[i],act_dattest2$d_t,99999)
  act_dattest2$minimum_n_t             = ifelse(act_dattest2$mortality == minimum_mort[i],act_dattest2$n_t,99999)
  act_dattest2$minimum_fed_t           = ifelse(act_dattest2$mortality == minimum_mort[i],act_dattest2$fed_t,99999)
  act_dattest2$minimum_time = ifelse(act_dattest2$mortality == minimum_mort[i],act_dattest2$time,99999)
  minimum_d_t[i]              = act_dattest2$minimum_d_t[act_dattest2$minimum_d_t < 9999]
  minimum_n_t[i]              = act_dattest2$minimum_n_t[act_dattest2$minimum_n_t < 9999]
  minimum_fed_t[i]            = act_dattest2$minimum_fed_t[act_dattest2$minimum_fed_t < 9999]
  times[i]                    = act_dattest2$minimum_time[act_dattest2$minimum_time < 9999]
}
for(i in 1:length(unique(act_dattest2$time))){ 
  act_dattest2$minimum_deterrence_IRS  = ifelse(act_dattest2$mortality == minimum_mort[i],act_dattest2$deterrence_IRS,99999)
  act_dattest2$minimum_deterrence_total= ifelse(act_dattest2$mortality == minimum_mort[i],act_dattest2$deterrence_total,99999)

  minimum_deterrence_IRS[i]   = act_dattest2$minimum_deterrence_IRS[act_dattest2$minimum_deterrence_IRS < 99999]
  minimum_deterrence_total[i] = act_dattest2$minimum_deterrence_total[act_dattest2$minimum_deterrence_total < 99999]
}
for(i in c(3,5:14)){
  act_dattest2$maximum_d_t             = ifelse(act_dattest2$mortality == maximum_mort[i],act_dattest2$d_t,99999)
  act_dattest2$maximum_n_t             = ifelse(act_dattest2$mortality == maximum_mort[i],act_dattest2$n_t,99999)
  act_dattest2$maximum_fed_t           = ifelse(act_dattest2$mortality == maximum_mort[i],act_dattest2$fed_t,99999)
  act_dattest2$maximum_time = ifelse(act_dattest2$mortality == maximum_mort[i],act_dattest2$time,99999)
  act_dattest2$maximum_deterrence_IRS  = ifelse(act_dattest2$mortality == maximum_mort[i],act_dattest2$deterrence_IRS,99999)
  act_dattest2$maximum_deterrence_total= ifelse(act_dattest2$mortality == maximum_mort[i],act_dattest2$deterrence_total,99999)
  maximum_d_t[i]              = act_dattest2$maximum_d_t[act_dattest2$maximum_d_t < 9999]
  maximum_n_t[i]              = act_dattest2$maximum_n_t[act_dattest2$maximum_n_t < 9999]
  maximum_fed_t[i]            = act_dattest2$maximum_fed_t[act_dattest2$maximum_fed_t < 9999]
  times[i]                    = act_dattest2$maximum_time[act_dattest2$maximum_time < 9999]
  maximum_deterrence_IRS[i]   = act_dattest2$maximum_deterrence_IRS[act_dattest2$maximum_deterrence_IRS < 99999]
  maximum_deterrence_total[i] = act_dattest2$maximum_deterrence_total[act_dattest2$maximum_deterrence_total < 99999]
  
}
for(i in c(1,2,4)){
  maximum_d_t[i]              = max(act_dattest2$d_t[which(act_dattest2$mortality == 1 & act_dattest2$time == unique(act_dattest2$time)[i])] )
  maximum_n_t[i]              = max(act_dattest2$n_t[which(act_dattest2$mortality == 1 & act_dattest2$time == unique(act_dattest2$time)[i])] )
  maximum_fed_t[i]            = max(act_dattest2$fed_t[which(act_dattest2$mortality == 1 & act_dattest2$time == unique(act_dattest2$time)[i])] )
  times[i]                    = max(act_dattest2$time[which(act_dattest2$mortality == 1 & act_dattest2$time == unique(act_dattest2$time)[i])] )
  maximum_deterrence_IRS[i]   = max(act_dattest2$deterrence_IRS[which(act_dattest2$mortality == 1 & act_dattest2$time == unique(act_dattest2$time)[i])] )
  maximum_deterrence_total[i] = max(act_dattest2$deterrence_total[which(act_dattest2$mortality == 1 & act_dattest2$time == unique(act_dattest2$time)[i])] )
}

## Select out these data and fit to these

dattest2_b=list(N=2*length(unique(act_dattest2$time)),
                n_t=c(minimum_n_t,maximum_n_t),
                d_t=c(minimum_d_t,maximum_d_t),
                fed_t=c(minimum_fed_t,maximum_fed_t),
                deterrence_IRS = c(minimum_deterrence_IRS,maximum_deterrence_IRS),
                deterrence_total = c(minimum_deterrence_total,maximum_deterrence_total),
                time=rep(unique(act_dattest2$time),2),
                N_IRS = 2,
                IRS = c(rep(1,length(unique(act_dattest2$time))),
                        rep(2,length(unique(act_dattest2$time)))) )


stan2_b <- stan(file="R code\\stan models\\random_effects_mean_fits_model.stan", 
              data=dattest2_b, 
              warmup=500,
              control = list(adapt_delta = 0.9,
                             max_treedepth = 20),
              iter=1000, chains=4)

## Useful to confirm model fitting well (diagnostics)
#library(shinystan)
#launch_shinystan(stan2_b)

base <- extract(stan2_b)

##Actellic minimum
mean(base$alpha1[,1]);mean(base$alpha2[,1]) ##  1.8322;-0.01385072
mean(base$beta1[,1]);mean(base$beta2[,1]) ##  -2.166622;0.01365688
mean(base$omega1[,1]);mean(base$omega2[,1]) ##  -4.471159;-0.01385072 ## using same deprecation as mortality


##Actellic maximum
mean(base$alpha1[,2]);mean(base$alpha2[,2]) ##  4.753691;-0.009835282
mean(base$beta1[,2]);mean(base$beta2[,2]) ##  -4.359774; 0.006331158
mean(base$omega1[,2]);mean(base$omega2[,2]) ## -1.219416;-0.001093292

DEAD = dattest2_b$d_t/dattest2_b$n_t
FED = dattest2_b$fed_t/dattest2_b$n_t

time = 1:365

plot(DEAD ~ dattest2_b$time,ylab="Mortlality range",ylim=c(0,1),col="gold4",pch=rep(c(15,17),each=14),
     main = "Actellic",cex.main=1.2,xlim=c(1,365),xaxt="n",
     xlab="",yaxt="n",cex.lab=1,cex.axis=1,cex=1.4)
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20),cex.lab=1.3,cex.axis=1)

mean_valssp_checker4Min = 1 / (1 + exp(-mean(base$alpha1[,1]) - mean(base$alpha2[,1])*time))
mean_valssp_checker4Max = 1 / (1 + exp(-mean(base$alpha1[,2]) - mean(base$alpha2[,2])*time))

polygon(c(time,rev(time)),c(mean_valssp_checker4Min,rev(mean_valssp_checker4Max)),col=transp("gold4",0.4),border=NA)

time=seq(1,365,length=365)

mean_valssp_checker4Min = 1 / (1 + exp(-1.832483 - -0.01385783*time))
mean_valssp_checker4Max = 1 / (1 + exp(-4.7590289 - -0.009848811*time))


plot(FED ~ dattest2_b$time,ylab="Mortlality range",ylim=c(0,1),col="gold4",pch=rep(c(15,17),each=14),
     main = "Actellic",cex.main=1.2,xlim=c(1,365),xaxt="n",
     xlab="",yaxt="n",cex.lab=1,cex.axis=1,cex=1.4)
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20),cex.lab=1.3,cex.axis=1)

mean_valsfp_checker4Min = 1 / (1 + exp(-mean(base$beta1[,1]) - mean(base$beta2[,1])*time))
mean_valsfp_checker4Max = 1 / (1 + exp(-mean(base$beta1[,2]) - mean(base$beta2[,2])*time))

polygon(c(time,rev(time)),c(mean_valsfp_checker4Min,rev(mean_valsfp_checker4Max)),col=transp("gold4",0.4),border=NA)

time=seq(1,365,length=365)


DET = dattest2_b$deterrence_IRS/dattest2_b$deterrence_total

mean_valsdet_checker4Min = 1 / (1 + exp(-mean(base$omega1[,1]) - mean(base$alpha2[,1])*time))
mean_valsdet_checker4Max = 1 / (1 + exp(-mean(base$omega1[,2]) - mean(base$alpha2[,2])*time))

plot(DET ~ dattest2_b$time,ylim=c(0,1),col="gold4",pch=rep(c(15,17),each=14),
     xlab="",xaxt="n",yaxt="n",cex.lab=1,cex.axis=1,xaxt="n",cex=1.4)
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20),cex.lab=1.3,cex.axis=1)
axis(1,at=seq(0,365,120),labels=seq(0,365,120),cex.lab=1.3,cex.axis=1)
polygon(c(time,rev(time)),c(mean_valsdet_checker4Min,rev(mean_valsdet_checker4Max)),col=transp("gold4",0.4),border=NA)


###########################
##
## 1_estimating the uncertainty 
## around IRS impacts

## Take the minimum estimate for each time
## point through the year

pyr_dattest2$mortality = pyr_dattest2$d_t/pyr_dattest2$n_t 

unique(pyr_dattest2$time)

minimum_mort = minimum_d_t = minimum_n_t = minimum_fed_t = times = 
  minimum_deterrence_IRS = minimum_deterrence_total = minimum_time = 
  maximum_mort = maximum_d_t = maximum_n_t = maximum_fed_t = 
  maximum_deterrence_IRS = maximum_deterrence_total = maximum_time  = pyr_times = numeric(length(unique(pyr_dattest2$time))) 
for(i in 1:length(unique(pyr_dattest2$time))){
  minimum_mort[i] = min(pyr_dattest2$mortality[pyr_dattest2$time == unique(pyr_dattest2$time)[i]])
  maximum_mort[i] = max(pyr_dattest2$mortality[pyr_dattest2$time == unique(pyr_dattest2$time)[i]])
  
  
}
for(i in 1:length(unique(pyr_dattest2$time))){
  pyr_dattest2$minimum_d_t             = ifelse(pyr_dattest2$mortality == minimum_mort[i],pyr_dattest2$d_t,99999)
  pyr_dattest2$minimum_n_t             = ifelse(pyr_dattest2$mortality == minimum_mort[i],pyr_dattest2$n_t,99999)
  pyr_dattest2$minimum_fed_t           = ifelse(pyr_dattest2$mortality == minimum_mort[i],pyr_dattest2$fed_t,99999)
  pyr_dattest2$minimum_time            = ifelse(pyr_dattest2$mortality == minimum_mort[i],pyr_dattest2$time,99999)
  minimum_d_t[i]              = max(pyr_dattest2$minimum_d_t[pyr_dattest2$minimum_d_t < 9999])
  minimum_n_t[i]              = max(pyr_dattest2$minimum_n_t[pyr_dattest2$minimum_n_t < 9999])
  minimum_fed_t[i]            = max(pyr_dattest2$minimum_fed_t[pyr_dattest2$minimum_fed_t < 9999])
  pyr_times[i]                = max(pyr_dattest2$minimum_time[pyr_dattest2$minimum_time < 9999])
}
for(i in 1:length(unique(pyr_dattest2$time))){ 
  pyr_dattest2$minimum_deterrence_IRS  = ifelse(pyr_dattest2$mortality == minimum_mort[i],pyr_dattest2$deterrence_IRS,99999)
  pyr_dattest2$minimum_deterrence_total= ifelse(pyr_dattest2$mortality == minimum_mort[i],pyr_dattest2$deterrence_total,99999)
  
  minimum_deterrence_IRS[i]   = max(pyr_dattest2$minimum_deterrence_IRS[pyr_dattest2$minimum_deterrence_IRS < 99999])
  minimum_deterrence_total[i] = max(pyr_dattest2$minimum_deterrence_total[pyr_dattest2$minimum_deterrence_total < 99999])
}
for(i in c(1:length(unique(pyr_dattest2$time)))){
  pyr_dattest2$maximum_d_t             = ifelse(pyr_dattest2$mortality == maximum_mort[i],pyr_dattest2$d_t,99999)
  pyr_dattest2$maximum_n_t             = ifelse(pyr_dattest2$mortality == maximum_mort[i],pyr_dattest2$n_t,99999)
  pyr_dattest2$maximum_fed_t           = ifelse(pyr_dattest2$mortality == maximum_mort[i],pyr_dattest2$fed_t,99999)
  pyr_dattest2$maximum_time = ifelse(pyr_dattest2$mortality == maximum_mort[i],pyr_dattest2$time,99999)
  pyr_dattest2$maximum_deterrence_IRS  = ifelse(pyr_dattest2$mortality == maximum_mort[i],pyr_dattest2$deterrence_IRS,99999)
  pyr_dattest2$maximum_deterrence_total= ifelse(pyr_dattest2$mortality == maximum_mort[i],pyr_dattest2$deterrence_total,99999)
  maximum_d_t[i]              = max(pyr_dattest2$maximum_d_t[pyr_dattest2$maximum_d_t < 9999])
  maximum_n_t[i]              = max(pyr_dattest2$maximum_n_t[pyr_dattest2$maximum_n_t < 9999])
  maximum_fed_t[i]            = max(pyr_dattest2$maximum_fed_t[pyr_dattest2$maximum_fed_t < 9999])
  times[i]                    = max(pyr_dattest2$maximum_time[pyr_dattest2$maximum_time < 9999])
  maximum_deterrence_IRS[i]   = max(pyr_dattest2$maximum_deterrence_IRS[pyr_dattest2$maximum_deterrence_IRS < 99999])
  maximum_deterrence_total[i] = max(pyr_dattest2$maximum_deterrence_total[pyr_dattest2$maximum_deterrence_total < 99999])
  
}

## Select out these data and fit to these

dattest2_b=list(N=2*length(unique(pyr_dattest2$time)),
                n_t=c(minimum_n_t,maximum_n_t),
                d_t=c(minimum_d_t,maximum_d_t),
                fed_t=c(minimum_fed_t,maximum_fed_t),
                deterrence_IRS = c(minimum_deterrence_IRS,maximum_deterrence_IRS),
                deterrence_total = c(minimum_deterrence_total,maximum_deterrence_total),
                time=rep(unique(pyr_dattest2$time),2),
                N_IRS = 2,
                IRS = c(rep(1,length(unique(pyr_dattest2$time))),
                        rep(2,length(unique(pyr_dattest2$time)))) )


stan2_b <- stan(file="R code\\stan models\\random_effects_mean_fits_model.stan", 
                data=dattest2_b, 
                warmup=500,
                control = list(adapt_delta = 0.9,
                               max_treedepth = 20),
                iter=1000, chains=4)

#library(shinystan)
#launch_shinystan(stan2_b)

base <- extract(stan2_b)

##pyrethroid minimum
mean(base$alpha1[,1]);mean(base$alpha2[,1]) ##  -2.754575;-0.000681545
mean(base$beta1[,1]);mean(base$beta2[,1]) ##   0.4413352;0.0003144774
mean(base$omega1[,1]);mean(base$omega2[,1]) ## -2.451382; 0.001873774


##pyrethroid maximum
mean(base$alpha1[,2]);mean(base$alpha2[,2]) ##  0.8431575;-0.009555796
mean(base$beta1[,2]);mean(base$beta2[,2]) ## -1.405392; 0.005824437
mean(base$omega1[,2]);mean(base$omega2[,2]) ## -0.636807;-0.003157326

DEAD = dattest2_b$d_t/dattest2_b$n_t
FED = dattest2_b$fed_t/dattest2_b$n_t

plot(DEAD ~ dattest2_b$time,ylab="Mortlality range",ylim=c(0,1),col="blue",pch=rep(c(15,17),each=42),
     main = "bendiocarb",cex.main=1.2,xlim=c(1,365),xaxt="n",
     xlab="",yaxt="n",cex.lab=1,cex.axis=1,cex=1.4)
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20),cex.lab=1.3,cex.axis=1)

mean_valssp_checker4Min = 1 / (1 + exp(-mean(base$alpha1[,1]) - mean(base$alpha2[,1])*time))
mean_valssp_checker4Max = 1 / (1 + exp(-mean(base$alpha1[,2]) - mean(base$alpha2[,2])*time))

polygon(c(time,rev(time)),c(mean_valssp_checker4Min,rev(mean_valssp_checker4Max)),col=transp("gold4",0.4),border=NA)

time=seq(1,365,length=365)

mean_valssp_checker4Min = 1 / (1 + exp(--2.754575 --0.000681545*time))
mean_valssp_checker4Max = 1 / (1 + exp(-0.8431575 - -0.009555796*time))
lines(mean_valssp_checker4Min ~ time)

plot(FED ~ dattest2_b$time,ylab="Mortlality range",ylim=c(0,1),col="gold4",pch=rep(c(15,17),each=14),
     main = "Actellic",cex.main=1.2,xlim=c(1,365),xaxt="n",
     xlab="",yaxt="n",cex.lab=1,cex.axis=1,cex=1.4)
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20),cex.lab=1.3,cex.axis=1)

mean_valsfp_checker4Min = 1 / (1 + exp(-mean(base$beta1[,1]) - mean(base$beta2[,2])*time))
mean_valsfp_checker4Max = 1 / (1 + exp(-mean(base$beta1[,2]) - mean(base$beta2[,1])*time))

polygon(c(time,rev(time)),c(mean_valsfp_checker4Min,rev(mean_valsfp_checker4Max)),col=transp("gold4",0.4),border=NA)

time=seq(1,365,length=365)

mean_valsfp_checker4Min = 1 / (1 + exp(-0.4413352 - 0.005824437*time))
mean_valsfp_checker4Max = 1 / (1 + exp(- -1.405392 - 0.0003144774*time))
lines(mean_valsfp_checker4Min ~ time)

DET = dattest2_b$deterrence_IRS/dattest2_b$deterrence_total

mean_valsdet_checker4Min = 1 / (1 + exp(-mean(base$omega1[,1]) - mean(base$alpha2[,2])*time))
mean_valsdet_checker4Max = 1 / (1 + exp(-mean(base$omega1[,2]) - mean(base$alpha2[,2])*time))

plot(DET ~ dattest2_b$time,ylim=c(0,1),col="gold4",pch=rep(c(15,17),each=14),
     xlab="",xaxt="n",yaxt="n",cex.lab=1,cex.axis=1,xaxt="n",cex=1.4)
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20),cex.lab=1.3,cex.axis=1)
axis(1,at=seq(0,365,120),labels=seq(0,365,120),cex.lab=1.3,cex.axis=1)
polygon(c(time,rev(time)),c(mean_valsdet_checker4Min,rev(mean_valsdet_checker4Max)),col=transp("gold4",0.4),border=NA)

mean_valsdet_checker4Min = 1 / (1 + exp(--2.451382 - -0.009555796*time))
mean_valsdet_checker4Max = 1 / (1 + exp(--0.636807- -0.009555796*time))



###########################
##
## 1_estimating the uncertainty 
## around IRS impacts

## Take the minimum estimate for each time
## point through the year

ben_dattest2$mortality = ben_dattest2$d_t/ben_dattest2$n_t 

unique(ben_dattest2$time)

minimum_mort = minimum_d_t = minimum_n_t = minimum_fed_t = 
  minimum_deterrence_IRS = minimum_deterrence_total = minimum_time = 
  maximum_mort = maximum_d_t = maximum_n_t = maximum_fed_t = 
  maximum_deterrence_IRS = maximum_deterrence_total = maximum_time  = ben_times = numeric(length(unique(ben_dattest2$time))) 
for(i in 1:length(unique(ben_dattest2$time))){
  minimum_mort[i] = min(ben_dattest2$mortality[ben_dattest2$time == unique(ben_dattest2$time)[i]])
  maximum_mort[i] = max(ben_dattest2$mortality[ben_dattest2$time == unique(ben_dattest2$time)[i]])
  
  
}
for(i in 1:length(unique(ben_dattest2$time))){
  ben_dattest2$minimum_d_t             = ifelse(ben_dattest2$mortality == minimum_mort[i],ben_dattest2$d_t,99999)
  ben_dattest2$minimum_n_t             = ifelse(ben_dattest2$mortality == minimum_mort[i],ben_dattest2$n_t,99999)
  ben_dattest2$minimum_fed_t           = ifelse(ben_dattest2$mortality == minimum_mort[i],ben_dattest2$fed_t,99999)
  ben_dattest2$minimum_time            = ifelse(ben_dattest2$mortality == minimum_mort[i],ben_dattest2$time,99999)
  minimum_d_t[i]              = max(ben_dattest2$minimum_d_t[ben_dattest2$minimum_d_t < 9999])
  minimum_n_t[i]              = max(ben_dattest2$minimum_n_t[ben_dattest2$minimum_n_t < 9999])
  minimum_fed_t[i]            = max(ben_dattest2$minimum_fed_t[ben_dattest2$minimum_fed_t < 9999])
  ben_times[i]                = max(ben_dattest2$minimum_time[ben_dattest2$minimum_time < 9999])
}
for(i in 1:length(unique(ben_dattest2$time))){ 
  ben_dattest2$minimum_deterrence_IRS  = ifelse(ben_dattest2$mortality == minimum_mort[i],ben_dattest2$deterrence_IRS,99999)
  ben_dattest2$minimum_deterrence_total= ifelse(ben_dattest2$mortality == minimum_mort[i],ben_dattest2$deterrence_total,99999)
  
  minimum_deterrence_IRS[i]   = max(ben_dattest2$minimum_deterrence_IRS[ben_dattest2$minimum_deterrence_IRS < 99999])
  minimum_deterrence_total[i] = max(ben_dattest2$minimum_deterrence_total[ben_dattest2$minimum_deterrence_total < 99999])
}
for(i in c(1:length(unique(ben_dattest2$time)))){
  ben_dattest2$maximum_d_t             = ifelse(ben_dattest2$mortality == maximum_mort[i],ben_dattest2$d_t,99999)
  ben_dattest2$maximum_n_t             = ifelse(ben_dattest2$mortality == maximum_mort[i],ben_dattest2$n_t,99999)
  ben_dattest2$maximum_fed_t           = ifelse(ben_dattest2$mortality == maximum_mort[i],ben_dattest2$fed_t,99999)
  ben_dattest2$maximum_time = ifelse(ben_dattest2$mortality == maximum_mort[i],ben_dattest2$time,99999)
  ben_dattest2$maximum_deterrence_IRS  = ifelse(ben_dattest2$mortality == maximum_mort[i],ben_dattest2$deterrence_IRS,99999)
  ben_dattest2$maximum_deterrence_total= ifelse(ben_dattest2$mortality == maximum_mort[i],ben_dattest2$deterrence_total,99999)
  maximum_d_t[i]              = max(ben_dattest2$maximum_d_t[ben_dattest2$maximum_d_t < 9999])
  maximum_n_t[i]              = max(ben_dattest2$maximum_n_t[ben_dattest2$maximum_n_t < 9999])
  maximum_fed_t[i]            = max(ben_dattest2$maximum_fed_t[ben_dattest2$maximum_fed_t < 9999])
  times[i]                    = max(ben_dattest2$maximum_time[ben_dattest2$maximum_time < 9999])
  maximum_deterrence_IRS[i]   = max(ben_dattest2$maximum_deterrence_IRS[ben_dattest2$maximum_deterrence_IRS < 99999])
  maximum_deterrence_total[i] = max(ben_dattest2$maximum_deterrence_total[ben_dattest2$maximum_deterrence_total < 99999])
  
}

## Select out these data and fit to these

dattest2_b=list(N=2*length(unique(ben_dattest2$time)),
                n_t=c(minimum_n_t,maximum_n_t),
                d_t=c(minimum_d_t,maximum_d_t),
                fed_t=c(minimum_fed_t,maximum_fed_t),
                deterrence_IRS = c(minimum_deterrence_IRS,maximum_deterrence_IRS),
                deterrence_total = c(minimum_deterrence_total,maximum_deterrence_total),
                time=rep(unique(ben_dattest2$time),2),
                N_IRS = 2,
                IRS = c(rep(1,length(unique(ben_dattest2$time))),
                        rep(2,length(unique(ben_dattest2$time)))) )


stan2_b <- stan(file="R code\\stan models\\random_effects_mean_fits_model.stan", 
                data=dattest2_b, 
                warmup=500,
                control = list(adapt_delta = 0.9,
                               max_treedepth = 20),
                iter=1000, chains=4)

#library(shinystan)
#launch_shinystan(stan2_b)

base <- extract(stan2_b)

##Bendiocarb minimum
mean(base$alpha1[,1]);mean(base$alpha2[,1]) ##  0.5348226;-0.02581418
mean(base$beta1[,1]);mean(base$beta2[,1]) ##   -0.4640009;0.01610091
mean(base$omega1[,1]);mean(base$omega2[,1]) ## -1.531136;0.008998386


##Bendiocarb maximum
mean(base$alpha1[,2]);mean(base$alpha2[,2]) ##  2.064936;-0.02998638
mean(base$beta1[,2]);mean(base$beta2[,2]) ##   -2.05132; 0.02299971
mean(base$omega1[,2]);mean(base$omega2[,2]) ## -1.341299; 0.007898068








###########################
##
## 1_estimating the uncertainty 
## around IRS impacts

## Take the minimum estimate for each time
## point through the year

sum_dattest2$mortality = sum_dattest2$d_t/sum_dattest2$n_t 

unique(sum_dattest2$time)

minimum_mort = minimum_d_t = minimum_n_t = minimum_fed_t = times = 
  minimum_deterrence_IRS = minimum_deterrence_total = minimum_time = 
  maximum_mort = maximum_d_t = maximum_n_t = maximum_fed_t = 
  maximum_deterrence_IRS = maximum_deterrence_total = maximum_time  = sum_times = numeric(length(unique(sum_dattest2$time))) 
for(i in 1:length(unique(sum_dattest2$time))){
  minimum_mort[i] = min(sum_dattest2$mortality[sum_dattest2$time == unique(sum_dattest2$time)[i]])
  maximum_mort[i] = max(sum_dattest2$mortality[sum_dattest2$time == unique(sum_dattest2$time)[i]])
  
  
}
for(i in 1:length(unique(sum_dattest2$time))){
  sum_dattest2$minimum_d_t             = ifelse(sum_dattest2$mortality == minimum_mort[i],sum_dattest2$d_t,99999)
  sum_dattest2$minimum_n_t             = ifelse(sum_dattest2$mortality == minimum_mort[i],sum_dattest2$n_t,99999)
  sum_dattest2$minimum_fed_t           = ifelse(sum_dattest2$mortality == minimum_mort[i],sum_dattest2$fed_t,99999)
  sum_dattest2$minimum_time            = ifelse(sum_dattest2$mortality == minimum_mort[i],sum_dattest2$time,99999)
  minimum_d_t[i]              = max(sum_dattest2$minimum_d_t[sum_dattest2$minimum_d_t < 9999])
  minimum_n_t[i]              = max(sum_dattest2$minimum_n_t[sum_dattest2$minimum_n_t < 9999])
  minimum_fed_t[i]            = max(sum_dattest2$minimum_fed_t[sum_dattest2$minimum_fed_t < 9999])
  sum_times[i]                = max(sum_dattest2$minimum_time[sum_dattest2$minimum_time < 9999])
}
for(i in 1:length(unique(sum_dattest2$time))){ 
  sum_dattest2$minimum_deterrence_IRS  = ifelse(sum_dattest2$mortality == minimum_mort[i],sum_dattest2$deterrence_IRS,99999)
  sum_dattest2$minimum_deterrence_total= ifelse(sum_dattest2$mortality == minimum_mort[i],sum_dattest2$deterrence_total,99999)
  
  minimum_deterrence_IRS[i]   = max(sum_dattest2$minimum_deterrence_IRS[sum_dattest2$minimum_deterrence_IRS < 99999])
  minimum_deterrence_total[i] = max(sum_dattest2$minimum_deterrence_total[sum_dattest2$minimum_deterrence_total < 99999])
}
for(i in c(1:length(unique(sum_dattest2$time)))){
  sum_dattest2$maximum_d_t             = ifelse(sum_dattest2$mortality == maximum_mort[i],sum_dattest2$d_t,99999)
  sum_dattest2$maximum_n_t             = ifelse(sum_dattest2$mortality == maximum_mort[i],sum_dattest2$n_t,99999)
  sum_dattest2$maximum_fed_t           = ifelse(sum_dattest2$mortality == maximum_mort[i],sum_dattest2$fed_t,99999)
  sum_dattest2$maximum_time = ifelse(sum_dattest2$mortality == maximum_mort[i],sum_dattest2$time,99999)
  sum_dattest2$maximum_deterrence_IRS  = ifelse(sum_dattest2$mortality == maximum_mort[i],sum_dattest2$deterrence_IRS,99999)
  sum_dattest2$maximum_deterrence_total= ifelse(sum_dattest2$mortality == maximum_mort[i],sum_dattest2$deterrence_total,99999)
  maximum_d_t[i]              = max(sum_dattest2$maximum_d_t[sum_dattest2$maximum_d_t < 9999])
  maximum_n_t[i]              = max(sum_dattest2$maximum_n_t[sum_dattest2$maximum_n_t < 9999])
  maximum_fed_t[i]            = max(sum_dattest2$maximum_fed_t[sum_dattest2$maximum_fed_t < 9999])
  times[i]                    = max(sum_dattest2$maximum_time[sum_dattest2$maximum_time < 9999])
  maximum_deterrence_IRS[i]   = max(sum_dattest2$maximum_deterrence_IRS[sum_dattest2$maximum_deterrence_IRS < 99999])
  maximum_deterrence_total[i] = max(sum_dattest2$maximum_deterrence_total[sum_dattest2$maximum_deterrence_total < 99999])
  
}

## Select out these data and fit to these

dattest2_b=list(N=2*length(unique(sum_dattest2$time)),
                n_t=c(minimum_n_t,maximum_n_t),
                d_t=c(minimum_d_t,maximum_d_t),
                fed_t=c(minimum_fed_t,maximum_fed_t),
                deterrence_IRS = c(minimum_deterrence_IRS,maximum_deterrence_IRS),
                deterrence_total = c(minimum_deterrence_total,maximum_deterrence_total),
                time=rep(unique(sum_dattest2$time),2),
                N_IRS = 2,
                IRS = c(rep(1,length(unique(sum_dattest2$time))),
                        rep(2,length(unique(sum_dattest2$time)))) )


stan2_b <- stan(file="R code\\stan models\\random_effects_mean_fits_model.stan", 
                data=dattest2_b, 
                warmup=500,
                control = list(adapt_delta = 0.9,
                               max_treedepth = 20),
                iter=1000, chains=4)

#library(shinystan)
#launch_shinystan(stan2_b)

base <- extract(stan2_b)

##Sumishield minimum
mean(base$alpha1[,1]);mean(base$alpha2[,1]) ##  -0.1274569;-0.005978999
mean(base$beta1[,1]);mean(base$beta2[,1]) ##   -0.5556669;0.007310865
mean(base$omega1[,1]);mean(base$omega2[,1]) ## -0.09236936; -0.005978999


##Sumishield maximum
mean(base$alpha1[,2]);mean(base$alpha2[,2]) ##  2.035248;-0.0104257
mean(base$beta1[,2]);mean(base$beta2[,2]) ##  -2.481011; 0.01119485
mean(base$omega1[,2]);mean(base$omega2[,2]) ## -1.444987;-0.001016367

DEAD = dattest2_b$d_t/dattest2_b$n_t
FED = dattest2_b$fed_t/dattest2_b$n_t

plot(DEAD ~ dattest2_b$time,ylab="Mortlality range",ylim=c(0,1),col="blue",pch=rep(c(15,17),each=14),
     main = "bendiocarb",cex.main=1.2,xlim=c(1,365),xaxt="n",
     xlab="",yaxt="n",cex.lab=1,cex.axis=1,cex=1.4)
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20),cex.lab=1.3,cex.axis=1)

mean_valssp_checker4Min = 1 / (1 + exp(-mean(base$alpha1[,1]) - mean(base$alpha2[,1])*time))
mean_valssp_checker4Max = 1 / (1 + exp(-mean(base$alpha1[,2]) - mean(base$alpha2[,2])*time))

polygon(c(time,rev(time)),c(mean_valssp_checker4Min,rev(mean_valssp_checker4Max)),col=transp("gold4",0.4),border=NA)

time=seq(1,365,length=365)

mean_valssp_checker4Min = 1 / (1 + exp(--0.1260343 --0.005987383*time))
mean_valssp_checker4Max = 1 / (1 + exp(-2.037771 -  -0.01043863*time))
lines(mean_valssp_checker4Min ~ time)

plot(FED ~ dattest2_b$time,ylab="Mortlality range",ylim=c(0,1),col="gold4",pch=rep(c(15,17),each=14),
     main = "Actellic",cex.main=1.2,xlim=c(1,365),xaxt="n",
     xlab="",yaxt="n",cex.lab=1,cex.axis=1,cex=1.4)
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20),cex.lab=1.3,cex.axis=1)

mean_valsfp_checker4Min = 1 / (1 + exp(-mean(base$beta1[,1]) - mean(base$beta2[,2])*time))
mean_valsfp_checker4Max = 1 / (1 + exp(-mean(base$beta1[,2]) - mean(base$beta2[,1])*time))

polygon(c(time,rev(time)),c(mean_valsfp_checker4Min,rev(mean_valsfp_checker4Max)),col=transp("gold4",0.4),border=NA)

time=seq(1,365,length=365)

mean_valsfp_checker4Min = 1 / (1 + exp(--0.9862012 - 0.008359408*time))
mean_valsfp_checker4Max = 1 / (1 + exp(--2.550309 - 0.004340422*time))
lines(mean_valsfp_checker4Min ~ time)

DET = dattest2_b$deterrence_IRS/dattest2_b$deterrence_total

mean_valsdet_checker4Min = 1 / (1 + exp(-mean(base$omega1[,1]) - mean(base$alpha2[,1])*time))
mean_valsdet_checker4Max = 1 / (1 + exp(-mean(base$omega1[,2]) - mean(base$alpha2[,2])*time))

plot(DET ~ dattest2_b$time,ylim=c(0,1),col="gold4",pch=rep(c(15,17),each=14),
     xlab="",xaxt="n",yaxt="n",cex.lab=1,cex.axis=1,xaxt="n",cex=1.4)
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20),cex.lab=1.3,cex.axis=1)
axis(1,at=seq(0,365,120),labels=seq(0,365,120),cex.lab=1.3,cex.axis=1)
polygon(c(time,rev(time)),c(mean_valsdet_checker4Min,rev(mean_valsdet_checker4Max)),col=transp("gold4",0.4),border=NA)

mean_valsdet_checker4Min = 1 / (1 + exp(--0.0923561 - -0.005987383*time))
mean_valsdet_checker4Max = 1 / (1 + exp(--1.44346 - -0.01043863*time))

