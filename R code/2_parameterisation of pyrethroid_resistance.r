################
##
## Probability adjusted estimates from West African huts
##
#################
## libraries used
library(rstan)
library(shinystan)
library(adegenet)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

pyr_dat = read.csv("data/data_summaryPyrethroids.csv",header=TRUE)
dim(pyr_dat)

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

pyr_dat3$timedays = pyr_dat3$time*30

pyr_dat2$GROUP = ifelse(pyr_dat2$GROUP == 22,19,
                        ifelse(pyr_dat2$GROUP == 23,20,
                               ifelse(pyr_dat2$GROUP ==24,21,pyr_dat2$GROUP)))
pyr_dat2 = subset(pyr_dat2,pyr_dat2$Ntotalfemalemosq_IRS > 5)

pyr_dat2$timedays = pyr_dat2$time*30

dfcheck = data.frame(pyr_dat2$Study,pyr_dat2$GROUP)
dfcheck[!duplicated(dfcheck), ]
k0 = 0.699 ##the prob of feeding in absence of interventions

pyr_dattest2 = list(N=nrow(pyr_dat2),
                    n_t=pyr_dat2$Ntotalfemalemosq_IRS,
                    d_t=pyr_dat2$Ntotaldied_IRS,
                    fed_t=round(pyr_dat2$Nbloodfed_IRS*(1 - pyr_dat2$Ntotaldied_IRS/pyr_dat2$Ntotalfemalemosq_IRS),0),
                    deterrence_IRS = ifelse(c(pyr_dat2$Ntotalfemalemosq_C-pyr_dat2$Ntotalfemalemosq_IRS)<0,0,
                                            c(pyr_dat2$Ntotalfemalemosq_C-pyr_dat2$Ntotalfemalemosq_IRS)),
                    deterrence_total = c(pyr_dat2$Ntotalfemalemosq_IRS+pyr_dat2$Ntotalfemalemosq_C),
                    time=(pyr_dat2$Months_since_IRS*30),
                    N_IRS = length(unique(pyr_dat2$GROUP)),
                    IRS = pyr_dat2$GROUP )

testall <- stan(file="R code/stan models/random_effects_fit_mortality_to_bioassay.stan", 
                data=pyr_dattest2, 
                warmup=500,
                control = list(adapt_delta = 0.8,
                               max_treedepth = 20),
                iter=1000, chains=4)

#################################################
##
## So, the logit function describes the hut_mortality, feeding and deterrence relationships
## as per those studies with a bioassay score
## These are taken directly from test_1


## Here are the parameters that describe
##  i) IRS_decay_mort1, IRS_decay_mort2 
## ii) IRS_decay_succ1, IRS_decay_succ2
##iii) IRS_decay_det1, IRS_decay_det2

##Check the model fit
my_shinystan <- as.shinystan(testall)
launch_shinystan(my_shinystan)

## This gives us the 
##    IRS_decay_MORT1a and IRS_decay_MORT2a  
##    IRS_decay_SUCC1a and IRS_decay_SUCC2a 
##    IRS_decay_DET1a and IRS_decay_DET2a  
## for a = 21 Pyrethroids with different starting 24 hour hut mortality (feeding and deterrence)
parm2 = extract(testall, permuted=TRUE);names(parm2)

alpha1 = mean(parm2$alpha1) ##
alpha2 = mean(parm2$alpha2) ##
beta1 = mean(parm2$beta1)   ##
beta2 = mean(parm2$beta2)   ##
omega1 = mean(parm2$omega1) ##
omega2 = mean(parm2$omega2) ##
#################################################
##
## So, the logit function describes the hut_mortality, feeding and deterrence relationships
## as per those studies with a bioassay score
## These are taken directly from test_1


## Here are the parameters that describe
##  i) IRS_decay_mort1, IRS_decay_mort2 
## ii) IRS_decay_succ1, IRS_decay_succ2
##iii) IRS_decay_det1, IRS_decay_det2

lowApyr <- uppApyr <- meanvalsApyr <- 
  lowBpyr <- uppBpyr <- meanvalsBpyr <- 
  lowCpyr <- uppCpyr <- meanvalsCpyr <- matrix(data=NA,nrow=365,ncol=21)
for(k in 1:365){
  for(n in 1:21){
    lowApyr[k,n] <- quantile(parm2$sp_ppc[,n,k],0.05)
    uppApyr[k,n] <- quantile(parm2$sp_ppc[,n,k],0.95)
    meanvalsApyr[k,n] <- mean(parm2$sp_ppc[,n,k])
    
    lowBpyr[k,n] <- quantile(parm2$fp_ppc[,n,k],0.05)
    uppBpyr[k,n] <- quantile(parm2$fp_ppc[,n,k],0.95)
    meanvalsBpyr[k,n] <- mean(parm2$fp_ppc[,n,k])
    
    lowCpyr[k,n] <- quantile(parm2$det_ppc[,n,k],0.05)
    uppCpyr[k,n] <- quantile(parm2$det_ppc[,n,k],0.95)
    meanvalsCpyr[k,n] <- mean(parm2$det_ppc[,n,k])
  }
}

func_plots = function(data1,labs){
  time =1:365
  plot(data1[,1] ~ time, ylim=c(0,1),
       xlim=c(1,365), ylab=labs,xlab="Time (days)")
  for(i in 1:21){
    lines(data1[,i] ~ time)
  }
  
}
par(mfrow=c(1,3))
func_plots(meanvalsApyr,"Hut mortality")
func_plots(meanvalsBpyr,"Hut feeding")
func_plots(meanvalsCpyr,"Hut deterrence")

LD50=LD20=mosqSurv=numeric(21)
for(i in 1:21){
  mosqSurv[i] = meanvalsApyr[1,i]
  LD50[i] = length(which(meanvalsApyr[,i] > 0.5))
  LD20[i] = length(which(meanvalsApyr[,i] > 0.2))
}
plot(LD50~mosqSurv,xlim=c(0,1))
fit1 = 1/(1+exp(5 +  8* mosqSurv ))
fit1 = lm(LD50 ~ I(mosqSurv^2) + mosqSurv) #fits a model with a quadratic term
fit2line = predict(fit1, data.frame(mosqSurv = seq(min(mosqSurv),1,length=9)),se=TRUE)
upp1 = fit2line$fit + 1.96*fit2line$se.fit
low1 = fit2line$fit - 1.96*fit2line$se.fit
lines(fit2line$fit~seq(min(mosqSurv),1,length=9))
polygon(c(seq(min(mosqSurv),1,length=9),rev(seq(min(mosqSurv),1,length=9))),c(upp1,rev(low1)),col=transp("red",0.2),border=NA)

plot(LD20~mosqSurv,xlim=c(0,1))
fit2 = 1/(1+exp(5 +  8* mosqSurv ))
fit2 = lm(LD20 ~ mosqSurv) #fits a model with a quadratic term
fit2line2 = predict(fit2, data.frame(mosqSurv = seq(min(mosqSurv),1,length=9)),se=TRUE)
upp2 = fit2line2$fit + 1.96*fit2line2$se.fit
low2 = fit2line2$fit - 1.96*fit2line2$se.fit
lines(fit2line2$fit~seq(min(mosqSurv),1,length=9))
polygon(c(seq(min(mosqSurv),1,length=9),rev(seq(min(mosqSurv),1,length=9))),c(upp2,rev(low2)),col=transp("red",0.2),border=NA)

res_mort_1 = c(mean(parm2$alpha1[,1]),mean(parm2$alpha1[,2]),mean(parm2$alpha1[,3]),mean(parm2$alpha1[,4]),
               mean(parm2$alpha1[,5]),mean(parm2$alpha1[,6]),mean(parm2$alpha1[,7]),mean(parm2$alpha1[,8]),
               mean(parm2$alpha1[,9]),mean(parm2$alpha1[,10]),mean(parm2$alpha1[,11]),mean(parm2$alpha1[,12]),
               mean(parm2$alpha1[,13]),mean(parm2$alpha1[,14]),mean(parm2$alpha1[,15]),mean(parm2$alpha1[,16]),
               mean(parm2$alpha1[,17]),mean(parm2$alpha1[,18]),mean(parm2$alpha1[,19]),mean(parm2$alpha1[,20]),
               mean(parm2$alpha1[,21]))

res_mort_2 = c(mean(parm2$alpha2[,1]),mean(parm2$alpha2[,2]),mean(parm2$alpha2[,3]),mean(parm2$alpha2[,4]),
               mean(parm2$alpha2[,5]),mean(parm2$alpha2[,6]),mean(parm2$alpha2[,7]),mean(parm2$alpha2[,8]),
               mean(parm2$alpha2[,9]),mean(parm2$alpha2[,10]),mean(parm2$alpha2[,11]),mean(parm2$alpha2[,12]),
               mean(parm2$alpha2[,13]),mean(parm2$alpha2[,14]),mean(parm2$alpha2[,15]),mean(parm2$alpha2[,16]),
               mean(parm2$alpha2[,17]),mean(parm2$alpha2[,18]),mean(parm2$alpha2[,19]),mean(parm2$alpha2[,20]),
               mean(parm2$alpha2[,21]))


res_succ_1 = c(mean(parm2$beta1[,1]),mean(parm2$beta1[,2]),mean(parm2$beta1[,3]),mean(parm2$beta1[,4]),
               mean(parm2$beta1[,5]),mean(parm2$beta1[,6]),mean(parm2$beta1[,7]),mean(parm2$beta1[,8]),
               mean(parm2$beta1[,9]),mean(parm2$beta1[,10]),mean(parm2$beta1[,11]),mean(parm2$beta1[,12]),
               mean(parm2$beta1[,13]),mean(parm2$beta1[,14]),mean(parm2$beta1[,15]),mean(parm2$beta1[,16]),
               mean(parm2$beta1[,17]),mean(parm2$beta1[,18]),mean(parm2$beta1[,19]),mean(parm2$beta1[,20]),
               mean(parm2$beta1[,21]))

res_succ_2 = c(mean(parm2$beta2[,1]),mean(parm2$beta2[,2]),mean(parm2$beta2[,3]),mean(parm2$beta2[,4]),
               mean(parm2$beta2[,5]),mean(parm2$beta2[,6]),mean(parm2$beta2[,7]),mean(parm2$beta2[,8]),
               mean(parm2$beta2[,9]),mean(parm2$beta2[,10]),mean(parm2$beta2[,11]),mean(parm2$beta2[,12]),
               mean(parm2$beta2[,13]),mean(parm2$beta2[,14]),mean(parm2$beta2[,15]),mean(parm2$beta2[,16]),
               mean(parm2$beta2[,17]),mean(parm2$beta2[,18]),mean(parm2$beta2[,19]),mean(parm2$beta2[,20]),
               mean(parm2$beta1[,21]))


res_det_1 = c(mean(parm2$omega1[,1]),mean(parm2$omega1[,2]),mean(parm2$omega1[,3]),mean(parm2$omega1[,4]),
              mean(parm2$omega1[,5]),mean(parm2$omega1[,6]),mean(parm2$omega1[,7]),mean(parm2$omega1[,8]),
              mean(parm2$omega1[,9]),mean(parm2$omega1[,10]),mean(parm2$omega1[,11]),mean(parm2$omega1[,12]),
              mean(parm2$omega1[,13]),mean(parm2$omega1[,14]),mean(parm2$omega1[,15]),mean(parm2$omega1[,16]),
              mean(parm2$omega1[,17]),mean(parm2$omega1[,18]),mean(parm2$omega1[,19]),mean(parm2$omega1[,20]),
              mean(parm2$omega1[,21]))

##*matching decay of the deterrence effect to mortality in absense of other data 
res_det_2 = c(mean(parm2$alpha2[,1]),mean(parm2$alpha2[,2]),mean(parm2$alpha2[,3]),mean(parm2$alpha2[,4]),
              mean(parm2$alpha2[,5]),mean(parm2$alpha2[,6]),mean(parm2$alpha2[,7]),mean(parm2$alpha2[,8]),
              mean(parm2$alpha2[,9]),mean(parm2$alpha2[,10]),mean(parm2$alpha2[,11]),mean(parm2$alpha2[,12]),
              mean(parm2$alpha2[,13]),mean(parm2$alpha2[,14]),mean(parm2$alpha2[,15]),mean(parm2$alpha2[,16]),
              mean(parm2$alpha2[,17]),mean(parm2$alpha2[,18]),mean(parm2$alpha2[,19]),mean(parm2$alpha2[,20]),
              mean(parm2$alpha2[,21]))


max_hut_mortality = pyr_dat2$Ntotaldied_IRS[pyr_dat2$time == 0] / pyr_dat2$Ntotalfemalemosq_IRS[pyr_dat2$time == 0]
max_hut_feeding = pyr_dat2$Nbloodfed_IRS[pyr_dat2$time == 0]/pyr_dat2$Ntotalfemalemosq_IRS[pyr_dat2$time == 0]
max_hut_feed_scaled = (pyr_dat2$Nbloodfed_IRS[pyr_dat2$time == 0] * 
                         (1 - pyr_dat2$Ntotaldied_IRS[pyr_dat2$time == 0]/pyr_dat2$Ntotalfemalemosq_IRS[pyr_dat2$time == 0]))/
  pyr_dat2$Ntotalfemalemosq_IRS[pyr_dat2$time == 0]
max_hut_deterrence = ifelse(pyr_dat2$Ntotalfemalemosq_C[pyr_dat2$time == 0]-pyr_dat2$Ntotalfemalemosq_IRS[pyr_dat2$time == 0] < 0, 0,
                            pyr_dat2$Ntotalfemalemosq_C[pyr_dat2$time == 0]-pyr_dat2$Ntotalfemalemosq_IRS[pyr_dat2$time == 0])/
  pyr_dat2$Ntotalfemalemosq_C[pyr_dat2$time == 0]
LD50_prediction = fit2line$fit
LD50_prediction_upper = upp1
LD50_prediction_lower = low1
LD20_prediction = fit2line2$fit
LD20_prediction_upper = upp2
LD20_prediction_lower = low2

step3data <- data.frame(max_hut_mortality,res_mort_1,res_mort_2,
                        max_hut_feeding,max_hut_feed_scaled,res_succ_1,res_succ_2,
                        max_hut_deterrence,res_det_1,res_det_2,LD20,LD50,mosqSurv,
                        LD50_prediction,LD50_prediction_upper,LD50_prediction_lower,
                        LD20_prediction,LD20_prediction_upper,LD20_prediction_lower)
write.csv(step3data,"H:\\Ellie\\IRS and resistance\\IRSresistance_paper25102016\\step3data_v2resubmission.csv")



#################################
##
## Plot out the relationships for parameters and t = 1 measures
##
##
##################################
step3data <- read.csv("data\\step3data_v2resubmission.csv",header=TRUE)

#par(mfrow=c(1,1))
par(mar=c(8,6,2,9))
plot(step3data$res_mort_1 ~ step3data$max_hut_mortality,pch=15,col="purple",frame=FALSE,
     ylab=expression(paste("Parameter: ",italic(IRS_decay_MORT1[1]))),xlim=c(0,1),
     cex.lab=1.6,cex=1.2,cex.axis=1.6,xlab="24 Hut mortality (t = 1)")
LM1 = lm(step3data$res_mort_1 ~ step3data$max_hut_mortality)
p_conf4 <- predict(LM1,interval="confidence")
dat=data.frame(p_conf4[,1:3],step3data$max_hut_mortality)
abline(lm(step3data$res_mort_1 ~ step3data$max_hut_mortality),col=transp("purple",0.4),lwd=2)
newdata <- dat[order(step3data$max_hut_mortality),] 
polygon(c(newdata$step3data.max_hut_mortality,rev(newdata$step3data.max_hut_mortality)),
        c(newdata$lwr,rev(newdata$upr)),col=transp("purple",0.2),border=NA)
summary(LM1)
LM1$coefficients


par(new=TRUE)

##dropping anomaly at -0.1


plot(step3data$res_mort_2[2:21] ~ step3data$max_hut_mortality[2:21],bty="n",cex.lab=1.8,pch=15,col=transp("darkgreen",0.5),cex=1.4,
     yaxt="n",ylim=c(-0.1,0.02),xaxt="n",ylab="",
     xlab="",xlim=c(0,1))
axis(4,las=2,at=seq(-0.1,0.02,0.02),
     labels=seq(-0.1,0.02,0.02),cex.axis=1.6,cex=1.6)
mtext(expression(paste("Parameter:  ", italic(IRS_decay_MORT2[1]))), side = 4, line = 5, cex = 1.2,col="darkgreen")
LM2 = lm(step3data$res_mort_2[2:21] ~ step3data$max_hut_mortality[2:21])
p_conf4 <- predict(LM2,interval="confidence")
dat=data.frame(p_conf4[,1:3],step3data$max_hut_mortality[2:21])
abline(lm(step3data$res_mort_2[2:21] ~ step3data$max_hut_mortality[2:21]),col=transp("darkgreen",0.4),lwd=2)
newdata <- dat[order(step3data$max_hut_mortality[2:21]),] 
polygon(c(newdata$step3data.max_hut_mortality,rev(newdata$step3data.max_hut_mortality)),
        c(newdata$lwr,rev(newdata$upr)),col=transp("darkgreen",0.2),border=NA)
summary(LM2)
LM2$coefficients

plot(step3data$res_succ_1 ~ step3data$max_hut_feed_scaled,pch=15,col="purple",frame=FALSE,
     ylab=expression(paste("Parameter: ",italic(IRS_decay_SUCC1[1]))),xlim=c(0,1),
     cex.lab=1.6,cex=1.2,cex.axis=1.6,xlab="24 Hut Successful feeding (t = 1)")
LM3 = lm(step3data$res_succ_1 ~ step3data$max_hut_feed_scaled)
p_conf4 <- predict(LM3,interval="confidence")
dat=data.frame(p_conf4[,1:3],step3data$max_hut_feed_scaled)
abline(lm(step3data$res_succ_1 ~ step3data$max_hut_feed_scaled),col=transp("purple",0.4),lwd=2)
newdata <- dat[order(step3data$max_hut_feed_scaled),] 
polygon(c(newdata$step3data.max_hut_feed_scaled,rev(newdata$step3data.max_hut_feed_scaled_inhibition)),
        c(newdata$lwr,rev(newdata$upr)),col=transp("purple",0.2),border=NA)
summary(LM3)
LM3$coefficients

par(new=TRUE)
plot(step3data$res_succ_2[1:20] ~ step3data$max_hut_feed_scaled[1:20],bty="n",cex.lab=1.8,pch=15,
     col=transp("darkgreen",0.5),cex=1.4,
     yaxt="n",xaxt="n",ylab="",ylim=c(min(step3data$res_succ_2),max(step3data$res_succ_2)),
     xlab="",xlim=c(0,1))
axis(4,las=2,at=seq(min(step3data$res_succ_2),max(step3data$res_succ_2),length=4),
     labels=round(seq(min(step3data$res_succ_2),max(step3data$res_succ_2),length=4),4),cex.axis=1.6,cex=1.6)
mtext(expression(paste("Parameter:  ", italic(IRS_decay_SUCC2[1]))), side = 4, line = 5, cex = 1.2,col="darkgreen")
LM4 = lm(step3data$res_succ_2[1:20] ~ step3data$max_hut_feed_scaled[1:20])
p_conf4 <- predict(LM4,interval="confidence")
dat=data.frame(p_conf4[,1:3],step3data$max_hut_feed_scaled[1:20])
abline(lm(step3data$res_succ_2[1:20] ~ step3data$max_hut_feed_scaled[1:20]),col=transp("darkgreen",0.4),lwd=2)
newdata <- dat[order(step3data$max_hut_feed_scaled),] 
polygon(c(newdata$step3data.max_hut_feed_scaled,rev(newdata$step3data.max_hut_feed_scaled)),
        c(newdata$lwr,rev(newdata$upr)),col=transp("darkgreen",0.2),border=NA)
summary(LM4)
LM4$coefficients

step3data2 <- subset(step3data,step3data$res_det_1 > -10)
plot(step3data2$res_det_1 ~ step3data2$max_hut_deterrence,pch=15,col="purple",frame=FALSE,
     ylab=expression(paste("Parameter: ",italic(IRS_decay_DET1[1]))),xlim=c(0,1),ylim=c(-15,5),
     cex.lab=1.6,cex=1.2,cex.axis=1.6,xlab="24 Hut Deterrence (t = 1)")
LM5 = lm(step3data2$res_det_1 ~ step3data2$max_hut_deterrence)
p_conf4 <- predict(LM5,interval="confidence")
dat=data.frame(p_conf4[,1:3],step3data2$max_hut_deterrence)
abline(lm(step3data2$res_det_1 ~ step3data2$max_hut_deterrence),col=transp("purple",0.4),lwd=2)
newdata <- dat[order(step3data2$max_hut_deterrence),] 
polygon(c(newdata$step3data2.max_hut_deterrence,rev(newdata$step3data2.max_hut_deterrence)),
        c(newdata$lwr,rev(newdata$upr)),col=transp("purple",0.2),border=NA)
summary(LM5)
LM5$coefficients

par(new=TRUE)
plot(step3data2$res_det_2 ~ step3data2$max_hut_deterrence,bty="n",cex.lab=1.8,pch=15,col=transp("darkgreen",0.5),cex=1.4,
     yaxt="n",xaxt="n",ylab="",ylim=c(min(step3data2$res_det_2),max(step3data2$res_det_2)),
     xlab="",xlim=c(0,1))
axis(4,las=2,at=seq(min(step3data2$res_det_2),max(step3data2$res_det_2),length=4),
     labels=round(seq(min(step3data2$res_det_2),max(step3data2$res_det_2),length=4),2),cex.axis=1.6,cex=1.6)
range(step3data2$res_det_2)
mtext(expression(paste("Parameter:  ", italic(IRS_decay_DET2[1]))), side = 4, line = 5, cex = 1.2,col="darkgreen")
LM6 = lm(step3data2$res_det_2 ~ step3data2$max_hut_deterrence)
p_conf4 <- predict(LM6,interval="confidence")
dat=data.frame(p_conf4[,1:3],step3data2$max_hut_deterrence)
abline(lm(step3data2$res_det_2 ~ step3data2$max_hut_deterrence),col=transp("darkgreen",0.4),lwd=2)
newdata <- dat[order(step3data2$max_hut_deterrence),] 
polygon(c(newdata$step3data2.max_hut_deterrence,rev(newdata$step3data2.max_hut_deterrence)),
        c(newdata$lwr,rev(newdata$upr)),col=transp("darkgreen",0.2),border=NA)
summary(LM6)
LM6$coefficients

##If AandB              
int1 =  LM1$coefficients[1]  #new -2.587528    #old -1.416835
int2 =  LM2$coefficients[1]  #new  -0.002989876     #old -0.009357982
int3 = LM3$coefficients[1]   #new -2.95455     #old -2.076946 #-1.2557
int4 =  LM4$coefficients[1]  #new 0.01200367    #old 0.002047915 #0.009720328
int5 =  LM5$coefficients[1]  #new    #old -2.960124  
int6 =  LM6$coefficients[1]  #new    #old 0.0003173997  

grad1 = LM1$coefficients[2] #new 5.777369    #old 2.972465 
grad2 = LM2$coefficients[2] #new -0.01362388     #old 0.001040958 
grad3 = LM3$coefficients[2] #new 5.231271   #old 2.352004 #4.2495
grad4 = LM4$coefficients[2] #new -0.004870386   #old 0.009000071 #-0.009046770
grad5 = LM5$coefficients[2] #new  8.542869   #old 4.771978
grad6 = LM6$coefficients[2] #new   0.005316776   #old -0.0077501779 


#####################################################
impb = read.csv("data/pyrethroid_resistance.csv",header = TRUE)
#impb_temp = subset(pyr_dat2, pyr_dat2$time == 0)
#impb_temp$BIOASSAY_MORTALITY = c()
## Prepare the data for pyrethroids:
imp_list_cleaner_f = function(data1){
  
  ##GIVEN THE INFORMATION IN ROWLAND ET AL 2013 RAW DATA AND MULLER UNPUBLISHED DATA
  ##CONTROLS 
  rW_adj = 0.05799048#0.07815834 ##REPELLED WITHOUT FEEDING
  dW_adj = 0.01831188#0.3343705  ##KILLED WITHOUT FEEDING
  dF_adj = 0.04560538#0.6656295  ##KILLED AND FED
  
  
  ## The probability that mosquitoes survive (Ps), feed (Pf) or exit (Pe)
  ##CONTROLS
  rW = ((data1$Ntotalfemalemosq_C - data1$Ntotaldied_C)/ data1$Ntotalfemalemosq_C)*rW_adj ## repelled without feeding
  dW = (data1$Ntotaldied_C/data1$Ntotalfemalemosq_C)*dW_adj                                 ## killed without feeding
  dF = (data1$Ntotaldied_C/data1$Ntotalfemalemosq_C)*dF_adj                                 ## killed with feeding
  
  sO      = 1 - rW - dW - dF  ## FED AT BASELINE
  
  rWI_adj = 0.02930113
  dWI_adj = 0.28
  dFI_adj = 0.72
  
  ## The probability that mosquitoes survive (Ps), feed (Pf) or exit (Pe)
  ##SPRAYED
  rWI = ((data1$Ntotalfemalemosq_IRS - data1$Ntotaldied_IRS)/ data1$Ntotalfemalemosq_IRS)*rWI_adj ## repelled without feeding
  dWI = (data1$Ntotaldied_IRS/data1$Ntotalfemalemosq_IRS)*dWI_adj                                 ## killed without feeding
  dFI = (data1$Ntotaldied_IRS/data1$Ntotalfemalemosq_IRS)*dFI_adj                                 ## killed with feeding
  
  sI      = 1 - rWI - dWI   ## FED AT BASELINE
  
  rP = mean((data1$Ntotalfemalemosq_C - data1$Ntotalfemalemosq_IRS) / data1$Ntotalfemalemosq_C) ##constant deterrence
  
  P_fed_and_alive = sI * sO * (1 - dFI)
  P_fed_and_dead =  sI * (dF + sO * dFI)
  P_unfed_and_dead =dWI + sI * dW
  P_unfed_and_alive=rWI + sI * (rW)
  
  P_killed2 = P_fed_and_dead + P_unfed_and_dead
  P_bloodfed2 = P_fed_and_alive
  P_other2 = 1 - P_killed2 - P_bloodfed2
  
  return(list(N=nrow(data1),
              x =   round(impb$BIOASSAY_MORTALITY),
              n_t=(c(data1$Ntotalfemalemosq_IRS)),
              d_t=round(c(data1$Ntotalfemalemosq_IRS)*P_killed2,0),
              n_det = ifelse(impb$Ntotalfemalemosq_C-impb$Ntotalfemalemosq_IRS < 0, 0, impb$Ntotalfemalemosq_C-impb$Ntotalfemalemosq_IRS),
              n_c = impb$Ntotalfemalemosq_C,
              f_t=round(c(data1$Ntotalfemalemosq_IRS)*P_bloodfed2,0),
              N_IRS = 1,
              IRS = rep(1,nrow(data1))) )
}


#mpb2 = impb[c(1:18),]
#imbp = data.frame(impb2)
#imp_list_1_real  <- imp_list_cleaner_f(impb)   

imp_list_1_real = list(N=nrow(impb),
                       x = round(impb$BIOASSAY_MORTALITY),
                       n_t=impb$Ntotalfemalemosq_IRS,
                       d_t=impb$Ntotaldied_IRS,
                       n_det = ifelse(impb$Ntotalfemalemosq_C-impb$Ntotalfemalemosq_IRS < 0, 0, impb$Ntotalfemalemosq_C-impb$Ntotalfemalemosq_IRS),
                       n_c = impb$Ntotalfemalemosq_C,
                       f_t=round(impb$Nbloodfed_IRS * (1 - impb$Ntotaldied_IRS/impb$Ntotalfemalemosq_IRS)),
                       N_IRS = 1,
                       IRS = rep(1,nrow(impb))) 
##Run the model to fit hut mortality to bioassay mortality for pyrethroids
test_1 <- stan(file="R code\\stan models\\fitting hut mortality and feeding success to bioassay_1.stan", 
               data=imp_list_1_real,
               warmup=500,
               control = list(adapt_delta = 0.9,max_treedepth = 20),#
               #sample_file = "C:\\Users\\Ellie\\Documents\\IRS resistance\\Pyrethroid_1_mort_fed_det.csv",
               iter=1000, chains=1)

##Confirm the fit of the model 
my_shinystan <- as.shinystan(test_1)
launch_shinystan(my_shinystan)

##Either pull out the parameters of interest
parm = extract(test_1, permuted=TRUE);names(parm)
## 
mean(parm$alpha1[,1]) ##-1.881695    ##NEW -1.20088
mean(parm$alpha2[,1]) ## 0.0265811   ##NEW 0.01878267
mean(parm$beta1[,1])  ## 1.435545    ##NEW 0.8682404
mean(parm$beta2[,1])  ## -0.02364654 ##NEW -0.01735203
mean(parm$theta1[,1]) ##-2.080929    ##NEW -1.672587
mean(parm$theta2[,1]) ##0.02         ##NEW -0.0007890287
#
quantile(parm$alpha1[,1],c(0.05,0.95)) ##-2.116358 ##NEW -1.285526 -1.098071 
quantile(parm$alpha2[,1],c(0.05,0.95)) ##0.02547475##NEW 0.01705237 0.02030704
quantile(parm$beta1[,1],c(0.05,0.95))  ##1.593385  ##NEW 0.7894646 0.9470872 
quantile(parm$beta2[,1],c(0.05,0.95))  ##-0.0218743##NEW -0.01889573 -0.01598732 
quantile(parm$theta1[,1],c(0.05,0.95)) ##-2.07176  ##NEW -1.777039 -1.566337
quantile(parm$theta2[,1],c(0.05,0.95)) ##0.02003083##NEW -0.002947599  0.001528303 



lowA <- uppA <- meanvalsA <- 
  lowB <- uppB <- meanvalsB <- 
  lowC <- uppC <- meanvalsC <- matrix(data=NA,nrow=101,ncol=1)
for(k in 0:100){
  for(n in 1:1){
    lowA[k+1,n] <- quantile(parm$y_ppc[,n,k],0.05)
    uppA[k+1,n] <- quantile(parm$y_ppc[,n,k],0.95)
    meanvalsA[k+1,n] <- mean(parm$y_ppc[,n,k])
    
    lowB[k+1,n] <- quantile(parm$k_ppc[,n,k],0.05)
    uppB[k+1,n] <- quantile(parm$k_ppc[,n,k],0.95)
    meanvalsB[k+1,n] <- mean(parm$k_ppc[,n,k])
    
    lowC[k+1,n] <- quantile(parm$det_ppc[,n,k],0.05)
    uppC[k+1,n] <- quantile(parm$det_ppc[,n,k],0.95)
    meanvalsC[k+1,n] <- mean(parm$det_ppc[,n,k])
  }
}


#########################
##
## Can now start from here
##
step3data <- read.csv("data\\step3data_v2resubmission.csv",header=TRUE)

## Bioassay mortality relationships
##Using the relationships defined as other chemistries and not accounting for Controls
alpha1 = -1.059923 ##THIS IS IF DOING DATA REARRANGING AS FUNCTION imp_list_cleaner_f -1.20088# -1.942755
alpha2 = 0.02387883##0.01878267# 0.03361892
beta1 =  0.7674771 ##0.8682404# 1.849397
beta2 =  -0.03166185 ##-0.01735203# -0.04921294
theta1 = -1.673563 ##-1.672587# -2.071873
theta2 = -0.0007511497 ##-0.0007890287# 0.02004906

alpha1qu = -1.1387719 ## -1.285526 # -2.112937
alpha2qu =0.02226354 ##  0.01705237#0.03114543 
beta1qu = 0.6885264  ##0.7894646  # 1.706689 
beta2qu = -0.03346408 ##-0.01889573  # -0.05236427 
theta1qu =-1.794143 ##-1.777039 # -2.225711
theta2qu = -0.003069259  ##-0.002947599   #0.01760429

alpha1ql = -0.9785583##  -1.098071 #-1.793964 
alpha2ql =0.02556298##   0.02030704#0.03624486
beta1ql = 0.8468829 ## 0.9470872 # 1.996350
beta2ql =  -0.02982156 ##  -0.01598732 #-0.04645465
theta1ql = -1.566623 ##  -1.566337# -1.908551 
theta2ql = 0.001473022 ##0.001528303 # 0.02252894
x_bioassay_mort = rev(seq(0,100,length = 101))

lowA <- uppA <- meanvalsA <- 
  lowB <- uppB <- meanvalsB <- 
  lowC <- uppC <- meanvalsC <- matrix(data=NA,nrow=101,ncol=1)
for(k in 1:101){
  for(n in 1:1){
    lowA[k,n] <- 1/ (1 + exp (-alpha1qu - alpha2qu * x_bioassay_mort[k]))
    uppA[k,n] <- 1/ (1 + exp (-alpha1ql - alpha2ql * x_bioassay_mort[k]))
    meanvalsA[k,n] <- 1/ (1 + exp (-alpha1 - alpha2 * x_bioassay_mort[k]))
    
    lowB[k,n] <- 1/ (1 + exp (-beta1qu - beta2qu * x_bioassay_mort[k]))
    uppB[k,n] <- 1/ (1 + exp (-beta1ql - beta2ql * x_bioassay_mort[k]))
    meanvalsB[k,n] <- 1/ (1 + exp (-beta1 - beta2 * x_bioassay_mort[k]))
    
    lowC[k,n] <- 1/ (1 + exp (-theta1qu - theta2qu * x_bioassay_mort[k]))
    uppC[k,n] <- 1/ (1 + exp (-theta1ql - theta2ql * x_bioassay_mort[k]))
    meanvalsC[k,n] <- 1/ (1 + exp (-theta1 - theta2 * x_bioassay_mort[k]))
  }
}

##If AandB              UPDATE THESE...***** ONCE PYRETHROID INDIVIDUAL STUDIES HAVE RUN
int1 = -2.587528 
int2 =  -0.002989876
int3 =  -2.95455 
int4 =   0.01200367 
int5 =  -3.917539 
int6 = 0.002780567 ##same as mortality decay...

grad1 =5.777369 
grad2 = -0.01362388 
grad3 = 5.231271 
grad4 =-0.004870386
grad5 =8.542869 
grad6 =  0.005316776 



## Estimating parameters where there is pyrethroid resistance
## For any given level of mortality in a bioassay (x) where x=1 indicates all mosquitoes die,
## these are linked by the relationship:
##  e.g.

x = 100 - x_bioassay_mort
mort1_temp = 1 - (1 / (1 + exp(-alpha1 - alpha2 * x_bioassay_mort)) )
succ1_temp = 1 - (1 / (1 + exp(-beta1 - beta2 * x_bioassay_mort)) )
det1_temp = 1 - (1 / (1 + exp(-theta1 - theta2 * x_bioassay_mort)) )

#############There are the best case scenario with no resistance
#PYRETHROIDS
#IRS_decay_MORT1	= 0.9766
#IRS_decay_MORT2	= -0.0085
#IRS_decay_SUCC1	= -1.9728
#IRS_decay_SUCC2	= 0.0024
#IRS_decay_DET1	= -0.6540
#IRS_decay_DET2	= -0.0085

time = 1:365
#y_mort_hut_orig = 1 / (1 + exp(-(IRS_decay_MORT1 + IRS_decay_MORT2 * time)))
#y_succ_hut_orig = 1 / (1 + exp(-(IRS_decay_SUCC1 + IRS_decay_SUCC2 * time)))
#y_det_hut_orig = 1 / (1 + exp(-(IRS_decay_DET1 + IRS_decay_DET2 * time)))
#y_rep_hut_orig = 1 - y_mort_hut_orig - y_succ_hut_orig
#plot(y_succ_hut_orig,ylim=c(0,1))
#points(y_mort_hut_orig,col="black")
#points(y_det_hut_orig,col="green")
#points(y_rep_hut_orig,col="orange")
#points(y_succ_hut_orig,col="darkred")

x=rev(seq(0.01,1,length=101))##when x is 1 all mosquitoes die
#x = 1 
death = feed = rep = deter = array(dim=c(365,length(x)))
death2 = feed2 = rep2 = deter2 = array(dim=c(365,length(x)))
params_estimated = array(dim=c(length(x),6))
for(i in 1:length(x)){
  temp_mort_1 = grad1 * (1 / (1 + exp(-alpha1 - alpha2 * 100 * x[i]))) + int1
  temp_mort_2 = grad2 * (1 / (1 + exp(-alpha1 - alpha2 * 100 * x[i]))) + int2
  
  #  temp_mort_1_adj = temp_mort_1/(grad1 * (1 / (1 + exp(-alpha1 - alpha2 * 100 * x[1]))) + int1)*IRS_decay_MORT1
  #  temp_mort_2_adj = temp_mort_2/(grad2 * (1 / (1 + exp(-alpha1 - alpha2 * 100 * x[1]))) + int2)*IRS_decay_MORT2
  
  death[,i] = 1 / (1 + exp(-(temp_mort_1+ temp_mort_2 * time)))
  
  temp_succ_1 = grad3 * (1 / (1 + exp(-beta1 - beta2 * 100 * x[i]))) + int3
  temp_succ_2 = grad4 * (1 / (1 + exp(-beta1 - beta2 * 100 * x[i]))) + int4
  
  #  temp_succ_1_adj = temp_succ_1/(grad3 * (1 / (1 + exp(-beta1 - beta2 * 100 * x[1]))) + int3)*IRS_decay_SUCC1
  #  temp_succ_2_adj = temp_succ_2/(grad4 * (1 / (1 + exp(-beta1 - beta2 * 100 * x[1]))) + int4)*IRS_decay_SUCC2
  #  
  feed[,i] = 1 / (1 + exp(-(temp_succ_1 + temp_succ_2 * time)))
  
  temp_det_1 = grad5 * (1 / (1 + exp(-theta1 - theta2 * 100 * x[i]))) + int5
  temp_det_2 = grad6 * (1 / (1 + exp(-theta1 - theta2 * 100 * x[i]))) + int6
  
  # temp_det_1_adj = temp_det_1/(grad5 * (1 / (1 + exp(-theta1 - theta2 * 100 * x[1]))) + int5)*IRS_decay_DET1
  # temp_det_2_adj = 0#temp_det_2/(grad6 * (1 / (1 + exp(-theta1 - theta2 * 100 * x[1]))) + int6)*IRS_decay_DET2
  
  deter[,i] = 1 / (1 + exp(-(temp_det_1 + temp_mort_2 * time)))
  
  params_estimated[i,1] = temp_mort_1
  params_estimated[i,2] = temp_mort_2
  params_estimated[i,3] = temp_succ_1
  params_estimated[i,4] = temp_succ_2
  params_estimated[i,5] = temp_det_1
  params_estimated[i,6] = temp_det_2
}
rep = 1 - death - feed

demo1 = seq(1,101,10)
##nOW CHECK HOW THIS LOOKS
par(mfrow=c(1,1))
par(mar=c(5,6,2,2))
plot(death[,length(x)]~time,ylim=c(0,1),pch="",yaxt="n",bty="n",
     ylab = "Probability of mosquito behavioural outcome",xlim=c(0,365),
     xlab="Time (days since spray)",cex.lab=1.2,cex.axis=1.2)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1.0),labels=c(0,0.2,0.4,0.6,0.8,1.0),cex.axis=1.2,cex.lab=1.2)

#points(feed_temp[,length(x)]~time,ylim=c(0,1),col="red",pch="")
#points(deter_temp[,length(x)]~time,ylim=c(0,1),col="darkgreen",pch="")
#points(rep_temp[,length(x)]~time,ylim=c(0,1),col="orange",pch="")
for(i in 1:length(demo1)){
 lines(death[,demo1[i]]~time)
  lines(feed[,demo1[i]]~time,col="red")
  lines(deter[,demo1[i]]~time,col="darkgreen")
  lines(rep[,demo1[i]]~time,col="orange")
}
points(y_mort_hut_orig,col="black")
points(y_det_hut_orig,col="darkgreen")
points(y_rep_hut_orig,col="orange")
points(y_succ_hut_orig,col="darkred")

lines(death[,1]~time,lwd=2)
max_p
params_estimated[1,]##best case scenario when all mosquitoes die and there is no resistance
params_estimated[40,]##best case scenario when 50% mosquitoes die - resistance
params_estimated[101,]##best case scenario when all mosquitoes survive - resistance

lines(1 / (1 + exp(-params_estimated[1,1] - params_estimated[1,2] * time)),lwd=4,col="purple")
lines(1 / (1 + exp(-params_estimated[1,3] - params_estimated[1,4] * time)),lwd=4,col="purple")
lines(1 / (1 + exp(-params_estimated[1,5] - params_estimated[1,6] * time)),lwd=4,col="purple")

params = rev(seq(0,1,length=21))
vals_vec = rev(seq(1,101,length=21))
params_store = array(dim=c(length(params),7))
params_store[,1] = params
for(i in 1:length(params)){
  temp_mort_1 = grad1 * (1 / (1 + exp(-alpha1 - alpha2 * vals_vec[i]))) + int1
  temp_mort_2 = grad2 * (1 / (1 + exp(-alpha1 - alpha2 * vals_vec[i]))) + int2
  
  temp_succ_1 = grad3 * (1 / (1 + exp(-beta1 - beta2 * vals_vec[i]))) + int3
  temp_succ_2 = grad4 * (1 / (1 + exp(-beta1 - beta2 * vals_vec[i]))) + int4
  
  temp_det_1 = grad5 * (1 / (1 + exp(-theta1 - theta2 * vals_vec[i]))) + int5
  temp_det_2 = grad6 * (1 / (1 + exp(-theta1 - theta2 * vals_vec[i]))) + int6
  
  params_store[i,2] = temp_mort_1
  params_store[i,3] = temp_mort_2
  params_store[i,4] = temp_succ_1
  params_store[i,5] = temp_succ_2
  params_store[i,6] = temp_det_1
  params_store[i,7] = temp_mort_2
  
}
colnames(params_store) = c("Prop_dying_at_bioassay",
                           "irs_decay_mort1_1",	"irs_decay_mort2_1",
                           "irs_decay_succ1_1",	"irs_decay_succ2_1",
                           "irs_decay_det1_1",	"irs_decay_det2_1")

write.csv(params_store,"estimates/params_IRS_pyrethroid_resistance_v1_each 5 percent.csv")

susceptible = 1 / (1 + exp(-2.0021878 - -0.013813096 * time)) 
resistant = 1 / (1 + exp(--1.0743542 --0.006558161 * time)) 
lines(resistant ~ time,col="aquamarine3",lwd=3)
lines(susceptible ~ time,col="yellow",lwd=3)

###################################
##
## Figure 2
TOTS = feed + rep + death + deter

first_line = feed / TOTS
second_line = (feed + rep) / TOTS
third_line = (feed + rep + deter ) / TOTS

## 1. PYRETHROIDS
Time = 1:365## Year by days

Time2 = rev(Time)
minimal = rep(0,length(Time2))
maximal = rep(1,length(Time2))



hut_mort1 = pyr_dat2$lp0[pyr_dat2$time == 0]  ##check this is there!
LD50 = step3data$LD50 # old c(39,   0,   0,   19,   0,   17,   0, 218,  160,  112, 172,  88,  87,   30)
mosqSurv =step3data$mosqSurv# old c(0.9232822, 0.4376658, 
             #0.2791288, 0.5282247, 0.3873151, 0.5669918, 
             #0.2489425, 0.8065589, 0.6462521, 0.8344493,
             #0.7526795, 0.6507726, 0.6183342, 0.5474685)
LD20 = step3data$LD20# old c(61, 85, 58, 241, 161, 112,  56, 357, 365, 223, 365, 288, 341, 238)

t0_hut_feeding_temp = t0_hut_total_temp = t0_hut_killed_temp = t0_hut_deterred_temp = t0_hut_deterTOT_temp = numeric(length(unique(pyr_dattest2$IRS)))
for(i in 1:length(unique(pyr_dattest2$IRS))){
  t0_hut_feeding_temp[i] = pyr_dattest2$fed_t[pyr_dattest2$IRS==unique(pyr_dattest2$IRS)[i]][1]
  t0_hut_total_temp[i] = pyr_dattest2$n_t[pyr_dattest2$IRS==unique(pyr_dattest2$IRS)[i]][1]
  t0_hut_killed_temp[i] = pyr_dattest2$d_t[pyr_dattest2$IRS==unique(pyr_dattest2$IRS)[i]][1]
  t0_hut_deterred_temp[i] = pyr_dattest2$deterrence_IRS[pyr_dattest2$IRS==unique(pyr_dattest2$IRS)[i]][1]
  t0_hut_deterTOT_temp[i] = pyr_dattest2$deterrence_total[pyr_dattest2$IRS==unique(pyr_dattest2$IRS)[i]][1]
}
t0_hut_mort = t0_hut_killed_temp/t0_hut_total_temp
t0_hut_feeding = t0_hut_feeding_temp/t0_hut_total_temp
t0_hut_deterred = t0_hut_deterred_temp/t0_hut_deterTOT_temp


#0.72 not sure what this is...

t0_hut_feedTre = (0.72*pyr_dat2$Nbloodfed_IRS[pyr_dat2$time == 0])/pyr_dat2$Ntotalfemalemosq_IRS[pyr_dat2$time == 0]
t0_hut_feedCon = pyr_dat2$Nbloodfed_C[pyr_dat2$time == 0]/pyr_dat2$Ntotalfemalemosq_C[pyr_dat2$time == 0]
t0_hut_feeding = t0_hut_feedCon * t0_hut_feedTre

par(mfrow=c(2,3))
par(mar=c(5,5,5,5))
nc = 100:0
impb$survival1 = (100 - impb$BIOASSAY_MORTALITY)
impb$bioassay = impb$BIOASSAY_MORTALITY / 100
plot(c(1 - (impb$Ntotaldied_IRS/impb$Ntotalfemalemosq_IRS)) ~ (1 - impb$bioassay),
     ylab="Mosquito survival in IRS hut trial (%)",cex.lab=1.8,
     xlab="Pyrethroid resistance test (% survival)",cex.axis=1.6,pch="",
     ylim=c(0,1),xlim=c(0,100),yaxt="n",bty="n",las=0)
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20),cex.axis=1.4)
cols = c("purple")
for(j in 1:ncol(lowA)){
  polygon(c(nc,rev(nc)),c(rev(1-lowA[,j]),(1-uppA[,j])),col=transp(cols[j],alpha=0.2),border=NA)
}
lines(rev(1-meanvalsA[,1]) ~ nc,col=cols[1],lwd=2)

par(xpd=NA)
di <- dev.size("in")
x <- grconvertX(c(0, di[1]), from="in", to="user")
y <- grconvertY(c(0, di[2]), from="in", to="user")
fig <- par("fig")
x <- x[1] + (x[2] - x[1]) * fig[1:2]
y <- y[1] + (y[2] - y[1]) * fig[3:4]
txt <- "(A)"
x <- x[1] + strwidth(txt, cex=3) / 2
y <- y[2] - strheight(txt, cex=3) / 2
text(x, y, txt, cex=1.2)
par(xpd=TRUE)
##**Change to matCh up with the studies from fIG 2 AND Supplementary figure 6
  morts = impb$Ntotaldied_IRS  /impb$Ntotalfemalemosq_IRS
points(c(1 - morts) ~ 
         impb$survival1, 
       col=c("red","red","purple","red","purple",
             "seagreen3","seagreen3","purple","purple","seagreen3",
             "seagreen3", "red","purple","seagreen3", "purple", 
             "purple", "purple", "purple"),
             pch=c(15,16,2,3,0,13,11,13,13,13,13,3,17,13,9,4,5,11),cex=2)
vec_type = c(8,9,14,15,17)
col_vec_type= c(rep("blue",4),"orange")
for(i in 1:length(vec_type)){
  points(c(1 - morts[vec_type[i]]) ~ 
           impb$survival1[vec_type[i]],pch=1,cex=4,lwd=2,col=col_vec_type[i])
  
}
#legend(0.05,0.95,legend=c("Matched data","Mis-matched chemistry","Data collated from 2 studies"),
#       ncol=1,pch=c(15,15,15),cex=1.6,col=c("purple","red","seagreen3",NA,"blue","orange"),bty="n")

#legend(0.05,0.95,legend=c(expression(italic("An. gambiae s.l.")),expression(italic("An. arabiensis")),
#                         expression(italic("An. funestus s.l."))),
#       ncol=1,pch=c(NA,1,1),cex=1.6,pt.cex=1.6,col=c(NA,"blue","orange"),bty="n")
#legend(0.05,0.95,legend=c(expression(italic("An. gambiae s.l.")),expression(italic("An. arabiensis")),
#                          expression(italic("An. funestus s.l."))),
#       ncol=1,pch=c(NA,1,1),cex=1.6,pt.cex=2,col=c(NA,"blue","orange"),bty="n")


#text(1,1,"a",cex=2)

feds = c(round(impb$Nbloodfed_IRS,0),round(pyr_dat2$Nbloodfed_IRS[pyr_dat2$time == 0],0))
tots = c(impb$Ntotalfemalemosq_IRS,pyr_dat2$Ntotalfemalemosq_IRS[pyr_dat2$time == 0])
morts2 = c(morts,t0_hut_mort)
prop_feds =feds/tots * morts2

df = data.frame(prop_feds,morts2)
plot(prop_feds ~ morts2,
     ylab="Successful blood-feeding (%)",cex.lab=1.8,
     xlab="Mosquito survival in IRS hut trial (%)",cex.axis=1.6,pch="",
     ylim=c(0,1),las=0,xlim=c(0,1),yaxt="n",bty="n")
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20),cex.axis=1.4)
#text(0,1,"b",cex=2)

par(xpd=NA)
di <- dev.size("in")
x <- grconvertX(c(0, di[1]), from="in", to="user")
y <- grconvertY(c(0, di[2]), from="in", to="user")
fig <- par("fig")
x <- x[1] + (x[2] - x[1]) * fig[1:2]
y <- y[1] + (y[2] - y[1]) * fig[3:4]
txt <- "(B)"
x <- x[1] + strwidth(txt, cex=3) / 2
y <- y[2] - strheight(txt, cex=3) / 2
text(x, y, txt, cex=1.2)
par(xpd=TRUE)

##**change these up to however many we have
points(df$prop_feds[1:18] ~ c(df$morts2[1:18]), 
       col=transp("purple",0.8), pch=c(15,16,2,3,0,13,11,13,13,13,13,3,17,13,9,4,5,11),cex=2)
points(df$prop_feds[19:39] ~ c(df$morts2[19:39]), 
       col=transp("blue",0.6), pch=c(15,1,16,19,18,2,3,rep(19,7),17,8,5,6,0,11,4),cex=2)
########################################################
#pyr_feed_die <- list(N=(nrow(pyr3)+length(morts2)),
#                     n_t=c(pyr3$Ntotalfemalemosq_IRS,imp_list_1_real$n_t),
#                     d_t=c(pyr3$kp_n,imp_list_1_real$f_t),
#                     x=c(pyr3$lp*100,morts2),
#                     N_TREAT = 1,
#                     TREAT = rep(1,(nrow(pyr3)+length(morts2))))

#pyr_feed_die <- list(N=25,
#                     n_t=c(imp_list_1_real$n_t,t0_hut_total_temp),#
#                     d_t=c(imp_list_1_real$f_t,t0_hut_feeding_temp),#
#                     x=c(100-morts2),
#                     N_TREAT = 1,
#                     TREAT = rep(1,25))

pyr_feed_die <- list(N=100,
                     n_t=c(tots,rep(100,61)),
                     d_t=c(round(feds * morts2),rep(0,61)),
                     x=c(round(c(morts2*100)),rep(0,61)),
                     N_TREAT = 1,
                     TREAT = rep(1,100))


test_3 <- stan(file="Q:\\RProjects\\Parameterising_IRS\\logistic_function_density_Treat.stan", 
               data=pyr_feed_die, 
               warmup=500,
               control = list(adapt_delta = 0.8,
                              max_treedepth = 20),
               iter=1000, chains=4)

parm3 = extract(test_3, permuted=TRUE);names(parm3)

lpyr <- upyr <- meanpyr <- 
  matrix(data=NA,nrow=100,ncol=1)
for(k in 1:100){
  for(n in 1:1){
    lpyr[k,n] <- quantile(parm3$y_ppc[,n,k],0.05)
    upyr[k,n] <- quantile(parm3$y_ppc[,n,k],0.95)
    meanpyr[k,n] <- mean(parm3$y_ppc[,n,k])
  }
}
x=seq(0,1,length=100)

for(j in 1:1){
  polygon(c(x,rev(x)),c(lpyr[,j],rev(upyr[,j])),col=transp(cols[j],alpha=0.2),border=NA)
}
lines(meanpyr[,1]~x,col=cols[1],lwd=2)

legend(0.05,0.95,title = "References (Table 2)",legend=c(1,11,2,5,4,14,15,3,7,17,18,10,6,16),
       ncol=3, cex=1.6,pch=c(15,1,16,19,18,2,3,17,8,5,6,0,11,4),col = "purple",bty="n")
#legend(0.05,0.6,legend="Additional data (Supplementary data 1)",
#       pch=13,col = "purple")

############################################
############################################
plot((impb$Ntotalfemalemosq_C-impb$Ntotalfemalemosq_IRS)/
       impb$Ntotalfemalemosq_C ~ impb$BIOASSAY_MORTALITY,
     ylab="Deterrence from entering a hut with IRS (%)",cex.lab=1.8,col="purple",cex=2,
     xlab="Mosquito survival in IRS hut trial (%)",cex.axis=1.6,pch="",
     ylim=c(0,1),las=0,xlim=c(0,100),yaxt="n",bty="n")
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20),cex.axis=1.4)
#text(0,1,"c",cex=2)

par(xpd=NA)
di <- dev.size("in")
x <- grconvertX(c(0, di[1]), from="in", to="user")
y <- grconvertY(c(0, di[2]), from="in", to="user")
fig <- par("fig")
x <- x[1] + (x[2] - x[1]) * fig[1:2]
y <- y[1] + (y[2] - y[1]) * fig[3:4]
txt <- "(C)"
x <- x[1] + strwidth(txt, cex=3) / 2
y <- y[2] - strheight(txt, cex=3) / 2
text(x, y, txt, cex=1.2)
par(xpd=TRUE)

detty = ifelse((impb$Ntotalfemalemosq_C-impb$Ntotalfemalemosq_IRS)/
                 impb$Ntotalfemalemosq_C > 0,
               (impb$Ntotalfemalemosq_C-impb$Ntotalfemalemosq_IRS)/
                 impb$Ntotalfemalemosq_C, 0)
points(detty ~ impb$BIOASSAY_MORTALITY, 
       col=transp("purple",0.8), pch=c(15,16,2,3,0,13,11,13,13,13,13,3,17,13,9,4,5,11),cex=2)

points(t0_hut_deterred ~ c(100*(1-t0_hut_mort)), 
       col=transp("blue",0.6), pch=c(15,1,16,19,18,2,3,rep(19,7),17,8,5,6,0,11,4),cex=2)       

########################################################
#pyr_DET_die <- list(N=nrow(pyr3),
#                    n_t=pyr3$deterrence_total,
#                    d_t=pyr3$deterrence_IRS,
#                    x=pyr3$lp*100,
#                    N_TREAT = 1,
#                    TREAT = rep(1,nrow(pyr3)))
pyr_DET_die <- list(N=39,
                    n_t=(c(impb$Ntotalfemalemosq_C+impb$Ntotalfemalemosq_IRS,
                           t0_hut_deterTOT_temp)),
                    d_t=c(round(detty*(impb$Ntotalfemalemosq_C+impb$Ntotalfemalemosq_IRS)),
                          t0_hut_deterred_temp),
                    x=c((100 - impb$BIOASSAY_MORTALITY),(1-t0_hut_mort)*100),
                    N_TREAT = 1,
                    TREAT = rep(1,39))

test_4 <- stan(file="H:\\Ellie\\IRS and resistance\\IRSresistance_paper25102016\\eLife v1\\logistic_function_density_Treat.stan", 
               data=pyr_DET_die, 
               warmup=500,
               control = list(adapt_delta = 0.8,
                              max_treedepth = 20),
               iter=1000, chains=4)

parm4 = extract(test_4, permuted=TRUE);names(parm4)
mean(parm4$alpha1);quantile(parm4$alpha1,c(0.1,0.9))
#mean(parm4$alpha2)
#xxx = seq(1,100,length(100))
#meanpyr = (1 - 1/(1 + exp(-mean(parm4$alpha1) - mean(parm4$alpha2)*xxx)))
#lpyr  = 1/(1 + exp(-quantile(parm4$alpha1,0.05) - quantile(parm4$alpha2,0.05)*xxx))
#upyr  = 1/(1 + exp(-quantile(parm4$alpha1,0.95) - quantile(parm4$alpha2,0.95)*xxx))

morta = 1:100
meanpyr <- 1/(1 + exp(-mean(parm4$alpha1) - -alpha2*morta))
lpyr <- 1/(1 + exp(-quantile(parm4$alpha1,0.99) - -alpha2qu*morta))
upyr <- 1/(1 + exp(-quantile(parm4$alpha1,0.01) - -alpha2ql*morta))
x=1:100
lines(meanpyr~x,col=cols[1],lwd=2)
polygon(c(x,rev(x)),c(lpyr,rev(upyr)),col=transp("purple",alpha=0.2),border=NA)

#for(k in 1:100){
#  for(n in 1:1){
#    lpyr[k,n] <- quantile(parm4$y_ppc[,n,k],0.05)
#    upyr[k,n] <- quantile(parm4$y_ppc[,n,k],0.95)
#    meanpyr[k,n] <- mean(parm4$y_ppc[,n,k])
#  }
#}
x=1:100
#detup=pyr3$deterrence_IRS/(pyr3$Ntotalfemalemosq_C + pyr3$Ntotalfemalemosq_IRS)
#dedup=pyr3$lp*100
#points(detup ~ dedup,col=transp("purple",0.6),pch=17,cex=1.5)
#for(j in 1:1){
#  polygon(c(x,rev(x)),c(lpyr[,j],rev(upyr[,j])),col=transp(cols[j],alpha=0.2),border=NA)
#}
#lines(meanpyr[,1]~x,col=cols[1],lwd=2)

#for(j in 1:1){
#  polygon(c(x,rev(x)),c(lpyr[,j],rev(upyr[,j])),col=transp(cols[j],alpha=0.2),border=NA)
#}
#lines(meanpyr[,1]~x,col=cols[1],lwd=2)


plot(step3data$LD50~hut_mort1,col=c("red","red","purple","red","purple",
                         "seagreen3","seagreen3","purple","purple","seagreen3",
                         "seagreen3", "red","purple","seagreen3", "purple", 
                         "purple", "purple", "purple"),
     pch=c(15,16,2,3,0,13,11,13,13,13,13,3,17,13,9,4,5,11),cex=2,
    ylab=expression("Number of days"),
     xlab="Susceptibility bioassay test (% mortality)",
     yaxt="n",cex.axis=1.6,cex.lab=1.8,bty="n",
     xaxt="n",xlim=c(0,1),ylim=c(0,365))
axis(2,las=2,at=seq(0,365,60),labels=seq(0,365,60),cex.axis=1.6)
axis(1,las=0,at=seq(0,1,0.2),labels=seq(0,100,20),cex.axis=1.6)
#text(0,360,"d",cex=2)


vec_type = c(8,9,14,15,17)
col_vec_type= c(rep("blue",4),"orange")
for(i in 1:length(vec_type)){
  points(c(1 - morts[vec_type[i]]) ~ 
           impb$survival1[vec_type[i]],pch=1,cex=3,lwd=2,col=col_vec_type[i])
  
}

par(xpd=FALSE)
lines(step3data$LD50_prediction~seq(min(hut_mort1),1,length=21))
polygon(c(seq(min(hut_mort1),1,length=21),rev(seq(min(hut_mort1),1,length=21))),c(step3data$LD50_prediction_upper,rev(step3data$LD50_prediction_lower)),col=transp("grey",0.2),border=NA)

par(xpd=NA)
di <- dev.size("in")
x <- grconvertX(c(0, di[1]), from="in", to="user")
y <- grconvertY(c(0, di[2]), from="in", to="user")
fig <- par("fig")
x <- x[1] + (x[2] - x[1]) * fig[1:2]
y <- y[1] + (y[2] - y[1]) * fig[3:4]
txt <- "(D)"
x <- x[1] + strwidth(txt, cex=3) / 2
y <- y[2] - strheight(txt, cex=3) / 2
text(x, y, txt, cex=1.2)
par(xpd=TRUE)

######################################
x=rev(seq(0,1,0.01))
minimal = rep(0,length(x))
maximal = rep(1,length(x))

plot(rev(first_line[1,]) ~ x,ylim=c(0,1),yaxt="n",
     ylab="Probabilities (%)",xlab="Pyrethroid resistance test (% survival)",
     main="Immediately after spraying",
     cex.axis=1.6,las=0,cex.main=1.6,cex.lab=1.8,bty="n",pch="")
axis(2,las=2,at=seq(0,1,0.2),labels = seq(0,100,20),cex.axis = 1.6)
polygon(c(x,rev(x)),c(rev(first_line[1,]),rev(minimal)),col=transp("red",0.4),border = NA)
polygon(c(x,rev(x)),c(rev(second_line[1,]),first_line[1,]),col=transp("orange",0.6),border = NA)
polygon(c(x,rev(x)),c(rev(third_line[1,]),second_line[1,]),col=transp("darkgreen",0.4),border = NA)
polygon(c(x,rev(x)),c(maximal,third_line[1,]),col=transp("royalblue",0.6),border = NA)
#text(0.03,0.96,"e",cex=2)

par(xpd=NA)
di <- dev.size("in")
x <- grconvertX(c(0, di[1]), from="in", to="user")
y <- grconvertY(c(0, di[2]), from="in", to="user")
fig <- par("fig")
x <- x[1] + (x[2] - x[1]) * fig[1:2]
y <- y[1] + (y[2] - y[1]) * fig[3:4]
txt <- "(E)"
x <- x[1] + strwidth(txt, cex=3) / 2
y <- y[2] - strheight(txt, cex=3) / 2
text(x, y, txt, cex=1.2)
par(xpd=TRUE)

text(0.3,0.9,"Killed",cex=1.5,col="blue")
text(0.5,0.4,"Deterred",cex=1.5,col="darkgreen")
text(0.55,0.3,"Exited no feeding",cex=1.5,col="darkorange3")
text(0.7,0.15,"Successfully blood fed",cex=1.5,col="darkred")


Time = 1:365
Time2 = rev(Time)
minimal = rep(0,length(Time))
maximal = rep(1,length(Time))



probs_plotB = function(res,main_title){
  plot(rev(first_line[,res]) ~ Time2,ylim=c(0,1),yaxt="n",
       ylab="Probabilities (%)",xlab="Time since spraying (days)",
       main=main_title,
       cex.axis=1.6,las=0,cex.main=1.6,cex.lab=1.8,bty="n",pch="")
  axis(2,las=2,at=seq(0,1,0.2),labels = seq(0,100,20),cex.axis = 1.6)
  polygon(c(Time2,rev(Time2)),c(rev(first_line[,res]),rev(minimal)),col=transp("red",0.4),border = NA)
  polygon(c(Time2,rev(Time2)),c(rev(second_line[,res]),first_line[,res]),col=transp("orange",0.6),border = NA)
  polygon(c(Time2,rev(Time2)),c(rev(third_line[,res]),second_line[,res]),col=transp("darkgreen",0.4),border = NA)
  polygon(c(Time2,rev(Time2)),c(maximal,third_line[,res]),col=transp("royalblue",0.6),border = NA)
  
}
probs_plotB(51,"50% mortality at WHO bioassay test")
#text(12,0.96,"f",cex=2)

par(xpd=NA)
di <- dev.size("in")
x <- grconvertX(c(0, di[1]), from="in", to="user")
y <- grconvertY(c(0, di[2]), from="in", to="user")
fig <- par("fig")
x <- x[1] + (x[2] - x[1]) * fig[1:2]
y <- y[1] + (y[2] - y[1]) * fig[3:4]
txt <- "(F)"
x <- x[1] + strwidth(txt, cex=3) / 2
y <- y[2] - strheight(txt, cex=3) / 2
text(x, y, txt, cex=1.2)
par(xpd=TRUE)

##Change in ESTIMATE at 50% resistance
##From xls transformation sheet
##No resistance (all mosquitoes die), 50% and 100%
## Year by days


#######
######
#####

probs_plot = function(res,main_title){
  plot(rev(first_line[,res]) ~ Time2,ylim=c(0,1),yaxt="n",
       ylab="Probabilities (%)",xlab="Time since spraying (days)",
       main=main_title,
       cex.axis=1.6,cex.main=1.6,cex.lab=1.8,bty="n",pch="")
  axis(2,las=2,at=seq(0,1,0.2),labels = seq(0,100,20),cex.axis = 1.6)
  polygon(c(Time2,rev(Time2)),c(rev(first_line[,res]),rev(minimal)),col=transp("red",0.4),border = NA)
  polygon(c(Time2,rev(Time2)),c(rev(second_line[,res]),first_line[,res]),col=transp("orange",0.6),border = NA)
  polygon(c(Time2,rev(Time2)),c(rev(third_line[,res]),second_line[,res]),col=transp("darkgreen",0.4),border = NA)
  polygon(c(Time2,rev(Time2)),c(maximal,third_line[,res]),col=transp("royalblue",0.6),border = NA)
  
  text(50,0.9,"Killed",cex=1.5,col="blue")
  text(70,0.63,"Deterred",cex=1.5,col="darkgreen")
  text(90,0.42,"Exited",cex=1.5,col="darkorange3")
  text(300,0.15,"Blood fed",cex=1.5,col="darkred")
  
}
#probs_plot(1,"100% mortality at WHO bioassay test")

probs_plot(21,"80% mortality at WHO bioassay test")
probs_plot(81,"20% mortality at WHO bioassay test")
probs_plot(101,"0% mortality at WHO bioassay test")

