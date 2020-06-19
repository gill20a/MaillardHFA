#Open necessary libraries
library(mvtnorm)
library(truncdist)
library(Runuran)
library(truncnorm)
library(fitdistrplus)
library(zoo)
library(sm)
library(vioplot)

#Download data file. 
urlfile = "https://raw.githubusercontent.com/gill20a/MaillardHFA/master/MaillardHFA_Exp2_masslossdata.csv"
data<-read.csv(url(urlfile))
head(data)

#Convert days to years
data$time_years<-data$time/365

#Convert percent mass remaining to fraction
data$prop_mass_remaining<-data$mass_remaining/100

#Break into treatments
Champ<-data[data$Site=="Champ",]
ChampChamp<-Champ[Champ$Origin=="Champ",]
ChampBreu<-Champ[Champ$Origin=="Breu",]

Breu<-data[data$Site=="Breu",]
BreuBreu<-Breu[Breu$Origin=="Breu",]
BreuChamp<-Breu[Breu$Origin=="Champ",]

#Make output matrix
output<-as.data.frame(array(NA, dim=c(4,23)))
colnames(output)<-c("Exp",	"Site",	"LitterOrigin",	"RMSE_single",	"RMSE_weibull",
                    "AIC_single",		"AIC_weibull",	"single.k",	"single.k.0.025",
                    "single.k.0.975",	"single.k.se","weibull.alpha",	"weibull.alpha.0.025",
                    "weibull.alpha.0.975", "weibull.alpha.se",	"weibull.beta",	"weibull.beta.0.025",	
                    "weibull.beta.0.975",	"weibull.beta.se",	"weibull.mrt",	"weibull.half.life",
                    "weibull.quarter.life",	"weibull.tenth.life")

half.life.calc=function(nls.mod){
  pars= coef(nls.mod)
  hl=pars[1] * (log(2))^(1/pars[2])
  names(hl)="half.life"
  return(hl)
}

quarter.life.calc=function(nls.mod){
  pars= coef(nls.mod)
  hl=pars[1] * (log(1/(3/4)))^(1/pars[2])
  names(hl)="half.life"
  return(hl)
}

tenth.life.calc=function(nls.mod){
  pars= coef(nls.mod)
  hl=pars[1] * (log(1/(9/10)))^(1/pars[2])
  names(hl)="half.life"
  return(hl)
}

mrt.calc=function(nls.mod){
  pars= coef(nls.mod)
  mrt=pars[1] * gamma(1+(1/pars[2]))
  names(mrt)="mrt"
  return(mrt)
}

LLS = function(y,k){
  Mhat=1*exp(-k*xNA$t) #creates a vector of=length to obs of preds
  ssq = sum((Mhat - xNA$Mt)^2)
  sigma = sqrt(ssq/length(xNA$Mt))
  return(-sum(dnorm(xNA$Mt,mean=Mhat,sd=sigma,log = TRUE)))
}

LLA = function(y,k,A){
  Mhat=A+((1-A)*exp(-k*xNA$t)) #creates a vector of=length to obs of preds
  ssq = sum((Mhat - xNA$Mt)^2)
  sigma = sqrt(ssq/length(xNA$Mt))
  return(-sum(dnorm(xNA$Mt,mean=Mhat,sd=sigma,log = TRUE)))
}

LLD = function(y,ks,k,A){
  Mhat=A*exp(-ks*xNA$t)+(1-A)*exp(-k*xNA$t) #creates a vector of=length to obs of preds
  ssq = sum((Mhat - xNA$Mt)^2)
  #browser()
  sigma = sqrt(ssq/length(xNA$Mt))
  return(-sum(dnorm(xNA$Mt,mean=Mhat,sd=sigma,log = TRUE)))
}


#################################################################################################
##Exp1; Champ Site/Champ Litter
i<-1
output[i,"Exp"]<-2
output[i,"Site"]<-"Champ"
output[i,"LitterOrigin"]<-"Champ"
s1=subset(Champ,Champ$Origin=="Champ")
t=(as.numeric(s1$time_years))
t=(as.numeric(s1$time_years))
class(t)
t1<-seq(0,3,length.out=10000)
Mt=s1$prop_mass_remaining
class(Mt)
x <- data.frame(t, Mt)
xNA<-na.exclude(x)
xNA
Mt<-xNA$Mt
t<-xNA$t

#Weibull function
Champ.fit<- nls((Mt) ~ exp(-(t/beta)^alpha), data=xNA, start =list(beta=1, alpha = .81), algorithm="port", lower=c(0.0001,0.0001))
summary(Champ.fit)

output[i,"weibull.alpha"]=coef(Champ.fit)[2]
output[i,"weibull.beta"]=coef(Champ.fit)[1]
output[i,"weibull.mrt"]<-mrt.calc(Champ.fit)
output[i,"weibull.half.life"]<-half.life.calc(Champ.fit)
output[i,"weibull.quarter.life"]<-quarter.life.calc(Champ.fit)
output[i,"weibull.tenth.life"]<-tenth.life.calc(Champ.fit)

weibull.fit<-exp(-(t1/output[i,"weibull.beta"])^output[i,"weibull.alpha"])
output[i,"AIC_weibull"]<-AIC(Champ.fit)
output[i,"weibull.beta.0.025"]<-confint(Champ.fit)[1]
output[i,"weibull.beta.0.975"]<-confint(Champ.fit)[3]
output[i,"weibull.beta.se"]<-summary(Champ.fit)$parameters[1,2]

output[i,"weibull.alpha.0.025"]<-confint(Champ.fit)[2]
output[i,"weibull.alpha.0.975"]<-confint(Champ.fit)[4]
output[i,"weibull.alpha.se"]<-summary(Champ.fit)$parameters[2,2]

#Single pool 
#Method 1: Liklihood function	
Champ.singleLL = mle2(minuslogl = LLS, start = list(k = 0.5), data = list(y=Mt),method="L-BFGS-B", lower=c(0),upper=c(1000))
output[i,"single.k"]=coef(Champ.singleLL)[1]
single.exp.fit<-exp(-output[i,"single.k"]*t1)
output[i,"AIC_single"]<-AIC(Champ.singleLL)

summary(Champ.singleLL)
confint(Champ.singleLL)
output[i,"single.k.0.025"]<-confint(Champ.singleLL)[1]
output[i,"single.k.0.975"]<-confint(Champ.singleLL)[2]
output[i,"single.k.se"]<-summary(Champ.fit)$parameters[2,2]

par(mfrow=c(2,2), oma=c(0,0,0,0), mar=c(4,4,1,1), mgp=c(2.4,0.8,0), cex.lab=1.25, cex.main=1.25, cex.axis=1.15)
plot(t,Mt, pch=16, col="dark grey", xlab="Yrs_Since_Deployment",ylab="Prop. C Remaining", main="Champ/Champ", xlim=c(0,2.5), ylim=c(0,1))
lines(t1,exp(-output[i,"single.k"]*t1), col="grey")
lines(t1, weibull.fit, col="dodgerblue1", lty=1)
legend("topright", c("single","weibull"), col=c("grey","dodgerblue1"), lty=c(1,1), bty="n")

xNA$single.predict<-exp(-output[i,"single.k"]*xNA$t)
output[i,"RMSE_single"]<-rmse(xNA$single.predict, Mt)
xNA$weibull.predict<-exp(-(xNA$t/output[i,"weibull.beta"])^output[i,"weibull.alpha"])
output[i,"RMSE_weibull"]<-rmse(xNA$weibull.predict, Mt)

#################################################################################################
##Exp1; Champ Site/Breu Litter
i<-2
output[i,"Exp"]<-2
output[i,"Site"]<-"Champ"
output[i,"LitterOrigin"]<-"Breu"
s1=subset(Champ,Champ$Origin=="Breu")
t=(as.numeric(s1$time_years))
t=(as.numeric(s1$time_years))
class(t)
t1<-seq(0,3,length.out=10000)
Mt=s1$prop_mass_remaining
class(Mt)
x <- data.frame(t, Mt)
xNA<-na.exclude(x)
xNA
Mt<-xNA$Mt
t<-xNA$t

#Weibull function
Champ.fit<- nls((Mt) ~ exp(-(t/beta)^alpha), data=xNA, start =list(beta=1, alpha = .81), algorithm="port", lower=c(0.0001,0.0001))
summary(Champ.fit)

output[i,"weibull.alpha"]=coef(Champ.fit)[2]
output[i,"weibull.beta"]=coef(Champ.fit)[1]
output[i,"weibull.mrt"]<-mrt.calc(Champ.fit)
output[i,"weibull.half.life"]<-half.life.calc(Champ.fit)
output[i,"weibull.quarter.life"]<-quarter.life.calc(Champ.fit)
output[i,"weibull.tenth.life"]<-tenth.life.calc(Champ.fit)

weibull.fit<-exp(-(t1/output[i,"weibull.beta"])^output[i,"weibull.alpha"])
output[i,"AIC_weibull"]<-AIC(Champ.fit)
output[i,"weibull.beta.0.025"]<-confint(Champ.fit)[1]
output[i,"weibull.beta.0.975"]<-confint(Champ.fit)[3]
output[i,"weibull.beta.se"]<-summary(Champ.fit)$parameters[1,2]

output[i,"weibull.alpha.0.025"]<-confint(Champ.fit)[2]
output[i,"weibull.alpha.0.975"]<-confint(Champ.fit)[4]
output[i,"weibull.alpha.se"]<-summary(Champ.fit)$parameters[2,2]

#Single pool 
#Method 1: Liklihood function	
Champ.singleLL = mle2(minuslogl = LLS, start = list(k = 0.5), data = list(y=Mt),method="L-BFGS-B", lower=c(0),upper=c(1000))
output[i,"single.k"]=coef(Champ.singleLL)[1]
single.exp.fit<-exp(-output[i,"single.k"]*t1)
output[i,"AIC_single"]<-AIC(Champ.singleLL)

summary(Champ.singleLL)
confint(Champ.singleLL)
output[i,"single.k.0.025"]<-confint(Champ.singleLL)[1]
output[i,"single.k.0.975"]<-confint(Champ.singleLL)[2]
output[i,"single.k.se"]<-summary(Champ.fit)$parameters[2,2]

plot(t,Mt, pch=16, col="dark grey", xlab="Yrs_Since_Deployment",ylab="Prop. C Remaining", main="Champ/Breu", xlim=c(0,2.5), ylim=c(0,1))
lines(t1,exp(-output[i,"single.k"]*t1), col="grey")
lines(t1, weibull.fit, col="dodgerblue1", lty=1)
legend("topright", c("single","weibull"), col=c("grey","dodgerblue1"), lty=c(1,1), bty="n")

xNA$single.predict<-exp(-output[i,"single.k"]*xNA$t)
output[i,"RMSE_single"]<-rmse(xNA$single.predict, Mt)
xNA$weibull.predict<-exp(-(xNA$t/output[i,"weibull.beta"])^output[i,"weibull.alpha"])
output[i,"RMSE_weibull"]<-rmse(xNA$weibull.predict, Mt)

#################################################################################################
##Exp1; Breu Site/Breu Litter
i<-3
output[i,"Exp"]<-2
output[i,"Site"]<-"Breu"
output[i,"LitterOrigin"]<-"Breu"
s1=subset(Breu,Breu$Origin=="Breu")
t=(as.numeric(s1$time_years))
t=(as.numeric(s1$time_years))
class(t)
t1<-seq(0,3,length.out=10000)
Mt=s1$prop_mass_remaining
class(Mt)
x <- data.frame(t, Mt)
xNA<-na.exclude(x)
xNA
Mt<-xNA$Mt
t<-xNA$t

#Weibull function
Champ.fit<- nls((Mt) ~ exp(-(t/beta)^alpha), data=xNA, start =list(beta=1, alpha = .81), algorithm="port", lower=c(0.0001,0.0001))
summary(Champ.fit)

output[i,"weibull.alpha"]=coef(Champ.fit)[2]
output[i,"weibull.beta"]=coef(Champ.fit)[1]
output[i,"weibull.mrt"]<-mrt.calc(Champ.fit)
output[i,"weibull.half.life"]<-half.life.calc(Champ.fit)
output[i,"weibull.quarter.life"]<-quarter.life.calc(Champ.fit)
output[i,"weibull.tenth.life"]<-tenth.life.calc(Champ.fit)

weibull.fit<-exp(-(t1/output[i,"weibull.beta"])^output[i,"weibull.alpha"])
output[i,"AIC_weibull"]<-AIC(Champ.fit)
output[i,"weibull.beta.0.025"]<-confint(Champ.fit)[1]
output[i,"weibull.beta.0.975"]<-confint(Champ.fit)[3]
output[i,"weibull.beta.se"]<-summary(Champ.fit)$parameters[1,2]

output[i,"weibull.alpha.0.025"]<-confint(Champ.fit)[2]
output[i,"weibull.alpha.0.975"]<-confint(Champ.fit)[4]
output[i,"weibull.alpha.se"]<-summary(Champ.fit)$parameters[2,2]

#Single pool 
#Method 1: Liklihood function	
Champ.singleLL = mle2(minuslogl = LLS, start = list(k = 0.5), data = list(y=Mt),method="L-BFGS-B", lower=c(0),upper=c(1000))
output[i,"single.k"]=coef(Champ.singleLL)[1]
single.exp.fit<-exp(-output[i,"single.k"]*t1)
output[i,"AIC_single"]<-AIC(Champ.singleLL)

summary(Champ.singleLL)
confint(Champ.singleLL)
output[i,"single.k.0.025"]<-confint(Champ.singleLL)[1]
output[i,"single.k.0.975"]<-confint(Champ.singleLL)[2]
output[i,"single.k.se"]<-summary(Champ.fit)$parameters[2,2]

plot(t,Mt, pch=16, col="dark grey", xlab="Yrs_Since_Deployment",ylab="Prop. C Remaining", main="Breu/Breu", xlim=c(0,2.5), ylim=c(0,1))
lines(t1,exp(-output[i,"single.k"]*t1), col="grey")
lines(t1, weibull.fit, col="dodgerblue1", lty=1)
legend("topright", c("single","weibull"), col=c("grey","dodgerblue1"), lty=c(1,1), bty="n")

xNA$single.predict<-exp(-output[i,"single.k"]*xNA$t)
output[i,"RMSE_single"]<-rmse(xNA$single.predict, Mt)
xNA$weibull.predict<-exp(-(xNA$t/output[i,"weibull.beta"])^output[i,"weibull.alpha"])
output[i,"RMSE_weibull"]<-rmse(xNA$weibull.predict, Mt)

################################################
#################################################################################################
##Exp1; Breu Site/Champ Litter
i<-4
output[i,"Exp"]<-2
output[i,"Site"]<-"Breu"
output[i,"LitterOrigin"]<-"Champ"
s1=subset(Breu,Breu$Origin=="Champ")
t=(as.numeric(s1$time_years))
t=(as.numeric(s1$time_years))
class(t)
t1<-seq(0,3,length.out=10000)
Mt=s1$prop_mass_remaining
class(Mt)
x <- data.frame(t, Mt)
xNA<-na.exclude(x)
xNA
Mt<-xNA$Mt
t<-xNA$t

#Weibull function
Champ.fit<- nls((Mt) ~ exp(-(t/beta)^alpha), data=xNA, start =list(beta=1, alpha = .81), algorithm="port", lower=c(0.0001,0.0001))
summary(Champ.fit)

output[i,"weibull.alpha"]=coef(Champ.fit)[2]
output[i,"weibull.beta"]=coef(Champ.fit)[1]
output[i,"weibull.mrt"]<-mrt.calc(Champ.fit)
output[i,"weibull.half.life"]<-half.life.calc(Champ.fit)
output[i,"weibull.quarter.life"]<-quarter.life.calc(Champ.fit)
output[i,"weibull.tenth.life"]<-tenth.life.calc(Champ.fit)

weibull.fit<-exp(-(t1/output[i,"weibull.beta"])^output[i,"weibull.alpha"])
output[i,"AIC_weibull"]<-AIC(Champ.fit)
output[i,"weibull.beta.0.025"]<-confint(Champ.fit)[1]
output[i,"weibull.beta.0.975"]<-confint(Champ.fit)[3]
output[i,"weibull.beta.se"]<-summary(Champ.fit)$parameters[1,2]

output[i,"weibull.alpha.0.025"]<-confint(Champ.fit)[2]
output[i,"weibull.alpha.0.975"]<-confint(Champ.fit)[4]
output[i,"weibull.alpha.se"]<-summary(Champ.fit)$parameters[2,2]

#Single pool 
#Method 1: Liklihood function	
Champ.singleLL = mle2(minuslogl = LLS, start = list(k = 0.5), data = list(y=Mt),method="L-BFGS-B", lower=c(0),upper=c(1000))
output[i,"single.k"]=coef(Champ.singleLL)[1]
single.exp.fit<-exp(-output[i,"single.k"]*t1)
output[i,"AIC_single"]<-AIC(Champ.singleLL)

summary(Champ.singleLL)
confint(Champ.singleLL)
output[i,"single.k.0.025"]<-confint(Champ.singleLL)[1]
output[i,"single.k.0.975"]<-confint(Champ.singleLL)[2]
output[i,"single.k.se"]<-summary(Champ.fit)$parameters[2,2]

plot(t,Mt, pch=16, col="dark grey", xlab="Yrs_Since_Deployment",ylab="Prop. C Remaining", main="Breu/Champ", xlim=c(0,2.5), ylim=c(0,1))
lines(t1,exp(-output[i,"single.k"]*t1), col="grey")
lines(t1, weibull.fit, col="dodgerblue1", lty=1)
legend("topright", c("single","weibull"), col=c("grey","dodgerblue1"), lty=c(1,1), bty="n")

xNA$single.predict<-exp(-output[i,"single.k"]*xNA$t)
output[i,"RMSE_single"]<-rmse(xNA$single.predict, Mt)
xNA$weibull.predict<-exp(-(xNA$t/output[i,"weibull.beta"])^output[i,"weibull.alpha"])
output[i,"RMSE_weibull"]<-rmse(xNA$weibull.predict, Mt)

##Export model parameters for exp 2
write.csv(output,"20200619_ModelParameters.csv")


##Bootstrap Appraoch ######
###Sampling Distributions
###########
pick1<-ChampChamp$time==118
fit2<-fitdist(ChampChamp[pick1, "prop_mass_remaining"], "weibull")
mat<-matrix(0,1000,1)
for(i in 1:nrow(mat)) {
  
  mat[i,1]<-mean(rweibull(length(ChampChamp[pick1, "prop_mass_remaining"]), shape=fit2$estimate[1], scale=fit2$estimate[2]))
}
ChampChamp118<-(mat[,1])

pick1<-ChampChamp$time==209
plot(fitdist(ChampChamp[pick1, "prop_mass_remaining"], "weibull"))
fit2<-fitdist(ChampChamp[pick1, "prop_mass_remaining"], "weibull")
mat<-matrix(0,1000,1)
for(i in 1:nrow(mat)) {
  
  mat[i,1]<-mean(rweibull(length(ChampChamp[pick1, "prop_mass_remaining"]), shape=fit2$estimate[1], scale=fit2$estimate[2]))
}
ChampChamp209<-(mat[,1])

pick1<-ChampChamp$time==271
plot(fitdist(ChampChamp[pick1, "prop_mass_remaining"], "weibull"))
fit2<-fitdist(ChampChamp[pick1, "prop_mass_remaining"], "weibull")
mat<-matrix(0,1000,1)
for(i in 1:nrow(mat)) {
  mat[i,1]<-mean(rweibull(length(ChampChamp[pick1, "prop_mass_remaining"]), shape=fit2$estimate[1], scale=fit2$estimate[2]))
}
ChampChamp271<-(mat[,1])

pick1<-ChampChamp$time==648
plot(fitdist(ChampChamp[pick1, "prop_mass_remaining"], "lnorm"))
fit2<-fitdist(ChampChamp[pick1, "prop_mass_remaining"], "lnorm")
mat<-matrix(0,1000,1)
for(i in 1:nrow(mat)) {
  mat[i,1]<-mean(rtrunc(length(ChampChamp[pick1, "prop_mass_remaining"]), spec="lnorm", a=0, b=Inf, mean=fit2$estimate[1], sd=fit2$estimate[2]))
}
ChampChamp648<-(mat[,1])

pick1<-ChampChamp$time==893
plot(fitdist(ChampChamp[pick1, "prop_mass_remaining"], "weibull"))
fit2<-fitdist(ChampChamp[pick1, "prop_mass_remaining"], "weibull")
mat<-matrix(0,1000,1)
for(i in 1:nrow(mat)) {
  mat[i,1]<-mean(rweibull(length(ChampChamp[pick1, "prop_mass_remaining"]), shape=fit2$estimate[1], scale=fit2$estimate[2]))
}
ChampChamp893<-(mat[,1])
ChampChamp0<-matrix(1,1000,1)
ChampChampBoot<-as.data.frame(cbind(ChampChamp0,ChampChamp118,ChampChamp209,ChampChamp271,ChampChamp648,ChampChamp893))
colnames(ChampChampBoot)<-c("ChampChamp0","ChampChamp118","ChampChamp209","ChampChamp271","ChampChamp648","ChampChamp893")
newt<-as.vector(c(0,118,209,271,648,893)/365)

mat<-matrix(0,1000,6)
colnames(mat)<-c("alpha", "beta","tenth.life", "quarter.life", "half.life", "mrt")
list(t(ChampChampBoot[i,]))
for(i in 1: nrow(ChampChampBoot)) {
  dat<-as.data.frame(cbind(t(ChampChampBoot[i,]), newt))
  colnames(dat)<-c("PropLitter", "Years")
  mat[i,1]<-coef(nls(PropLitter ~ exp(-(newt/beta)^alpha), data=dat, start =list(beta=1, alpha = .81), algorithm="port", lower=c(0.0001,0.0001)))[2]
  mat[i,2]<-coef(nls(PropLitter ~ exp(-(newt/beta)^alpha), data=dat, start =list(beta=1, alpha = .81), algorithm="port", lower=c(0.0001,0.0001)))[1]
  mat[i,3]<-mat[i,2] * (log(1/(9/10)))^(1/mat[i,1])
  mat[i,4]<-mat[i,2] * (log(1/(3/4)))^(1/mat[i,1])
  mat[i,5]<-mat[i,2] * (log(2))^(1/mat[i,1])
  mat[i,6]<-mat[i,2] * gamma(1+(1/mat[i,1]))
}
ChampChampParams<-as.data.frame(mat)
ChampChampParams$Site<-"Champ"
ChampChampParams$Origin<-"Champ"
ChampChampParams$Trt<-"Home"

##########################
########################
unique(ChampBreu$time)
head(data)
pick1<-ChampBreu$time==118
plot(fitdist(ChampBreu[pick1, "prop_mass_remaining"], "weibull"))
fit2<-fitdist(ChampBreu[pick1, "prop_mass_remaining"], "weibull")
mat<-matrix(0,1000,1)
for(i in 1:nrow(mat)) {
  
  mat[i,1]<-mean(rweibull(length(ChampBreu[pick1, "prop_mass_remaining"]), shape=fit2$estimate[1], scale=fit2$estimate[2]))
}
ChampBreu118<-(mat[,1])

pick1<-ChampBreu$time==209
plot(fitdist(ChampBreu[pick1, "prop_mass_remaining"], "weibull"))
fit2<-fitdist(ChampBreu[pick1, "prop_mass_remaining"], "weibull")
mat<-matrix(0,1000,1)
for(i in 1:nrow(mat)) {
  
  mat[i,1]<-mean(rweibull(length(ChampBreu[pick1, "prop_mass_remaining"]), shape=fit2$estimate[1], scale=fit2$estimate[2]))
}
ChampBreu209<-(mat[,1])

pick1<-ChampBreu$time==271
plot(fitdist(ChampBreu[pick1, "prop_mass_remaining"], "norm"))
fit2<-fitdist(ChampBreu[pick1, "prop_mass_remaining"], "norm")
mat<-matrix(0,1000,1)
for(i in 1:nrow(mat)) {
  mat[i,1]<-mean(rtrunc(length(ChampBreu[pick1, "prop_mass_remaining"]), spec="norm", a=0, b=Inf, mean=fit2$estimate[1], sd=fit2$estimate[2]))
}
ChampBreu271<-(mat[,1])

pick1<-ChampBreu$time==648
plot(fitdist(ChampBreu[pick1, "prop_mass_remaining"], "norm"))
fit2<-fitdist(ChampBreu[pick1, "prop_mass_remaining"], "norm")
mat<-matrix(0,1000,1)
for(i in 1:nrow(mat)) {
  mat[i,1]<-mean(rtrunc(length(ChampBreu[pick1, "prop_mass_remaining"]), spec="norm", a=0, b=Inf, mean=fit2$estimate[1], sd=fit2$estimate[2]))
}
ChampBreu648<-(mat[,1])

pick1<-ChampBreu$time==893
plot(fitdist(ChampBreu[pick1, "prop_mass_remaining"], "lnorm"))
fit2<-fitdist(ChampBreu[pick1, "prop_mass_remaining"], "lnorm")
mat<-matrix(0,1000,1)
for(i in 1:nrow(mat)) {
  mat[i,1]<-mean(rtrunc(length(ChampBreu[pick1, "prop_mass_remaining"]), spec="lnorm", a=0, b=Inf, mean=fit2$estimate[1], sd=fit2$estimate[2]))
}
ChampBreu893<-(mat[,1])

ChampBreu0<-matrix(1,1000,1)
ChampBreuBoot<-as.data.frame(cbind(ChampBreu0,ChampBreu118,ChampBreu209,ChampBreu271,ChampBreu648,ChampBreu893))
colnames(ChampBreuBoot)<-c("ChampBreu0","ChampBreu118","ChampBreu209","ChampBreu271","ChampBreu648","ChampBreu893")
newt<-as.vector(c(0,118,209,271,648,893)/365)

mat<-matrix(0,1000,6)
colnames(mat)<-c("alpha", "beta","tenth.life", "quarter.life", "half.life", "mrt")
for(i in 1: nrow(ChampBreuBoot)) {
  dat<-as.data.frame(cbind(t(ChampBreuBoot[i,]), newt))
  colnames(dat)<-c("PropLitter", "Years")
  mat[i,1]<-coef(nls(PropLitter ~ exp(-(newt/beta)^alpha), data=dat, start =list(beta=1, alpha = .81), algorithm="port", lower=c(0.0001,0.0001)))[2]
  mat[i,2]<-coef(nls(PropLitter ~ exp(-(newt/beta)^alpha), data=dat, start =list(beta=1, alpha = .81), algorithm="port", lower=c(0.0001,0.0001)))[1]
  mat[i,3]<-mat[i,2] * (log(1/(9/10)))^(1/mat[i,1])
  mat[i,4]<-mat[i,2] * (log(1/(3/4)))^(1/mat[i,1])
  mat[i,5]<-mat[i,2] * (log(2))^(1/mat[i,1])
  mat[i,6]<-mat[i,2] * gamma(1+(1/mat[i,1]))
  
}
ChampBreuParams<-as.data.frame(mat)
ChampBreuParams$Site<-"Champ"
ChampBreuParams$Origin<-"Breu"
ChampBreuParams$Trt<-"Away"

#########################################################################
########################
unique(BreuBreu$time)
head(data)
pick1<-BreuBreu$time==77
plot(fitdist(BreuBreu[pick1, "prop_mass_remaining"], "weibull"))
fit2<-fitdist(BreuBreu[pick1, "prop_mass_remaining"], "weibull")
mat<-matrix(0,1000,1)
for(i in 1:nrow(mat)) {
  
  mat[i,1]<-mean(rweibull(length(BreuBreu[pick1, "prop_mass_remaining"]), shape=fit2$estimate[1], scale=fit2$estimate[2]))
}
BreuBreu77<-(mat[,1])

pick1<-BreuBreu$time==271
plot(fitdist(BreuBreu[pick1, "prop_mass_remaining"], "weibull"))
fit2<-fitdist(BreuBreu[pick1, "prop_mass_remaining"], "weibull")
mat<-matrix(0,1000,1)
for(i in 1:nrow(mat)) {
  
  mat[i,1]<-mean(rweibull(length(BreuBreu[pick1, "prop_mass_remaining"]), shape=fit2$estimate[1], scale=fit2$estimate[2]))
}
BreuBreu271<-(mat[,1])

pick1<-BreuBreu$time==544
plot(fitdist(BreuBreu[pick1, "prop_mass_remaining"], "lnorm"))
fit2<-fitdist(BreuBreu[pick1, "prop_mass_remaining"], "lnorm")
mat<-matrix(0,1000,1)
for(i in 1:nrow(mat)) {
  mat[i,1]<-mean(rtrunc(length(BreuBreu[pick1, "prop_mass_remaining"]), spec="lnorm", a=0, b=Inf, mean=fit2$estimate[1], sd=fit2$estimate[2]))
}
BreuBreu544<-(mat[,1])

pick1<-BreuBreu$time==903
plot(fitdist(BreuBreu[pick1, "prop_mass_remaining"], "weibull"))
fit2<-fitdist(BreuBreu[pick1, "prop_mass_remaining"], "weibull")
mat<-matrix(0,1000,1)
for(i in 1:nrow(mat)) {
  mat[i,1]<-mean(rweibull(length(BreuBreu[pick1, "prop_mass_remaining"]), shape=fit2$estimate[1], scale=fit2$estimate[2]))
}
BreuBreu903<-(mat[,1])

BreuBreu0<-matrix(1,1000,1)
BreuBreuBoot<-as.data.frame(cbind(BreuBreu0,BreuBreu77,BreuBreu271,BreuBreu544,BreuBreu903))
colnames(BreuBreuBoot)<-c("BreuBreu0","BreuBreu77","BreuBreu271","BreuBreu544","BreuBreu903")
newt<-as.vector(c(0,77,271,544,903)/365)

mat<-matrix(0,1000,6)
colnames(mat)<-c("alpha", "beta","tenth.life", "quarter.life", "half.life", "mrt")
for(i in 1: nrow(BreuBreuBoot)) {
  dat<-as.data.frame(cbind(t(BreuBreuBoot[i,]), newt))
  colnames(dat)<-c("PropLitter", "Years")
  mat[i,1]<-coef(nls(PropLitter ~ exp(-(newt/beta)^alpha), data=dat, start =list(beta=1, alpha = .81), algorithm="port", lower=c(0.0001,0.0001)))[2]
  mat[i,2]<-coef(nls(PropLitter ~ exp(-(newt/beta)^alpha), data=dat, start =list(beta=1, alpha = .81), algorithm="port", lower=c(0.0001,0.0001)))[1]
  mat[i,3]<-mat[i,2] * (log(1/(9/10)))^(1/mat[i,1])
  mat[i,4]<-mat[i,2] * (log(1/(3/4)))^(1/mat[i,1])
  mat[i,5]<-mat[i,2] * (log(2))^(1/mat[i,1])
  mat[i,6]<-mat[i,2] * gamma(1+(1/mat[i,1]))
}
BreuBreuParams<-as.data.frame(mat)
BreuBreuParams$Site<-"Breu"
BreuBreuParams$Origin<-"Breu"
BreuBreuParams$Trt<-"Home"

########################
unique(BreuChamp$time)
head(data)
pick1<-BreuChamp$time==77
plot(fitdist(BreuChamp[pick1, "prop_mass_remaining"], "norm"))
fit2<-fitdist(BreuChamp[pick1, "prop_mass_remaining"], "norm")
mat<-matrix(0,1000,1)
for(i in 1:nrow(mat)) {
  
  mat[i,1]<-mean(rtrunc(length(BreuChamp[pick1, "prop_mass_remaining"]), spec="norm", a=0, b=Inf, mean=fit2$estimate[1], sd=fit2$estimate[2]))
}
BreuChamp77<-(mat[,1])

pick1<-BreuChamp$time==271
plot(fitdist(BreuChamp[pick1, "prop_mass_remaining"], "weibull"))
fit2<-fitdist(BreuChamp[pick1, "prop_mass_remaining"], "weibull")
mat<-matrix(0,1000,1)
for(i in 1:nrow(mat)) {
  
  mat[i,1]<-mean(rweibull(length(BreuChamp[pick1, "prop_mass_remaining"]), shape=fit2$estimate[1], scale=fit2$estimate[2]))
}
BreuChamp271<-(mat[,1])

pick1<-BreuChamp$time==544
plot(fitdist(BreuChamp[pick1, "prop_mass_remaining"], "lnorm"))
fit2<-fitdist(BreuChamp[pick1, "prop_mass_remaining"], "lnorm")
mat<-matrix(0,1000,1)
for(i in 1:nrow(mat)) {
  mat[i,1]<-mean(rtrunc(length(BreuChamp[pick1, "prop_mass_remaining"]), spec="lnorm", a=0, b=Inf, mean=fit2$estimate[1], sd=fit2$estimate[2]))
}
BreuChamp544<-(mat[,1])

pick1<-BreuChamp$time==903
plot(fitdist(BreuChamp[pick1, "prop_mass_remaining"], "norm"))
fit2<-fitdist(BreuChamp[pick1, "prop_mass_remaining"], "norm")
mat<-matrix(0,1000,1)
for(i in 1:nrow(mat)) {
  mat[i,1]<-mean(rtrunc(length(BreuChamp[pick1, "prop_mass_remaining"]), spec="norm", a=0, b=Inf, mean=fit2$estimate[1], sd=fit2$estimate[2]))
}
BreuChamp903<-(mat[,1])

BreuChamp0<-matrix(1,1000,1)
BreuChampBoot<-as.data.frame(cbind(BreuChamp0,BreuChamp77,BreuChamp271,BreuChamp544,BreuChamp903))
colnames(BreuChampBoot)<-c("BreuChamp0","BreuChamp77","BreuChamp271","BreuChamp544","BreuChamp903")
newt<-as.vector(c(0,77,271,544,903)/365)

mat<-matrix(0,1000,6)
colnames(mat)<-c("alpha", "beta","tenth.life", "quarter.life", "half.life", "mrt")
for(i in 1: nrow(BreuChampBoot)) {
  dat<-as.data.frame(cbind(t(BreuChampBoot[i,]), newt))
  colnames(dat)<-c("PropLitter", "Years")
  mat[i,1]<-coef(nls(PropLitter ~ exp(-(newt/beta)^alpha), data=dat, start =list(beta=1, alpha = .81), algorithm="port", lower=c(0.0001,0.0001)))[2]
  mat[i,2]<-coef(nls(PropLitter ~ exp(-(newt/beta)^alpha), data=dat, start =list(beta=1, alpha = .81), algorithm="port", lower=c(0.0001,0.0001)))[1]
  mat[i,3]<-mat[i,2] * (log(1/(9/10)))^(1/mat[i,1])
  mat[i,4]<-mat[i,2] * (log(1/(3/4)))^(1/mat[i,1])
  mat[i,5]<-mat[i,2] * (log(2))^(1/mat[i,1])
  mat[i,6]<-mat[i,2] * gamma(1+(1/mat[i,1]))
}

BreuChampParams<-as.data.frame(mat)
BreuChampParams$Site<-"Breu"
BreuChampParams$Origin<-"Champ"
BreuChampParams$Trt<-"Away"

#Matrix of bootstrapped parameter fits
AllParams<-rbind(ChampChampParams, ChampBreuParams, BreuBreuParams, BreuChampParams)
write.csv(AllParams, "20200620_Exp2_BootstrapParameters.csv")


#######
dev.off()
par(mfrow=c(2,2), mar=c(1, 5.1, 2.1, 1),cex.axis=1, cex.lab=1)
position.text<-c(1,2,3,4)
plot(1, type="n", axes=F, xlab="", ylab="", xlim=c(0.5,4.35), ylim=c(0.5,1.75))
rect(1.5, -1, 2.5, 2.5, col= "light grey", border=NA)
rect(3.5, -1, 4.5, 2.5, col= "light grey", border=NA)
axis(side=2, tck=-0.02, labels=NA)
axis(side=2, lwd=0, line = -.4, las=1, cex.axis=1.25)
mtext(side=2, line=2.5, "Weibull alpha", cex=1)
vioplot(ChampChampParams$alpha, BreuChampParams$alpha,  BreuBreuParams$alpha,ChampBreuParams$alpha, 
        col=c(rgb(0.59,0.8,0.8), rgb(0.59,0.8,0.8),rgb(0.933,0.91,0.667),rgb(0.933,0.91,0.667)),
        main="Weibull alpha", border=NA, add=TRUE)
legend("topleft", c("Champ Litter", "Breu Litter"), pch=c(16,16), col=c(rgb(0.59,0.8,0.8),rgb(0.933,0.91,0.667)),
       bty="n")
legend("bottomleft", c("Home", "Away"), pch=c(0,15), col=c("black", "light grey"),
       bty="n")

# par(mfrow=c(1,1), mar=c(5, 5.1, 2.1, 3.1),cex.axis=1, cex.lab=1)
position.text<-c(1,2,3,4)
plot(1, type="n", axes=F, xlab="", ylab="", xlim=c(0.5,4.35), ylim=c(1.5,2.8))
rect(1.5, -1, 2.5, 3, col= "light grey", border=NA)
rect(3.5, -1, 4.5, 3, col= "light grey", border=NA)
mtext(side=2, line=2.5, "Weibull MRT", cex=1)
axis(side=2, tck=-0.02, labels=NA)
axis(side=2, lwd=0, line = -.4, las=1, cex.axis=1.25)
vioplot(ChampChampParams$mrt, BreuChampParams$mrt,  BreuBreuParams$mrt,ChampBreuParams$mrt, 
        col=c(rgb(0.59,0.8,0.8), rgb(0.59,0.8,0.8),rgb(0.933,0.91,0.667),rgb(0.933,0.91,0.667)),
        main="Weibull MRT", border=NA, add=T)
legend("topleft", c("Champ Litter", "Breu Litter"), pch=c(16,16), col=c(rgb(0.59,0.8,0.8),rgb(0.933,0.91,0.667)),
       bty="n")
legend("bottomleft", c("Home", "Away"), pch=c(0,15), col=c("black", "light grey"),
       bty="n")

plot(1, type="n", axes=F, xlab="", ylab="", xlim=c(0.5,4.35), ylim=c(0.1,0.6))
rect(1.5, -1, 2.5, 2.5, col= "light grey", border=NA)
rect(3.5, -1, 4.5, 2.5, col= "light grey", border=NA)
axis(side=2, tck=-0.02, labels=NA)
axis(side=2, lwd=0, line = -.4, las=1, cex.axis=1.25)
mtext(side=2, line=2.5, expression('Weibull '*italic(t[1/10])*''), cex=1)
vioplot(ChampChampParams$tenth.life, BreuChampParams$tenth.life,  BreuBreuParams$tenth.life,ChampBreuParams$tenth.life, 
        col=c(rgb(0.59,0.8,0.8), rgb(0.59,0.8,0.8),rgb(0.933,0.91,0.667),rgb(0.933,0.91,0.667)),
        main=expression('Weibull '*italic(t[1/10])*''), border=NA, add=T)
legend("topleft", c("Champ Litter", "Breu Litter"), pch=c(16,16), col=c(rgb(0.59,0.8,0.8),rgb(0.933,0.91,0.667)),
       bty="n")
legend("bottomleft", c("Home", "Away"), pch=c(0,15), col=c("black", "light grey"),
       bty="n")

# par(mfrow=c(1,1), mar=c(5, 5.1, 2.1, 3.1),cex.axis=1, cex.lab=1)
position.text<-c(1,2,3,4)
plot(1, type="n", axes=F, xlab="", ylab="", xlim=c(0.5,4.35), ylim=c(0.5,1.1))
rect(1.5, -1, 2.5, 3, col= "light grey", border=NA)
rect(3.5, -1, 4.5, 3, col= "light grey", border=NA)
mtext(side=2, line=2.5, expression('Weibull '*italic(t[1/4])*''), cex=1)
axis(side=2, tck=-0.02, labels=NA)
axis(side=2, lwd=0, line = -.4, las=1, cex.axis=1.25)
vioplot(ChampChampParams$quarter.life, BreuChampParams$quarter.life,  BreuBreuParams$quarter.life,ChampBreuParams$quarter.life, 
        col=c(rgb(0.59,0.8,0.8), rgb(0.59,0.8,0.8),rgb(0.933,0.91,0.667),rgb(0.933,0.91,0.667)),
        main=expression('Weibull '*italic(t[1/4])*''), border=NA, add=T)
legend("topleft", c("Champ Litter", "Breu Litter"), pch=c(16,16), col=c(rgb(0.59,0.8,0.8),rgb(0.933,0.91,0.667)),
       bty="n")
legend("bottomleft", c("Home", "Away"), pch=c(0,15), col=c("black", "light grey"),
       bty="n")


##################################
t1<-seq(0,3,length.out=1000)

ChampChamp_alpha_mean<- mean(ChampChampParams$alpha)
ChampChamp_alpha_upper<- quantile(ChampChampParams$alpha, 0.975)
ChampChamp_alpha_lower<- quantile(ChampChampParams$alpha, 0.025)
ChampChamp_beta_mean<- mean(ChampChampParams$beta)
ChampChamp_beta_upper<- quantile(ChampChampParams$beta, 0.975)
ChampChamp_beta_lower<- quantile(ChampChampParams$beta, 0.025)

ChampChamp_mean<- exp(-(t1/ChampChamp_beta_mean)^ChampChamp_alpha_mean)
ChampChamp_975<- exp(-(t1/ChampChamp_beta_upper)^ChampChamp_alpha_upper)
ChampChamp_025<- exp(-(t1/ChampChamp_beta_lower)^ChampChamp_alpha_lower)
#
ChampBreu_alpha_mean<- mean(ChampBreuParams$alpha)
ChampBreu_alpha_upper<- quantile(ChampBreuParams$alpha, 0.975)
ChampBreu_alpha_lower<- quantile(ChampBreuParams$alpha, 0.025)
ChampBreu_beta_mean<- mean(ChampBreuParams$beta)
ChampBreu_beta_upper<- quantile(ChampBreuParams$beta, 0.975)
ChampBreu_beta_lower<- quantile(ChampBreuParams$beta, 0.025)

ChampBreu_mean<- exp(-(t1/ChampBreu_beta_mean)^ChampBreu_alpha_mean)
ChampBreu_975<- exp(-(t1/ChampBreu_beta_upper)^ChampBreu_alpha_upper)
ChampBreu_025<- exp(-(t1/ChampBreu_beta_lower)^ChampBreu_alpha_lower)
#
BreuBreu_alpha_mean<- mean(BreuBreuParams$alpha)
BreuBreu_alpha_upper<- quantile(BreuBreuParams$alpha, 0.975)
BreuBreu_alpha_lower<- quantile(BreuBreuParams$alpha, 0.025)
BreuBreu_beta_mean<- mean(BreuBreuParams$beta)
BreuBreu_beta_upper<- quantile(BreuBreuParams$beta, 0.975)
BreuBreu_beta_lower<- quantile(BreuBreuParams$beta, 0.025)

BreuBreu_mean<- exp(-(t1/BreuBreu_beta_mean)^BreuBreu_alpha_mean)
BreuBreu_975<- exp(-(t1/BreuBreu_beta_upper)^BreuBreu_alpha_upper)
BreuBreu_025<- exp(-(t1/BreuBreu_beta_lower)^BreuBreu_alpha_lower)
#
BreuChamp_alpha_mean<- mean(BreuChampParams$alpha)
BreuChamp_alpha_upper<- quantile(BreuChampParams$alpha, 0.975)
BreuChamp_alpha_lower<- quantile(BreuChampParams$alpha, 0.025)
BreuChamp_beta_mean<- mean(BreuChampParams$beta)
BreuChamp_beta_upper<- quantile(BreuChampParams$beta, 0.975)
BreuChamp_beta_lower<- quantile(BreuChampParams$beta, 0.025)

BreuChamp_mean<- exp(-(t1/BreuChamp_beta_mean)^BreuChamp_alpha_mean)
BreuChamp_975<- exp(-(t1/BreuChamp_beta_upper)^BreuChamp_alpha_upper)
BreuChamp_025<- exp(-(t1/BreuChamp_beta_lower)^BreuChamp_alpha_lower)

ChampChamp=subset(Champ,Champ$Origin=="Champ")
ChampBreu=subset(Champ,Champ$Origin=="Breu")
BreuBreu=subset(Breu,Breu$Origin=="Breu")
BreuChamp=subset(Breu,Breu$Origin=="Champ")

par(mfrow=c(1,1), mgp=c(2.5,1,0), mar=c(4,4,3,1))
plot(t1, ChampChamp_mean, type="l", lwd=3, col="paleturquoise3", ylim=c(0,1.2), xlab="Years",
     ylab="Proportion Litter Mass Remaining", main="Exp. 2\n Weibull Model Bootstrap 95% CI")
points(ChampChamp$time_years, ChampChamp$prop_mass_remaining, col="paleturquoise3")
polygon(c(t1[2:1000],rev(t1[2:1000])), c((ChampChamp_975[2:1000]),rev((ChampChamp_025[2:1000]))),
        col= rgb(0.59,0.8,0.8, alpha=0.3), border=F)

points(t1, ChampBreu_mean, type="l", lwd=3, col="darkgoldenrod", lty=2)
points(ChampBreu$time_years, ChampBreu$prop_mass_remaining, col="darkgoldenrod", pch=2)
polygon(c(t1[2:1000],rev(t1[2:1000])), c((ChampBreu_975[2:1000]),rev((ChampBreu_025[2:1000]))),
        col= rgb(0.933,0.91,0.667, alpha=0.3), border=F)

points(t1, BreuChamp_mean, type="l", lwd=3, lty=2, col="paleturquoise3")
points(BreuChamp$time_years, BreuChamp$prop_mass_remaining, col="paleturquoise3", pch=2)
polygon(c(t1[2:1000],rev(t1[2:1000])), c((BreuChamp_975[2:1000]),rev((BreuChamp_025[2:1000]))),
        col= rgb(0.59,0.8,0.8, alpha=0.3), border=F)

points(t1, BreuBreu_mean, type="l", lwd=3, col="darkgoldenrod")
points(BreuBreu$time_years, BreuBreu$prop_mass_remaining, col="darkgoldenrod")
polygon(c(t1[2:1000],rev(t1[2:1000])), c((BreuBreu_975[2:1000]),rev((BreuBreu_025[2:1000]))),
        col= rgb(0.933,0.91,0.667, alpha=0.3), border=F)
legend("topright", c("Champ Litter", "Breu Litter"), col=c("paleturquoise3", "darkgoldenrod"), pch=16, bty="n")
legend("topleft", c("Home", "Away"), lty=c(1,2), pch=c(1,2), bty="n")



#######################################
## Output for bootstrap model summaries
new.output<-as.data.frame(array(NA, dim=c(4,18)))
colnames(new.output)<-c("alpha_mean", "alpha_0.975", "alpha_0.025", "beta_mean", "beta_0.975", "beta_0.025",
                        "tenth.life_mean", "tenth.life_0.975", "tenth.life_0.025", "quarter.life_mean", "quarter.life_0.975", "quarter.life_0.025",
                        "half.life_mean", "half.life_0.975", "half.life_0.025", "mrt_mean", "mrt_0.975", "mrt_0.025")
rownames(new.output)<-c("ChampChamp", "ChampBreu", "BreuBreu", "BreuChamp")
new.output[1,"alpha_mean"]<- mean(ChampChampParams$alpha)
new.output[1,"alpha_0.975"]<- quantile(ChampChampParams$alpha, 0.975)
new.output[1,"alpha_0.025"]<- quantile(ChampChampParams$alpha, 0.025)
new.output[1,"beta_mean"]<- mean(ChampChampParams$beta)
new.output[1,"beta_0.975"]<- quantile(ChampChampParams$beta, 0.975)
new.output[1,"beta_0.025"]<- quantile(ChampChampParams$beta, 0.025)
new.output[1,"tenth.life_mean"]<- mean(ChampChampParams$tenth.life)
new.output[1,"tenth.life_0.975"]<- quantile(ChampChampParams$tenth.life, 0.975)
new.output[1,"tenth.life_0.025"]<- quantile(ChampChampParams$tenth.life, 0.025)
new.output[1,"quarter.life_mean"]<- mean(ChampChampParams$quarter.life)
new.output[1,"quarter.life_0.975"]<- quantile(ChampChampParams$quarter.life, 0.975)
new.output[1,"quarter.life_0.025"]<- quantile(ChampChampParams$quarter.life, 0.025)
new.output[1,"half.life_mean"]<- mean(ChampChampParams$half.life)
new.output[1,"half.life_0.975"]<- quantile(ChampChampParams$half.life, 0.975)
new.output[1,"half.life_0.025"]<- quantile(ChampChampParams$half.life, 0.025)
new.output[1,"mrt_mean"]<- mean(ChampChampParams$mrt)
new.output[1,"mrt_0.975"]<- quantile(ChampChampParams$mrt, 0.975)
new.output[1,"mrt_0.025"]<- quantile(ChampChampParams$mrt, 0.025)

new.output[2,"alpha_mean"]<- mean(ChampBreuParams$alpha)
new.output[2,"alpha_0.975"]<- quantile(ChampBreuParams$alpha, 0.975)
new.output[2,"alpha_0.025"]<- quantile(ChampBreuParams$alpha, 0.025)
new.output[2,"beta_mean"]<- mean(ChampBreuParams$beta)
new.output[2,"beta_0.975"]<- quantile(ChampBreuParams$beta, 0.975)
new.output[2,"beta_0.025"]<- quantile(ChampBreuParams$beta, 0.025)
new.output[2,"tenth.life_mean"]<- mean(ChampBreuParams$tenth.life)
new.output[2,"tenth.life_0.975"]<- quantile(ChampBreuParams$tenth.life, 0.975)
new.output[2,"tenth.life_0.025"]<- quantile(ChampBreuParams$tenth.life, 0.025)
new.output[2,"quarter.life_mean"]<- mean(ChampBreuParams$quarter.life)
new.output[2,"quarter.life_0.975"]<- quantile(ChampBreuParams$quarter.life, 0.975)
new.output[2,"quarter.life_0.025"]<- quantile(ChampBreuParams$quarter.life, 0.025)
new.output[2,"half.life_mean"]<- mean(ChampBreuParams$half.life)
new.output[2,"half.life_0.975"]<- quantile(ChampBreuParams$half.life, 0.975)
new.output[2,"half.life_0.025"]<- quantile(ChampBreuParams$half.life, 0.025)
new.output[2,"mrt_mean"]<- mean(ChampBreuParams$mrt)
new.output[2,"mrt_0.975"]<- quantile(ChampBreuParams$mrt, 0.975)
new.output[2,"mrt_0.025"]<- quantile(ChampBreuParams$mrt, 0.025)

new.output[3,"alpha_mean"]<- mean(BreuBreuParams$alpha)
new.output[3,"alpha_0.975"]<- quantile(BreuBreuParams$alpha, 0.975)
new.output[3,"alpha_0.025"]<- quantile(BreuBreuParams$alpha, 0.025)
new.output[3,"beta_mean"]<- mean(BreuBreuParams$beta)
new.output[3,"beta_0.975"]<- quantile(BreuBreuParams$beta, 0.975)
new.output[3,"beta_0.025"]<- quantile(BreuBreuParams$beta, 0.025)
new.output[3,"tenth.life_mean"]<- mean(BreuBreuParams$tenth.life)
new.output[3,"tenth.life_0.975"]<- quantile(BreuBreuParams$tenth.life, 0.975)
new.output[3,"tenth.life_0.025"]<- quantile(BreuBreuParams$tenth.life, 0.025)
new.output[3,"quarter.life_mean"]<- mean(BreuBreuParams$quarter.life)
new.output[3,"quarter.life_0.975"]<- quantile(BreuBreuParams$quarter.life, 0.975)
new.output[3,"quarter.life_0.025"]<- quantile(BreuBreuParams$quarter.life, 0.025)
new.output[3,"half.life_mean"]<- mean(BreuBreuParams$half.life)
new.output[3,"half.life_0.975"]<- quantile(BreuBreuParams$half.life, 0.975)
new.output[3,"half.life_0.025"]<- quantile(BreuBreuParams$half.life, 0.025)
new.output[3,"mrt_mean"]<- mean(BreuBreuParams$mrt)
new.output[3,"mrt_0.975"]<- quantile(BreuBreuParams$mrt, 0.975)
new.output[3,"mrt_0.025"]<- quantile(BreuBreuParams$mrt, 0.025)

new.output[4,"alpha_mean"]<- mean(BreuChampParams$alpha)
new.output[4,"alpha_0.975"]<- quantile(BreuChampParams$alpha, 0.975)
new.output[4,"alpha_0.025"]<- quantile(BreuChampParams$alpha, 0.025)
new.output[4,"beta_mean"]<- mean(BreuChampParams$beta)
new.output[4,"beta_0.975"]<- quantile(BreuChampParams$beta, 0.975)
new.output[4,"beta_0.025"]<- quantile(BreuChampParams$beta, 0.025)
new.output[4,"tenth.life_mean"]<- mean(BreuChampParams$tenth.life)
new.output[4,"tenth.life_0.975"]<- quantile(BreuChampParams$tenth.life, 0.975)
new.output[4,"tenth.life_0.025"]<- quantile(BreuChampParams$tenth.life, 0.025)
new.output[4,"quarter.life_mean"]<- mean(BreuChampParams$quarter.life)
new.output[4,"quarter.life_0.975"]<- quantile(BreuChampParams$quarter.life, 0.975)
new.output[4,"quarter.life_0.025"]<- quantile(BreuChampParams$quarter.life, 0.025)
new.output[4,"half.life_mean"]<- mean(BreuChampParams$half.life)
new.output[4,"half.life_0.975"]<- quantile(BreuChampParams$half.life, 0.975)
new.output[4,"half.life_0.025"]<- quantile(BreuChampParams$half.life, 0.025)
new.output[4,"mrt_mean"]<- mean(BreuChampParams$mrt)
new.output[4,"mrt_0.975"]<- quantile(BreuChampParams$mrt, 0.975)
new.output[4,"mrt_0.025"]<- quantile(BreuChampParams$mrt, 0.025)

write.csv(new.output,"20200619_Exp2_Bootstrap_parameters.csv")

ChampChamp_beta_mean<- mean(ChampChampParams$beta)
ChampChamp_beta_upper<- quantile(ChampChampParams$beta, 0.975)
ChampChamp_beta_lower<- quantile(ChampChampParams$beta, 0.025)


