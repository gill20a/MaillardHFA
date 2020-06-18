# Code to fit weibull and single exponential decomposition models to Maillard HFA litter 
## decomposition experiment.

#Open necessary libraries
library(stats4)
library(bbmle)
library(readr)
library(mltools)
library(nlstools)
library(RCurl)

#Download data file. 

urlfile = "https://raw.githubusercontent.com/gill20a/MaillardHFA/master/MaillardHFA_Exp1_masslossdata.csv?token=ABRRWDQMY2RFLO7AFQ7NLYC65PKPU"
data<-read.csv(url(urlfile))
head(data)

#Convert days to years
data$time_years<-data$time/365

#Convert percent mass remaining to fraction
data$prop_mass_remaining<-data$mass_remaining/100

#divide dataset by site
Champ<-data[data$Site=="Champ",]

#Make output matrix
output<-as.data.frame(array(NA, dim=c(2,23)))
colnames(output)<-c("Exp",	"Site",	"LitterOrigin",	"RMSE_single",	"RMSE_weibull",
                    "AIC_single",		"AIC_weibull",	"single.k",	"single.k.0.025",
                    "single.k.0.975",	"single.k.se","weibull.alpha",	"weibull.alpha.0.025",
                    "weibull.alpha.0.975", "weibull.alpha.se",	"weibull.beta",	"weibull.beta.0.025",	
                    "weibull.beta.0.975",	"weibull.beta.se",	"weibull.mrt",	"weibull.half.life",
                    "weibull.quarter.life",	"weibull.tenth.life")


#Here, models are fit one at a time. You could make this a loop if you 
#are convinced all your models will converge using the same starting values.
#I never have such luck.


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


#################################################################################################
##Exp1; Champ Site/Champ Litter
i<-1
output[i,"Exp"]<-1
output[i,"Site"]<-"Champ"
output[i,"LitterOrigin"]<-"Champ"
s1=subset(Champ,Champ$Origin=="Champ")
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
output[i,"Exp"]<-1
output[i,"Site"]<-"Champ"
output[i,"LitterOrigin"]<-"Breu"
s1=subset(Champ,Champ$Origin=="Breu")
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

# par(mfrow=c(2,2), oma=c(0,0,0,0), mar=c(4,4,1,1), mgp=c(2.4,0.8,0), cex.lab=1.25, cex.main=1.25, cex.axis=1.15)
plot(t,Mt, pch=16, col="dark grey", xlab="Yrs_Since_Deployment",ylab="Prop. C Remaining", main="Champ/Breu", xlim=c(0,2.5), ylim=c(0,1))
lines(t1,exp(-output[i,"single.k"]*t1), col="grey")
lines(t1, weibull.fit, col="dodgerblue1", lty=1)
legend("topright", c("single","weibull"), col=c("grey","dodgerblue1"), lty=c(1,1), bty="n")

xNA$single.predict<-exp(-output[i,"single.k"]*xNA$t)
output[i,"RMSE_single"]<-rmse(xNA$single.predict, Mt)
xNA$weibull.predict<-exp(-(xNA$t/output[i,"weibull.beta"])^output[i,"weibull.alpha"])
output[i,"RMSE_weibull"]<-rmse(xNA$weibull.predict, Mt)

###This file contains summarized model parameters for the two litter origin treatments in 
#experiment one
output

#################################################################################################
##Need to run curve fitting and export output file on Exp 2 before running this part
para<-read.csv("20200506_ModelParameters.csv", header=T, sep=",")
head(para)

s1=subset(Champ,Champ$Origin=="Champ")
t=(as.numeric(s1$time_years))
class(t)
t1<-seq(0,3,length.out=1000)
Mt=s1$prop_mass_remaining
class(Mt)
x <- data.frame(t, Mt)
Champ.xNA<-na.exclude(x)


Champ.Mass.mean<- exp(-(t1/para[1,"weibull.beta"])^para[1,"weibull.alpha"])
Champ.Mass.975<- exp(-(t1/para[1,"weibull.beta.0.975"])^para[1,"weibull.alpha.0.975"])
Champ.Mass.025<- exp(-(t1/para[1,"weibull.beta.0.025"])^para[1,"weibull.alpha.0.025"])

Champ.Asymp.Mass.mean<- para[1,"asymp.A"]+((1-para[1,"asymp.A"])*exp(-para[1,"asymp.k"]*t1))
Champ.Asymp.Mass.975<- para[1,"asymp.A.0.975"]+((1-para[1,"asymp.A.0.975"])*exp(-para[1,"asymp.k.0.975"]*t1))
Champ.Asymp.Mass.025<- para[1,"asymp.A.0.025"]+((1-para[1,"asymp.A.0.025"])*exp(-para[1,"asymp.k.0.025"]*t1))

##
s1=subset(Champ,Champ$Origin=="Breu")
t=(as.numeric(s1$time_years))
class(t)
t1<-seq(0,3,length.out=1000)
Mt=s1$prop_mass_remaining
class(Mt)
x <- data.frame(t, Mt)
Breu.xNA<-na.exclude(x)


Breu.Mass.mean<- exp(-(t1/para[2,"weibull.beta"])^para[2,"weibull.alpha"])
Breu.Mass.975<- exp(-(t1/para[2,"weibull.beta.0.975"])^para[2,"weibull.alpha.0.975"])
Breu.Mass.025<- exp(-(t1/para[2,"weibull.beta.0.025"])^para[2,"weibull.alpha.0.025"])

Breu.Asymp.Mass.mean<- para[2,"asymp.A"]+((1-para[2,"asymp.A"])*exp(-para[2,"asymp.k"]*t1))
Breu.Asymp.Mass.975<- para[2,"asymp.A.0.975"]+((1-para[2,"asymp.A.0.975"])*exp(-para[2,"asymp.k.0.975"]*t1))
Breu.Asymp.Mass.025<- para[2,"asymp.A.0.025"]+((1-para[2,"asymp.A.0.025"])*exp(-para[2,"asymp.k.0.025"]*t1))

par(mfrow=c(1,2), mgp=c(2.5,1,0), mar=c(4,4,3,1))
##
plot(t1, Champ.Mass.mean, type="l", lwd=3, col="paleturquoise3", ylim=c(0,1.2), xlab="Years",
     ylab="Proportion Litter Mass Remaining", main="Exp. 1 Weibull Model \n95% CI")
points(Champ.xNA$t, Champ.xNA$Mt, col="paleturquoise3")
polygon(c(t1[2:1000],rev(t1[2:1000])), c((Champ.Mass.975[2:1000]),rev((Champ.Mass.025[2:1000]))),
        col= rgb(0.59,0.8,0.8, alpha=0.3), border=F)

points(t1, Breu.Mass.mean, type="l", lwd=3, col="darkgoldenrod")
points(Breu.xNA$t, Breu.xNA$Mt, col="darkgoldenrod")
polygon(c(t1[2:1000],rev(t1[2:1000])), c((Breu.Mass.975[2:1000]),rev((Breu.Mass.025[2:1000]))),
        col= rgb(0.933,0.91,0.667, alpha=0.3), border=F)
legend("topright", c("Champ", "Breu"), col=c("paleturquoise3", "darkgoldenrod"), pch=16, bty="n")

##
plot(t1, Champ.Asymp.Mass.mean, type="l", lwd=3, col="paleturquoise3", ylim=c(0,1.2), xlab="Years",
     ylab="Proportion Litter Mass Remaining", main="Exp. 1 Asymptotic Model \n95% CI")
points(Champ.xNA$t, Champ.xNA$Mt, col="paleturquoise3")
polygon(c(t1[2:1000],rev(t1[2:1000])), c((Champ.Asymp.Mass.975[2:1000]),rev((Champ.Asymp.Mass.025[2:1000]))),
        col= rgb(0.59,0.8,0.8, alpha=0.3), border=F)

points(t1, Breu.Asymp.Mass.mean, type="l", lwd=3, col="darkgoldenrod")
points(Breu.xNA$t, Breu.xNA$Mt, col="darkgoldenrod")
polygon(c(t1[2:1000],rev(t1[2:1000])), c((Breu.Asymp.Mass.975[2:1000]),rev((Breu.Asymp.Mass.025[2:1000]))),
        col= rgb(0.933,0.91,0.667, alpha=0.3), border=F)
legend("topright", c("Champ", "Breu"), col=c("paleturquoise3", "darkgoldenrod"), pch=16, bty="n")


