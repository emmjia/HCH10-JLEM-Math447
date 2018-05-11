####################################################
##     Estimation in multiple linear regression    #
####################################################
#viscosity example (eg.10.1,  pp. 454)
viscosity=read.table("data/viscosity.txt", header = TRUE)
viscosity.lm=lm(Viscosity ~ Temperature + CatalystFeedRate, viscosity)
summary(viscosity.lm)

#use your own code to estimate beta
n=nrow(viscosity)
X=cbind(rep(1,n),viscosity$Temperature,viscosity$CatalystFeedRate)
X=as.matrix(X)
Y=viscosity$Viscosity
solve(t(X)%*%X)%*%t(X)%*%Y #beta_hat
viscosity.lm

#use your own code to esimate sigma^2
k=2
p=k+1
H=X%*%solve(t(X)%*%X)%*%t(X)
sigma2_hat=t(Y)%*%(diag(n)-H)%*%Y/(n-p) #diag(n) is n-dim identity matrix
sqrt(sigma2_hat)

#use your own code to esimate standard error of beta
cova=as.numeric(sigma2_hat)*solve(t(X)%*%X)
cova
sqrt(diag(cova))
summary(viscosity.lm)


#regresssion analysis of a 2^3 factorial design with center points
#yield example (eg.10.2,  pp. 454)
Temperature=c(rep(c(-1,1),4),rep(0,4))
Pressure=c(rep(c(-1,-1, 1, 1), 2), rep(0,4))
Catalyst=c(rep(-1,4),rep(1,4),rep(0,4))
Yield=c(32, 46, 57, 65, 36, 48, 57, 68, 50, 44, 53, 56)
n=length(Yield);k=3;p=k+1
yield=data.frame(Temperature, Pressure, Catalyst, Yield)
yield.lm=lm(Yield ~ Temperature + Pressure + Catalyst, yield)
summary(yield.lm)
summary(yield.lm, correlation = TRUE)$correlation
X=cbind(rep(1,n),Temperature,Pressure, Catalyst)
Y=Yield
H=X%*%solve(t(X)%*%X)%*%t(X)
sigma2_hat=t(Y)%*%(diag(n)-H)%*%Y/(n-p)
cova=as.numeric(sigma2_hat)*solve(t(X)%*%X)
cova
t(X)%*%X # is diagonal matrix because of columns arare orthogonal

#same example, with run 8 missing
yield8.lm=lm(formula = Yield ~ Temperature + Pressure + Catalyst, data = yield[-8, ])
summary(yield8.lm)
vcov(yield8.lm)
summary(yield8.lm, correlation = TRUE)$correlation
t(X[-8,])%*%X[-8,]

########################################################
## hypothesis testing in multiple linear regression    #
########################################################
#test the significance of regression or a single coefficient
summary(viscosity.lm)

#test a group of coefficients (test CatalystFeedRate only, could be more complex group)
viscosity.lm1=lm(Viscosity ~ Temperature, viscosity) #reduced model
anova(viscosity.lm1,viscosity.lm) #anova(reduced model, full model)
3.52^2 #here F test stat = squared t test stat

#test a group of coefficients (test CatalystFeedRate*Temperature, CatelystFeedRate^2, Temperature^2)
viscosity.lm2=lm(Viscosity ~ Temperature + CatalystFeedRate, viscosity)
viscosity.lm3=lm(Viscosity ~ Temperature + CatalystFeedRate + I(Temperature^2) + I(CatalystFeedRate^2)
                 +I(CatalystFeedRate * Temperature), viscosity)
anova(viscosity.lm2,viscosity.lm3)

#################################################################
## Confidence interval and prediction interval in regression    #
#################################################################
#confidence interval of coefficients
confint(viscosity.lm)
confint(viscosity.lm,level=0.99)

#confidence interval of mean response
predict(viscosity.lm, newdata = data.frame(Temperature = 90, CatalystFeedRate = 10), interval = "confidence")
predict(viscosity.lm, newdata = data.frame(Temperature = 90, CatalystFeedRate = 10), interval = "confidence",level=0.99)

#prediction interval of a new observation
predict(viscosity.lm, newdata = data.frame(Temperature = 90, CatalystFeedRate = 10), interval = "prediction")
predict(viscosity.lm, newdata = data.frame(Temperature = 90, CatalystFeedRate = 10), interval = "prediction",level=0.99)

##############################
## regression diagnostics    #
##############################
#standard residual plots
Y=viscosity$Viscosity
res=Y-fitted(viscosity.lm)
plot(fitted(viscosity.lm),res)
plot(viscosity$Temperature,res)
plot(viscosity$CatalystFeedRate,res)
library(car)
qqPlot(res) #for point esimator normal assumption is not required, is require for testing and CI

##############################
## influence diagnostics     #
##############################
#prediction R2
RsqPred = function(lm) {
  infl = influence(lm)
  PRESS = sum((infl$wt.res / (1 - infl$hat))^2)
  rsq = summary(lm)$r.squared
  sst = sum(infl$wt.res^2) / (1 - rsq)
  return (1 - PRESS / sst)
}
RsqPred(viscosity.lm) #penalizes over-fitting more than adj.R2
summary(viscosity.lm)$r.squared
summary(viscosity.lm)$adj.r.squared

#cook's distance (one measure of influence)
# high influence point if D > 1.
cooks.distance(viscosity.lm)
max(cooks.distance(viscosity.lm))
which.max(cooks.distance(viscosity.lm))

#leverage (high leverage if h>2*p/n)
p=3;n=nrow(viscosity)
lev=influence(viscosity.lm)$hat
2*p/n
which(lev>2*p/n)

#S_(i)
influence(viscosity.lm)$sigma

#difference between beta.hat and beta.hat_(i)
influence(viscosity.lm)$coefficients

#influence() produces four measures
influence(viscosity.lm)

#Leverage and Cook's Distance are shown in the fourth residual plot
plot(viscosity.lm)

##################################
##     Box-Cox transformation    #
##################################
peakDischarge=read.table("data/peak-discharge.txt", header = TRUE)
peakDischargeLong=reshape(peakDischarge, varying = 2:7,
                          v.names = "Discharge", timevar = "Obs",
                          idvar = "Code", direction = "long")
#checking model adequacy
peak.lm=lm(Discharge ~ factor(Method), data = peakDischargeLong)
library(car)
# check normaility
qqPlot(residuals(peak.lm))
# check structureless(independence)
plot(fitted(peak.lm), residuals(peak.lm))

#after transformation
peak.lm1=lm(sqrt(Discharge) ~ factor(Method), data = peakDischargeLong)
# check normaility
qqPlot(residuals(peak.lm1))
# check structureless(independence)
plot(fitted(peak.lm1), residuals(peak.lm1))
#test homogeneity of variances
bartlett.test(sqrt(Discharge)~factor(Method),data=peakDischargeLong)
with(peakDischargeLong, leveneTest(sqrt(Discharge), factor(Method)))

#Box-cox transofrmation
peak.aov=aov(Discharge ~ factor(Method), data = peakDischargeLong)
library(MASS)
bc=boxcox(peak.aov,lambda=seq(-1,1.5,0.1))
idx=which.max(bc$y)
bc$x[idx] #the lambda maximizes the log-likelihood
    
    