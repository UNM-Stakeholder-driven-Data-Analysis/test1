#### read me ####

# the purpose of this script is to demonstrate statisitcal testing of trends by two methods

# Here, I use monthly data with no NAs and with seasonality removed via decomposition. 
# It's worth exploring what happens and how you need to adjust models if you use higher frequency data, retain NAs, retain seasonality, or use other methods to remove it.

#### libraries ####

library(tidyverse)
library(lubridate)
library(forecast)
library(MARSS)
library(nlme)
library(zoo)
library(beepr)
library(gridExtra)
library(lme4)
library(car)
library(visreg)

#### load data and format date/time ####

RG_abq = read.csv("Week 6/riogrande_ABQ.csv", header = T)

RG_abq$datetimeNM = as.POSIXct(RG_abq$datetime, "%m/%d/%y %H:%M", tz="America/Denver")

#### USGS metadata ####

# Some of the data that you have obtained from this U.S. Geological Survey database
# may not have received Director's approval. Any such data values are qualified
# as provisional and are subject to revision. Provisional data are released on the
# condition that neither the USGS nor the United States Government may be held liable
# for any damages resulting from its use.
#
# Additional info: https://help.waterdata.usgs.gov/policies/provisional-data-statement
#
# File-format description:  https://help.waterdata.usgs.gov/faq/about-tab-delimited-output
# Automated-retrieval info: https://help.waterdata.usgs.gov/faq/automated-retrievals
#
# Contact:   gs-w_support_nwisweb@usgs.gov
# retrieved: 2020-08-11 15:21:50 EDT       (nadww01)
#
# Data for the following 1 site(s) are contained in this file
#    USGS 08330000 RIO GRANDE AT ALBUQUERQUE, NM

#
# Data provided for site 08330000
#            TS   parameter     Description
#        101248       00060     Discharge, cubic feet per second
#        246224       00095     Specific conductance, water, unfiltered, microsiemens per centimeter at 25 degrees Celsius
#        246225       00010     Temperature, water, degrees Celsius
#        246226       00300     Dissolved oxygen, water, unfiltered, milligrams per liter
#        246227       63680     Turbidity, water, unfiltered, monochrome near infra-red LED light, 780-900 nm, detection angle 90 +-2.5 degrees, formazin nephelometric units (FNU)
#        246228       00400     pH, water, unfiltered, field, standard units
#
# Data-value qualification codes included in this output:
#        
#     A  Approved for publication -- Processing and review completed.
#     P  Provisional data subject to revision.
#     >  Actual value is known to be greater than reported value.
#     e  Value has been estimated.
#    91  Daily mean calculated from data on this day matches published daily mean within 1 percent
#    92  Daily mean calculated from data on this day matches published daily mean within 5 percent
#    93  Daily mean calculated from data on this day matches published daily mean within 10 percent


#### prep time series ####

# check quality of Q_csf 
table(RG_abq$Q_cfs_QAQC)
# all data is either approved for publication or has no QAQC code

sum(is.na(RG_abq$datetimeNM))
# there are 4 NA dates due to the time zone conversion. Remove these obs below.

# subset just discharge and time and remove obs with no date/time
RG_abq_Q = RG_abq %>% select(datetimeNM, Q_cfs) %>% drop_na(datetimeNM) %>% arrange(datetimeNM)

# round all date/time stamps to most common time step: 15 min
RG_abq_Q$datetimeNM = round_date(RG_abq_Q$datetimeNM, "15 minute")

# remove dup time stamps by taking mean of values
sum(duplicated(RG_abq_Q$datetimeNM))
RG_abq_Q[which(duplicated(RG_abq_Q$datetimeNM)),]
RG_abq_Q = RG_abq_Q %>%
  group_by(datetimeNM) %>%
  summarise(Q_cfs = mean(Q_cfs, na.rm = T))

# make dates at regular intervals
time <- data.frame(
  datetimeNM = seq.POSIXt(
    from = as.POSIXct(min(RG_abq_Q$datetimeNM, na.rm = T)),
    to = as.POSIXct(max(RG_abq_Q$datetimeNM, na.rm = T)),
    by = "15 min" ))
RG_abq_Q_reg = left_join(time, RG_abq_Q, by="datetimeNM")
# chec for duplicate date/time stamps
anyDuplicated(RG_abq_Q_reg$datetimeNM)
# check percentage of dataset with NAs - this is important to document!
sum(is.na(RG_abq_Q_reg))/nrow(RG_abq_Q_reg)*100

## fill gaps with spline interpolation ##
# for calculating long-term trends, you can be pretty liberal with how large of gaps you fill. However, if you go too big with a spline interpolation, you'll get wacky results (ALWAYS EXAMINE THE RESULTS OF GAP FILLING!!). To strike this balance, I'm filling gaps of up to five days here. 
par(mfrow=c(2,1)) # set up plotting window to comapare ts before and after gap filling
# Make univariate zoo time series #
ts.temp<-read.zoo(RG_abq_Q_reg, index.column=1, format="%Y-%m-%d %H:%M:%S", tz="America/Denver")
# ‘order.by’ are not unique warning suggests duplicate time stamps. I found that this is due to time zone changes, so nothing to worry about for regular time steps. 
plot(ts.temp)
# Apply NA interpolation method
RG_abq_filled = na.spline(ts.temp, na.rm = T, maxgap = 96*5)
plot(RG_abq_filled)
par(mfrow=c(1,1)) # reset plotting window
# revert back to df
RG_abq_filled = as.data.frame(RG_abq_filled)
RG_abq_filled$datetimeNM = RG_abq_Q_reg$datetimeNM
names(RG_abq_filled) = c(colnames(RG_abq_Q_reg)[2],colnames(RG_abq_Q_reg)[1])
RG_abq_filled = RG_abq_filled %>% select(datetimeNM, Q_cfs)
# check NAs that are left
sum(is.na(RG_abq_filled$Q_cfs))
sum(is.na(RG_abq_filled$Q_cfs))/length(RG_abq_filled$Q_cfs) * 100
# I've filled ~1% of dataset gaps with spline interpolation - report this in results!

# for calculating long-term trends, you want to aggregate to a fairly large time step. Goal is to balance how much data you retain (you want as much as possible) with smoothing over gaps and temporal trends unrelated to the long-term trend. I am going with monthly here to smooth over gaps but retain a lot of data. This will still leave me with a seasonal pattern that will complicate my trend analysis, so I'll remove this below.

# aggregate to monthly
RG_abq_mo = 
  RG_abq_filled %>% 
  mutate(yr = lubridate::year(datetimeNM)) %>%
  mutate(mo = lubridate::month(datetimeNM)) %>%
  select(yr, mo, Q_cfs) %>%
  group_by(yr, mo) %>%
  summarise(Q_cfs = mean(Q_cfs, na.rm = T))
RG_abq_mo$date = as.Date(paste(RG_abq_mo$yr, RG_abq_mo$mo, "15", sep="-"))
# check NAs that are left
sum(is.na(RG_abq_mo$Q_cfs))
sum(is.na(RG_abq_mo$Q_cfs))/length(RG_abq_mo$Q_cfs) * 100

# removing all NAs through some combination of gap filling and aggregating will make your life much easier, and it's pretty easy to do without messing up patterns for long-term trends. 

#### create time series ####

# need to do this to prep for removing seasonality

RG_abq_mo_ts = ts(RG_abq_mo$Q_cfs, start = c(1989, 10), frequency = 12)
head(RG_abq_mo_ts)

par(mfrow=c(1,1))
plot(RG_abq_mo_ts)

#### remove seasonality ####

# examine seasonality
par(mfrow=c(3,1))
plot(RG_abq_mo_ts)
Acf(RG_abq_mo_ts)
Pacf(RG_abq_mo_ts)

# decompose into additive components
plot(decompose(RG_abq_mo_ts))
# decompose into multiplicative components
plot(decompose(RG_abq_mo_ts, type="multiplicative"))
# extract components from multiplicative
RG_abq_mo_ts_decomp = decompose(RG_abq_mo_ts, type="multiplicative")
RG_abq_mo_ts_trend = RG_abq_mo_ts_decomp$trend
RG_abq_mo_ts_remainder = RG_abq_mo_ts_decomp$random
# save de-seasoned ts
RG_abq_mo_ts_DEs = RG_abq_mo_ts_trend * RG_abq_mo_ts_remainder

# compare original to de-seasoned ts
par(mfrow=c(3,2))
plot(RG_abq_mo_ts)
plot(RG_abq_mo_ts_DEs)
Acf(RG_abq_mo_ts)
Acf(RG_abq_mo_ts_DEs)
Pacf(RG_abq_mo_ts)
Pacf(RG_abq_mo_ts_DEs)

# revert back to df
RG_abq_mo_DEs = as.data.frame(RG_abq_mo_ts_DEs)
RG_abq_mo_DEs$date = RG_abq_mo$date
names(RG_abq_mo_DEs) = c("Q_cfs","date")
RG_abq_mo_DEs = RG_abq_mo_DEs %>% select(date, Q_cfs) %>% arrange(date)
RG_abq_mo_DEs = na.trim(RG_abq_mo_DEs, "both")

ggplot(RG_abq_mo_DEs, aes(x=date, y=Q_cfs))+
  geom_path() + geom_point() + theme_bw()

#### linear trends ####

# add simple time steps to df
RG_abq_mo_DEs$t = c(1:nrow(RG_abq_mo_DEs))

mod = lm(Q_cfs ~ t, RG_abq_mo_DEs)

summary(mod)

visreg(mod,"t")

confint(mod, 't', level=0.95)

## diagnostics ##
Acf(resid(mod))
forecast::checkresiduals(mod)

#### test & calculate trends - nlme::gls ####

# see package manual: https://cran.r-project.org/web/packages/nlme/nlme.pdf

# ask auto.arima what it thinks the autocorrelation structure is
auto.arima(RG_abq_mo_DEs$Q_cfs)

# fit AR(1) regression model with time as a predictor
mod_Ar1 = gls(Q_cfs ~ t, data=RG_abq_mo_DEs, correlation=corAR1(), method="ML")

# fit some other candidate structures
mod_AMRAp1q1 = gls(Q_cfs ~ t, data=RG_abq_mo_DEs, correlation=corARMA(p=1,q=1), method="ML")
mod_AMRAp2 = gls(Q_cfs ~ t, data=RG_abq_mo_DEs, correlation=corARMA(p=2), method="ML")
mod_AMRAp3 = gls(Q_cfs ~ t, data=RG_abq_mo_DEs, correlation=corARMA(p=3), method="ML")

# compare models with AIC, AICc, and BIC
# For small data, use AICc – the small sample correction which provides greater penalty for each parameter but approaches AIC as n becomes large. If it makes a difference, you should use it. 
# For large data and especially time series data, consider BIC. BIC is better in situations where a false positive is more misleading than a false negative. Remember that false positives are more common with time series. 
bbmle::AICtab(mod_Ar1,mod_AMRAp1q1,mod_AMRAp2,mod_AMRAp3)
bbmle::AICctab(mod_Ar1,mod_AMRAp1q1,mod_AMRAp2,mod_AMRAp3)
bbmle::BICtab(mod_Ar1,mod_AMRAp1q1,mod_AMRAp2,mod_AMRAp3)

summary(mod_Ar1)
# intervals() for nlme is equivelant to confint() for lm
intervals(mod_Ar1)

# notice that p value 0.0002 for mod_Ar1 is much higher than in linear model
# and 95% CIs are much broader for mod_Ar1
# this demonstrates the type 1 error of lm with ts data!

par(mfrow=c(1,1))
visreg(mod_Ar1,"t")

### Important notes about extracting residuals from model fits!! ###

# It's important to understand that many extraction fxns in R, such as residuals(modelfit) (same as resid(modelfit)), will detect the object type and call on methods from that package appropriate for that object. So, residuals(modelfit) is using different methods for different model types when the model package requires it, and you need to look up the options for these different methods.
# E.g., residuals(nlme model) calls residuals.lme(nlme model), which has different options than if you call residuals(model fit) on a different kind of model. 
# see ?residuals.gls for the methods avaiable for this model type
# type ?residuals. into your console and scroll through the options for other residuals methods for loaded packages

# For gls, you want to assess assumptions on normalized residuals, which is not an option for standard linear models.
# normalized residuals = standardized residuals pre-multiplied by the inverse square-root factor of the estimated error correlation matrix
# see https://stats.stackexchange.com/questions/80823/do-autocorrelated-residual-patterns-remain-even-in-models-with-appropriate-corre

Acf(resid(mod_Ar1))

# extract and assess residuals
par(mfrow=c(1,3))
Acf(resid(mod_Ar1, type = "normalized"), main="GLS AR(1) model residuals")
plot(resid(mod_Ar1, type = "normalized")~c(1:length(RG_abq_mo_DEs$t)), main="GLS AR(1) model residuals"); abline(h=0)
qqnorm(resid(mod_Ar1, type = "normalized"), main="GLS AR(1) model residuals", pch=16, 
       xlab=paste("shapiro test: ", round(shapiro.test(resid(mod_Ar1, type = "normalized"))$statistic,2))); qqline(resid(mod_Ar1, type = "normalized"))

# exctract parameter estimates for comparison with MARSS
mod_Ar1.phi = coef(mod_Ar1$modelStruct[[1]], unconstrained=FALSE)
ests.gls = c(b=mod_Ar1.phi, alpha=coef(mod_Ar1)[1],
             time=coef(mod_Ar1)[2],
             logLik=logLik(mod_Ar1))

#### test & calculate trends - UARSS ####

## MARSS ##
# User;s guide: https://cran.r-project.org/web/packages/MARSS/vignettes/UserGuide.pdf
# Package manual: https://cran.r-project.org/web/packages/MARSS/MARSS.pdf
# Lectures: https://nwfsc-timeseries.github.io/atsa/
# Lab book: https://nwfsc-timeseries.github.io/atsa-labs/
# Quick start guide: https://cran.r-project.org/web/packages/MARSS/vignettes/Quick_Start.pdf
# MARSS is a extremely flexible ts modeling framework with a steep learning curve requiring lots of matrix algebra for many applications, so wade into these documents with caution and prepare to get overwhelmed quickly! Also the user's guide is written, in my opinion, pretty poorly for the average user. This class will give you a few common "recipes", but know that almost anything is possible if you take time to learn the modeling framework in full.


## MARSS equivelent to gls corAR1 model ##

## format response var y as a vector
dat = as.vector(RG_abq_mo_DEs$Q_cfs)
# remove leading and trailing NAs #
dat = na.trim(dat, is.na="any")

## set model parameters - see powerpoint for 
mod.list.AR1 = list(
  B=matrix("b"),          # state process model
  Q=matrix("q"),          # state process model
  U="zero",               # state process model
  C="zero",               # state process model
  Z="identity",           # observation process model
  R=matrix("r"),          # observation process model
  A=matrix("intercept"),  # observation process model
  D=matrix("time"),       # observation process model
  d=matrix(c(1:length(dat)), nrow=1), # observation process model
  x0=matrix(dat[1]), 
  tinitx=0
)

## fit model
# mod.AR1 <- MARSS(dat, model=mod.list.AR1, method="BFGS")
mod.AR1 <- MARSS(dat, model=mod.list.AR1, control=list(maxit=10000))
beep(2)
# est.AR1 <- MARSSparamCIs(mod.AR1, method = "parametric", alpha = 0.05, nboot = 2000, silent=F)
est.AR1 <- MARSSparamCIs(mod.AR1, method = "hessian", alpha = 0.05)

## extract parameter estimates for comparison to gls
ests.marss = c(b=coef(mod.AR1)$B, alpha=coef(mod.AR1)$A,
               time=coef(mod.AR1)$D[1],
               logLik=logLik(mod.AR1))

## compare UARSS and gls results
# parameter estimates
ests.marss
ests.gls
# 95% CIs on trend estimate
intervals(mod_Ar1)
est.AR1
# notice that UARSS provides a slightly narrower confidence interval

## test residuals for ac
# extract residuals
resids.1 <- residuals(mod.AR1) # see ?residuals.marssMLE
# plot residuals
par(mfrow=c(2,3))
Acf(resids.1$model.residuals[1,], main="Observation process model residuals")
plot(resids.1$model.residuals[1,]~c(1:length(dat)), main="Observation process model residuals"); abline(h=0)
qqnorm(resids.1$model.residuals[1,], main="Observation process model residuals", pch=16, 
       xlab=paste("shapiro test: ", round(shapiro.test(resids.1$model.residuals[1,])$statistic,2))); qqline(resids.1$model.residuals[1,])
#
Acf(resids.1$state.residuals[1,], main="State process model residuals", na.action = na.pass)
plot(resids.1$state.residuals[1,]~c(1:length(dat)), main="State process model residuals"); abline(h=0)
qqnorm(resids.1$state.residuals[1,], main="Observation process model residuals", pch=16, 
       xlab=paste("shapiro test: ", round(shapiro.test(resids.1$state.residuals[1,])$statistic,2))); qqline(resids.1$state.residuals[1,])


# MARSS is focused on explaining temporal dynamics and is not "willing" to shunt all remaining autocorrelation to error, unlike nlme options. If it can't be explained by an autocorrelated process or covars, it will retain the remaining autocorrelation. This is good information! But not often the most practical option if you can't include the right covars. 
# Here, this tells us that a trend over time in observation model is inadequate to capture all the systematic variation in Rio Grande discharge, suggests additional covars are needed.  


## Plot fitted values over observations ###
# extract MARSS results
kf=print(mod.AR1, what="kfs") # Kalman filter and smoother output
# plot observed data (y)
par(mfrow=c(1,1),oma = c(0, 0, 2, 0))
plot(as.vector(dat) ~ RG_abq_mo_DEs$date, type="p", pch=19,
     main = "UARSS model predictions conditioned on all y",
     ylab = "Discharge (cfs)", xlab="")
# calc and plot predicted values
predicts = as.vector(kf$xtT) + 
  coef(mod.AR1)$A[1] + 
  (as.vector(mod.AR1[["model"]][["fixed"]][["d"]])* coef(mod.AR1)$D[1])
lines(predicts ~ RG_abq_mo_DEs$date, col="blue",lwd=2) 
lines(RG_abq_mo_DEs$date, predicts-1.96*mod.AR1$states.se,
      type="l",lwd=1,lty=2,col="blue")
lines(RG_abq_mo_DEs$date, predicts+1.96*mod.AR1$states.se,
      type="l",lwd=1,lty=2,col="blue")
# calc and plot predicted values without trend
predicts.trendless = as.vector(kf$xtT) + 
  coef(mod.AR1)$A[1] 
lines(predicts.trendless ~ RG_abq_mo_DEs$date, col="red",lwd=2) 
lines(RG_abq_mo_DEs$date, predicts.trendless-1.96*mod.AR1$states.se,
      type="l",lwd=1,lty=2,col="red")
lines(RG_abq_mo_DEs$date, predicts.trendless+1.96*mod.AR1$states.se,
      type="l",lwd=1,lty=2,col="red")
mtext("Fitted values over observations", outer = TRUE, cex = 1.5)
# dashed lines are 95% CIs, calculated from the standard error. The value of 1.96 is based on the fact that 95% of the area of a normal distribution is within 1.96 standard deviations of the mean.

# trend model predictions and non-trend model predictions start diverging ~1996
# red line indicates Rio Grande discharge if there weren't a long-term decline. 

# you can make a similar plot from nlme::gls results! go try and figure that out :)

