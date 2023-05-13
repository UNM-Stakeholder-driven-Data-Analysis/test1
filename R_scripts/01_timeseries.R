#### read me ####

# the purpose of this script is to demonstrate time series exploration, simple ARIMA model fitting, forecasting, and wavelet decomposition. There is also some time series data wrangling. 

#### libraries ####

library(tidyverse)
library(lubridate)
library(forecast)
library(zoo)
library(xts)
library(imputeTS)
library(tseries)
library(astsa)
library(WaveletComp)

#### load and filter data ####

# this is stream chemistry, discharge, and precipitation in 4 watersheds in Alaska from summer 2017. Data is collected every 15 min for all variables.
# note that this means that there are 96 observations per day
setwd("C:/Users/marag.DESKTOP-NEKBNJE/OneDrive/Desktop/Stakeholder data analysis/test1/Data")

dat = read.csv("CPCRW.2017.csv", row.names = 1, header=T)

# format date/time
dat$date_timeAK = as.POSIXct(dat$date_timeAK, tz="America/Anchorage")
# check that date/time formatted correctly
class(dat$date_timeAK); tz(dat$date_timeAK)

# site = C2, variable = nitrate (uM)
C2_no3 = dat %>% 
  filter(site.ID=="C2") %>% 
  select(date_timeAK, nitrate_uM_c_bc)

#### explore filling gaps ####

# general guidelines:
# - if you're interested in diel patters, be very cautious about filling gaps > 6 h
# - if you're interested in day to day patters, case by case basis - examine whether gaps are during events like storms
# - in general, examine time series before and after you fill gaps and ask whether you've influenced patterns of interest in a defensible way or not
# note that your data must be evenly spaced in time before doing this!!

### 15 min data ###

## fill with linear interpolation
# Make univariate zoo time series#
ts.temp<-read.zoo(C2_no3, index.column=1, format="%Y-%m-%d %H:%M:%S", tz="America/Anchorage")
# Apply NA interpolation method
C2_no3_filled_linearinterp = na.approx(ts.temp, na.rm = T, maxgap = 24*4)
# revert back to df
C2_no3_filled_linearinterp = as.data.frame(C2_no3_filled_linearinterp)
C2_no3_filled_linearinterp$date_timeAK = as.POSIXct(row.names(C2_no3_filled_linearinterp), tz="America/Anchorage")
names(C2_no3_filled_linearinterp) = c(colnames(C2_no3)[2],colnames(C2_no3)[1])

# fill with spline interpolation
# Make univariate zoo time series #
ts.temp<-read.zoo(C2_no3, index.column=1, format="%Y-%m-%d %H:%M:%S", tz="America/Anchorage")
# Apply NA interpolation method
C2_no3_filled_splineinterp = na.spline(ts.temp, na.rm = T, maxgap = 24*4)
# revert back to df
C2_no3_filled_splineinterp = as.data.frame(C2_no3_filled_splineinterp)
C2_no3_filled_splineinterp$date_timeAK = as.POSIXct(row.names(C2_no3_filled_splineinterp), tz="America/Anchorage")
names(C2_no3_filled_splineinterp) = c(colnames(C2_no3)[2],colnames(C2_no3)[1])

# explore ?na.approx in zoo package for  more options

## fill with arima model + kalman smoother
# Document gaps  
largegap.num = 24*4 # 4*6 = 6 hours
is.na.rle <- rle(is.na(C2_no3$nitrate_uM_c_bc))
is.na.rle$values <- is.na.rle$values & is.na.rle$lengths >= (largegap.num)
biggaps = C2_no3[inverse.rle(is.na.rle), ]
tz(biggaps$date_timeAK) = "America/Anchorage"
biggaps = subset(biggaps, select = "date_timeAK")
# Make univariate time series: covert to zoo, then to ts #
xts<-read.zoo(C2_no3, index.column=1, format="%Y-%m-%d %H:%M:%S", tz="America/Anchorage")
xts<-as.xts(xts)
# remove leading and trailing NAs #
xts = na.trim(xts, is.na="any")
# Apply auto.arima and kalman filter to impute missing values #
fit = auto.arima(xts) 
kal = KalmanSmooth(xts, fit$model)
id.na<-which(is.na(xts))
for(i in id.na) {
  xts[i]<-fit$model$Z %*% kal$smooth[i,]}
# revert to dataframe #
C2_no3_filled_ARIMA = as.data.frame((xts))
C2_no3_filled_ARIMA$date_timeAK = as.POSIXct(row.names(C2_no3_filled_ARIMA), tz="America/Anchorage")
names(C2_no3_filled_ARIMA) = c("nitrate_uM_c_bc", "date_timeAK")
# remove large gaps # 
C2_no3_filled_ARIMA$nitrate_uM_c_bc[C2_no3_filled_ARIMA$date_timeAK %in% as.POSIXct(biggaps$date_timeAK)] = NA


par(mfrow=c(4,1))
xmin = min(C2_no3$date_timeAK); xmax = max(C2_no3$date_timeAK)
plot(C2_no3$nitrate_uM_c_bc ~ C2_no3$date_timeAK, type="l", main="original", xlim=c(xmin, xmax))
plot(C2_no3_filled_linearinterp$nitrate_uM_c_bc ~ C2_no3_filled_linearinterp$date_timeAK, type="l", main="linear interpolation", xlim=c(xmin, xmax))
plot(C2_no3_filled_splineinterp$nitrate_uM_c_bc ~ C2_no3_filled_splineinterp$date_timeAK, type="l", main="spline interpolation", xlim=c(xmin, xmax))
plot(C2_no3_filled_ARIMA$nitrate_uM_c_bc ~ C2_no3_filled_ARIMA$date_timeAK, type="l", main="ARIMA interpolation", xlim=c(xmin, xmax))

par(mfrow=c(4,1))
xmin = min(C2_no3$date_timeAK); xmax = as.POSIXct("2017-06-05")
ymin= 39; ymax=42
plot(C2_no3$nitrate_uM_c_bc ~ C2_no3$date_timeAK, type="l", main="original", xlim=c(xmin, xmax), ylim=c(ymin, ymax))
plot(C2_no3_filled_linearinterp$nitrate_uM_c_bc ~ C2_no3_filled_linearinterp$date_timeAK, type="l", main="linear interpolation", xlim=c(xmin, xmax), ylim=c(ymin, ymax))
plot(C2_no3_filled_splineinterp$nitrate_uM_c_bc ~ C2_no3_filled_splineinterp$date_timeAK, type="l", main="spline interpolation", xlim=c(xmin, xmax), ylim=c(ymin, ymax))
plot(C2_no3_filled_ARIMA$nitrate_uM_c_bc ~ C2_no3_filled_ARIMA$date_timeAK, type="l", main="ARIMA interpolation", xlim=c(xmin, xmax), ylim=c(ymin, ymax))

# pause here and play around with gap sizes!
# note that best practice is to fill gaps before aggregating data to daily, weekly, etc.

# I will go forward with the spline interploation and max map size of 2 days

#### create time series ####

# create time series object
# note that dates must be sequential before doing this!! xts and other ts-making functions will assume that they are. 
C2_no3_xts = xts(C2_no3_filled_splineinterp$nitrate_uM_c_bc, order.by = C2_no3_filled_splineinterp$date_timeAK)
class(C2_no3_xts)
plot(C2_no3_xts)
summary(C2_no3_xts)
# 1201 NAs
# decompose() and stl() require a ts() or mtst() object
# I might use mtst() if this were multi-annual data with a seasonal signal to capture both diel and seasonal frequencies
C2_no3_ts = ts(C2_no3_filled_splineinterp$nitrate_uM_c_bc, 
               frequency = 96, 
               start=hour(min(C2_no3_filled_splineinterp$date_timeAK)))
print(C2_no3_ts, calendar = T)
# create one for the contiguous time series
C2_no3_ts_contig = na.contiguous(C2_no3_filled_splineinterp)
C2_no3_ts_contig = ts(C2_no3_ts_contig$nitrate_uM_c_bc, 
               frequency = 96, 
               start=hour(min(C2_no3_ts_contig$date_timeAK)))

# aggregate to daily and create time series object
C2_no3_daily = C2_no3_filled_splineinterp %>% 
  mutate(day = as.Date(date_timeAK)) %>%
  select(day, nitrate_uM_c_bc) %>%
  group_by(day) %>%
  summarise(nitrate_uM = mean(nitrate_uM_c_bc, na.rm = T))
C2_no3_daily_xts = xts(C2_no3_daily$nitrate_uM, order.by = C2_no3_daily$day)
plot(C2_no3_daily_xts)
summary(C2_no3_daily_xts)
# 12 NAs
C2_no3_daily_ts = ts(C2_no3_daily$nitrate_uM,  
                     frequency = 1, 
                     start=yday(min(C2_no3_daily$day)))
print(C2_no3_daily_ts, calendar = T)
# create one for the contiguous time series
C2_no3_daily_ts_contig = na.contiguous(C2_no3_daily)
C2_no3_daily_ts_contig = ts(C2_no3_daily_ts_contig$nitrate_uM, 
                      frequency = 1, 
                      start=yday(min(C2_no3_daily_ts_contig$day)))

# This is only one summer of data, so I won't aggregate to monthly or yearly, but I would if it were longer


#### explore autocorrelation structure ####

# The autocorrelation function ACF measures the correlation of a time series against a time-shifted version of itself
# The partial autocorrelation function PACF measures the correlation between a time series against a time-shifted version of itself while accounting for how autocorrelationperpetuates through lags.  

# reset plotting window
par(mfrow=c(1,1))

# examine acf and pacf for 15 min data
Acf(C2_no3_xts, na.action = na.interp, lag.max = 96*2)
Pacf(C2_no3_xts, na.action = na.interp, lag.max = 96*2)
Acf(C2_no3_xts, na.action = na.contiguous, lag.max = 96*2)
Pacf(C2_no3_xts, na.action = na.contiguous, lag.max = 96*2)
# strong autocorrelation at lag 1, 2, 3

# examine acf and pacf for daily data
Acf(C2_no3_daily_xts, na.action = na.pass, lag.max = 10)
Pacf(C2_no3_daily_xts, na.action = na.pass, lag.max = 10)
Acf(C2_no3_daily_xts, na.action = na.interp, lag.max = 10)
Pacf(C2_no3_daily_xts, na.action = na.interp, lag.max = 10)
Acf(C2_no3_daily_xts, na.action = na.contiguous, lag.max = 10)
Pacf(C2_no3_daily_xts, na.action = na.contiguous, lag.max = 10)
# strong autocorrelation at lag 1, possibly (but unlikely) 6


#### classic decomposition ####

# note that trend is not a consistent up to down trend, it's just the time series minus seaosnal pattern + noise

## 15 min data ##
# decompose into additive components
plot(decompose(C2_no3_ts_contig))
plot(decompose(na.interp(C2_no3_ts)))
plot(stl(C2_no3_ts, s.window = 96, na.action = na.contiguous))
plot(stl(C2_no3_ts, s.window = 96, na.action = na.interp))
# decompose into multiplicative components
plot(decompose(na.contiguous(C2_no3_ts), type="multiplicative"))
plot(decompose(na.interp(C2_no3_ts), type="multiplicative"))
# extract components from additive
C2_no3_ts_decomp = stl(C2_no3_ts, s.window = 96, na.action = na.contiguous)
C2_no3_ts_diel = C2_no3_ts_decomp$time.series[,"seasonal"]
plot(C2_no3_ts_diel)
C2_no3_ts_trend = C2_no3_ts_decomp$time.series[,"trend"]
plot(C2_no3_ts_trend)
C2_no3_ts_noise = C2_no3_ts_decomp$time.series[,"remainder"]
plot(C2_no3_ts_noise)

plot(C2_no3_ts_trend + C2_no3_ts_noise)

## daily data ##
# decompose
decompose(C2_no3_daily_ts)
plot(stl(C2_no3_daily_ts, s.window = 1, na.action = na.contiguous))
plot(stl(C2_no3_daily_ts, s.window = 1, na.action = na.interp))
# this does not work because there is no seasonal pattern. You can't decompose data that doesn't have a seasonal component! 

#### differencing ####

# An alternative to decomposition for removing trends or seasonality is differencing. 
# this method does not require having a ts object or seasonaility

# first-differencing a time series will remove a linear trend (i.e., differences = 1)
# twice-differencing will remove a quadratic trend (i.e., differences = 2). 
# first-differencing a time series at a lag equal to the period will remove a seasonal trend (e.g., set lag = 12 for monthly data).

## 15 min data ##
# 1st-difference removes linear trend 
C2_no3_ts_d1 <- diff(C2_no3_filled_splineinterp$nitrate_uM_c_bc, differences = 1, lag=1)
# 2nd-difference removes quadradic trend
C2_no3_ts_d2 <- diff(C2_no3_filled_splineinterp$nitrate_uM_c_bc, differences = 2, lag=1)
# 1st-difference + lag = frequency/period removes linear trend  + seasonality
C2_no3_ts_d1_lag96 <- diff(C2_no3_filled_splineinterp$nitrate_uM_c_bc, differences = 1, lag=96)
# 2nd-difference + lag = frequency/period removes quadradic trend  + differences
C2_no3_ts_d2_lag96 <- diff(C2_no3_filled_splineinterp$nitrate_uM_c_bc, differences = 2, lag=96)
## plot the differenced data
par(mfrow=c(2,2))
plot(C2_no3_ts_d1, type="l", main="1st differenced")
plot(C2_no3_ts_d2, type="l", main="2nd differenced")
plot(C2_no3_ts_d1_lag96, type="l", main="1st differenced with lag")
plot(C2_no3_ts_d2_lag96, type="l", main="2nd differenced with lag")

## daily data ##
## once-difference 
C2_no3_daily_ts_d1 <- diff(C2_no3_daily$nitrate_uM, differences = 1, lag=1)
## twice-difference 
C2_no3_daily_ts_d2 <- diff(C2_no3_daily$nitrate_uM,, differences = 2, lag=1)
## plot the differenced data
par(mfrow=c(2,1))
plot(C2_no3_daily_ts_d1, type="l")
plot(C2_no3_daily_ts_d2, type="l")

#### forecasting (no covariates) - 15 min data ####

## Box-Jenkins method
# A. Model form selection
# - Evaluate stationarity
# - Selection of the differencing level (d) – to fix stationarity problems
# - Selection of the AR level (p)
# - Selection of the MA level (q)
# B. Parameter estimation
# C. Model checking
# D. Forecast

### 15 min data ###

## A. Model form selection - Evaluate stationarity and select differencing level (if needed)
# The basic stationarity diagnostics are the following:
# 1. Plot your data. 
par(mfrow=c(1,1))
plot(C2_no3_ts)
# Look for
#  - An increasing or decreasing trend: Possibly decreasing
#  - A non-zero mean (if no trend): Yes
#  - Does it look like it might be stationary around a trend?: No
#  - Strange shocks or steps in your data: No


# 2. Apply stationarity tests
#  - adf.test() p-value less than 0.05 (reject null) suggests ts is stationary

adf.test(C2_no3_ts) # does not allow NAs
adf.test(na.interp(C2_no3_ts)) 
adf.test(na.contiguous(C2_no3_ts)) 
adf.test(na.interp(C2_no3_ts), k=0) 
adf.test(na.contiguous(C2_no3_ts), k=0) 

#  - kpss.test() p-value greater than 0.05 (do not reject null) suggests ts is stationary

kpss.test(C2_no3_ts_contig, null="Level") 
kpss.test(C2_no3_ts_contig, null="Trend") 

# 3. If stationarity tests are failed, then try differencing to correct
#  - Try ndiffs() in the forecast package or manually try different differences.

ndiffs(C2_no3_ts_contig)

kpss.test(diff(C2_no3_ts_contig, differences = 1), null="Level") 
kpss.test(diff(C2_no3_ts_contig, differences = 1, lag=96), null="Level") 
kpss.test(diff(C2_no3_ts_contig, differences = 1), null="Trend") 
kpss.test(diff(C2_no3_ts_contig, differences = 1, lag=96), null="Trend") 

# note that differencing just once is sufficient to achieve stationarity... suggests that diel pattern may be somewhat irregular. 

## A. Model form selection - Selection of the AR level (p) and AR level (p)

fit.auto.arima = forecast::auto.arima(C2_no3_ts_contig, start.p = 0, max.p = 3, start.q = 0, max.q = 3)
fit.auto.arima # see best fit parameters

fit.sarima = sarima(C2_no3_ts_contig, p=3, d=1, q =0)

fit.arima = arima(C2_no3_ts_contig, order = c(3, 1, 0))
par(mfrow=c(2,1)); Acf(resid(fit.arima)); Pacf(resid(fit.arima))

fit.Arima = forecast::Arima(C2_no3_ts_contig, order = c(3, 1, 0), include.constant = F)
par(mfrow=c(2,1)); Acf(resid(fit.Arima)); Pacf(resid(fit.Arima))

# this fit is pretty marginal (lots of autocorrelation remaining, lots of non-normal residuals), which suggests that there are things happening in this time series that autoregressive processes alone cannot account for.

## Forecast!
par(mfrow=c(3,1))
fore.auto.arima <- forecast::forecast(fit.auto.arima, h = 96*10)
plot(fore.auto.arima)
fore.arima <- forecast::forecast(fit.arima, h = 96*10)
plot(fore.arima)
fore.Arima <- forecast::forecast(fit.Arima, h = 96*10)
plot(fore.Arima)

# ask auto.arima to detect seasonal pattern
fit.auto.arima = forecast::auto.arima(C2_no3_ts_contig, seasonal = TRUE)
fit.auto.arima # see best fit parameters
# auto.arima is not detecting a seasonal pattern, suggesting that it is not regular enough to improve model fit

fit.Arima.seasonal = forecast::Arima(C2_no3_ts_contig, 
                                     order = c(3, 1, 0), seasonal = c(0,1,1),
                                     include.constant = F)
fit.Arima.seasonal <- forecast::forecast(fit.Arima.seasonal, h = 96*10)
par(mfrow=c(1,1))
plot(fit.Arima.seasonal)
ggtsdisplay(residuals(fit.Arima.seasonal))

# residuals are still terrible - this model is not meeting assumptions. 
# pause and play around with seasonal specification, but note that it will take a long time to fit each

# forecast the seaosnally adjusted data?
C2_no3_ts_NS = C2_no3_ts_trend + C2_no3_ts_noise
fit.auto.arima = forecast::auto.arima(C2_no3_ts_NS, start.p = 0, max.p = 3, start.q = 0, max.q = 3)
fit.auto.arima 
fore.auto.arima <- forecast::forecast(fit.auto.arima, h = 96*10)
par(mfrow=c(1,1))
plot(fore.auto.arima)
ggtsdisplay(residuals(fit.auto.arima))

#
#### forecasting (no covariates) - daily data ####

## Box-Jenkins method
# A. Model form selection
# - Evaluate stationarity
# - Selection of the differencing level (d) – to fix stationarity problems
# - Selection of the AR level (p)
# - Selection of the MA level (q)
# B. Parameter estimation
# C. Model checking
# D. Forecast

### daily data ###

## A. Model form selection - Evaluate stationarity and select differencing level (if needed)
# The basic stationarity diagnostics are the following:
# 1. Plot your data. 
plot(C2_no3_daily_ts)
# Look for
#  - An increasing or decreasing trend: Possibly decreasing
#  - A non-zero mean (if no trend): Yes
#  - Does it look like it might be stationary around a trend?: No
#  - Strange shocks or steps in your data: No


# 2. Apply stationarity tests
#  - adf.test() p-value less than 0.05 (reject null) suggests ts is stationary

adf.test(C2_no3_daily_ts_contig) # does not allow NAs
adf.test(na.interp(C2_no3_daily_ts)) 
adf.test(C2_no3_daily_ts_contig) 
adf.test(na.interp(C2_no3_daily_ts), k=0) 
adf.test(C2_no3_daily_ts_contig, k=0) 

#  - kpss.test() p-value greater than 0.05 (do not reject null) suggests ts is stationary

kpss.test(C2_no3_daily_ts_contig, null="Level") 
kpss.test(C2_no3_daily_ts_contig, null="Trend") 

# 3. If stationarity tests are failed, then try differencing to correct
#  - Try ndiffs() in the forecast package or manually try different differences.

kpss.test(diff(C2_no3_daily_ts_contig, differences = 1), null="Level") 
kpss.test(diff(C2_no3_daily_ts_contig, differences = 1), null="Trend") 

# suggests theres no trend

## A. Model form selection - Selection of the AR level (p) and AR level (p)

fit.auto.arima = forecast::auto.arima(na.contiguous(C2_no3_daily_ts), start.p = 0, max.p = 3, start.q = 0, max.q = 3)
fit.auto.arima # see best fit parameters

fit.sarima = sarima(na.contiguous(C2_no3_daily_ts), p=0, d=1, q =0) # run without saving object for diagnostic plots

fit.arima = arima(na.contiguous(C2_no3_daily_ts), order = c(0, 1, 0))
par(mfrow=c(2,1)); Acf(resid(fit.arima)); Pacf(resid(fit.arima))

fit.Arima = forecast::Arima(na.contiguous(C2_no3_daily_ts), order = c(0, 1, 0), include.constant = F)
par(mfrow=c(2,1)); Acf(resid(fit.Arima)); Pacf(resid(fit.Arima))

# great fit!

## Forecast!
par(mfrow=c(1,1))
fore.auto.arima <- forecast::forecast(fit.auto.arima, h = 10)
plot(fore.auto.arima)
fore.arima <- forecast::forecast(fit.arima, h = 10)
plot(fore.arima)
fore.Arima <- forecast::forecast(fit.Arima, h = 10)
plot(fore.Arima)


#### wavelet decomposition ####

# time series should be an ordered dataframe with timesteps in "date" column, data in 2nd column, and no NAs
for.wavelet = C2_no3_filled_splineinterp %>%
  na.contiguous() %>%
  select(date_timeAK, nitrate_uM_c_bc)%>%
  rename(date = date_timeAK, x =nitrate_uM_c_bc) %>%
  arrange(date)
# check for NAs
any(is.na(for.wavelet))

# check out options for generating the null hypotheses time series
?analyze.wavelet # look at 'method' argument

C2_no3_15min_wavelet = 
  analyze.wavelet(for.wavelet, 
                my.series = 2, 
                method="AR", # null hypothesis is red noise
                params = list(3),
                dt=1/96, # change dt to match time steps per period. 1/96 for 15 min data.
                dj=1/20, # changes the resolution of the analysis across the frequency domain. default usually works well.
                lowerPeriod=1/4, #probably no need for less than 0.25 d
                make.pval=T, n.sim=100,# this should be something like n.sim=100
                date.format = "%b-%d", date.tz = "America/Anchorage") 

## Plot of wavelet power spectrum ##
?wt.image

par(mfrow=c(1,1))
wt.image(C2_no3_15min_wavelet, 
         color.key="quantile", 
         label.time.axis=TRUE, 
         plot.ridge=FALSE, 
         show.date = TRUE, 
         date.format = "%b-%d", 
         legend.params=list(lab="wavelet power levels"), 
         main="Nitrate (15 min) in C2", 
         graphics.reset = T)

## Plot of global wavelet power: ##
?wt.avg 

wt.avg(C2_no3_15min_wavelet,  show.siglvl = TRUE)

## Plot together
par(mfrow=c(1,2))
wt.image(C2_no3_15min_wavelet, 
         color.key="quantile", 
         label.time.axis=TRUE, 
         plot.ridge=FALSE, 
         show.date = TRUE, 
         date.format = "%b-%d", 
         legend.params=list(lab="wavelet power levels"), 
         main="Nitrate (15 min) in C2", 
         graphics.reset = F)
wt.avg(C2_no3_15min_wavelet,  show.siglvl = TRUE)

