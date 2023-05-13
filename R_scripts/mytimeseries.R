####my time series####
#### read me ####

# the purpose of this script is to demonstrate time series exploration, simple ARIMA model fitting, forecasting, and wavelet decomposition. There is also some time series data wrangling. 
# this code works very time its redited depending on the varibles of interest. 
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

setwd("C:/Users/marag.DESKTOP-NEKBNJE/OneDrive/Desktop/Stakeholder data analysis/test1/Data")
dat = read.csv("anionsf.csv")
dat <- read.csv("anionsf.csv", na.strings = c("."," ","NA"))

#dat<- na.omit(dat)

str(dat)
# reduce to 3 sites for the purposes of most of this tutorial
dat = dat[dat$Samplingtype =="Center"| dat$Samplingtype =="River"| dat$Samplingtype =="Ditch",]
dat = dat[dat$SiteName =="Alameda"| dat$SiteName =="Badger"| dat$SiteName =="Savannah" | dat$SiteName =="Bobcat",]

#remove time sampling column
dat <- select(dat, - c("Timesampling", "Timeofsitearrival", "Bromide"))

# format date/time: this code worked well to format date, i have dateF and Date
dat$DateF = as.POSIXct(dat$Date, format="%m/%d/%Y %H:%M", tz="US/Mountain")
dat$Date = as.POSIXct(dat$Date, format="%m/%d/%Y %H:%M", tz="US/Mountain")
dat$Datecombined = as.POSIXct(dat$Datecombined, format="%m/%d/%Y %H:%M", tz="US/Mountain")
str(dat)
#dat$datetime_NM = as.POSIXct(dat$Sample_DateTime, format="%m/%d/%y %H:%M", tz="US/Mountain")

# convert characters that should be factors (categories) to factors
dat$SiteName = as.factor(dat$SiteName)
dat$Samplingtype = as.factor(dat$Samplingtype)

# convert water quality data to numeric
#dat$Bromide = as.numeric(dat$Bromide)
dat$Fluoride = as.numeric(dat$Fluoride)
dat$NH4N  = as.numeric(dat$NH4N )
dat$Sulfate = as.numeric(dat$Sulfate)
dat$Chloride = as.numeric(dat$Chloride)
dat$PhosphateP = as.numeric(dat$PhosphateP)
dat$Turbidity = as.numeric(dat$Turbidity)
dat$condTemp = as.numeric(dat$condTemp)
dat$DO  = as.numeric(dat$DO)
dat$NitrateN  = as.numeric(dat$NitrateN)
dat$pH1   = as.numeric(dat$pH1)
dat$pHTemp    = as.numeric(dat$pHTemp)
dat$Conductivity    = as.numeric(dat$Conductivity )
dat$DOTemp   = as.numeric(dat$DOTemp)
dat$DoPercent  = as.numeric(dat$DoPercent)
dat$Airtemp   = as.numeric(dat$Airtemp)
dat$Day   = as.numeric(dat$Day)
dat$Nitrite  = as.numeric(dat$Nitrite)

#datS <- select(dat, c("SiteName","Airtemp", "DoPercent", "DOTemp", "Conductivity", "condTemp" ,"pH1" , "Chloride", "Sulfate"))
#str(datS)

#Means by year for variables

#dat_means <- dat %>%
 # group_by(Year, SiteName) %>%
  #summarize(across
   #         Meantemp = mean(Airtemp, na.rm=TRUE))

#str(dat_means)

# check that date/time formatted correctly
class(dat$Date); tz(dat$Date)

# sitename = Savannah, variable = phosphate 
Savr_phos = dat %>% 
  filter(SiteName=="Savannah") %>% 
  select(DateF, PhosphateP)

#Savr_phos<- na.omit(Savr_phos)

#### explore filling gaps ####

# general guidelines:
# - if you're interested in diel patters, be very cautious about filling gaps > 6 h
# - if you're interested in day to day patters, case by case basis - examine whether gaps are during events like storms
# - in general, examine time series before and after you fill gaps and ask whether you've influenced patterns of interest in a defensible way or not
# note that your data must be evenly spaced in time before doing this!!

### 1 month data ###

## fill with linear interpolation
# Make univariate zoo time series#
ts.temp<-read.zoo(Savr_phos, index.column=1, format="%Y-%m-%d %H:%M:%S", tz="US/Mountain")
#plot
plot.zoo(ts.temp)
# Apply NA interpolation method
Savr_phos_filled_linearinterp = na.approx(ts.temp, na.rm = T, maxgap = 1)
# revert back to df
Savr_phos_filled_linearinterp = as.data.frame(Savr_phos_filled_linearinterp)
Savr_phos_filled_linearinterp$Date = as.POSIXct(row.names(Savr_phos_filled_linearinterp), tz="US/Mountain")
names(Savr_phos_filled_linearinterp) = c(colnames(Savr_phos)[2],colnames(Savr_phos)[1])


###No need for interpolation because the data is not rich at all, very few observations

plot.zoo(Savr_phos_filled_linearinterp)

# fill with spline interpolation
# Make univariate zoo time series #
ts.temp<-read.zoo(Savr_phos, index.column=1, format="%Y-%m-%d %H:%M:%S", tz="US/Mountain")
# Apply NA interpolation method
Center_Sul_filled_splineinterp = na.spline(ts.temp, na.rm = T, maxgap = 1)
# revert back to df
Center_Sul_filled_splineinterp = as.data.frame(Center_Sul_filled_splineinterp)
Center_Sul_filled_splineinterp$DateF = as.POSIXct(row.names(Center_Sul_filled_splineinterp), tz="US/Mountain")
names(Center_Sul_filled_splineinterp) = c(colnames(Savr_phos)[2],colnames(Center_Sul)[1])


# explore ?na.approx in zoo package for  more options

## fill with arima model + kalman smoother
# Document gaps  
largegap.num = 1 # its one day
is.na.rle <- rle(is.na(Center_Sul$Center))
is.na.rle$Sulfate <- is.na.rle$Sulfate & is.na.rle$lengths >= (largegap.num)
biggaps = Center_Sul[inverse.rle(is.na.rle), ]
tz(biggaps$Date) = "US/Mountain"
biggaps = subset(biggaps, select = "Date")
# Make univariate time series: covert to zoo, then to ts #
xts<-read.zoo(Center_Sul, index.column=1, format="%Y-%m-%d %H:%M:%S", tz="US/Mountain")
xts<-as.xts(xts)
# remove leading and trailing NAs #
xts = na.trim(xts, is.na="any")
# Apply auto.arima and kalman filter to impute missing Sulfate #
fit = auto.arima(xts) 
kal = KalmanSmooth(xts, fit$model)
id.na<-which(is.na(xts))
for(i in id.na) {
  xts[i]<-fit$model$Z %*% kal$smooth[i,]}
# revert to dataframe #
Center_Sul_filled_ARIMA = as.data.frame((xts))
Center_Sul_filled_ARIMA$DateF = as.POSIXct(row.names(Center_Sul_filled_ARIMA), tz="US/Mountain")
names(Center_Sul_filled_ARIMA) = c("Sulfate", "DateF")
# remove large gaps # 
Center_Sul_filled_ARIMA$Sulfate[Center_Sul_filled_ARIMA$Date %in% as.POSIXct(biggaps$Date)] = NA


par(mfrow=c(4,1))
xmin = min(Center_Sul$DateF); xmax = max(Center_Sul$Date)
plot(Center_Sul$Center ~ Center_Sul$Date, type="l", main="original", xlim=c(xmin, xmax))
plot(Center_Sul_filled_linearinterp$Center ~ Center_Sul_filled_linearinterp$DateF, type="l", main="linear interpolation", xlim=c(xmin, xmax))
plot(Center_Sul_filled_splineinterp$Center ~ Center_Sul_filled_splineinterp$DateF, type="l", main="spline interpolation", xlim=c(xmin, xmax))
plot(Center_Sul_filled_ARIMA$Center ~ Center_Sul_filled_ARIMA$DateF, type="l", main="ARIMA interpolation", xlim=c(xmin, xmax))


par(mfrow=c(4,1))
xmin = min(Center_Sul$DateF); xmax = as.POSIXct("2011-08-01")
ymin= 1; ymax=50
plot(Center_Sul$Sulfate ~ Center_Sul$Date, type="l", main="original", xlim=c(xmin, xmax), ylim=c(ymin, ymax))
plot(Center_Sul_filled_linearinterp$Sulfate ~ Center_Sul_filled_linearinterp$Date, type="l", main="linear interpolation", xlim=c(xmin, xmax), ylim=c(ymin, ymax))
plot(Center_Sul_filled_splineinterp$Sulfate ~ Center_Sul_filled_splineinterp$Date, type="l", main="spline interpolation", xlim=c(xmin, xmax), ylim=c(ymin, ymax))
plot(Center_Sul_filled_ARIMA$Sulfate ~ Center_Sul_filled_ARIMA$Date, type="l", main="ARIMA interpolation", xlim=c(xmin, xmax), ylim=c(ymin, ymax))


# create time series object
# note that dates must be sequential before doing this!! xts and other ts-making functions will assume that they are. 
Center_Sul_xts = xts(Center_Sul_filled_splineinterp$Sulfate, order.by = Center_Sul_filled_splineinterp$DateF)
class(Center_Sul_xts)
plot(Center_Sul_xts)
summary(Center_Sul_xts)


# decompose() and stl() require a ts() or mtst() object
# I might use mtst() if this were multi-annual data with a seasonal signal to capture both diel and seasonal frequencies
Center_Sul_ts = ts(Center_Sul_filled_splineinterp$Sulfate, 
               frequency = 1, 
               start=hour(min(Center_Sul_filled_splineinterp$DateF)))
print(Center_Sul_ts, calendar = T)
# create one for the contiguous time series
Center_Sul_ts_contig = na.contiguous(Center_Sul_filled_splineinterp)
Center_Sul_ts_contig = ts(Center_Sul_ts_contig$Sulfate, 
                      frequency = 1, 
                      start=hour(min(Center_Sul_ts_contig$DateF)))


# aggregate to daily and create time series object
Center_Sul_daily = Center_Sul_filled_splineinterp %>% 
  mutate(day = as.Date(DateF)) %>%
  select(day, Sulfate) %>%
  group_by(day) %>%
  summarise(Sulfateum = mean(Sulfate, na.rm = T))
Center_Sul_daily_xts = xts(Center_Sul_daily$Sulfate, order.by = Center_Sul_daily$day)
plot(Center_Sul_daily_xts)
summary(Center_Sul_daily_xts)

# 12 NAs
Center_Sul_daily_ts = ts(Center_Sul_daily$Sulfate,  
                     frequency = 1, 
                     start=yday(min(Center_Sul_daily$day)))
print(Center_Sul_daily_ts, calendar = T)
# create one for the contiguous time series
Center_Sul_daily_ts_contig = na.contiguous(Center_Sul_daily)
Center_Sul_daily_ts_contig = ts(Center_Sul_daily_ts_contig$nitrate_uM, 
                            frequency = 1, 
                            start=yday(min(Center_Sul_daily_ts_contig$day)))


#### explore autocorrelation structure ####

# reset plotting window
par(mfrow=c(1,1))

# examine acf and pacf for Monthly  data
Acf(Center_Sul_xts, na.action = na.interp, lag.max = 7 )
Pacf(Center_Sul_xts, na.action = na.interp, lag.max = 7)
Acf(Center_Sul_xts, na.action = na.contiguous, lag.max = 7)
Pacf(Center_Sul_xts, na.action = na.contiguous, lag.max = 7)
# strong autocorrelation at lag 1, 2, 3


#### classic decomposition ####


## 1 month data ##
# decompose into additive components
plot(decompose(Center_Sul_ts_contig))
plot(decompose(na.interp(Center_Sul_ts)))
plot(stl(Center_Sul_ts, s.window = 1, na.action = na.contiguous))
plot(stl(Center_Sul_ts, s.window = 1, na.action = na.interp))
# decompose into multiplicative components
plot(decompose(na.contiguous(Center_Sul_ts), type="multiplicative"))
plot(decompose(na.interp(Center_Sul_ts), type="multiplicative"))
# extract components from additive
Center_Sul_ts_decomp = stl(Center_Sul_ts, s.window = 1, na.action = na.contiguous)
Center_Sul_ts_diel = Center_Sul_ts_decomp$time.series[,"seasonal"]
plot(Center_Sul_ts_diel)
Center_Sul_ts_trend = Center_Sul_ts_decomp$time.series[,"trend"]
plot(Center_Sul_ts_trend)
Center_Sul_ts_noise = Center_Sul_ts_decomp$time.series[,"remainder"]
plot(Center_Sul_ts_noise)

plot(Center_Sul_ts_trend + Center_Sul_ts_noise)

## daily data ##
# decompose
decompose(Center_Sul_daily_ts)
plot(stl(Center_Sul_daily_ts, s.window = 1, na.action = na.contiguous))
plot(stl(Center_Sul_daily_ts, s.window = 1, na.action = na.interp))
# this does not work because there is no seasonal pattern. You can't decompose data that doesn't have a seasonal component! 



