####my time series####
#### read me ####

# the purpose of this script is to demonstrate time series exploration, simple ARIMA model fitting, forecasting, and wavelet decomposition. There is also some time series data wrangling. 

#### libraries ####

library(tidyverse)
library(lubridate)
library(forecast)data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAWElEQVR42mNgGPTAxsZmJsVqQApgmGw1yApwKcQiT7phRBuCzzCSDSHGMKINIeDNmWQlA2IigKJwIssQkHdINgxfmBBtGDEBS3KCxBc7pMQgMYE5c/AXPwAwSX4lV3pTWwAAAABJRU5ErkJggg==
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

dat = read.csv("anions.csv")

str(dat)



# test date formating
head(dat$Date)
head(as.POSIXct(dat$Date, format="%m/%d/%y %H:%M:%S", tz="MST"))
tail(dat$Date)
tail(as.POSIXct(dat$Date, format="%m/%d/%y %H:%M:%S", tz="MST"))

# format date/time
dat$Date = as.POSIXct(dat$Date, format="%m/%d/%y %H:%M:%S", tz="MST")

# create new date/time col with correct format (always keep the original!)
dat$Date = as.POSIXct(data$DateF, format="%m/%d/%y %H:%M:%S", tz="MST")


# check that date/time formatted correctly
class(dat$Date); tz(dat$Date)

# convert characters that should be factors (categories) to factors
dat$Site.Name = as.factor(dat$Site.Name)
dat$Sampling.type = as.factor(dat$Sampling.type)
#dat$Is_Nondetect = as.factor(dat$Is_Nondetect)

# convert water quality data to numeric
dat$Bromide..mg.L. = as.numeric(dat$Bromide..mg.L.)
dat$Fluoride..mg.L. = as.numeric(dat$Fluoride..mg.L.)
dat$NH4.N..mg.L. = as.numeric(dat$NH4.N..mg.L.)
dat$Sulfate..mg.L. = as.numeric(dat$Sulfate..mg.L.)
dat$Chloride..mg.L. = as.numeric(dat$Chloride..mg.L.)
dat$Phosphate.P..mg.L. = as.numeric(dat$Phosphate.P..mg.L.)
dat$Turbidity..NTU. = as.numeric(dat$Turbidity..NTU.)
dat$conductivity.Temp...C. = as.numeric(dat$conductivity.Temp...C.)
dat$Do....  = as.numeric(dat$Do.... )


# reduce to 3 sites for the purposes of most of this tutorial
dat = dat[dat$Site.Name =="Savannah"| dat$Site.Name =="Calabacillas"| dat$Site.Name =="Minnow",]

str(dat)

# site = Savannah, variable = sulphate (mg.L.)
Savnah_Sul = dat %>% 
  filter(Site.Name=="Savannah") %>% 
  select(Date, Sulfate..mg.L.)

#### explore filling gaps ####

# general guidelines:
# - if you're interested in diel patters, be very cautious about filling gaps > 6 h
# - if you're interested in day to day patters, case by case basis - examine whether gaps are during events like storms
# - in general, examine time series before and after you fill gaps and ask whether you've influenced patterns of interest in a defensible way or not
# note that your data must be evenly spaced in time before doing this!!

### 15 min data ###

## fill with linear interpolation
# Make univariate zoo time series#
ts.temp<-read.zoo(Savnah_Sul, index.column=1, format="%Y-%m-%d %H:%M:%S", tz="MST")
# Apply NA interpolation method
Savnah_Sul_filled_linearinterp = na.approx(ts.temp, na.rm = T, maxgap = 24*0.25)
# revert back to df
Savnah_Sul_filled_linearinterp = as.data.frame(Savnah_Sul_filled_linearinterp)
Savnah_Sul_filled_linearinterp$date_timeAK = as.POSIXct(row.names(Savnah_Sul_filled_linearinterp), tz="MST")
names(Savnah_Sul_filled_linearinterp) = c(colnames(Savnah_Sul)[2],colnames(Savnah_Sul)[1])
