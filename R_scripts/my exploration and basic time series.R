#### read me ####

# The purpose of this script is to explore my data, including
# - describing dataset size (# variables & # observations)
# - describing data types
# - checking data distributions
# - checking for spatial autocorrelation
# - checking for temporal autocorrelation
# - checking for correlation between variables


#### libraries ####

library(tidyverse)
library(lubridate)
library(psych) # to plot pair-wise correlations
library(car) # I like their qq plot fxn
library(tsibble) # useful for creating time series objects
library(forecast) # I like their Acf fxn
library(ape) # for spatial autocorrelation
library(ade4)# for spatial autocorrelation
library(rgdal) # for mapping
library(ggplot2) # for plotting
library(dplyr)
library(magrittr)


#### load and tidy data ####

setwd("C:/Users/marag.DESKTOP-NEKBNJE/OneDrive/Desktop/Stakeholder data analysis/test1/Data")
dat = read.csv("anionsf.csv")
dat <- read.csv("anionsf.csv", na.strings = c("."))
#dat<- na.omit(dat)

str(dat)
# reduce to 3 sites for the purposes of most of this tutorial
#dat = dat[dat$Samplingtype =="Center"| dat$Samplingtype =="River"| dat$Samplingtype =="Ditch",]
dat = dat[dat$Samplingtype =="River"| dat$Samplingtype =="Ditch"| dat$Samplingtype =="Center",]
#dat = dat[dat$SiteName =="Alameda"| dat$SiteName =="Badger"| dat$SiteName =="Harrison",]


#remove time sampling column
dat <- select(dat, - c("Timesampling", "Timeofsitearrival"))

str(dat)


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
dat$Bromide = as.numeric(dat$Bromide)
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

str(dat)

#try filtering for one site
dat = dat[dat$SiteName =="Bobcat"| dat$SiteName =="Badger"| dat$SiteName =="Savannah" |dat$SiteName =="Harrison" | dat$SiteName =="Belen",]


#### simple plotting basic time series ####
#sulfate
library(ggplot2)
ggplot(dat = dat, aes(x=Year, y= Sulfate, color = SiteName))+
  geom_point(position = position_jitter(w = 0.5, h = 0.2)) +
  geom_smooth()+
  xlab("year") +
  ylab ("Sulphate")+
  labs(title = "Sulphate over years")
print()


#chloride
library(ggplot2)
ggplot(dat = dat, aes(x=Year, y= Chloride, color = SiteName))+
  geom_point(position = position_jitter(w = 0.5, h = 0.2)) +
  geom_smooth()+
  xlab("year") +
  ylab ("chloride")+
  labs(title = "chloride over years")
print()

#phosphate
library(ggplot2)
ggplot(dat = dat, aes(x=Year, y= PhosphateP, color = SiteName))+
  geom_point(position = position_jitter(w = 0.5, h = 0.2)) +
  geom_smooth()+
  xlab("year") +
  ylab ("phosphate")+
  labs(title = "phosphate over years")
print()


#Nh4n
library(ggplot2)
ggplot(dat = dat, aes(x=Year, y= NH4N, color = SiteName))+
  geom_point(position = position_jitter(w = 0.5, h = 0.2)) +
  geom_smooth()+
  xlab("year") +
  ylab ("NH4N")+
  labs(title = "NH4N over years")
print()

#### describe dataset size and structure ####

head(dat)
str(dat)

with(dat, table(Samplingtype, SiteName ))

### check timesteps by looking and time series of most frequently collected parameters
# make dataset of one of the most frequently collected parameters
dat_lts = 
  dat %>% 
  group_by(Samplingtype) %>% 
  filter(n() > 1) %>% 
  arrange(DateF)



dat_lts_center = 
  dat %>% 
  filter(Samplingtype =="Center") %>% 
  arrange(DateF)
# add year and day of year for plotting
dat_lts_center$year = lubridate::year(dat_lts_center$DateF)
dat_lts_center$day = lubridate::day(dat_lts_center$DateF)
# plot
ggplot(data=dat_lts_center, aes(x=Day, y=PhosphateP, color=SiteName))+
  geom_point() + geom_path()+ 
  facet_wrap(~year, scales="free_y")+
  theme(legend.title = element_blank()) +
  theme_bw()
# can also look at single years in detail - Not sure where am messing, i excluded x limits
ggplot(data=dat_lts_center, aes(x=Year, y=PhosphateP, color=SiteName))+
  geom_point() + geom_path()+
  #xlim(c(as.POSIXct("2006-07-01"), as.POSIXct("2011-07-011")))+
  theme(legend.title = element_blank()) +
  theme_bw()
# timesteps are all over the place from year to year and site to site, what we would call "irregular"
# timesteps are not sub-daily, at most frequent are approximately monthly

# ........ etc. for each parameter I'm interested in using in analysis .........


### How many variables are in your dataset?
str(dat)
# 25 parameters

### How many observations are in your dataset?
nrow(dat)
# 261 total
with(dat, table(Samplingtype , SiteName))# Select sites with enough sampling types
range(with(dat, table(Samplingtype , SiteName)))
# there are a variable # of observations for each water quality parameter in each site, from 0 to 21 total

### Are the data nested in time or space?
# Yes in time - observations were collected repeatedly on an irregular schedule
# Yes in space - observations were collected in three different sites, need more research/exploration to find out if sites are connected in any way

#### describe data types ####

str(dat)
summary(dat$Samplingtype)
# most water quality parameters are numerical continous ratios

#### check distributions ####

# I'm only going to check the distributions of data with at least 100 obs in each site, as I am unlikely to analyze less frequently gathered data
dat_r = 
  dat %>% 
  group_by(Samplingtype, SiteName) %>% 
  filter(n() > 10) %>% 
  arrange(DateF)
summary(dat_r$Samplingtype)
library(ggplot2)
temp = dat_r[dat_r$Samplingtype == "Center",]
qqPlot(temp$Sulfate); shapiro.test(temp$Sulfate) # normal
qqPlot(temp$Sulfate[temp$SiteName=='Bobcat']); shapiro.test(temp$Sulfate[temp$SiteName=='Bobcat']) # normal,few outliers 10 and 7
qqPlot(temp$Sulfate[temp$SiteName=='Badger']); shapiro.test(temp$Sulfate[temp$SiteName=='Badger']) # normal
qqPlot(temp$Sulfate[temp$SiteName=='Harrison']); shapiro.test(temp$Sulfate[temp$SiteName=='Harrison']) # normal

temp = dat_r[dat_r$Samplingtype == "Ditch",]
qqPlot(temp$Sulfate); shapiro.test(temp$Sulfate) # normal
qqPlot(temp$Sulfate[temp$SiteName=='Bobcat']); shapiro.test(temp$Sulfate[temp$SiteName=='Bobcat']) # normal 
qqPlot(temp$Sulfate[temp$SiteName=='Badger']); shapiro.test(temp$Sulfate[temp$SiteName=='Badger']) # normal
qqPlot(temp$Sulfate[temp$SiteName=='Harrison']); shapiro.test(temp$Sulfate[temp$SiteName=='Harrison']) # normal


temp = dat_r[dat_r$Samplingtype == "River",]
qqPlot(temp$Sulfate); shapiro.test(temp$Sulfate) # normal
qqPlot(temp$Sulfate[temp$SiteName=='Bobcat']); shapiro.test(temp$Sulfate[temp$SiteName=='Bobcat']) # not showing, few outliers 10 and 7
qqPlot(temp$Sulfate[temp$SiteName=='Badger']); shapiro.test(temp$Sulfate[temp$SiteName=='Badger']) # normal
qqPlot(temp$Sulfate[temp$SiteName=='Harrison']); shapiro.test(temp$Sulfate[temp$SiteName=='Harrison']) # normal




# etc........ for the rest of the parameters that I think I'll use in this analysis

### Examine non-normal data closely ###
# ask:
# are outliers making it non-normal?
# can I justify removing outliers based on my knowledge of the data?
# if data is still non-normal, what distribution is it?


temp = dat_r[dat_r$Samplingtype == "Center",]
summary(temp$Bromide)
hist(temp$Bromide)
plot(density(temp$Bromide)) # not showing !


# this data has 1 an extreme negative outlier. Center Bromide cannot be negative, so this is an error. I will remove it in the main datasets an re-check the data's normality
dat$Bromide[dat$Samplingtype=="Center" & dat$Bromide<0] = NA # rplace it in main dataset
dat_r$Bromide[dat_r$Samplingtype =="Center" & dat_r$Bromide<0] = NA # replace it in reduced dataset
temp = dat_r[dat_r$Samplingtype == "Center",]
qqPlot(temp$Bromide); shapiro.test(temp$Bromide)
# there is now a high outlier to examine
temp = dat_r[dat_r$Samplingtype == "Center",]
summary(temp$Bromide)
hist(temp$Bromide)
plot(density(temp$Bromide, na.rm = T))
# this data has 1 an extreme positive outlier. Center Bromides do not get this high in natural conditions. This is coal data so maybe it isn't natural, but even still, we'd expect to see more than one point if this were not an error. I will remove it in the main datasets an re-check the data's normality
dat$Bromide[dat$Sampling.type=="Center" & dat$Bromide>11] = NA # rplace it in main dataset
dat_r$Bromide[dat_r$Sampling.type=="Center" & dat_r$Bromide>11] = NA # replace it in reduced dataset
temp = dat_r[dat_r$Samplingtype == "Center",]
qqPlot(temp$Bromide); shapiro.test(temp$Bromide)
hist(temp$Bromide)
plot(density(temp$Bromide, na.rm = T))
range(temp$Bromide, na.rm = T)

# still not normal!
# this looks like a lognormal, Gamma, or Weibull distribution
# it is bounded above zero and is right-skewed
# what happens if I log-transform it?
temp = dat_r[dat_r$Samplingtype == "Center",]
qqPlot(log10(temp$Bromide)); shapiro.test(log10(temp$Bromide))


#### check for temporal autocorrelation ####

# I'm going to check these one site at a time and only of data with at least 100 (i used 10)obs in each site, as I am unlikely to analyze less frequently gathered data
dat_r = 
  dat %>% 
  group_by(Samplingtype, SiteName) %>% 
  filter(n() > 10) %>% 
  arrange(DateF)
summary(dat_r$Samplingtype) # note which parameters are left after filtering
summary(dat_r$SiteName) # note that VR-1 no longer has any observations, so I will focus on the other two sites

# checking for temporal autocorrelation requires the data to be a time series object (read ?ts for details on this)
# To achieve this, I need regularly spaced data. This data is irregularly spaced, approximately monthly, but sometimes there are more than one observations per month or fewer
# I will start by averaging observations within the same month:

# this code not running well, review afterwards

dat_monthly = dat_r

dat_monthly = 
  dat_r %>%
  mutate(yr = lubridate::year(Datecombined)) %>%
  mutate(mo = lubridate::month(Datecombined)) %>%
  dplyr::select(SiteName, Samplingtype, yr, mo, Bromide) %>%
  group_by(Bromide.mn = mean(Bromide, na.rm = T)) %>%
  mutate(date = paste(yr, mo, "7", Jul="-")) %>%
  mutate(date = as.Date(date))


#### Sulfate in Bobcat center
### subset data to be one site and one sampling ty[e]
temp = dat_monthly[dat_monthly$Samplingtype == "Center" & dat_monthly$SiteName=="Bobcat" ,]
### make this a time series object
## first, make doubly sure that the data is arranged by time before converting to ts object!
temp = temp %>% arrange(DateF) 
## second, make the spacing of dates consistent and fill in missing obs with NA. This is a handy fxn. You can also create a df of evenly spaced dates and left_join the data to this.
temp_ts =
  temp %>% 
  complete(DateF = seq(min(date), max(DateF), by = "1 month"), 
           fill = list(Sulfate = NA)) %>%
  as_tsibble(index = DateF)


## finally, convert to a ts object
# a ts object is a vector of data taken sequentially through time. Required arguments are:
# - the data vector
# - the frequency, which is the number of observations per unit of time. Lots of ways to specify this. For monthly data, you can put in 12 and it will assume that's 12 obs in a year. Google for help for other frequencies.
# - the start, which specifies when the first obs occured. Lots of ways to specify this. For monthly data, you can put in c(year, month) and it will know what you mean. 
head (temp_ts)
temp_ts = ts(temp_ts$Sulfate.mn, frequency=12, start=c(2005, 10)) 
# check that you specified the ts correctly
print(temp_ts, calendar = T) 
### now we're ready to check for temporal autocorrelation in this ts!
# I prefer the forecast pkg's Acf fxn over base R acf() because Acf() doesn't include 0 (which is always 1) and shows month #s by default instead of decimal years. Note the different options for dealing with NAs and how this changes the results (see ?na.fail and ?Acf for details). 
forecast::Acf(temp_ts, na.action = na.pass) 
forecast::Acf(temp_ts, na.action = na.contiguous) 
forecast::Acf(temp_ts, na.action = na.interp)

forecast::Pacf(temp_ts, na.action = na.pass)
forecast::Pacf(temp_ts, na.action = na.contiguous)
forecast::Pacf(temp_ts, na.action = na.interp)

# acf tells me that there is temporal autocorrelation. The sin-wave-like pattern is typical of a ts impacted by seasonality
# pcaf tells me that strongest source of autocorrelation is at lag 1, which indicates a random walk/AR1 process. There is possibly ac at other lags, depending on how NAs are handled. 


#### Sulphate in VR-3
### subset data to be one site and one parameter
temp = dat_monthly[dat_monthly$Samplingtype == "River" & dat_monthly$SiteName=="Bobcat" ,]
### make this a time series object
## first, make doubly sure that the data is arranged by time before converting to ts object!
temp = temp %>% arrange(date) 
## second, make the spacing of dates consistent and fill in missing obs with NA. This is a handy fxn. You can also create a df of evenly spaced dates and left_join the data to this.
temp_ts =
  temp %>% 
  complete(date = seq(min(date), max(date), by = "1 month"), 
           fill = list(Sulfate = NA)) %>%
  as_tsibble(index = date)
## finally, convert to a ts object
# a ts object is a vector of data taken sequentially through time. Required arguments are:
# - the data vector
# - the frequency, which is the number of observations per unit of time. Lots of ways to specify this. For monthly data, you can put in 12 and it will assume that's 12 obs in a year. Google for help for other frequencies.
# - the start, which specifies when the first obs occured. Lots of ways to specify this. For monthly data, you can put in c(year, month) and it will know what you mean. 
head (temp_ts)
temp_ts = ts(temp_ts$Value.mn, frequency=12, start=c(2005, 1)) 
# check that you specified the ts correctly
print(temp_ts, calendar = T) 
### now we're ready to check for temporal autocorrelation in this ts!
# I prefer the forecast pkg's Acf fxn over base R acf() because Acf() doesn't include 0 (which is always 1) and shows month #s by default instead of decimal years. Note the different options for dealing with NAs and how this changes the results (see ?na.fail and ?Acf for details).
forecast::Acf(temp_ts, na.action = na.pass) 
forecast::Acf(temp_ts, na.action = na.contiguous) 
forecast::Acf(temp_ts, na.action = na.interp) 
forecast::Pacf(temp_ts, na.action = na.pass)
forecast::Pacf(temp_ts, na.action = na.contiguous)
forecast::Pacf(temp_ts, na.action = na.interp)

# acf tells me that there is temporal autocorrelation. The sin-wave-like pattern is typical of a ts impacted by seasonality
# pcaf tells me that strongest source of autocorrelation is at lag 1, which indicates a random walk/AR1 process. There is possibly ac at other lags, depending on how NAs are handled. 


# ....... ect. for each parameter and site combination I might include in the analysis .......



#### check for spatial autocorrelation ####

# I'm interested in spatial and not temporal autocorrelation, so I am going to look at just a few observations across all sites

# reload and format data with all sites
dat_all <- dat

#dat_all = read.csv("anionsf")
dat_all$Datecombined = as.POSIXct(dat_all$Datecombined, format="%m/%d/%Y", tz="US/Mountain")
dat_all$Sampling.type = as.factor(dat_all$Samplingtype)
dat_all$Site.Name = as.factor(dat_all$SiteName)
dat_all$Chloride = as.factor(dat_all$Chloride)
dat_all$Chloride = as.numeric(dat_all$Chloride)
# how many sites are there?
length(unique(dat_all$SiteName))
# 3

# what parameters were collected across all sites in Nov 2005? shld s
dat_Nov2005 = dat_all[dat_all$Datecombined >= as.POSIXct("2005-11-07") &
                         dat_all$Datecombined < as.POSIXct("2005-11-07"),]
tb = as.data.frame( with(dat_Nov2005, table(SiteName, Samplingtype)) )
tb = tb[tb$Freq>0,]
tb2 = tb %>% group_by(Samplingtype) %>% summarise(n = n()) %>% arrange(desc(n))
head(tb2, 15) # I have 0 observations for 2 variables

#I may not test spatial autocorrelation

# ^ these are not good options for testing for spatial autocorrelation

### Savannah in Nov 2005
dat_Nov2005 = dat_all[dat_all$Datecombined >= as.POSIXct("2005-11-07") &
                         dat_all$Datecombined < as.POSIXct("2005-11-07"),]
temp = dat_Nov2005 %>%  filter(Sampling.type =="Savannah")
# randomly generate lat/lon for demo
set.seed(42)
temp$lat = runif(nrow(temp),35.090956,35.634117)
temp$lon = runif(nrow(temp),-107.65829,-106.65829)
## Moran.I
# generate an inverse distance matrix 
dists = as.matrix(dist(cbind(temp$lon, temp$lat)))
dists.inv = 1/dists
diag(dists.inv) = 0
# calculate Moran.I
Moran.I(temp$Chloride, dists.inv)
# we can NOT reject the null hypothesis that there is zero spatial autocorrelation present. In other words, there doesn't seem to be a lot of spatial autocorrelation. 
## Mantel test
# generate spatial distance matrix
site_dists = dist(cbind(temp$lon, temp$lat))
# generate response distance matrix 
resp_dists = dist(temp$Chloride)
# run Mantel test
mantel.rtest(site_dists, resp_dists, nrepet = 9999)
# 'observation' is the correlation between the distance matrices
# p value suggests that they are NOT correlated
# So, based on this test, there is no detectable correlation
## Map
proj = CRS("+proj=longlat +datum=WGS84")
temp_spatial  <- SpatialPointsDataFrame(coords= cbind(temp$lon, temp$lat),
                                        data = as.data.frame(cbind(temp$Site.Name, temp$Chloride)),
                                        proj4string = proj)
plot(temp_spatial)

# ect.......... for other parameters of interest and for a few other time points, depending on how your data is structured ...........

#
#### check correlation between variables ####

# first, returning to the dataset of just 3 sites and more than 100 obs per parameter (dat_r), reformat data to make it wider, such that parameters get their own columns. 

dat_r_long = dat_r %>% 
  select(c(SiteName, Datecombined, Samplingtype, Bromide))%>%
  group_by(SiteName, Datecombined, Samplingtype) %>%
  summarise(Bromide = mean(Bromide, na.rm = T)) %>%
  pivot_wider(names_from = Samplingtype, 
              values_from = Bromide,
              values_fill = list(Bromide = NA))

# reduce data down to one site - Savannah
temp = dat_r_long %>% filter(SiteName=="Alameda") 
# plot correlations (of data columns only)
pairs.panels(temp[,3:24], scale=T)
pairs.panels(temp[,3:24], scale=F)
# make table of correlations (I am rounding and replacing low values with text so that it is easier to see results)
tab = round(as.data.frame(cor(cov(temp[,3:24], use="na.or.complete"))), 2)
tab[abs(tab)<0.4] = "no_corr"

# reduce data down to one site - Calabacillas
temp = dat_r_long %>% filter(SiteName=="Harrison")
# plot correlations (of data columns only)
pairs.panels(temp[,3:24], scale=T)
pairs.panels(temp[,3:24], scale=F)
# make table of correlations (I am rounding and replacing low values with text so that it is easier to see results)
tab = round(as.data.frame(cor(cov(temp[,3:24], use="na.or.complete"))), 2)
tab[abs(tab)<0.4] = "no_corr"

#for linear model and 


datfiltered <-
  dat %>%
  filter(SiteName %in% c("Badger","Belen","Bobcat","Harrison","Savannah")) %>%
  mutate(SiteName = factor(SiteName)) %>%
  filter(Chloride > 0) %>%
  filter(Sulfate>0 ) %>%
  filter(PhosphateP >0) %>% 
  filter(NH4N >0)



library(ggplot2)
p <- ggplot(datfiltered %>% filter(Sulfate < 150), aes(x = Date, y = Sulfate))
p <- p + theme_bw()
#p <- p + geom_point(alpha = 1/10)
p <- p + geom_jitter(width = 0.25, height = 5, alpha = 1/10)
p <- p + geom_smooth()
p <- p + stat_smooth(method = lm, colour = "red")
p <- p + labs(
  title = "Sulfate over time"
  #, x = "Date"
  , y = "Sulfate_mg/l "
  , caption = "Sulfate over years"
)

print(p)


str(datfiltered)

dat = read.csv("anionsf.csv")

datfiltered %>%
  select(SiteName , Sulfate) %>%
  drop_na() %>%
  ggplot(aes(x = SiteName, y = Sulfate)) + 
  geom_boxplot()


# create the linear model
Siteanions <- lm(Sulfate ~ SiteName + Samplingtype, data = datfiltered, na.action=na.omit)
# check assumptions
plot(Siteanions)


