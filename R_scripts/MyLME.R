# packages to load
library(tidyverse) 
library(palmerpenguins) # for first example
library(lme4) # for creating mixed models
library(car) # for Anova(), vif()
library(MuMIn) # for AICc
library(emmeans) # for emmeans, emtrends, all the post hoc tests and plotting
library(dplyr)

## here is a cheatsheet for emmeans. Super helpful!!
# https://timmastny.rbind.io/blog/tests-pairwise-categorical-mean-emmeans-contrast/



#### load and tidy data ####

setwd("C:/Users/marag.DESKTOP-NEKBNJE/OneDrive/Desktop/Stakeholder data analysis/test1/Data")
dat = read.csv("anionsf.csv")
dat <- read.csv("anionsf.csv", na.strings = c("."))
#dat<- na.omit(dat)

str(dat)
# reduce to 3 sites for the purposes of most of this tutorial
dat = dat[dat$Samplingtype =="River"| dat$Samplingtype =="Ditch"| dat$Samplingtype =="Center",]

#dat = dat[dat$SiteName =="Bobcat"| dat$SiteName =="Badger"| dat$SiteName =="Harrison"| dat$SiteName =="Savannah"| dat$SiteName =="Belen",]
dat = dat[dat$SiteName =="Badger"| dat$SiteName =="Bobcat"| dat$SiteName =="Savannah"| dat$SiteName =="Harrison"| dat$SiteName =="Belen",]

#arrange from north to south
dat %>% arrange("Badger","Bobcat","Savannah","Harrison","Belen")

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


#plot after reordering site and sampling type :sulfate
dat %>%
  mutate(SiteName = fct_relevel(SiteName, 
                                "Badger", "Bobcat", "Savannah", 
                                "Harrison", "Belen")) %>%
  mutate(Samplingtype = fct_relevel(Samplingtype, 
                                "River", "Ditch", "Center")) %>%
  ggplot( aes(x=SiteName, y=Sulfate, color = Samplingtype)) +
  geom_boxplot()


#plot after reordering site and sampling type :Chloride
dat %>%
  mutate(SiteName = fct_relevel(SiteName, 
                                "Badger", "Bobcat", "Savannah", 
                                "Harrison", "Belen")) %>%
  mutate(Samplingtype = fct_relevel(Samplingtype, 
                                    "River", "Ditch", "Center")) %>%
  ggplot( aes(x=SiteName, y=Chloride, color = Samplingtype)) +
  geom_boxplot()
dat %>%
  ggplot(aes(x = SiteName, y = Chloride, color = Samplingtype)) + 
  geom_boxplot()


#plot after reordering site and sampling type :PhosphateP
dat %>%
  mutate(SiteName = fct_relevel(SiteName, 
                                "Badger", "Bobcat", "Savannah", 
                                "Harrison", "Belen")) %>%
  mutate(Samplingtype = fct_relevel(Samplingtype, 
                                    "River", "Ditch", "Center")) %>%
  ggplot( aes(x=SiteName, y=PhosphateP, color = Samplingtype)) +
  geom_boxplot()
dat %>%
  ggplot(aes(x = SiteName, y = Chloride, color = Samplingtype)) + 
  geom_boxplot()


#plot after reordering site and sampling type :NH4N
dat %>%
  mutate(SiteName = fct_relevel(SiteName, 
                                "Badger", "Bobcat", "Savannah", 
                                "Harrison", "Belen")) %>%
  mutate(Samplingtype = fct_relevel(Samplingtype, 
                                    "River", "Ditch", "Center")) %>%
  ggplot( aes(x=SiteName, y=PhosphateP, color = Samplingtype)) +
  geom_boxplot()
dat %>%
  ggplot(aes(x = SiteName, y = NH4N, color = Samplingtype)) + 
  geom_boxplot()


# create the linear model sulfate
Sul_m1 <- lm(Sulfate ~ SiteName + Samplingtype, data = dat, na.action=na.omit)
# check assumptions
plot(Sul_m1)
# run type 3 ANOVA
Anova(Sul_m1, type = 3)



# create the linear model phsophate
Phos_m1 <- lm(PhosphateP ~ SiteName + Samplingtype, data = dat, na.action=na.omit)
# check assumptions
plot(Phos_m1)
# run type 3 ANOVA

Anova(Phos_m1, type = 3)



# create the linear model chloride
Chlo_m1 <- lm(Chloride ~ SiteName + Samplingtype, data = dat, na.action=na.omit)
# check assumptions
plot(Chlo_m1)
# run type 3 ANOVA
Anova(Chlo_m1, type = 3)


# create the linear model NH4N
NH4n_m1 <- lm(NH4N ~ SiteName + Samplingtype, data = dat, na.action=na.omit)
# check assumptions
plot(NH4n_m1)
# run type 3 ANOVA
Anova(NH4n_m1, type = 3)


# create the linear model river predicted either by  well and ditch
RDG_m1 <- lm(Samplingtype$River ~ Samplingtype$Center + Samplingtype$Ditch, data = dat, na.action=na.omit)
# check assumptions
plot(RDG_m1)
# run type 3 ANOVA
Anova(RDG_m1, type = 3)



# run type 3 ANOVA
Anova(Sul_m1, type = 3)

#Model for sulfate

Sul_m1 <- lm(Sulfate ~ SiteName * Samplingtype, data = dat, na.action=na.omit)
Sul_additive <- lm(Sulfate ~ SiteName + Samplingtype, data = dat, na.action=na.omit)
Sul_m1_site <- lm(Sulfate ~ SiteName, data = dat, na.action=na.omit)
Sul_m1_Samplingtype <- lm(Sulfate ~ Samplingtype, data = dat, na.action=na.omit)
Sul_m1_null <- lm(Sulfate ~  1, data = dat, na.action=na.omit)
Sul_m1year <- lm(Sulfate ~ SiteName * Year, data = dat, na.action=na.omit)


AICc(Sul_m1,Sul_additive,Sul_m1_site,Sul_m1_null,Sul_m1_Samplingtype,Sul_m1year)

anova(Sul_m1, Sul_additive,Sul_m1_site)


#Model for phosphate

Pho_m1 <- lm(PhosphateP ~ SiteName * Samplingtype, data = dat, na.action=na.omit)
Pho_additive <- lm(PhosphateP ~ SiteName + Samplingtype, data = dat, na.action=na.omit)
Pho_m1_site <- lm(PhosphateP ~ SiteName, data = dat, na.action=na.omit)
Pho_m1_Samplingtype <- lm(PhosphateP ~ Samplingtype, data = dat, na.action=na.omit)
Pho_m1_null <- lm(PhosphateP ~  1, data = dat, na.action=na.omit)
Pho_m1year <- lm(PhosphateP ~ SiteName * Year, data = dat, na.action=na.omit)


AICc(Pho_m1,Pho_additive,Pho_m1_site,Pho_m1_null,Pho_m1_Samplingtype,Pho_m1year)

anova(Pho_m1, Pho_additive,Pho_m1_site)


#Model for chloride

Chl_m1 <- lm(Chloride ~ SiteName * Samplingtype, data = dat, na.action=na.omit)
Chl_additive <- lm(Chloride ~ SiteName + Samplingtype, data = dat, na.action=na.omit)
Chl_m1_site <- lm(Chloride ~ SiteName, data = dat, na.action=na.omit)
Chl_m1_Samplingtype <- lm(Sulfate ~ Samplingtype, data = dat, na.action=na.omit)
Chl_m1_null <- lm(Chloride ~  1, data = dat, na.action=na.omit)
Chl_m1year <- lm(Chloride ~ SiteName * Year, data = dat, na.action=na.omit)


AICc(Chl_m1,Chl_additive,Chl_m1_site,Chl_m1_null,Chl_m1_Samplingtype,Chl_m1year)

anova(Chl_m1, Chl_additive,Sul_m1_site,)


#post-hoc test
dat %>%
  select(Samplingtype, Sulfate, SiteName) %>%
  drop_na() %>%
  ggplot(aes(x = Sulfate, y = Sulfate)) + 
  geom_boxplot() + facet_grid(~Samplingtype)


# tukey test comparing species for females and for males
emmeans(Sul_m1, pairwise ~ SiteName | Samplingtype)

# tukey test comparing males vs famales for each species
emmeans(Sul_m1, pairwise ~ Samplingtype | SiteName)

post_hoc_dat <- emmeans(Sul_m1, pairwise ~ SiteName | Samplingtype)

post_hoc_dat$emmeans %>%
  as.data.frame() %>%
  ggplot(aes(x = SiteName, y = emmean, color = Samplingtype)) + 
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), size = 2) + 
  ylab("Sulfate_mg") + theme_bw()





##### plotting my data by anions over years #####
dat %>%
  ggplot(aes(x = Year, y = Sulfate, color = Samplingtype)) + geom_point() + facet_wrap(~SiteName)

#sulfate
dat %>%
  ggplot(aes(x = Year, y = Sulfate, color = SiteName)) + geom_point() + facet_wrap(~Samplingtype)

#phosphate
dat %>%
  ggplot(aes(x = Year, y = PhosphateP, color = Samplingtype)) + geom_smooth() + facet_wrap(~SiteName)

#chlorides
dat %>%
  ggplot(aes(x = Year, y = Chloride, color = Samplingtype)) + geom_smooth() + facet_wrap(~SiteName)


#chlorides
dat %>%
  ggplot(aes(x = Year, y = NH4N, color = Samplingtype)) + geom_smooth() + facet_wrap(~SiteName)




