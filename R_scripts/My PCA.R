#### read me ####

# The purpose of this script is to demonstrate calculating and exploring a PCA on a stream chemistry dataset

# Steps:
# 1. Check that data meets assumptions of PCAs
# 2. Compute PCA
# 3. Explore PCA
# 4. Plot PCA
# 5. Check that results meet assumptions - outliers
# 6. (actually 1!) Decide which axes to retain
# 7. Interpret axes
# 8. Extract scores for further analysis (if applicable)

#### libraries ####

library(tidyverse)
library(ade4)
library(psych)
library(FactoMineR)
library(factoextra)
library(corrplot)

#### load and tidy data ####
setwd("C:/Users/marag.DESKTOP-NEKBNJE/OneDrive/Desktop/Stakeholder data analysis/test1/Data")
dat = read.csv("anionsf.csv")
dat <- read.csv("anionsf.csv", na.strings = c("."," ","NA"))
str(dat)
#dat<- na.omit(dat)


# reduce to sampling types and sites for the purposes of most of this tutorial
dat = dat[dat$Samplingtype=="Center"| dat$Samplingtype =="River"| dat$Samplingtype =="Ditch"| dat$Samplingtype =="East"| dat$Samplingtype =="South"| dat$Samplingtype =="North"| dat$Samplingtype =="West"| dat$Samplingtype =="East well"| dat$Samplingtype =="South well"| dat$Samplingtype =="North well"| dat$Samplingtype =="West well",]
#other sites | dat$Samplingtype =="East well"| dat$Samplingtype =="South well"| dat$Samplingtype =="North well"| dat$Samplingtype =="West well"

#dat = dat[dat$SiteName =="Badger"| dat$SiteName =="Harrison"| dat$SiteName =="Savannah"| dat$SiteName =="Bobcat"| dat$SiteName =="Belen",]


# check data classes 
str(dat)

# format date/time

# format date/time: this code worked well to format date
dat$DateF = as.POSIXct(dat$Date, format="%m/%d/%Y %H:%M", tz="US/Mountain")
str(dat)
#dat$datetime_NM = as.POSIXct(dat$Sample_DateTime, format="%m/%d/%y %H:%M", tz="US/Mountain")

# convert characters that should be factors (categories) to factors
dat$SiteName = as.factor(dat$SiteName)
dat$Samplingtype = as.factor(dat$Samplingtype)
#dat$Is_Nondetect = as.factor(dat$Is_Nondetect)

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

#datF <- 
 # dat %>%  subset(dat, select = -c(Year,Month,Day,Sampling.type, Do,Fluoride,Nitrite,Bromide,NitrateN,PhosphateP,NH4.N,DateF, Timesampling, Timeofsitearrival, Date, Datecombined ))

# only keep select columns
#datS <-
 # dat %>%  subset(dat, select = -c("Airtemp", "Turbidity", "DoPercent", "DOTemp", "Conductivity", "condTemp" ,"pH1" , "Chloride", "Sulfate"))

#datS <- select(dat, c("SiteName","Airtemp", "DoPercent", "DOTemp", "Conductivity", "condTemp" ,"pH1" , "Chloride", "Sulfate"))
str(datS)


# trying to find the mean
dat_mean <- dat %>%
  group_by(Samplingtype) %>%
  summarise(across(everything(), mean, na.rm=TRUE))

str(dat_mean)

#select columns to keep
#all variables with anions

#dat_mean <- select(dat_mean, c("Samplingtype","Airtemp", "Turbidity", "DoPercent", "DOTemp", "Conductivity", "condTemp" ,"pH1", "Fluoride","NitrateN","PhosphateP","NH4N","Chloride", "Sulfate"))
# select necessary anions
dat_mean <- select(dat_mean, c("Samplingtype","Fluoride","NitrateN","PhosphateP","NH4N","Chloride", "Sulfate"))



#dat = dat %>% 
 # rownames_to_column("Site.Name") %>% 
 # separate(Site, c("Year","Day","Month","Sampling.type","Time","DO","Fluoride","Nitrite","Bromide","NitrateN", "PhosphateP", "NH4.N ", "DateF "), remove = F)
# turn columns and variables to numeric and factor


#exclude variables and work on non zero inflated variables (Will have ~ 9 variables not zero inflated)


#### 1. check assumption ####

dat_mean<- na.omit(dat_mean)

pairs.panels(dat_mean[,2:7], scale=TRUE)

range(dat_mean$Sulfate)
range(dat_mean$Chloride)
range(dat_mean$PhosphateP)

#### 2. compute PCA ####

d2 = dat_mean[,2:7]
rownames(d2) = dat_mean$Samplingtype

pca = prcomp(d2, scale = T)

summary(pca)

str(pca)

# The value of each observation on each PC axis
# also known as "scores"
pca$x 

#### 3. explore PCA ####

### view eigenvalues: ##
get_eigenvalue(pca) 
fviz_eig(pca, addlabels = TRUE, choice = c("eigenvalue")) # as scree plot with eigenvalues
fviz_eig(pca, addlabels = TRUE, choice = c("variance")) # as scree plot with % variance explained

#### 4. plot pca ####

# biplot of sites + variables #
pca.p.1.2 = 
  fviz_pca_biplot(pca, 
                  axes = c(1, 2), # specify axes
                  repel = TRUE,
                  col.var = "blue", # Variables color
                  col.ind = "#696969"  # Individuals color
  )

pca.p.1.3 =
  fviz_pca_biplot(pca, 
                  axes = c(1, 3), # specify axes
                  repel = TRUE,
                  col.var = "#2E9FDF", # Variables color
                  col.ind = "#696969"  # Individuals color
  )


# plot of sites + 95% CI groupings #

fviz_pca_ind(pca, habillage=dat_mean$Samplingtype, 
             axes = c(1, 2),
             label="none", 
             addEllipses=TRUE, ellipse.level=0.95)

fviz_pca_ind(pca, habillage=dat_mean$Samplingtype, 
             axes = c(1, 3),
             label="none", 
             addEllipses=TRUE, ellipse.level=0.95)

fviz_pca_ind(pca, habillage=dat_mean$Samplingtype, 
             axes = c(2, 3),
             label="none", 
             addEllipses=TRUE, ellipse.level=0.95)

# fviz_pca_ind(pca, habillage=dat_mean$Samplingtype, 
             # axes = c(1, 4),
             # label="none", 
             # addEllipses=TRUE, ellipse.level=0.95)

# fviz_pca_ind(pca, habillage=dat_mean$Samplingtype,
             # axes = c(1, 5),
            # label="none", 
            # addEllipses=TRUE, ellipse.level=0.95)



#### 5. check weight of outliers ####

# specify outliers #
#outliers = c("au_2")

# remove outliers #
d3 = d2 #[!(row.names(d2) %in% outliers),]

# re-compute pca #
pca2 = prcomp(d3, scale = T)

# re-plot pca #
fviz_pca_biplot(pca2, 
                axes = c(1, 2), # specify axes
                repel = TRUE,
                col.var = "blue", # Variables color
                col.ind = "#696969"  # Individuals color
)

fviz_pca_biplot(pca2, 
                axes = c(1, 3), # specify axes
                repel = TRUE,
                col.var = "blue", # Variables color
                col.ind = "#696969"  # Individuals color
)

fviz_pca_biplot(pca2, 
                axes = c(2, 3), # specify axes
                repel = TRUE,
                col.var = "blue", # Variables color
                col.ind = "#696969"  # Individuals color
)
# looks similar, so can use unaltered d #

#### 6. (actually 1!) decide which axes to retain ####

### YOU SHOULD CHOOSE WHICH RULE TO FOLLOW A-PRIORI!!! ###
# Possible a priori rules:
# a)  1 eigenvalue rule: values of 1 indicate that those axes account for more variance than any one original variable (standardized data only). Rule is to retain all axes with an eigenvalue of 1 or more
# b) % total variance rule: retain all axes for which % variance explained sums to > X% (e.g., 90%)
# c) % variance rule: retain all axes for which at least X% variance is explained (e.g., 10%)
# d) Broken stick rule: retain all axes before break in scree plot



fviz_eig(pca, addlabels = TRUE, choice = c("eigenvalue")) # as scree plot with eigenvalues
fviz_eig(pca, addlabels = TRUE, choice = c("variance")) # as scree plot with % variance explained

#### 7. interpret axes ####

# extract results
var = get_pca_var(pca)

# extract contribution of each variable to each PC (unit = %) #
contrib = (var$contrib)
contrib

# plot contribution #

corrplot(contrib, is.corr = F)

# axis 1:

# axis 2:

# axis 3:


#### 8. extract scores #### at this point we have 3 axes we consider latent variables(orthogonal) that can explain the variability for the variables, we can use them as predictor varibles in a model

pca.scrs = as.data.frame(pca$x[,1:3])
pca.scrs

