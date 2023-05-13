#### read me ####


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
library(psych)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(dplyr)

#### load and tidy data ####
setwd("C:/Users/marag.DESKTOP-NEKBNJE/OneDrive/Desktop/Stakeholder data analysis/test1/Data")
dat = read.csv("anionsp.csv")
dat <- read.csv("anionsp.csv", na.strings = c("."))
str(dat)
#dat<- na.omit(dat)



#### load and tidy data ####


#dat$datetime_NM = as.POSIXct(dat$Sample_DateTime, format="%m/%d/%y %H:%M", tz="US/Mountain")

# convert characters that should be factors (categories) to factors
dat$site_group = as.factor(dat$site_group)
dat$siteID = as.factor(dat$siteID)
#dat$Is_Nondetect = as.factor(dat$Is_Nondetect)

# convert water quality data to numeric
dat$Bromide = as.numeric(dat$Bromide)
dat$Fluoride = as.numeric(dat$Fluoride)
dat$NH4N  = as.numeric(dat$NH4N )
dat$Sulfate = as.numeric(dat$Sulfate)
dat$Chloride = as.numeric(dat$Chloride)
dat$PhosphateP = as.numeric(dat$PhosphateP)
dat$NitrateN  = as.numeric(dat$NitrateN)
dat$site_num   = as.numeric(dat$site_num)


#exclude variables and work on non zero inflated variables (Will have ~ 9 variables not zero inflated)

dat <- select(dat, -c("Fluoride","NitrateN","Bromide","Nitrite"))
str(dat)
# reduce to sampling types and sites for the purposes of most of this tutorial
#dat = dat[dat$Samplingtype=="Center"| dat$Samplingtype =="River"| dat$Samplingtype =="Ditch",]
dat = dat[dat$site_group=="Al"| dat$site_group =="Sv"| dat$site_group =="Bc"| dat$site_group =="Bp"| dat$site_group =="Bg",]



#### 1. check assumption ####

dat <- na.omit(dat)

pairs.panels(dat[,4:7], scale=TRUE)

range(dat$Sulfate)
range(dat$Chloride)
range(dat$PhosphateP)

#### 2. compute PCA ####

d2 = dat[,4:7]
rownames(d2) = dat$siteID

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

pca.p.2.3 =
  fviz_pca_biplot(pca, 
                  axes = c(2, 3), # specify axes
                  repel = TRUE,
                  col.var = "#2E9FDF", # Variables color
                  col.ind = "#696969"  # Individuals color
  )

pca.p.1.4 =
  fviz_pca_biplot(pca, 
                  axes = c(1, 4), # specify axes
                  repel = TRUE,
                  col.var = "#2E9FDF", # Variables color
                  col.ind = "#696969"  # Individuals color
  )

# plot of sites + 95% CI groupings #

fviz_pca_ind(pca, habillage=dat$site_group, 
             axes = c(1, 2),
             label="none", 
             addEllipses=TRUE, ellipse.level=0.95)

fviz_pca_ind(pca, habillage=dat$site_group, 
             axes = c(1, 3),
             label="none", 
             addEllipses=TRUE, ellipse.level=0.95)

fviz_pca_ind(pca, habillage=dat$site_group, 
             axes = c(2, 3),
             label="none", 
             addEllipses=TRUE, ellipse.level=0.95)

fviz_pca_ind(pca, habillage=dat$site_group, 
             axes = c(1, 4),
             label="none", 
             addEllipses=TRUE, ellipse.level=0.95)
fviz_pca_ind(pca, habillage=dat$site_group, 
             axes = c(2, 4),
             label="none", 
             addEllipses=TRUE, ellipse.level=0.95)
fviz_pca_ind(pca, habillage=dat$site_group, 
             axes = c(4, 3),
             label="none", 
             addEllipses=TRUE, ellipse.level=0.95)

# fviz_pca_ind(pca, habillage=dat$siteID, 
# axes = c(1, 4),
# label="none", 
# addEllipses=TRUE, ellipse.level=0.95)

# fviz_pca_ind(pca, habillage=dat$SiteName,
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

fviz_pca_biplot(pca2, 
                axes = c(1, 4), # specify axes
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

# axis 1: chloride and sulfate explain a lot of variability in this dimension 1.

# axis 2: phosphate explain a lot of variability in this dimension 2.

# axis 3: Ammonia nitrogen combination variability in this dimension 3.



#### 8. extract scores #### at this point we have 3 axes we consider latent variables(orthogonal) that can explain the variability for the variables, we can use them as predictor varibles in a model

pca.scrs = as.data.frame(pca$x[,1:4])
pca.scrs

