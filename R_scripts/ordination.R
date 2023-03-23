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

data(meaudret)

d = meaudret$env
d = d %>% 
  rownames_to_column("siteID") %>% 
  separate(siteID, c("site_group","site_num"), remove = F)

#### 1. check assumption ####

pairs.panels(d[,4:12], scale=TRUE)

range(d$Bdo5)
range(d$Flow)
range(d$Ammo)

#### 2. compute PCA ####

d2 = d[,4:12]
rownames(d2) = d$siteID

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

pca.p.1.4 =
  fviz_pca_biplot(pca, 
                  axes = c(1, 4), # specify axes
                  repel = TRUE,
                  col.var = "#2E9FDF", # Variables color
                  col.ind = "#696969"  # Individuals color
  )

# plot of sites + 95% CI groupings #

fviz_pca_ind(pca, habillage=d$site_group, 
             axes = c(1, 2),
             label="none", 
             addEllipses=TRUE, ellipse.level=0.95)

fviz_pca_ind(pca, habillage=d$site_group, 
             axes = c(1, 3),
             label="none", 
             addEllipses=TRUE, ellipse.level=0.95)

fviz_pca_ind(pca, habillage=d$site_group, 
             axes = c(2, 3),
             label="none", 
             addEllipses=TRUE, ellipse.level=0.95)

fviz_pca_ind(pca, habillage=d$site_group, 
             axes = c(1, 4),
             label="none", 
             addEllipses=TRUE, ellipse.level=0.95)


#### 5. check weight of outliers ####

# specify outliers #
outliers = c("au_2")

# remove outliers #
d3 = d2[!(row.names(d2) %in% outliers),]

# re-compute pca #
pca2 = prcomp(d3, scale = T)

# re-plot pca #
fviz_pca_biplot(pca2, 
                axes = c(1, 2), # specify axes
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


#### 8. extract scores ####

pca.scrs = as.data.frame(pca$x[,1:3])
pca.scrs
