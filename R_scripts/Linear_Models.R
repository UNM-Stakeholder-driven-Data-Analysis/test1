# packages to load
library(tidyverse) 
library(palmerpenguins) # for first example
library(lme4) # for creating mixed models
library(car) # for Anova(), vif()
library(MuMIn) # for AICc
library(emmeans) # for emmeans, emtrends, all the post hoc tests and plotting


## here is a cheatsheet for emmeans. Super helpful!!
# https://timmastny.rbind.io/blog/tests-pairwise-categorical-mean-emmeans-contrast/


#### penguin examples with linear model (no random effects) ####

# going to compare body mass across the species #

## plot the data
penguins %>%
  ggplot(aes(x = species, y = body_mass_g)) + 
  geom_boxplot()

penguins %>%
  ggplot(aes(x = species, y = body_mass_g, color = sex)) + 
  geom_boxplot()

penguins %>%
  select(species, body_mass_g, sex) %>%
  drop_na() %>%
  ggplot(aes(x = species, y = body_mass_g, color = sex)) + 
  geom_boxplot()


# create the linear model
mass_species_m1 <- lm(body_mass_g ~ species * sex, data = penguins, na.action=na.omit)
# check assumptions
plot(mass_species_m1)


# run type 3 ANOVA
Anova(mass_species_m1, type = 3)

mass_species_m1 <- lm(body_mass_g ~ species * sex, data = penguins, na.action=na.omit)
mass_species_additive <- lm(body_mass_g ~ species + sex, data = penguins, na.action=na.omit)
mass_species_sex <- lm(body_mass_g ~ sex, data = penguins, na.action=na.omit)
mass_species_species <- lm(body_mass_g ~ species, data = penguins, na.action=na.omit)
mass_species_null <- lm(body_mass_g ~  1, data = penguins, na.action=na.omit)

AICc(mass_species_m1,mass_species_additive,mass_species_sex,  mass_species_species, mass_species_null)

anova(mass_species_m1, mass_species_additive)

#post-hoc test
penguins %>%
  select(species, body_mass_g, sex) %>%
  drop_na() %>%
  ggplot(aes(x = species, y = body_mass_g)) + 
  geom_boxplot() + facet_grid(~sex)


# tukey test comparing species for females and for males
emmeans(mass_species_m1, pairwise ~ species | sex)

# tukey test comparing males vs famales for each species
emmeans(mass_species_m1, pairwise ~ sex | species)

post_hoc_penguins <- emmeans(mass_species_m1, pairwise ~ species | sex)

post_hoc_penguins$emmeans %>%
  as.data.frame() %>%
  ggplot(aes(x = species, y = emmean, color = sex)) + 
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), size = 2) + 
  ylab("Body Mass") + theme_bw()




#### snowpack example #####
## using linear mixed models

# get snow sites and their elevation
snow_sites <- 
  read_csv("data/snowdepthsites.csv") %>%
  separate(site_name, into = c("Site", "Station"), sep = " \\(") %>% #clean site and statino so they match the dataset
  mutate(Station = str_replace(Station, pattern = "\\)", "")) %>%
  rename(elevation = `elev (ft)`) %>%  # clean up elevation column
  select(Site, Station, elevation) %>%
  mutate(Site = str_replace(Site, "PanchueLa", "Panchuela"))



snowpack <- 
  read_csv("data/historicmonthlysnowdepth.csv") %>%
  rename(Year = `Water Year`) %>%
  select(Site, Station, Year, Feb_2) %>%
  drop_na()




### data clean up ####
# dont worry about what I did here. It gets technical

snowpack_clean <- 
  snowpack %>%
  left_join(snow_sites, by = c("Site", "Station")) %>%
  mutate(Feb_2 = as.numeric(Feb_2))

# so some of the stations aren't present in the metadata. 
# I decided to just average the sites that lacked a matching station. 
# If you only merge by site, there are multiple stations at some of the sites and it will inflate your samples!
# this is just a quick workaround for that problem. The models do not change much even with it improperly inflated.


# make a vector of the sites with missing matching stations
missing_sites <- 
  snowpack_clean %>%
  filter(is.na(elevation)) %>%
  select(Site) %>%
  unique() %>% 
  pull()

# extract elevation for missing sites. average if multiple stations are present at that site
extra_sites <- 
  snow_sites %>%
  filter(Site %in% missing_sites) %>%
  group_by(Site) %>% 
  summarise(elevation = mean(elevation))

# join it again. This will make elevation now be elevation.x and elevation.y. This ok because there is a function
# called coalesce that will merge these two so the NAs are replaced with real values.
snowpack_clean_elevations <- 
  snowpack_clean %>%
  left_join(extra_sites, by = "Site")


# rewrite elevation with the newly merged elevation data
snowpack_clean <- 
  snowpack_clean %>%
  mutate(elevation = coalesce(snowpack_clean_elevations$elevation.x, snowpack_clean_elevations$elevation.y)) 



##### plotting snowpack data #####
snowpack_clean %>%
  ggplot(aes(x = Year, y = Feb_2, color = elevation)) + geom_point() + facet_wrap(~Site)





#### building and checking assumptions ####
# build model
lme4::lmer(Feb_2 ~ Year*elevation + (1|Site), data = snowpack_clean, REML = FALSE)

# rescale the predictors
snow_model <- lme4::lmer(Feb_2 ~ scale(Year)*scale(elevation) + (1|Site), data = snowpack_clean,  REML = FALSE)

# fails homoscedasticity predictor so we will need to transform it
plot(snow_model)

# no evidence of multicollinearity 
car::vif(snow_model)


# not the best fit but it should work
qqnorm(resid(snow_model))
qqline(resid(snow_model))

# transform model
snow_model <- lme4::lmer(sqrt(Feb_2) ~ scale(Year)*scale(elevation) + (1|Site), data = snowpack_clean, REML = FALSE)

# much better!!
plot(snow_model)

# not perfect but definitely an improvement
qqnorm(resid(snow_model))
qqline(resid(snow_model))

# building simpler models to compare to the full model

# additive model. No interaction
snow_model_add <- lme4::lmer(sqrt(Feb_2) ~ scale(Year) + scale(elevation) + (1|Site), data = snowpack_clean, REML = FALSE)

# removed year so only elevation
snow_model_elev <- lme4::lmer(sqrt(Feb_2) ~  scale(elevation) + (1|Site), data = snowpack_clean, REML = FALSE)

# removed  elevation and only kept year
snow_model_year <- lme4::lmer(sqrt(Feb_2) ~ scale(Year)  + (1|Site), data = snowpack_clean, REML = FALSE)

# removed all predictors to create a null model
snow_model_null <- lme4::lmer(sqrt(Feb_2 ) ~ 1 + (1|Site), data = snowpack_clean, REML = FALSE)





#### evaluating the snowpack model ####

# way #1. Just run type 3 ANOVA.
car::Anova(snow_model, type = 3)



# way #2: compare models using log-likelihood
anova(snow_model, snow_model_add) # need to simplify model. Throw out interaction
anova(snow_model_add,snow_model_elev)
anova(snow_model_add,snow_model_year)




# way #3: AICc approach
# the best model is snow_model_add
AICc(snow_model, snow_model_add, snow_model_elev,snow_model_year,  snow_model_null)

AICc(snow_model_add)

# clean it up and show delta AICc
AICc(snow_model, snow_model_add, snow_model_elev,snow_model_year,  snow_model_null) %>% 
  rownames_to_column("Model") %>% # makes row names become a column called "Model"
  mutate(delta_AICc = AICc - AICc(snow_model_add)) %>% # compute new column that substracts the AICc score for the best model
  arrange((delta_AICc)) # arrange so the best model is on top




# plot best fit lines

emmip(snow_model_add, Year ~ elevation,  
      cov.reduce = FALSE, # need this to plot lines
      transform = "response", # back-transforms response 
      at = list(Year = c(1940, 1960, 1980, 2000, 2020))) + # years to plot the lines
  ylab("Snow Water Equivalent (in) Start of February") +  
  ggtitle("Predicted Snow Pack in February") + 
  xlab("Elevation (ft)") + 
  theme_bw()





# you can add confidence invtervals as well
# i reduced it to only predict for 1945 and 2020
emmip(snow_model_add, Year ~ elevation,  
               cov.reduce = FALSE, transform = "response", 
               CI = TRUE, # add 95% confidence intervals
               at = list(Year = c(1945, 2020))) + 
  ylab("Snow Water Equivalent (in) Start of February") + 
  ggtitle("Predicted Snow Pack in February") + 
  xlab("Elevation (ft)") + 
  theme_bw()


               
               