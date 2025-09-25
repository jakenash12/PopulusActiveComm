
#This script requires you to have run 16S_Analysis_SILVA.R,
#ITS_Analysis.R, SpatialDistanceMatrix.R and have the resulting data objects in memory

library(lme4)
library(car)
library(dplyr)
library(tidyr)
library(purrr)
library(broom.mixed)
library(ggplot2)


mod_fixed=
  AUE2021_ITS_rarefied_class_sensor_t_RNA %>%
  lmer(Dothideomycetes~meansoiltemp*Time + (1|Site),.)

mod_random=
  AUE2021_ITS_rarefied_class_sensor_t_RNA %>%
  lmer(Dothideomycetes~meansoiltemp*Time + (meansoiltemp|Site),.)

mod_fixed=
  AUE2021_16S_rarefied_phylum_t_DNA %>%
  lmer(Crenarchaeota   ~CECConc_CaConc_TOCConc*Time + (1|Site),.)

mod_random=
  AUE2021_16S_rarefied_phylum_t_DNA %>%
  lmer(Crenarchaeota   ~CECConc_CaConc_TOCConc*Time + (Time*CECConc_CaConc_TOCConc|Site),.)

anova(mod_fixed, mod_random)

Anova(mod_fixed)

ggplot(AUE2021_16S_rarefied_phylum_t_DNA,
       aes(x=CECConc_CaConc_TOCConc, y =Crenarchaeota , color=Site)) +
  geom_point(size=4) +
  geom_smooth(method="lm", se=FALSE) +
  facet_grid(.~Time)

# Define variables
responses <- DomClass_ITS
predictors <- names(select(AUE2021_ITS_rarefied_class_t, NO3Flux:NConc))

# Initialize results list
results <- list()

# Loop over response and predictor combinations
for (resp in responses) {
  for (pred in predictors) {
    df <- AUE2021_ITS_rarefied_class_t
    
    # Construct formulae
    formula_fixed <- as.formula(paste0(resp, " ~ ", pred, "*Time + source + (1 | Site)"))
    formula_random <- as.formula(paste0(resp, " ~ ", pred, "*Time + source + (", pred, " | Site)"))
    
    # Fit models
    model_fixed <- tryCatch(lmer(formula_fixed, data = df), error = function(e) NULL)
    model_random <- tryCatch(lmer(formula_random, data = df), error = function(e) NULL)
    
    if (!is.null(model_fixed) && !is.null(model_random)) {
      # Type II ANOVA p-values for each model
      anova_fixed <- tryCatch(Anova(model_fixed, type = 2), error = function(e) NULL)
      anova_random <- tryCatch(Anova(model_random, type = 2), error = function(e) NULL)
      
      # Likelihood ratio test
      model_comp <- tryCatch(anova(model_fixed, model_random), error = function(e) NULL)
      
      if (!is.null(anova_fixed) && !is.null(anova_random) && !is.null(model_comp)) {
        results[[length(results) + 1]] <- tibble(
          response = resp,
          predictor = pred,
          p_fixed = anova_fixed[match(pred, rownames(anova_fixed)), "Pr(>Chisq)"],
          p_random = anova_random[match(pred, rownames(anova_random)), "Pr(>Chisq)"],
          p_compare = model_comp$`Pr(>Chisq)`[2]
        )
      }
    }
  }
}

# Combine into long-format dataframe
results_df <- bind_rows(results)
nrow(results_df)
results_df %>%
  filter((p_random<0.05 | p_fixed <0.05),p_compare<0.05) %>% View
  filter(p_compare<0.05, p_random<0.05) %>% nrow
###################################################

mod_Fixedslope=
  rf_df_16S_sensor_DNA %>%
  left_join(samples_df2_16S) %>%
  lmer(MDS1~meansoiltemp + (1 | Site),.)

mod_Randomslope=
  rf_df_16S_sensor_DNA %>%
  left_join(samples_df2_16S) %>%
  lmer(MDS1~meansoiltemp + (meansoiltemp | Site),.)

anova(mod_Fixedslope, mod_Randomslope)

AUE2021_ITS_rarefied_class_sensor_t_DNA %>%
  filter(Time=="October") %>%
  lmer(Eurotiomycetes~meansoiltemp + (1|Site),.) %>%
  Anova

mod_Fixedslope=
  AUE2021_ITS_rarefied_class_sensor_t_RNA %>%
  lmer(Dothideomycetes~meansoiltemp*Time + (1 | Site),.)

lmer(response~predictor*Time + source+ (predictor | Site),
     AUE2021_ITS_rarefied_class_t)

DomClass_ITS
NO3Flux:NConc

anova(mod_Fixedslope, mod_Randomslope)

AUE2021_ITS_rarefied_class_sensor_t_RNA %>%
  ggplot(aes(x=meansoiltemp, y =Dothideomycetes, group=Site)) +
  geom_point(size=3, aes(color=Site)) +
  geom_smooth(se=FALSE, method="lm")


mod_Fixedslope=
  rf_df_ITS_sensor_DNA %>%
  left_join(samples_df2) %>%
  lmer(MDS1~meansoiltemp + (1 | Site),.)

mod_Randomslope=
  rf_df_ITS_sensor_DNA %>%
  left_join(samples_df2) %>%
  lmer(MDS1~meansoiltemp + (meansoiltemp | Site),.)

anova(mod_Fixedslope, mod_Randomslope)


rf_df_16S_sensor_DNA %>%
  left_join(samples_df2_16S) %>%
  ggplot(aes(x=meansoiltemp,y=MDS1)) +
  geom_point(size=5) +
  facet_wrap(SiteType~Site, nrow=3, scales="free")


AUE2021_ITS_rarefied_class_sensor_t_DNA %>%
  filter(Time == "June") %>%
  lmer(Eurotiomycetes ~ meansoiltemp + (meansoiltemp | Site), .) %>%
  Anova()
