#This script tests for the scale dependence of relationships between environmental factors
#and the sequence abundance of microbial lineages.
#This is done by comparing models that have a random intercept with those that have a fixed
#intercept. If the random intercept model is significantly better than the fixed intercept
#model, then it would be evidence of scale dependence.

#This script requires you to have run 16S_Analysis_SILVA.R,
#ITS_Analysis.R, RandomForests_Taxa.R and SoilVarDimensionalReduction.R 
#and have the resulting data objects in memory

library(lme4)
library(car)
library(dplyr)
library(tidyr)
library(purrr)
library(broom.mixed)
library(ggplot2)

####################ITS2 scale dependence#######################################


# Define variables
predictors_ITS <- names(select(AUE2021_ITS_rarefied_class_t, NO3Flux:NConc))

# Initialize results list
results_ITS <- list()

# Loop over response and predictor combinations
for (resp in DomClass_ITS) {
  for (pred in predictors_ITS) {
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
        results_ITS[[length(results_ITS) + 1]] <- tibble(
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
results_ITS_df <- bind_rows(results_ITS)

#filter for relationships where either the fixed or random slope model was significant
#and there was a significant difference between the two
results_ITS_df %>%
  filter((p_random<0.05 | p_fixed <0.05),p_compare<0.05) 


####################16S scale dependence#######################################
  
  # Define variables
  predictors_16S <- names(select(AUE2021_16S_rarefied_phylum_t, NO3Flux:NConc))
  
  # Initialize results list
  results_16S <- list()
  
  # Loop over response and predictor combinations
  for (resp in DomPhyla_16S) {
    for (pred in predictors_16S) {
      df <- AUE2021_16S_rarefied_phylum_t
      
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
          results_16S[[length(results_16S) + 1]] <- tibble(
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
results_16S_df <- bind_rows(results_16S)
  
#filter for relationships where either the fixed or random slope model was significant
#and there was a significant difference between the two
results_16S_df %>%
  filter((p_random<0.05 | p_fixed <0.05),p_compare<0.05) %>% View

###################################################

