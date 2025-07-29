#this script imports vegetation and soils data
#and performs statistical tests

library(ggplot2)
library(lme4)
library(emmeans)
library(car)
library(tidyr)
library(tidyverse)
library(multcomp)
library(readxl)

#read in vegetation survey data
Vegsurvey=read.delim("./InputFiles/VegSurveyFormattedData.txt")

#read in soil test results from UGA testing lab
UGA_Results=read.delim("./InputFiles/20210616_UGASoilTestResults.txt")


###############Tests for effect site type on ########
#####################soil properties#################
#models down with lmer using random effect of site

#creates list of soil variables to test
UGA_vars=c("LBC", "pH", "BaseSat", "CEC", "Ca", "Cd", "Fe", "K", "Mg", "Mn",
           "Na", "Ni", "P", "Zn", "N", "TOC")

# Placeholder for significance results
results_list <- list()

# Loop through each response variable
for (var in UGA_vars) {
  
  # Dynamically construct the formula
  formula <- as.formula(paste(var, "~ Habitat + (1|Site)"))
  
  # Fit the lmer model
  mod <- lmer(formula, data = UGA_Results)
  
  # Perform pairwise comparisons
  pairwise_emmeans <- emmeans(mod, ~ Habitat)
  pairwise_results <- pairs(pairwise_emmeans)
  
  # Generate compact letter display
  cld <- cld(pairwise_emmeans, Letters = letters)
  
  # Format the results into a dataframe for easy merging
  cld_df <- as.data.frame(cld)
  cld_df$ResponseVariable <- var
  
  # Store the results
  results_list[[var]] <- cld_df
}

# Combine all results into a single dataframe
final_results <- do.call(rbind, results_list)

# Print the final results
print(final_results)


#################################################################################
# Pivot data longer for ggplot2
UGA_Results_long <- UGA_Results %>%
  pivot_longer(
    cols = all_of(UGA_vars),
    names_to = "ResponseVariable",
    values_to = "Value"
  )


# Calculate additional space for each variable
y_limits <- UGA_Results_long %>%
  group_by(ResponseVariable) %>%
  summarize(y_min = min(Value, na.rm = TRUE),
            y_max = max(Value, na.rm = TRUE)) %>%
  mutate(y_upper_limit = y_max + (y_max - y_min) * 0.2,  # Add 20% extra space on top
         y_lower_limit = y_min - (y_max - y_min) * 0.1)  # Add 10% extra space at bottom

# Merge y-axis limits back to the dataset
UGA_Results_long <- UGA_Results_long %>%
  left_join(y_limits, by = "ResponseVariable")

# Create multipanel plot with custom y-axis limits
UGA_plot <- ggplot(UGA_Results_long, aes(x = Habitat, y = Value, color = Habitat)) +
  geom_boxplot(size = 1) +
  theme_test() +
  scale_color_brewer(palette = "Dark2") +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    panel.spacing = unit(1, "lines")
  ) +
  facet_wrap(~ ResponseVariable, scales = "free_y") +
  labs(
    color = "Habitat",
    x = NULL,
    y = "Value"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)))  # Add space at top (20%) and bottom (10%)

#saves plot of soil properties as pdf
pdf(file="./OutputFiles/UGA_plot.pdf", 
    width = 8.5, height = 10)
UGA_plot
dev.off()

#########################################################################
