#This script performs a followup to visualize and model the presence of Glomeromycetes in response
#to nitrate flux that was identified through random forest modelling.
#
#Based on the distribution of the data, Glomeromycetes sequence abundance is recoded as
#presence/absence using a sequence abundance threshold.
#A logistic regression model is used to moodel Glomeromycetes presence/absence as a
#function of NO3 flux and an accompanying regression is made

library(tidyverse)
library(car)

#50 sequences is used as the cutoff for determing Glomeromycetes presence/absence
AMF_presence_df <- AUE2021_ITS_rarefied_class_t_RNA %>%
  mutate(Glom_present = ifelse(Glomeromycetes > 50, 1, 0))

#creates a logistic regression plot
ggplot(filter(AMF_presence_df, NO3Flux<200), 
       aes(x = NO3Flux, y = Glom_present)) +
  geom_jitter(height = 0.04, 
              width = 0, 
              size = 5, 
              alpha = 0.2,
              color="darkblue") +  # adds jittered binary points
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = TRUE,
              color="black") +  # logistic curve
  labs(
    y = "Probability of Glomeromycetes Presence",
  ) +
  theme_test() +
  theme(aspect.ratio=1) + 
  scale_x_continuous(breaks = c(0, 65, 130))

#runs logistic regression
model <- glm(Glom_present ~ NO3Flux, 
             data = filter(AMF_presence_df, NO3Flux < 200), 
             family = "binomial")
Anova(model)
