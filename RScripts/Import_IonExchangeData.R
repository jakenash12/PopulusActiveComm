library(car)
library(ggplot2)
#library(dplyr)
#library(tibble)
library(cowplot)
#library(grid)
#library(gridExtra) 
library(lme4)
#library(tidyr)
library(car)

#imports ion exchange data
PRS_Results = read.delim("./InputFiles/IonExchangeData_Formatted.txt") %>%
  mutate(PlotDate=paste(.$Plot,.$Time,sep="_"),
         Time=factor(Time, levels=c("June", "August", "October")),
         SiteType=dplyr::recode(SiteType, H="High Elevation",
                                M="Mid Elevation",
                                R="Riparian")) %>%
  mutate(SiteType=factor(SiteType, levels=c("Riparian", "Mid Elevation", "High Elevation")))

#recodes site type variable for plotting
PRS_Results_ed=
  PRS_Results %>%
  mutate(SiteType=dplyr::recode(SiteType, "Riparian"="Rip.", 
                                "Mid Elevation"="Mid El.",
                                "High Elevation"="High El."))

#generates plots by site type and season for each ion flux
NO3_plot=
  ggplot(PRS_Results_ed, aes(x=SiteType, y=NO3, colour=Time)) +
  geom_boxplot() +
  theme_test() +
  theme(legend.position="none", text=element_text(color="black"), axis.text=element_text(color="black"), axis.title.x=element_blank())

Ca_plot=
  ggplot(PRS_Results_ed, aes(x=SiteType, y=Ca, colour=Time)) +
  geom_boxplot() +
  theme_test() +
  theme(legend.position="none", text=element_text(color="black"), axis.text=element_text(color="black"), axis.title.x=element_blank())

Mg_plot=
  ggplot(PRS_Results_ed, aes(x=SiteType, y=Mg, colour=Time)) +
  geom_boxplot() +
  theme_test() +
  theme(legend.position="none", text=element_text(color="black"), axis.text=element_text(color="black"), axis.title.x=element_blank())

K_plot=
  ggplot(PRS_Results_ed, aes(x=SiteType, y=K, colour=Time)) +
  geom_boxplot() +
  theme_test() +
  theme(legend.position="none", text=element_text(color="black"), axis.text=element_text(color="black"), axis.title.x=element_blank())

P_plot=
  ggplot(PRS_Results_ed, aes(x=SiteType, y=P, colour=Time)) +
  geom_boxplot() +
  theme_test() +
  theme(legend.position="none", text=element_text(color="black"), axis.text=element_text(color="black"), axis.title.x=element_blank())

Fe_plot=
  ggplot(PRS_Results_ed, aes(x=SiteType, y=Fe, colour=Time)) +
  geom_boxplot() +
  theme_test() +
  theme(legend.position="none", text=element_text(color="black"), axis.text=element_text(color="black"), axis.title.x=element_blank())

Mn_plot=
  ggplot(PRS_Results_ed, aes(x=SiteType, y=Mn, colour=Time)) +
  geom_boxplot() +
  theme_test() +
  theme(legend.position="none", text=element_text(color="black"), axis.text=element_text(color="black"), axis.title.x=element_blank())

Zn_plot=
  ggplot(PRS_Results_ed, aes(x=SiteType, y=Zn, colour=Time)) +
  geom_boxplot() +
  theme_test() +
  theme(legend.position="none", text=element_text(color="black"), axis.text=element_text(color="black"), axis.title.x=element_blank())

B_plot=
  ggplot(PRS_Results_ed, aes(x=SiteType, y=B, colour=Time)) +
  geom_boxplot() +
  theme_test() +
  theme(legend.position="none", text=element_text(color="black"), axis.text=element_text(color="black"), axis.title.x=element_blank())

S_plot=
  ggplot(PRS_Results_ed, aes(x=SiteType, y=S, colour=Time)) +
  geom_boxplot() +
  theme_test() +
  theme(legend.position="none", text=element_text(color="black"), axis.text=element_text(color="black"), axis.title.x=element_blank())

Al_plot=
  ggplot(PRS_Results_ed, aes(x=SiteType, y=Al, colour=Time)) +
  geom_boxplot() +
  theme_test() +
  theme(legend.position="none",text=element_text(color="black"), axis.text=element_text(color="black"), axis.title.x=element_blank())

#function to plot legend
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 

#makes legend
PRS_legend=
  ggplot(mutate(PRS_Results_ed, Season=Time), aes(x=SiteType, y=Al, colour=Season)) +
  geom_boxplot() +
  theme_test() 

#formats legend
PRS_legend = g_legend(PRS_legend)

#makes multipanel plot of all ion fluxes
PRS_plot=plot_grid(PRS_legend,NO3_plot, Ca_plot, Mg_plot, K_plot,
          P_plot, Fe_plot, Mn_plot, S_plot,
          B_plot, Zn_plot, Al_plot, 
          align = "v", axis = 'l')

#conducts mixed modelling of all PRS variables
#of interest and formats into tidy df
###############################################
#make list of response vars
response_vars=c("NO3", "Ca", "Mg", "K",
                "P", "Fe", "Mn", "S",
                "B", "Zn", "Al")

# Initialize an empty list to store results
results_list <- list()

# Loop through each response variable
for (response in response_vars) {
  # Build the formula for the mixed model
  formula <- as.formula(paste(response, "~ SiteType * Time + (1|Site)"))
  
  # Fit the model
  model <- tryCatch({
    lmer(formula, data = PRS_Results_ed)
  }, error = function(e) {
    message("Error in fitting model for ", response, ": ", e$message)
    return(NULL)  # If the model fails, return NULL
  })
  
  # If model fitting was successful
  if (!is.null(model)) {
    # Use the Anova function from the car package for type II or type III ANOVA
    anova_results <- Anova(model, type = 2)  # You can change to type = 3 if needed
    
    # Extract the p-values for SiteType, Time, and Interaction
    p_values <- anova_results["SiteType", "Pr(>Chisq)"]
    p_values_time <- anova_results["Time", "Pr(>Chisq)"]
    p_values_interaction <- anova_results["SiteType:Time", "Pr(>Chisq)"]
    
    # Create a tidy result dataframe
    result <- data.frame(
      Response = response,
      SiteType_p_value = p_values,
      Time_p_value = p_values_time,
      Interaction_p_value = p_values_interaction
    )
    
    # Append to results list
    results_list[[response]] <- result
  }
}

# Combine all results into a single dataframe
results_df <- bind_rows(results_list)

# Optional: reorder columns for better readability
results_df <- results_df %>%
  select(Response, SiteType_p_value, Time_p_value, Interaction_p_value)

# View the final results
print(results_df)

#writes model results to tsv file
write_tsv(results_df, "./OutputFiles/PRS_lmer_results.tsv")

