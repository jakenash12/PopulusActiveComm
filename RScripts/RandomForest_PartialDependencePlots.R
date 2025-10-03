#This script creates a partial dependence plot of mean soil tempearture
#from random forest models that will be used in a publication figure.
#Partial dependence plots are meant to display the effect of one variable
#while accounting for the covariance of other predictor variables.
#Thus, this plot shows the effect of soil temperature on microbial community
#composition (NMDS axis 1) while accounting for effects of all other soil variables
#that were used in the random forest models

#This script requires that you run the scripts 16S_Analysis_SILVA.R,
#ITS2_Analysis.R, PERMANOVA.R, SoilVarDimensionalReduction.R, and 
#RandomForests_PCoA_Axes.R and have the resulting data objects loaded in memory


library(randomForest)
library(pdp)
library(patchwork)

# Calculate partial dependence for each model
pdp_ITS_DNA <- partial(rf_ITS_sensor_DNA_MDS1, pred.var = "meansoiltemp", train = rf_df_ITS_sensor_DNA)
pdp_ITS_RNA <- partial(rf_ITS_sensor_RNA_MDS1, pred.var = "meansoiltemp", train = rf_df_ITS_sensor_RNA)
pdp_16S_DNA <- partial(rf_16S_sensor_DNA_MDS1, pred.var = "meansoiltemp", train = rf_df_16S_sensor_DNA)
pdp_16S_RNA <- partial(rf_16S_sensor_RNA_MDS1, pred.var = "meansoiltemp", train = rf_df_16S_sensor_RNA)

# Combine all the partial dependence data into one data frame
pdp_combined <- rbind(
  data.frame(meansoiltemp = pdp_ITS_DNA$meansoiltemp, yhat = pdp_ITS_DNA$yhat, model = "ITS_DNA"),
  data.frame(meansoiltemp = pdp_ITS_RNA$meansoiltemp, yhat = pdp_ITS_RNA$yhat, model = "ITS_RNA"),
  data.frame(meansoiltemp = pdp_16S_DNA$meansoiltemp, yhat = pdp_16S_DNA$yhat, model = "16S_DNA"),
  data.frame(meansoiltemp = pdp_16S_RNA$meansoiltemp, yhat = pdp_16S_RNA$yhat, model = "16S_RNA")
)

##########################
#creates histogram to plot above partial dependence plot
# 1. Histogram of meansoiltemp
density_plot <- ggplot(rf_df_ITS_sensor_DNA, aes(x = meansoiltemp)) +
  geom_density(fill = "grey60", color = "black", alpha = 0.7, size=0.6) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    panel.grid = element_blank(),
    axis.line.y = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    plot.margin = ggplot2::margin(0, 10, 0, 10)
  ) +
  labs(y = "Density") +
  scale_x_continuous(expand = c(0, 0)) +
  theme(aspect.ratio = 0.25)

pdp_colors=c(
  "ITS_DNA"="#1b9e77",
  "ITS_RNA"="#a6dcd1",
  "16S_DNA" = "#7570b3",
  "16S_RNA"="#b3b3cc"
)

pdp_plot <- ggplot(pdp_combined, aes(x = meansoiltemp, y = yhat, color = model)) +
  geom_rect(aes(xmin = 13.25, xmax = 14, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE,
            fill = "red", alpha = 0.003) +
  geom_smooth(se = FALSE, method = "loess", span = 0.2, linewidth = 1.5) +
  labs(
    x = "Mean Soil Temperature",
    y = "Partial Dependence (MDS1)",
    color = "Model"
  ) +
  theme_test() +
  theme(
    legend.position = "bottom",
    plot.margin = ggplot2::margin(0, 10, 10, 10)
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_color_manual(values = pdp_colors) +
  theme(aspect.ratio = 1)

# Combine histogram and PDP
combined_density_pdp_plot=
  density_plot / pdp_plot + plot_layout(heights = c(0.75, 3))

#saves partial dependence plot as pdf to be used in figure
pdf(file="combined_density_pdp_plot.pdf", 
    width = 4, height = 4.3)
combined_density_pdp_plot
dev.off()


