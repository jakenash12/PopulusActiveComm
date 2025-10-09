#Makes some plots displaying response of microbial lineages to soil factors
#that were identified as significant by random forest modelling
#1st plot: response of Alphaproteobacteria to mean soil temp
#2nd plot: response of Dothideomycetes to mean soil temp

#This script requires that you run the scripts 16S_Analysis_SILVA.R,
#ITS2_Analysis.R, PlotBacterialPhyla.R, RandomForests_Taxa.R and 
#SoilVarDimensionalReduction.R and have the resulting data objects loaded in memory

library(tidyverse)
library(lme4)

#make the plot 
Alphaproteobacteria_meansoiltemp_plot=
  AUE2021_16S_rarefied_phylum_sensor_t %>%
  mutate(source=factor(source,levels=c("DNA", "cDNA"))) %>%
  ggplot(aes(x=meansoiltemp, 
           y =Alphaproteobacteria/13407)) +
  geom_jitter(size = 4, width = 0.1, aes(color = meansoiltemp)) +  # <- jitter added here
  geom_smooth(method="lm", se =FALSE, color="black") +
  theme_test() +
  scale_color_gradientn(
    colors = c("#66a866", "#b2df8a", "#fc8d62"),  # light green → medium green → orange
    values = scales::rescale(c(11.5, 12.5, 15.5))     # adjust as needed based on your data range
  ) +
  theme(legend.position="none", aspect.ratio=2) +
  facet_grid(.~source, scales="free") + 
  scale_y_continuous(breaks = c(0.11, 0.15, 0.19, 0.23))

#do linear mixed models of Alphaproteobacteria relative abundances vs soil temp
#in active and total community
AUE2021_16S_rarefied_phylum_sensor_t %>%
  filter(source=="cDNA") %>%
  lmer(Alphaproteobacteria ~ meansoiltemp + (1|Site),.) %>% 
  Anova

#make another plot
Dothideomycetes_meansoiltemp_plot=
  AUE2021_ITS_rarefied_class_sensor_t %>%
  filter(source=="cDNA") %>%
  ggplot(aes(x=meansoiltemp, 
             y =Dothideomycetes/17083 )) +
  geom_jitter(size = 5, width = 0.1, aes(color = meansoiltemp)) +  # <- jitter added here
  geom_smooth(method="lm", se =FALSE, color="black") +
  scale_color_gradientn(
    colors = c("#66a866", "#b2df8a", "#fc8d62"),  # light green → medium green → orange
    values = scales::rescale(c(11.5, 12.5, 15.5))     # adjust as needed based on your data range
  ) +
  theme_test() +
  theme(legend.position = "none", aspect.ratio = 1)
  facet_grid(source~Time, scales="free")

pdf(file="Alphaproteobacteria_meansoiltemp_plot.pdf", 
    width = 4.8, height = 4.8)
Alphaproteobacteria_meansoiltemp_plot
dev.off()

pdf(file="Dothideomycetes_meansoiltemp_plot.pdf", 
    width = 4.8, height = 4.8)
Dothideomycetes_meansoiltemp_plot
dev.off()

