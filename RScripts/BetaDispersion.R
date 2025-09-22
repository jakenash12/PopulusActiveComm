#script to test for differences in within-site beta dispersion
#across site types, seasons, and between RNA/DNA

#This script calculates the distance to centroids for each
#unique combination of site/season/source and then uses a 
#mixed model to test whether distance to centroids varies

library(vegan)
library(lme4)
library(tidyverse)
library(magrittr)
library(car)

##############ITS Beta Dispersion#############################

#use the variable SiteTimeSource 
#(combination of Site, Time and source variables)
#as grouping variable for beta dispersion function from vegan
betadisper_ITS_SiteTimeSource=
  betadisper(AUE2021_ITS_rarefied_dm_bc, 
             AUE2021_ITS_rarefied_dm_metadata$SiteTimeSource)

#formats distance to centroids as a dataframe and joins
#with metadata
betadisper_ITS_SiteTimeSource_df=
  betadisper_ITS_SiteTimeSource$distances %>%
  as.data.frame %>% 
  set_colnames("DistCentr") %>%
  mutate(sample=rownames(.)) %>%
  left_join(AUE2021_ITS_rarefied_dm_metadata)

#makes mixed model predicting distance to centroid
#with random effect of plot
mod=lmer(DistCentr~SiteType*Time*source+(1|Plot), 
         betadisper_ITS_SiteTimeSource_df)
Anova(mod)

betadisper_ITS_SiteTimeSource_df %>%
  ggplot(aes(x=Time, y=DistCentr, color=SiteType)) +
  geom_boxplot() +
  facet_grid(.~source)


Betadisp_palette=c("High Elevation" = "#60bba0", 
                   "Mid Elevation" = "#e48f4e",
                   "Riparian"="#9f9bca")
betadisp_ITS_plot=
  betadisper_ITS_SiteTimeSource_df %>%
  group_by(SiteType, Time, source) %>%
  summarise(mean_betadisp=mean(DistCentr),
            se_betadisp=sd(DistCentr)/sqrt(n())) %>%
  mutate(source=factor(source, levels=c("DNA", "cDNA"))) %>%
  ggplot(aes(x=Time, y=mean_betadisp, color=SiteType, group=SiteType)) +
  geom_line(size=1.5,position=position_dodge(width=Dodge)) +
  geom_errorbar(aes(ymin=mean_betadisp-se_betadisp, ymax=mean_betadisp+se_betadisp), 
                width=0, size=1.5, position=position_dodge(width=Dodge)) +
  geom_point(stroke=1.5, size=3, position=position_dodge(width=Dodge), 
             shape=21, color="gray20", aes(fill=SiteType)) +
  theme_test() +
  theme(legend.position="none") +
  theme(panel.border = element_rect(linetype = "solid",
                                    colour = "gray20", linewidth = 1.5),
        strip.background = element_blank(),
        text=element_blank(),
        axis.ticks=element_line(size=1.5, colour = "gray20")) +
  scale_color_manual(values=Betadisp_palette) +
  scale_fill_manual(values=Betadisp_palette) +
  facet_wrap(.~source)

betadisper_ITS_SiteTimeSource_df %>%
  select(sample, DistCentr) %>%
  left_join(AUE2021_alphadiv) %>%
  filter(!is.na(source)) %>%
  ggplot(aes(x=Shannon, y=DistCentr)) +
  geom_point() +
  facet_wrap(.~source, scales="free") +
  geom_smooth(method="lm", se=FALSE)

betadisper_ITS_SiteTimeSource_df %>%
  select(sample, DistCentr) %>%
  left_join(AUE2021_alphadiv) %>%
  lm(DistCentr~Shannon+source, .) %>%
  Anova

betadisper_ITS_SiteTimeSource_df %>%
  select(sample, DistCentr) %>%
  left_join(AUE2021_alphadiv) %>%
  filter(source=="RNA") %>%
  lm(DistCentr~Shannon, .) %>%
  summary

betadisper_ITS_SiteTimeSource_df %>%
  select(sample, DistCentr) %>%
  left_join(AUE2021_alphadiv) %>%
  filter(source=="RNA") %>%
  lm(DistCentr~Shannon, .) %>%
  summary
##############16S Beta Dispersion#############################

#use the variable SiteTimeSource 
#(combination of Site, Time and source variables)
#as grouping variable for beta dispersion function from vegan
betadisper_16S_SiteTimeSource=
  betadisper(AUE2021_16S_rarefied_dm_bc, 
             AUE2021_16S_rarefied_dm_metadata$SiteTimeSource)

#formats distance to centroids as a dataframe and joins
#with metadata
betadisper_16S_SiteTimeSource_df=
  betadisper_16S_SiteTimeSource$distances %>%
  as.data.frame %>% 
  set_colnames("DistCentr") %>%
  mutate(sample=rownames(.)) %>%
  left_join(AUE2021_16S_rarefied_dm_metadata)

#makes mixed model predicting distance to centroid
#with random effect of plot
mod=lmer(DistCentr~SiteType*Time*source+(1|Plot), 
         betadisper_16S_SiteTimeSource_df)
Anova(mod)

betadisp_16S_plot=
  betadisper_16S_SiteTimeSource_df %>%
  group_by(SiteType, Time, source) %>%
  summarise(mean_betadisp=mean(DistCentr),
            se_betadisp=sd(DistCentr)/sqrt(n())) %>%
  mutate(source=factor(source, levels=c("DNA", "cDNA"))) %>%
  ggplot(aes(x=Time, y=mean_betadisp, color=SiteType, group=SiteType)) +
  geom_line(size=1.5,position=position_dodge(width=Dodge)) +
  geom_errorbar(aes(ymin=mean_betadisp-se_betadisp, ymax=mean_betadisp+se_betadisp), 
                width=0, size=1.5, position=position_dodge(width=Dodge)) +
  geom_point(stroke=1.5, size=3, position=position_dodge(width=Dodge), 
             shape=21, color="gray20", aes(fill=SiteType)) +
  theme_test() +
  theme(legend.position="none") +
  theme(panel.border = element_rect(linetype = "solid",
                                    colour = "gray20", linewidth = 1.5),
        strip.background = element_blank(),
        text=element_blank(),
        axis.ticks=element_line(size=1.5, colour = "gray20")) +
  scale_color_manual(values=Betadisp_palette) +
  scale_fill_manual(values=Betadisp_palette) +
  facet_wrap(.~source)

pdf(file="C:/Users/akeja/OneDrive - Duke University/Documents/Duke PhD/Projects/PMI/AUE_2021/MetabarcodingManuscript/Misc Figure, Tables, Etc/betadisp_16S_plot.pdf", 
    width = 5.11811, height = 4.33071)
betadisp_16S_plot
dev.off()

pdf(file="C:/Users/akeja/OneDrive - Duke University/Documents/Duke PhD/Projects/PMI/AUE_2021/MetabarcodingManuscript/Misc Figure, Tables, Etc/betadisp_ITS_plot.pdf", 
    width = 5.11811, height = 4.33071)
betadisp_ITS_plot
dev.off()

betadisper_16S_SiteTimeSource_df %>%
  select(sample, DistCentr) %>%
  left_join(AUE2021_16S_alphadiv) %>%
  ggplot(aes(x=Shannon, y=DistCentr)) +
  geom_point() +
  facet_wrap(.~source, scales="free")


betadisper_16S_SiteTimeSource_df %>%
  select(sample, DistCentr) %>%
  left_join(AUE2021_16S_alphadiv) %>%
  filter(source=="RNA") %>%
  lm(DistCentr~Shannon, .) %>%
  Anova

betadisper_16S_SiteTimeSource_df %>%
  select(sample, DistCentr) %>%
  left_join(AUE2021_16S_alphadiv) %>%
  filter(source=="RNA") %>%
  lm(DistCentr~Shannon, .) %>%
  summary
