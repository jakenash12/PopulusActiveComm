#Conducts PERMANOVA and NMDS to test for effects
#of site type and season
#
#This script relies on objects created in the 
#ITS2_Analysis.R and 16S_Analysis.R scripts.
#Make sure to run those first

#library(cowplot)
library(vegan)
library(phyloseq)
library(tidyverse)
#library(RRPP)
library(magrittr)
##########################ITS PERMANOVA##################################
#creates otu table df from phyloseq object - full dataset
AUE2021_ITS_rarefied_df=
  AUE2021_ITS_rarefied %>%
  otu_table %>%
  t %>%
  as.data.frame

#creates otu table df from phyloseq object - DNA only
AUE2021_ITS_rarefied_df_DNA=
  AUE2021_ITS_rarefied %>%
  subset_samples(source=="DNA") %>%
  otu_table %>%
  t %>%
  as.data.frame

#creates otu table df from phyloseq object - RNA only
AUE2021_ITS_rarefied_df_RNA=
  AUE2021_ITS_rarefied %>%
  subset_samples(source=="cDNA") %>%
  otu_table %>%
  t %>%
  as.data.frame

#creates bray-curtis dist matrix - full dataset
AUE2021_ITS_rarefied_dm_bc=
  AUE2021_ITS_rarefied %>%
  otu_table %>%
  t %>%
  as.data.frame %>%
  metaMDSdist(distance="bray")

#creates bray-curtis dist matrix - RNA-only
AUE2021_ITS_rarefied_RNA_dm_bc=
  AUE2021_ITS_rarefied %>%
  subset_samples(source=="cDNA") %>%
  otu_table %>%
  t %>%
  as.data.frame %>%
  metaMDSdist(distance="bray")

#creates bray-curtis dist matrix - DNA-only
AUE2021_ITS_rarefied_DNA_dm_bc=
  AUE2021_ITS_rarefied %>%
  subset_samples(source=="DNA") %>%
  otu_table %>%
  t %>%
  as.data.frame %>%
  metaMDSdist(distance="bray")

#formats metadata to match order in dist mat - full dataset
AUE2021_ITS_rarefied_dm_metadata =
  AUE2021_ITS_rarefied_dm_bc %>%
  as.matrix %>%
  as.data.frame %>% 
  mutate(sample=rownames(.)) %>%
  select(sample) %>%
  left_join(samples_df2) %>%
  mutate(SiteTimeSource=paste(Site, Time, source, sep="_"),
         PlotSource=paste(Plot, source,sep="_"))

#formats metadata to match order in dist mat - RNA-only
AUE2021_ITS_rarefied_RNA_dm_metadata =
  AUE2021_ITS_rarefied_RNA_dm_bc %>%
  as.matrix %>%
  as.data.frame %>% 
  mutate(sample=rownames(.)) %>%
  select(sample) %>%
  left_join(samples_df2)

#formats metadata to match order in dist mat - DNA-only
AUE2021_ITS_rarefied_DNA_dm_metadata =
  AUE2021_ITS_rarefied_DNA_dm_bc %>%
  as.matrix %>%
  as.data.frame %>% 
  mutate(sample=rownames(.)) %>%
  select(sample) %>%
  left_join(samples_df2)

#sets seed for PERMANOVA
set.seed(72)
#runs permanova-full dataset
ITS_adonis=
  adonis2(AUE2021_ITS_rarefied_dm_bc ~ 
            AUE2021_ITS_rarefied_dm_metadata$SiteType+
            AUE2021_ITS_rarefied_dm_metadata$Site+
            AUE2021_ITS_rarefied_dm_metadata$Plot+
            AUE2021_ITS_rarefied_dm_metadata$Time+
            AUE2021_ITS_rarefied_dm_metadata$source)

#sets seed
set.seed(72)
#runs permanova - DNA-only
ITS_DNA_adonis=
  adonis2(AUE2021_ITS_rarefied_DNA_dm_bc ~ 
            AUE2021_ITS_rarefied_DNA_dm_metadata$SiteType+
            AUE2021_ITS_rarefied_DNA_dm_metadata$Site+
            AUE2021_ITS_rarefied_DNA_dm_metadata$Plot+
            AUE2021_ITS_rarefied_DNA_dm_metadata$Time)
ITS_DNA_adonis

#sets seed
set.seed(72)
#runs permanova - RNA-only
ITS_RNA_adonis=
  adonis2(AUE2021_ITS_rarefied_RNA_dm_bc ~ 
            AUE2021_ITS_rarefied_RNA_dm_metadata$SiteType+
            AUE2021_ITS_rarefied_RNA_dm_metadata$Site+
            AUE2021_ITS_rarefied_RNA_dm_metadata$Plot+
            AUE2021_ITS_rarefied_RNA_dm_metadata$Time)
ITS_RNA_adonis

#conducts test of beta dispersion by site type and
#season to validate permanova results
ITS_betadisper_SiteType=
  betadisper(AUE2021_ITS_rarefied_dm_bc, 
             AUE2021_ITS_rarefied_dm_metadata$SiteType)
ITS_betadisper_Time=
  betadisper(AUE2021_ITS_rarefied_dm_bc, 
             AUE2021_ITS_rarefied_dm_metadata$Time)

ITS_betadisper_SiteTimeSource=
  betadisper(AUE2021_ITS_rarefied_dm_bc, 
             AUE2021_ITS_rarefied_dm_metadata$SiteTimeSource)

anova(ITS_betadisper_SiteType)
anova(ITS_betadisper_Time)

#conducts tests of beta disp in DNA-only dataset
#by season
betadisper_Time_DNA_ITS=
  betadisper(AUE2021_ITS_rarefied_DNA_dm_bc, 
             AUE2021_ITS_rarefied_DNA_dm_metadata$Time)
anova(betadisper_Time_DNA_ITS)

#conducts tests of beta disp in RNA-only dataset
#by season
betadisper_Time_RNA_ITS=
  betadisper(AUE2021_ITS_rarefied_RNA_dm_bc, 
             AUE2021_ITS_rarefied_RNA_dm_metadata$Time)
anova(betadisper_Time_RNA_ITS)

#formats distance to season centroids as a dataframe
BetaDisp_ITS=
  ITS_betadisper_Time$distances %>% 
  as.data.frame %>%
  set_colnames("Distance") %>%
  mutate(sample=rownames(.), amplicon="ITS")%>% 
  left_join(select(AUE2021_ITS_rarefied_dm_metadata, sample, Time, SiteType, source))

#conducts NMDS on ITS DNA and ITS RNA dataframes
set.seed(72)
ITS_DNA_PCoA=metaMDS(AUE2021_ITS_rarefied_df_DNA, k=10, trymax=100)
set.seed(72)
ITS_RNA_PCoA=metaMDS(AUE2021_ITS_rarefied_df_RNA, k=10, trymax=100)

#formats ITS DNA NMDS as a dataframe
ITS_DNA_PCoA_Points=
  ITS_DNA_PCoA$points %>%
  as.data.frame %>%
  mutate(sample=rownames(.)) %>%
  left_join(samples_df2) %>%
  mutate(SourceSiteType=paste(SiteType, source),
         SourceSiteTypeTime=paste(SiteType, source, Time))

#formats ITS RNA NMDS as a dataframe
ITS_RNA_PCoA_Points=
  ITS_RNA_PCoA$points %>%
  as.data.frame %>%
  mutate(sample=rownames(.)) %>%
  left_join(samples_df2) %>%
  mutate(SourceSiteType=paste(SiteType, source),
         SourceSiteTypeTime=paste(SiteType, source, Time))

#conducts NMDS on full ITS dataset
set.seed(72)
ITS_PCoA=metaMDS(AUE2021_ITS_rarefied_df, k=10, trymax=100)
ITS_PCoA_Points=
  ITS_PCoA$points %>%
  as.data.frame %>%
  mutate(sample=rownames(.)) %>%
  left_join(samples_df2) %>%
  mutate(source=factor(source, levels=c("DNA", "cDNA")),
         SourceSiteType=paste(SiteType, source),
         SourceSiteTypeTime=paste(SiteType, source, Time)) %>%
  mutate(Time_source=paste(Time, source, sep="_"))

#plots ITS NMDS
alpha_nmds=0.7
pointscale=0.7
cDNA_penalty=0.8
ITS_NMDS_open=
  mutate(ITS_PCoA_Points, Time_source=paste(Time, source, sep="_")) %>%
  ggplot(aes(x = -MDS1, y = MDS2)) +
  geom_point(data=filter(ITS_PCoA_Points, source=="DNA", Time!="June"), 
             aes(shape = Time_source, color = SiteType, fill=SiteType), 
             size = 5*pointscale, stroke = 1, alpha = alpha_nmds, color="grey20") +
  geom_point(data=filter(ITS_PCoA_Points, source=="DNA", Time=="June"), 
             aes(shape = Time_source, color = SiteType, fill=SiteType), 
             size = 6*pointscale, stroke = 1, alpha = alpha_nmds, color="grey20") +
  geom_point(data=filter(ITS_PCoA_Points, source=="cDNA", Time!="June"), 
             aes(shape = Time_source, color = SiteType, fill=SiteType), 
             size = 3*pointscale, stroke = 2.5, alpha = alpha_nmds) +
  geom_point(data=filter(ITS_PCoA_Points, source=="cDNA", Time=="June"), 
             aes(shape = Time_source, color = SiteType, fill=SiteType), 
             size = 4*pointscale, stroke = 2.5, alpha = alpha_nmds) +
  geom_point(data=filter(ITS_PCoA_Points, source=="cDNA", Time!="June"), 
             aes(shape = Time_source), 
             size = (3*pointscale)+2.3, stroke = 1, alpha = alpha_nmds, color="grey20") +
  geom_point(data=filter(ITS_PCoA_Points, source=="cDNA", Time=="June"), 
             aes(shape = Time_source, color = SiteType, fill=SiteType), 
             size = (4*pointscale)+2.3, stroke = 1, alpha = alpha_nmds, color="grey20") +
  theme_test() +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  scale_shape_manual(values = c("June_DNA" = 21, 
                                "June_cDNA" = 1,
                                "August_DNA"=23,
                                "August_cDNA"=5,
                                "October_DNA"=24,
                                "October_cDNA"=2)) +
  theme(legend.position = "none", text = element_blank()) +
  theme(panel.border = element_rect(linetype = "solid",
                                    colour = "gray20", linewidth = 1.5),
        strip.background = element_blank(),
        axis.ticks = element_line(size = 1.5, colour = "gray20"),
        aspect.ratio = 1)  +
  stat_ellipse(aes(linetype=source, group=source), 
               level=0.95, linewidth=1.5, colour = "gray20")
  
ITS_NMDS_open

###################################16S PERMANOVA#################
#converts phyloseq otu table to a dataframe and transposes
AUE2021_16S_rarefied_df=
  AUE2021_16S_rarefied %>%
  otu_table %>%
  t %>%
  as.data.frame

#makes a DNA-only dataframe from phyloseq otu table and transposes
AUE2021_16S_rarefied_df_DNA=
  AUE2021_16S_rarefied %>%
  subset_samples(source=="DNA") %>%
  otu_table %>%
  t %>%
  as.data.frame

#makes an RNA-only dataframe from phyloseq otu table and transposes
AUE2021_16S_rarefied_df_RNA=
  AUE2021_16S_rarefied %>%
  subset_samples(source=="cDNA") %>%
  otu_table %>%
  t %>%
  as.data.frame

#greates an RNA-only bray-curtis dist matrix
AUE2021_16S_rarefied_RNA_dm_bc=
  AUE2021_16S_rarefied %>%
  subset_samples(source=="cDNA") %>%
  otu_table %>%
  t %>%
  as.data.frame %>% 
  metaMDSdist(distance="bray")

#greates a DNA-only bray-curtis dist matrix
AUE2021_16S_rarefied_DNA_dm_bc=
  AUE2021_16S_rarefied %>%
  subset_samples(source=="DNA") %>%
  otu_table %>%
  t %>%
  as.data.frame %>%
  metaMDSdist(distance="bray")

#creates a bray-curtis dist matrix from full dataset
AUE2021_16S_rarefied_dm_bc=
  AUE2021_16S_rarefied_df %>%
  metaMDSdist(distance="bray")

#joins metadata w/ dist matrix to ensure rows are in same order
AUE2021_16S_rarefied_dm_metadata =
  AUE2021_16S_rarefied_dm_bc %>%
  as.matrix %>%
  as.data.frame %>% 
  mutate(sample=rownames(.)) %>%
  select(sample) %>%
  left_join(samples_df2_16S) %>%
  mutate(SiteTimeSource=paste(Site, Time, source, sep="_"),
         PlotSource=paste(Plot, source, sep="_"))

#joins metadata w/ RNA-only dist matrix to ensure rows are in same order
AUE2021_16S_rarefied_RNA_dm_metadata =
  AUE2021_16S_rarefied_RNA_dm_bc %>%
  as.matrix %>%
  as.data.frame %>% 
  mutate(sample=rownames(.)) %>%
  select(sample) %>%
  left_join(samples_df2_16S)

#joins metadata w/ DNA-only dist matrix to ensure rows are in same order
AUE2021_16S_rarefied_DNA_dm_metadata =
  AUE2021_16S_rarefied_DNA_dm_bc %>%
  as.matrix %>%
  as.data.frame %>% 
  mutate(sample=rownames(.)) %>%
  select(sample) %>%
  left_join(samples_df2_16S)

#sets random seed to ensure reproducibiltiy
set.seed(72)
#conducts PERMANOVA on full dataset
adonis_16S=
  adonis2(AUE2021_16S_rarefied_dm_bc ~ 
            AUE2021_16S_rarefied_dm_metadata$SiteType+
            AUE2021_16S_rarefied_dm_metadata$Site+
            AUE2021_16S_rarefied_dm_metadata$Plot+
            AUE2021_16S_rarefied_dm_metadata$Time+
            AUE2021_16S_rarefied_dm_metadata$source)
adonis_16S

#sets random seed to ensure reproducibility
set.seed(72)
#conducts PERMANOVA on DNA-only dataset
adonis_DNA_16S=
  adonis2(AUE2021_16S_rarefied_DNA_dm_bc ~ 
            AUE2021_16S_rarefied_DNA_dm_metadata$SiteType+
            AUE2021_16S_rarefied_DNA_dm_metadata$Site+
            AUE2021_16S_rarefied_DNA_dm_metadata$Plot+
            AUE2021_16S_rarefied_DNA_dm_metadata$Time)
adonis_DNA_16S

#sets random seed to ensure reproducibility
set.seed(72)
#conducts PERMANOVA on RNA-only dataset
adonis_RNA_16S=
  adonis2(AUE2021_16S_rarefied_RNA_dm_bc ~ 
            AUE2021_16S_rarefied_RNA_dm_metadata$SiteType+
            AUE2021_16S_rarefied_RNA_dm_metadata$Site+
            AUE2021_16S_rarefied_RNA_dm_metadata$Plot+
            AUE2021_16S_rarefied_RNA_dm_metadata$Time)
adonis_RNA_16S

betadisper_SiteType_16S=
  betadisper(AUE2021_16S_rarefied_dm_bc, 
             AUE2021_16S_rarefied_dm_metadata$SiteType)
anova(betadisper_SiteType_16S)

betadisper_Time_DNA_16S=
  betadisper(AUE2021_16S_rarefied_DNA_dm_bc, 
             AUE2021_16S_rarefied_DNA_dm_metadata$Time)
anova(betadisper_Time_DNA_16S)

betadisper_Time_RNA_16S=
  betadisper(AUE2021_16S_rarefied_RNA_dm_bc, 
             AUE2021_16S_rarefied_RNA_dm_metadata$Time)
anova(betadisper_Time_RNA_16S)


betadisper_SiteType_DNA_16S=
  betadisper(AUE2021_16S_rarefied_DNA_dm_bc, 
             AUE2021_16S_rarefied_DNA_dm_metadata$SiteType)
anova(betadisper_SiteType_DNA_16S)

betadisper_SiteType_RNA_16S=
  betadisper(AUE2021_16S_rarefied_RNA_dm_bc, 
             AUE2021_16S_rarefied_RNA_dm_metadata$SiteType)
anova(betadisper_SiteType_RNA_16S)

set.seed(72)
PCoA_DNA_16S=metaMDS(AUE2021_16S_rarefied_DNA_dm_bc, k=10, trymax=100)
set.seed(72)
PCoA_RNA_16S=metaMDS(AUE2021_16S_rarefied_RNA_dm_bc, k=10, trymax=100)

PCoA_16S_DNA_Points=
  PCoA_DNA_16S$points %>%
  as.data.frame %>%
  mutate(sample=rownames(.)) %>%
  left_join(samples_df2_16S) %>%
  mutate(SourceSiteType=paste(SiteType, source),
         SourceSiteTypeTime=paste(SiteType, source, Time))

PCoA_16S_RNA_Points=
  PCoA_RNA_16S$points %>%
  as.data.frame %>%
  mutate(sample=rownames(.)) %>%
  left_join(samples_df2_16S) %>%
  mutate(SourceSiteType=paste(SiteType, source),
         SourceSiteTypeTime=paste(SiteType, source, Time))

ggplot(PCoA_16S_DNA_Points, aes(x=-MDS1, y=MDS2))+
  geom_point(aes(shape=Time, fill=Time), 
             size=3, stroke=1,alpha=0.8, color="black") +
  theme_test()+
  scale_fill_brewer(palette = "Dark2") +
  scale_shape_manual(values = c("June" = 21, "August" = 23, "October" = 24)) +
  #theme(legend.position="none", text=element_blank()) +
  theme(panel.border = element_rect(linetype = "solid",
                                    colour = "gray20", linewidth = 1.5),
        strip.background = element_blank(),
        axis.ticks=element_line(size=1.5, colour = "gray20"),
        aspect.ratio=1) +
  facet_wrap(.~Site, scales="free")

set.seed(72)
PCoA_16S=metaMDS(AUE2021_16S_rarefied_df, k=10, trymax=100)
PCoA_16S_Points=
  PCoA_16S$points %>%
  as.data.frame %>%
  mutate(sample=rownames(.)) %>%
  left_join(samples_df2_16S) %>%
  mutate(source=factor(source, levels=c("DNA", "cDNA")),
         SourceSiteType=paste(SiteType, source),
         SourceSiteTypeTime=paste(SiteType, source, Time),
         Time_source=paste(Time, source, sep="_")) 
  
alpha_nmds=0.7
pointscale=0.7
NMDS_16S_open_plot=
  mutate(PCoA_16S_Points, Time_source=paste(Time, source, sep="_")) %>%
  ggplot(aes(x = -MDS1, y = MDS2)) +
  geom_point(data=filter(PCoA_16S_Points, source=="DNA", Time!="June"), 
             aes(shape = Time_source, color = SiteType, fill=SiteType), 
             size = 5*pointscale, stroke = 1, alpha = alpha_nmds, color="grey20") +
  geom_point(data=filter(PCoA_16S_Points, source=="DNA", Time=="June"), 
             aes(shape = Time_source, color = SiteType, fill=SiteType), 
             size = 6*pointscale, stroke = 1, alpha = alpha_nmds, color="grey20") +
  geom_point(data=filter(PCoA_16S_Points, source=="cDNA", Time!="June"), 
             aes(shape = Time_source, color = SiteType, fill=SiteType), 
             size = 3*pointscale, stroke = 2.5, alpha = alpha_nmds) +
  geom_point(data=filter(PCoA_16S_Points, source=="cDNA", Time=="June"), 
             aes(shape = Time_source, color = SiteType, fill=SiteType), 
             size = 4*pointscale, stroke = 2.5, alpha = alpha_nmds) +
  geom_point(data=filter(PCoA_16S_Points, source=="cDNA", Time!="June"), 
             aes(shape = Time_source), 
             size = (3*pointscale)+2.3, stroke = 1, alpha = alpha_nmds, color="grey20") +
  geom_point(data=filter(PCoA_16S_Points, source=="cDNA", Time=="June"), 
             aes(shape = Time_source, color = SiteType, fill=SiteType), 
             size = (4*pointscale)+2.3, stroke = 1, alpha = alpha_nmds, color="grey20") +
  theme_test() +
  #scale_color_manual(values=Troye_palette) +
  #scale_fill_manual(values=Troye_palette) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  scale_shape_manual(values = c("June_DNA" = 21, 
                                "June_cDNA" = 1,
                                "August_DNA"=23,
                                "August_cDNA"=5,
                                "October_DNA"=24,
                                "October_cDNA"=2)) +
  theme(legend.position = "none", text = element_blank()) +
  theme(panel.border = element_rect(linetype = "solid",
                                    colour = "gray20", linewidth = 1.5),
        strip.background = element_blank(),
        axis.ticks = element_line(size = 1.5, colour = "gray20"),
        aspect.ratio = 1)  +
  stat_ellipse(aes(linetype=source, group=source), 
               level=0.95, linewidth=1.5, colour = "gray20")

NMDS_16S_open_plot

#conducts test of beta dispersion by site type and
#season to validate permanova results
betadisper_SiteType_16S=
  betadisper(AUE2021_16S_rarefied_dm_bc, 
             AUE2021_16S_rarefied_dm_metadata$SiteType)
betadisper_Time_16S=
  betadisper(AUE2021_16S_rarefied_dm_bc, 
             AUE2021_16S_rarefied_dm_metadata$Time)

betadisper_SiteTimeSource_16S=
  betadisper(AUE2021_16S_rarefied_dm_bc, 
             AUE2021_16S_rarefied_dm_metadata$SiteTimeSource)

anova(betadisper_SiteType_16S)
anova(betadisper_Time_16S)

##############################
BetaDisp_16S_ITS=
  rbind(BetaDisp_16S, BetaDisp_ITS)

ggplot(BetaDisp_16S_ITS, aes(x=Time, y=Distance)) +
  geom_boxplot(aes(color=SiteType), size=1) +
  theme_test() +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position="none") +
  facet_grid(amplicon~source, scales="free_y")

plot_grid(PCoA_16S_plot, PCoA_16S_plot)
plot_grid(BetaDisp_ITS, BetaDisp_16S)
  