library(cowplot)
library(vegan)
library(phyloseq)
library(tidyverse)
library(RRPP)
library(magrittr)
##########################ITS PERMANOVA##################################
AUE2021_ITS_rarefied_df=
  AUE2021_ITS_rarefied %>%
  otu_table %>%
  t %>%
  as.data.frame

AUE2021_ITS_rarefied_df_DNA=
  AUE2021_ITS_rarefied %>%
  subset_samples(source=="DNA") %>%
  otu_table %>%
  t %>%
  as.data.frame

AUE2021_ITS_rarefied_df_RNA=
  AUE2021_ITS_rarefied %>%
  subset_samples(source=="cDNA") %>%
  otu_table %>%
  t %>%
  as.data.frame

AUE2021_ITS_rarefied_RNA_dm_bc=
  AUE2021_ITS_rarefied %>%
  subset_samples(source=="cDNA") %>%
  otu_table %>%
  t %>%
  as.data.frame %>%
  metaMDSdist(distance="bray")

AUE2021_ITS_rarefied_DNA_dm_bc=
  AUE2021_ITS_rarefied %>%
  subset_samples(source=="DNA") %>%
  otu_table %>%
  t %>%
  as.data.frame %>%
  metaMDSdist(distance="bray")

AUE2021_ITS_rarefied_dm_bc=
  AUE2021_ITS_rarefied %>%
  otu_table %>%
  t %>%
  as.data.frame %>%
  metaMDSdist(distance="bray")

AUE2021_ITS_rarefied_dm_metadata =
  AUE2021_ITS_rarefied_dm_bc %>%
  as.matrix %>%
  as.data.frame %>% 
  mutate(sample=rownames(.)) %>%
  select(sample) %>%
  left_join(samples_df2) %>%
  mutate(SiteTimeSource=paste(Site, Time, source, sep="_"),
         PlotSource=paste(Plot, source,sep="_"))


AUE2021_ITS_rarefied_RNA_dm_metadata =
  AUE2021_ITS_rarefied_RNA_dm_bc %>%
  as.matrix %>%
  as.data.frame %>% 
  mutate(sample=rownames(.)) %>%
  select(sample) %>%
  left_join(samples_df2)

AUE2021_ITS_rarefied_DNA_dm_metadata =
  AUE2021_ITS_rarefied_DNA_dm_bc %>%
  as.matrix %>%
  as.data.frame %>% 
  mutate(sample=rownames(.)) %>%
  select(sample) %>%
  left_join(samples_df2)

set.seed(72)

ITS_adonis=
  adonis2(AUE2021_ITS_rarefied_dm_bc ~ 
            AUE2021_ITS_rarefied_dm_metadata$SiteType+
            AUE2021_ITS_rarefied_dm_metadata$Site+
            AUE2021_ITS_rarefied_dm_metadata$Plot+
            AUE2021_ITS_rarefied_dm_metadata$Time+
            AUE2021_ITS_rarefied_dm_metadata$source)

set.seed(72)
ITS_DNA_adonis=
  adonis2(AUE2021_ITS_rarefied_DNA_dm_bc ~ 
            AUE2021_ITS_rarefied_DNA_dm_metadata$SiteType+
            AUE2021_ITS_rarefied_DNA_dm_metadata$Site+
            AUE2021_ITS_rarefied_DNA_dm_metadata$Plot+
            AUE2021_ITS_rarefied_DNA_dm_metadata$Time)
ITS_DNA_adonis

set.seed(72)
ITS_RNA_adonis=
  adonis2(AUE2021_ITS_rarefied_RNA_dm_bc ~ 
            AUE2021_ITS_rarefied_RNA_dm_metadata$SiteType+
            AUE2021_ITS_rarefied_RNA_dm_metadata$Site+
            AUE2021_ITS_rarefied_RNA_dm_metadata$Plot+
            AUE2021_ITS_rarefied_RNA_dm_metadata$Time)
ITS_RNA_adonis

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

betadisper_Time_DNA_ITS=
  betadisper(AUE2021_ITS_rarefied_DNA_dm_bc, 
             AUE2021_ITS_rarefied_DNA_dm_metadata$Time)
anova(betadisper_Time_DNA_ITS)

betadisper_Time_RNA_ITS=
  betadisper(AUE2021_ITS_rarefied_RNA_dm_bc, 
             AUE2021_ITS_rarefied_RNA_dm_metadata$Time)
anova(betadisper_Time_RNA_ITS)

BetaDisp_ITS=
  betadisper_Time_ITS$distances %>% 
  as.data.frame %>%
  set_colnames("Distance") %>%
  mutate(sample=rownames(.), amplicon="ITS")%>% 
  left_join(select(AUE2021_ITS_rarefied_dm_metadata, sample, Time, SiteType, source))

set.seed(72)
ITS_DNA_PCoA=metaMDS(AUE2021_ITS_rarefied_df_DNA, k=10, trymax=100)
set.seed(72)
ITS_RNA_PCoA=metaMDS(AUE2021_ITS_rarefied_df_RNA, k=10, trymax=100)

ITS_DNA_PCoA_Points=
  ITS_DNA_PCoA$points %>%
  as.data.frame %>%
  mutate(sample=rownames(.)) %>%
  left_join(samples_df2) %>%
  mutate(SourceSiteType=paste(SiteType, source),
         SourceSiteTypeTime=paste(SiteType, source, Time))
ITS_RNA_PCoA_Points=
  ITS_RNA_PCoA$points %>%
  as.data.frame %>%
  mutate(sample=rownames(.)) %>%
  left_join(samples_df2) %>%
  mutate(SourceSiteType=paste(SiteType, source),
         SourceSiteTypeTime=paste(SiteType, source, Time))

ggplot(ITS_RNA_PCoA_Points, aes(x=-MDS1, y=MDS2))+
  geom_point(aes(shape=Time, fill=Time), 
             size=4, stroke=1,alpha=0.8, color="black") +
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

ITS_NMDS=
  ggplot(ITS_PCoA_Points, aes(x=-MDS1, y=MDS2))+
  geom_point(aes(shape=Time, fill=SiteType), 
             size=3, stroke=1,alpha=0.8, color="black") +
  geom_point(data=filter(ITS_PCoA_Points, source=="cDNA"), size=0.75) +
  theme_test()+
  scale_fill_brewer(palette = "Dark2") +
  scale_shape_manual(values = c("June" = 21, "August" = 23, "October" = 24)) +
  theme(legend.position="none", text=element_blank()) +
  theme(panel.border = element_rect(linetype = "solid",
                                    colour = "gray20", linewidth = 1.5),
        strip.background = element_blank(),
        axis.ticks=element_line(size=1.5, colour = "gray20"),
        aspect.ratio=1) +
  stat_ellipse(aes(linetype=source, group=source), 
               level=0.95, linewidth=1.5, colour = "gray20")

Troye_palette=c("#c03a2e", "#0397b9", "#5f7b3e")

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
  
ITS_NMDS_open

tiff("C:/Users/akeja/OneDrive - Duke University/Documents/Duke PhD/Projects/PMI/AUE_2021/MetabarcodingManuscript/Misc Figure, Tables, Etc/ITS_NMDS.tiff", 
     width = 13, height = 13, units = "cm", res=4000)
ITS_NMDS_open
dev.off()

pdf(file="C:/Users/akeja/OneDrive - Duke University/Documents/Duke PhD/Projects/PMI/AUE_2021/MetabarcodingManuscript/Misc Figure, Tables, Etc/ITS_NMDS.pdf", 
    width = 5.11811, height = 5.11811)
ITS_NMDS
dev.off()

ITS_PCoA_plot=
  ggplot(ITS_PCoA_Points, aes(x=MDS1, y=MDS2, color=source))+
  geom_point(size=2) +
  theme_test()+
  scale_color_brewer(palette = "Dark2") +
  stat_ellipse(size=1) +
  coord_fixed() +
  theme(legend.position="none",
      panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.background = element_rect(fill='transparent'),
      legend.box.background = element_rect(fill='transparent'),
      panel.border=element_rect(size=2)
    )

AUE2021_ITS_rarefied_rel_phylum %>%
  otu_table %>%
  as.data.frame %>%
  mutate(otu=rownames(.)) %>% 
  left_join(select(tax_mat_df, otu, Phylum)) %>% 
  set_rownames(.$Phylum) %>% 
  select(-otu,-Phylum) %>% 
  t %>%
  as.data.frame %>%
  mutate(sample=rownames(.)) %>%
  left_join(samples_df2) %>%
  mutate(source=dplyr::recode(source, cDNA="RNA")) %>%
  mutate(source=factor(source, levels=c("DNA", "RNA"))) %>%
  ggplot((aes(x=source, y=Glomeromycota, fill=source))) +
  geom_boxplot(size=1) +
  scale_y_continuous(trans='log10') +
  scale_fill_brewer(palette = "Dark2") +
  coord_fixed() +
  theme_test() +
  theme(legend.position="none")

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

betadisper_Time_RNA_16S=
  betadisper(AUE2021_16S_rarefied_RNA_dm_bc, 
             AUE2021_16S_rarefied_RNA_dm_metadata$Time)
anova(betadisper_Time_RNA_16S)

betadisper_Time_DNA_16S=
  betadisper(AUE2021_16S_rarefied_DNA_dm_bc, 
             AUE2021_16S_rarefied_DNA_dm_metadata$Time)
anova(betadisper_Time_DNA_16S)

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
  

NMDS_16S_plot=
  ggplot(PCoA_16S_Points, aes(x=-MDS1, y=MDS2))+
  geom_point(aes(shape=Time, fill=SiteType), 
             size=3, stroke=1,alpha=0.8, color="black") +
  geom_point(data=filter(PCoA_16S_Points, source=="cDNA"), size=0.75) +
  theme_test()+
  scale_fill_brewer(palette = "Dark2") +
  scale_shape_manual(values = c("June" = 21, "August" = 23, "October" = 24)) +
  theme(legend.position="none", text=element_blank()) +
  theme(panel.border = element_rect(linetype = "solid",
                                    colour = "gray20", linewidth = 1.5),
        strip.background = element_blank(),
        axis.ticks=element_line(size=1.5, colour = "gray20"),
        aspect.ratio=1) +
  stat_ellipse(aes(linetype=source, group=source), 
               level=0.95, linewidth=1.5, colour = "gray20")


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

tiff("C:/Users/akeja/OneDrive - Duke University/Documents/Duke PhD/Projects/PMI/AUE_2021/MetabarcodingManuscript/Misc Figure, Tables, Etc/NMDS_16S_plot.tiff", 
     width = 13, height = 13, units = "cm", res=4000)
NMDS_16S_open_plot
dev.off()

tiff("C:/Users/akeja/OneDrive - Duke University/Documents/Duke PhD/Projects/PMI/AUE_2021/MetabarcodingManuscript/Misc Figure, Tables, Etc/NMDS_16S_plot.tiff", 
     width = 13, height = 13, units = "cm", res=4000)
NMDS_16S_plot
dev.off()

pdf(file="C:/Users/akeja/OneDrive - Duke University/Documents/Duke PhD/Projects/PMI/AUE_2021/MetabarcodingManuscript/Misc Figure, Tables, Etc/NMDS_16S_plot.pdf", 
    width = 5.11811, height = 5.11811)
NMDS_16S_plot
dev.off()

BetaDisp_16S=  
  betadisper_Time_16S$distances %>% 
  as.data.frame %>%
  set_colnames("Distance") %>%
  mutate(sample=rownames(.), amplicon="16S")%>% 
  left_join(select(AUE2021_16S_rarefied_dm_metadata, sample, Time, SiteType, source))
  
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
  