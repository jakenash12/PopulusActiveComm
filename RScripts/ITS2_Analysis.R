library(fungarium)
library(lme4)
library(tidyverse)
library(vegan)
library(phyloseq)
library(magrittr)
library(tibble)
library(car)

#reads in base metadata, then joins with comprehensive 
#sample metadata (soils, sensors, etc)
samples_df2 <- read.delim("./InputFiles/Metadata_ITS.txt") %>% 
  left_join((dplyr::select(PRS_Results,-Time, -Site, -Plot, -SiteType)%>%rename_with(~ paste0(., "_PRS"), .cols = -PlotDate)),by="PlotDate") %>% 
  left_join(dplyr::select(Vegsurvey, -Site), by="Plot") %>% 
  left_join((dplyr::select(UGA_Results, -Habitat, -Site)%>%rename_with(~ paste0(., "_UGA"), .cols = -Plot)),by="Plot")  %>%
  left_join(dplyr::select(PlotAvgSoilMoist, Plot, plotavg_moist)) %>%
  left_join(dplyr::select(PlotAvgSoilTemp, Plot, plotavg_soiltemp)) %>%
  mutate(SiteType=dplyr::recode(.$SiteType, "H"="High Elevation", "R"="Riparian", "M"="Mid Elevation"),
         Time=factor(Time, levels=c("June","August","October")),
         TimeSiteType=paste(Time,SiteType, sep="_"),
         TimeSite=paste(Time,Site,sep="_"),
         alias=as.character(alias),
         TimeSite=factor(TimeSite, levels=c(
           "June_HA","June_HB","June_HC","June_MA","June_MB","June_MC","June_RA","June_RB","June_RC",
           "August_HA","August_HB","August_HC","August_MA","August_MB","August_MC","August_RA","August_RB","August_RC",
           "October_HA","October_HB","October_HC","October_MA","October_MB","October_MC","October_RA","October_RB","October_RC")))

#reads in otu table
otu_mat <- read.delim("./InputFiles/OTU_table_ITS.txt")

#creates version of otu table only with samples for which
#there is available soil sensor data
otu_mat_temp_moist <- 
  read.delim("./InputFiles/OTU_table_ITS.txt") %>%
  dplyr::select(filter(samples_df2, !is.na(plotavg_soiltemp))$sample,otu)

#second copy of the sample metadata that will be formatted for phyloseq
samples_df <- read.delim("./InputFiles/Metadata_ITS.txt") %>% 
  left_join((dplyr::select(PRS_Results,-Time, -Site, -Plot, -SiteType)%>%rename_with(~ paste0(., "_PRS"), .cols = -PlotDate)),by="PlotDate") %>% 
  left_join(dplyr::select(Vegsurvey, -Site), by="Plot") %>% 
  left_join((dplyr::select(UGA_Results, -Habitat, -Site)%>%rename_with(~ paste0(., "_UGA"), .cols = -Plot)),by="Plot")  %>%
  left_join(dplyr::select(PlotAvgSoilMoist, Plot, plotavg_moist)) %>%
  left_join(dplyr::select(PlotAvgSoilTemp, Plot, plotavg_soiltemp)) %>%
  mutate(SiteType=dplyr::recode(.$SiteType, "H"="High Elevation", "R"="Riparian", "M"="Mid Elevation"),
         Time=factor(Time, levels=c("June","August","October")),
         TimeSiteType=paste(Time,SiteType, sep="_"),
         TimeSite=paste(Time,Site,sep="_"),
         alias=as.character(alias),
         TimeSite=factor(TimeSite, levels=c(
           "June_HA","June_HB","June_HC","June_MA","June_MB","June_MC","June_RA","June_RB","June_RC",
           "August_HA","August_HB","August_HC","August_MA","August_MB","August_MC","August_RA","August_RB","August_RC",
           "October_HA","October_HB","October_HC","October_MA","October_MB","October_MC","October_RA","October_RB","October_RC")))

#reads in UNITE taxonomy from QIIME2
tax_mat <- read.delim("./InputFiles/taxonomy_UNITE.txt")
tax_mat_df <- read.delim("./InputFiles/taxonomy_UNITE.txt") %>%
  mutate(TaxClass=ifelse(!is.na(Species),Species,ifelse(!is.na(Genus),Genus,ifelse(!is.na(Family),Family,ifelse(!is.na(Order),Order,ifelse(!is.na(Class),Class,ifelse(!is.na(Phylum),Phylum,Domain)))))))

#reads in custom list of dark septate endophyte taxa
DSE_List=read.delim("./InputFiles/DSE_List.txt")

#runs funguild on taxonomy matrix
tax_mat_fg=
  tax_mat %>%
  as.data.frame %>% 
  fungarium::fg_assign(.,tax_cols=c("Domain", "Phylum", "Class", "Order", "Family", "Genus")) %>%
  as.data.frame


tax_mat_fg_edited=
  tax_mat_fg %>%
  mutate(Genus=dplyr::recode(.$Genus, "Pithomyces"="Aquilomyces")) %>%
  mutate(primary_guild = str_extract(guild, "\\|([^|]+)\\|") %>%  # Extract text between two |
           str_replace_all("\\|", "") %>%                         # Remove the remaining |
           str_trim()) %>%                                        # Trim leading and trailing spaces
  mutate(guild_edited=primary_guild,
         guild_edited = ifelse(str_detect(guild, "Saprotroph"), "Sap", guild_edited)) %>%
  mutate(guild_edited=ifelse(.$Genus %in% DSE_List$Genus, "DSE", .$guild_edited)) %>%
  mutate(guild_edited=case_when(
           Phylum=="Glomeromycota" ~ "Arbuscular Mycorrhizal",
           .default=guild_edited))
  

otu_mat_temp_moist %<>% column_to_rownames("otu")
otu_mat %<>% column_to_rownames("otu")
tax_mat %<>% column_to_rownames("otu")
samples_df %<>% column_to_rownames("sample")

otu_mat_temp_moist <-as.matrix(otu_mat_temp_moist )
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

OTU_moist = otu_table(otu_mat_temp_moist, taxa_are_rows = TRUE )
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

AUE2021_ITS <- phyloseq(OTU, TAX, samples)
AUE2021_ITS_moist <- phyloseq(OTU_moist, TAX, samples)
seq_depth <-
  sample_sums(AUE2021_ITS) %>%
  as.data.frame %>%
  rename(SeqCount=".") %>%
  mutate(sample=rownames(.)) %>%
  left_join(read_excel(PhyloPath, sheet = "MappingFile"))

ITS_rarecurve=
  rarecurve(t(otu_table(AUE2021_ITS)), step=100, cex=0.5,label=FALSE,tidy=TRUE)

rare <- lapply(ITS_rarecurve, function(x){
  b <- as.data.frame(x)
  b <- data.frame(OTU = b[,1], raw.read = rownames(b))
  b$raw.read <- as.numeric(gsub("N", "",  b$raw.read))
  return(b)
})
rare <- map_dfr(rare, function(x){
  z <- data.frame(x)
  return(z)
}, .id = "sample")

rarecurve_ITS_plot=
  ggplot(rare, aes(x=raw.read, y=OTU, group=sample)) +
  geom_line() +
  geom_vline(xintercept=17182, color="red") +
  xlab("Number of Reads") +
  ylab("Number of OTUs") +
  ggtitle("ITS")

plot_grid(rarecurve_ITS_plot,rarecurve_16S_plot)

AUE2021_ITS_rarefied = rarefy_even_depth(AUE2021_ITS,17083, rngseed=TRUE)
AUE2021_ITS_rarefied_rel = transform_sample_counts(AUE2021_ITS_rarefied, function(x) x / sum(x) )

AUE2021_ITS_moist_rarefied = rarefy_even_depth(AUE2021_ITS_moist,17083, rngseed=TRUE)
AUE2021_ITS__moistrarefied_rel = transform_sample_counts(AUE2021_ITS_moist_rarefied, function(x) x / sum(x) )

AUE2021_alphadiv =
  estimate_richness(AUE2021_ITS_rarefied) %>%
  mutate(sample=rownames(.)) %>%
  left_join(samples_df2, by="sample") %>%
  mutate(Site=factor(Site, levels=c("RA","RB","RC","MA","MB","MC","HA","HB","HC")),
         Time=factor(Time, levels=c("June","August","October")),
         SiteType=factor(.$SiteType, levels=c("High Elevation","Mid Elevation","Riparian")),
         source=dplyr::recode(source, "cDNA"="RNA")) %>%
  mutate(source=factor(.$source, levels=c("DNA", "RNA")))




AUE2021_alphadiv_DNA=
  AUE2021_alphadiv %>%
  filter(source=="DNA")

AUE2021_alphadiv_RNA=
  AUE2021_alphadiv %>%
  filter(source=="RNA")

mod=lmer(Shannon~SiteType*Time+(1|Plot), AUE2021_alphadiv_RNA)
Anova(mod)
summary(mod)

mod=lmer(Shannon~SiteType*Time+(1|Plot), AUE2021_alphadiv_DNA)
Anova(mod)
summary(mod)

AUE2021_alphadiv_summary=
  AUE2021_alphadiv %>%
  group_by(Time, source, SiteType) %>%
  summarise(mean_Shannon=mean(Shannon),
            se_Shannon=sd(Shannon)/sqrt(n()),
            mean_InvSimpson=mean(InvSimpson),
            se_InvSimpson=sd(InvSimpson)/sqrt(n()),
            mean_Richness=mean(Observed),
            se_Richness=sd(Observed)/sqrt(n()))

Dodge=0.4

#color palette to match the slightly transparent points w/ dark2 palette
#used in the ordination
Alphadivpalette=c("High Elevation" = "#72b599", 
                  "Mid Elevation" = "#e48c5a",
                  "Riparian"="#9791c5")

AlphaDivPlot_ITS=
  ggplot(AUE2021_alphadiv_summary, aes(x=Time, y=mean_Shannon, color=SiteType, group=SiteType)) +
  geom_line(size=1.5,position=position_dodge(width=Dodge)) +
  geom_errorbar(aes(ymin=mean_Shannon-se_Shannon, ymax=mean_Shannon+se_Shannon), 
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
  scale_y_continuous(breaks=c(2.6, 3.1, 3.6)) +
  scale_color_manual(values=Alphadivpalette) +
  scale_fill_manual(values=Alphadivpalette) +
  facet_wrap(.~source)

RichnessPlot_ITS=
  ggplot(AUE2021_alphadiv_summary, aes(x=Time, y=mean_Richness, color=SiteType, group=SiteType)) +
  geom_line(size=1.5,position=position_dodge(width=Dodge)) +
  geom_errorbar(aes(ymin=mean_Richness-se_Richness, ymax=mean_Richness+se_Richness), 
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
  #scale_y_continuous(breaks=c(2.6, 3.1, 3.6)) +
  scale_color_manual(values=Alphadivpalette) +
  scale_fill_manual(values=Alphadivpalette) +
  facet_wrap(.~source)

InvSimpsonPlot_ITS=
  ggplot(AUE2021_alphadiv_summary, aes(x=Time, y=mean_InvSimpson, color=SiteType, group=SiteType)) +
  geom_line(size=1.5,position=position_dodge(width=Dodge)) +
  geom_errorbar(aes(ymin=mean_InvSimpson-se_InvSimpson, ymax=mean_InvSimpson+se_InvSimpson), 
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
  #scale_y_continuous(breaks=c(2.6, 3.1, 3.6)) +
  scale_color_manual(values=Alphadivpalette) +
  scale_fill_manual(values=Alphadivpalette) +
  facet_wrap(.~source)

tiff("C:/Users/akeja/OneDrive - Duke University/Documents/Duke PhD/Projects/PMI/AUE_2021/MetabarcodingManuscript/Misc Figure, Tables, Etc/AlphaDivPlot_ITS.tiff", 
     width = 13, height = 11, units = "cm", res=4000)
AlphaDivPlot_ITS
dev.off()

pdf(file="C:/Users/akeja/OneDrive - Duke University/Documents/Duke PhD/Projects/PMI/AUE_2021/MetabarcodingManuscript/Misc Figure, Tables, Etc/AlphaDivPlot_ITS.pdf", 
     width = 5.11811, height = 4.33071)
AlphaDivPlot_ITS
dev.off()


AlphaDivITS_Plot=
  ggplot(AUE2021_alphadiv, aes(x=Time, y=Shannon, color=SiteType)) +
  geom_boxplot(size=1) +
  theme_test() +
  #theme(legend.position="none") +
  scale_color_brewer(palette = "Dark2") +
  facet_wrap(.~source)


AlphaDivITS_DNA_Plot=
  ggplot((AUE2021_alphadiv_DNA %>%mutate(SiteType=factor(SiteType,levels=c("High Elevation","Mid Elevation","Riparian")))), 
         aes(x=Time, y=Shannon, color=SiteType)) +
  geom_boxplot(size=1) +
  theme_light() +
  theme(legend.position="none") +
  scale_color_brewer(palette = "Dark2") 

AlphaDivITS_RNA_Plot=
  ggplot((AUE2021_alphadiv_RNA %>%mutate(SiteType=factor(SiteType,levels=c("High Elevation","Mid Elevation","Riparian")))), 
         aes(x=Time, y=Shannon, color=SiteType)) +
  geom_boxplot(size=1) +
  theme_light() +
  theme(legend.position="none") +
  scale_color_brewer(palette = "Dark2") 

plot_grid(AlphaDivITS_DNA_Plot,AlphaDivITS_RNA_Plot,
          AlphaDiv16S_DNA_Plot, AlphaDiv16S_RNA_Plot)

AirtempITS_plot=ggplot(AUE2021_alphadiv_summary, aes(x=mean_airtemp,y=mean_Observed)) +
  geom_point(aes(color=SiteType), size=4) +
  #geom_errorbar(aes(ymin=mean_Observed-se_Observed, ymax=mean_Observed+se_Observed, color=SiteType), 
   #             size=1, width=0) +
  ylab("# of OTUs per sample") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  #theme(legend.position="none") +
  scale_color_brewer(palette = "Dark2") +
  geom_smooth(method="lm",se=F, color="gray40") +
  theme(aspect.ratio=1) +
  theme(axis.title=element_blank())

SoiltempITS_plot=ggplot(AUE2021_alphadiv_summary, aes(x=mean_soiltemp,y=mean_Observed)) +
  geom_point(aes(color=SiteType), size=4) +
  geom_errorbar(aes(ymin=mean_Observed-se_Observed, ymax=mean_Observed+se_Observed, color=SiteType), 
                size=1, width=0) +
  ylab("# of OTUs per sample") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="none") +
  scale_color_brewer(palette = "Dark2") +
  geom_smooth(method="lm",se=F, color="gray40") +
  theme(aspect.ratio=1) +
  theme(axis.title=element_blank())

ggplot(AUE2021_alphadiv_summary, aes(x=SiteType,y=mean_Shannon)) +
  geom_boxplot(aes(color=SiteType), size=1) +
  ylab("Shannon Diversity") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="none") +
  #scale_color_brewer(palette = "Dark2") +
  geom_smooth(method="lm",se=F, color="gray40") +
  theme(aspect.ratio=1)

plot_grid(AirtempITS_plot, Airtemp16S_plot, SoiltempITS_plot, Soiltemp16S_plot)

mod=lm(mean_Observed~mean_airtemp, AUE2021_alphadiv_summary)  
Anova(mod)
summary(mod)

ggplot(AUE2021_alphadiv, aes(x=Time,y=Shannon)) +
  geom_boxplot() +
  ylab("# OTUs")+
  facet_grid(.~source)

mod=lmer(Shannon~SiteType*Time*source +(1|Site), AUE2021_alphadiv)
Anova(mod, type=2)
summary(mod)


AUE2021_ITS_rarefied_rel_species=
  AUE2021_ITS_rarefied_rel %>%
  tax_glom("Species")
AUE2021_ITS_rarefied_rel_genus=
  AUE2021_ITS_rarefied_rel %>%
  tax_glom("Genus")
AUE2021_ITS_rarefied_rel_family=
  AUE2021_ITS_rarefied_rel %>%
  tax_glom("Family")
AUE2021_ITS_rarefied_rel_order=
  AUE2021_ITS_rarefied_rel %>%
  tax_glom("Order")
AUE2021_ITS_rarefied_rel_phylum=
  AUE2021_ITS_rarefied_rel %>%
  tax_glom("Phylum")

AUE2021_ITS_rarefied_rel_phylum_df=
  AUE2021_ITS_rarefied_rel_phylum %>%
  otu_table %>%
  as.data.frame %>%
  mutate(otu=rownames(.)) %>%
  left_join(select(tax_mat_df, otu,Phylum)) %>%
  filter(Phylum=="Glomeromycota") %>%
  select(-otu, -Phylum) %>% 
  set_rownames("Glomeromycota") %>%
  t %>%
  as.data.frame %>% 
  mutate(sample=rownames(.)) %>%
  left_join(samples_df2) 

mod=lmer(Glomeromycota~SiteType*Time+source +(1|Site), AUE2021_ITS_rarefied_rel_phylum_df)
Anova(mod)
summary(mod)

glomplot_pal=c("#8d4653","#265DAB", "#999999")

ggplot(AUE2021_ITS_rarefied_rel_phylum_df,
       aes(x=source, y=(Glomeromycota+0.0001)*100, color=source)) +
  geom_boxplot(size=1) +
  scale_y_continuous(trans='log10', 
                     breaks = c(0, 0.1, 1, 0.01)) +
  scale_colour_manual(values=glomplot_pal) +
  theme_classic() +
  theme(aspect.ratio=1) +
  theme(legend.position="none", axis.title.y=element_blank())

Anova(mod)
AUE2021_OTU_ord <- ordinate(AUE2021_ITS_rarefied, "PCoA", "bray")
AUE2021_Genus_ord <- ordinate(AUE2021_ITS_rarefied_rel_genus, "PCoA", "bray")
AUE2021_family_ord <- ordinate(AUE2021_ITS_rarefied_rel_family, "PCoA", "bray")
AUE2021_order_ord <- ordinate(AUE2021_ITS_rarefied_rel_order, "PCoA", "bray")

cbPalette <- c("#E69F00",  "#0072B2", "#D55E00", "#CC79A7")
plot_ordination(AUE2021_ITS_rarefied,AUE2021_OTU_ord) +
  geom_point(aes(color=source), size=3) +
  theme_test() +
  theme(aspect.ratio = 1) +
  facet_grid(.~source)
  scale_color_manual(values=cbPalette)

SampleDf_for_Ord=
  samples_df2 %>%
  filter(sample %in% sample_names(AUE2021_ITS_rarefied)) %>%
  select(source:TotalCover,avg_airtemp)

samples_df_betadisp=
  samples_df2 %>%
  filter(sample%in%row.names(sample_data(AUE2021_ITS_rarefied_rel))) %>%
  mutate(Treatment=paste(.$source,.$Time,.$Site,sep=""))

Site_betadisp=
  betadisper(bray_dist, samples_df_betadisp$Time)


samples_df_betadisp_site_level=
  samples_df_betadisp %>%
  group_by(Site,Time,source) %>%
  summarise(avg_temp=mean(avg_temp, na.rm=T),
            avg_moist=mean(avg_moist,na.rm=T),
            Treatment=Treatment,
            SiteType=SiteType,
            avg_airtemp=mean(avg_airtemp, na.rm=T))

Site_betadisp_df=
  read_excel("C:/Users/akeja/OneDrive - Duke University/Documents/Duke PhD/Projects/PMI/AUE_2021/AmpliconSeqFiles/ITS2/ITS2_BetaDispersion.xlsx") %>%
  left_join(samples_df_betadisp_site_level)

AirTempBD_plot=
  ggplot(Site_betadisp_df, aes(x=avg_airtemp, y=BetaDisp)) +
  geom_point(size=3, aes(color=SiteType)) +
  facet_grid(source~Time, scales="free_y") +
  ylab("Beta dispersion within site") +
  xlab("Average Air Temperature")+
  theme(legend.position="none")

SoilTempBD_plot=
  ggplot(Site_betadisp_df, aes(x=avg_temp, y=BetaDisp)) +
  geom_point(size=3, aes(color=SiteType)) +
  #facet_grid(source~Time, scales="free_y") +
  ylab("Beta dispersion within site") +
  xlab("Average Soil Temperature")+
  theme(legend.position="none")

SoilMoistBD_plot=
  ggplot(Site_betadisp_df, aes(x=avg_moist, y=BetaDisp)) +
  geom_point(size=3, aes(color=SiteType)) +
  facet_grid(source~Time, scales="free_y") +
  ylab("Beta dispersion within site") +
  xlab("Average Soil Moisture") +
  theme(legend.position="none")

ggplot(Site_betadisp_df, aes(x=Time, y=BetaDisp)) +
  geom_boxplot()  +
  facet_grid(source~., scales="free_y") +
  ylab("Beta dispersion within site") +
  xlab("Average Soil Temperature")+
  theme(legend.position="none")

plot_grid(AirTempBD_plot,SoilTempBD_plot,SoilMoistBD_plot,ncol=1)

mod=lm(BetaDisp~avg_temp+Time+source,Site_betadisp_df)
Anova(mod)
summary(mod)

Site_betadisp$call

samples_df2_filtered=
  samples_df2 %>%
  filter(sample %in% sample_names(AUE2021_ITS_rarefied))
Metadata_dist_mat=
  samples_df2_filtered %>%
  select(NO3:avg_airtemp, -c(Habitat, Plotnum, Cr, Cu.y, Mo, Pb.y)) %>%
  set_rownames(samples_df2_filtered$sample) %>%
  vegdist(method="bray")

AUE2021_ITS_mds <- wcmdscale(bray_dist)

AUE2021_ITS_mds_points=
  AUE2021_ITS_mds %>% 
  as.data.frame %>%
  mutate(sample=rownames(.)) %>%
  left_join(samples_df2)


AUE2021_OTU_ord_df =
  AUE2021_OTU_ord$vectors

AUE2021_OTU_ord_genus_envfit=
  envfit(AUE2021_ITS_mds, AUE2021_ITS_rarefied_rel_genus_df)

AUE2021_OTU_ord_genus_envfit_df=
  AUE2021_OTU_ord_genus_envfit %>%
  scores(display="vector") %>%
  as.data.frame %>%
  mutate(R2=AUE2021_OTU_ord_genus_envfit$vectors$r,
         p=AUE2021_OTU_ord_genus_envfit$vectors$pvals,
         Genus=rownames(.))

AUE2021_OTU_ord_genus_envfit_df_filtered=
  AUE2021_OTU_ord_genus_envfit_df %>%
  filter(R2>0.1&p<0.01&Genus!="Chaetothyriales_unidentified_1") %>%
  mutate(Dim1_scale=Dim1*0.5, Dim2_scale=Dim2*0.5)

AUE2021_ITS_mds_points_Format=
AUE2021_ITS_mds_points %>%
  group_by(Plot,Time) %>%
  summarise(V1=mean(V1),
            V2=mean(V2),
            avg_airtemp=mean(avg_airtemp)) %>%
  mutate(SiteType=as.factor(str_extract(Plot, "^."))) %>%
  ungroup

ggplot(AUE2021_ITS_mds_points, aes(x=V1, y=V2)) +
  geom_point(size=4, aes(color=source)) +
  #stat_ellipse(level=0.8) +
  xlab("PCoA 1") +
  ylab("PCoA 2") +
  #scale_color_gradient(low = 'blue', high = 'red') +
  theme_test() +
  geom_segment(data = AUE2021_OTU_ord_genus_envfit_df_filtered,
               aes(x = 0, xend = Dim1_scale, y = 0, yend = Dim2_scale),
               arrow = arrow(length = unit(0.2, "cm"),type="closed"), colour = "black",
               size=1) +
  geom_text(data = AUE2021_OTU_ord_genus_envfit_df_filtered, aes(x = Dim1_scale, y = Dim2_scale, label = Genus),
            size = 3) + 
  theme(aspect.ratio=1) +
  #theme(legend.position="None") +
  xlim(-.4, 0.23) +
  ylim(-.37, 0.31) +
  scale_color_brewer(palette = "Dark2")

AUE2021_ITS_rarefied_rel_genus_df =
  AUE2021_ITS_rarefied_rel_genus %>%
  otu_table %>%
  as.data.frame %>%
  mutate(otu=rownames(.)) %>%
  left_join(select(tax_mat_df, otu,Genus)) %>%
  set_rownames(.$Genus) %>%
  select(-otu,-Genus) %>%
  t %>%
  as.data.frame

AUE2021_ITS_rarefied_rel_genus_df2 =
  AUE2021_ITS_rarefied_rel_genus %>%
  otu_table %>%
  as.data.frame %>%
  mutate(otu=rownames(.)) %>%
  left_join(select(tax_mat_df, otu,Genus)) %>%
  set_rownames(.$Genus) %>%
  select(-otu,-Genus) %>%
  t %>%
  as.data.frame %>%
  mutate(sample=rownames(.)) %>%
  left_join(filter(samples_df2, sample %in% .$sample))

AUE2021_ITS_rarefied_rel_genus_df2_sitesummary=
  AUE2021_ITS_rarefied_rel_genus_df2 %>%
  group_by(Site.x) %>% 
  summarise(Aquilomyces=mean(Pithomyces),
            Inocybe=mean(Inocybe),
            Ilyonectria=mean(Ilyonectria),
            Mortierella=mean(Mortierella),
            mean_airtemp=mean(avg_airtemp)) %>%
  left_join(distinct(samples_df2, Site.x, .keep_all=T), by="Site.x")

ggplot(AUE2021_ITS_rarefied_rel_genus_df2_sitesummary, aes(x=mean_airtemp, y =Aquilomyces*100)) +
  geom_point(aes(color=SiteType), size=4) +
  ylab("Aquilomyces % Abundance") +
  xlab("Mean Site Air Temperature (C)") +
  theme_classic() +
  theme(legend.position="none") +
  scale_color_brewer(palette = "Dark2") +
  geom_smooth(method="lm",se=F, color="gray40") +
  theme(aspect.ratio=1)

mod=lm(Pithomyces~mean_airtemp, AUE2021_ITS_rarefied_rel_genus_df2_sitesummary)
Anova(mod)
summary(mod)
AUE2021_OTU_ord_genus_envfit$vectors$arrows %>% View

mod=lm(V1~avg_airtemp, AUE2021_ITS_mds_points)

plot_ordination(AUE2021_ITS_rarefied, AUE2021_OTU_ord, type="sample",
                color="SiteType") +
  geom_point(size=3) +
  stat_ellipse(aes(group=SiteType))+
  facet_grid(source~factor(Time,levels=c("June","August","October")))

plot_ordination(AUE2021_ITS_rarefied_rel_genus, AUE2021_Genus_ord, type="sample",
                color="SiteType") +
  geom_text(size=3) +
  facet_grid(source~Time)

plot_ordination(AUE2021_ITS_rarefied_rel_family, AUE2021_family_ord, type="sample",
                color="SiteType") +
  geom_point(size=3) +
  facet_grid(source~Time)

plot_ordination(AUE2021_ITS_rarefied_rel_order, AUE2021_order_ord, type="sample",
                color="SiteType") +
  geom_point(size=3) +
  facet_grid(source~Time)

bray_dist_moist = phyloseq::distance(AUE2021_ITS__moistrarefied_rel, method="bray", weighted=T)


adonis(bray_dist_moist ~sample_data(AUE2021_ITS__moistrarefied_rel)$avg_airtemp)


bray_dist = phyloseq::distance(AUE2021_ITS_rarefied, method="bray", weighted=T)
adonis2(bray_dist ~ sample_data(AUE2021_ITS_rarefied)$source)

bray_dist_genus = phyloseq::distance(AUE2021_ITS_rarefied_rel_genus, method="bray", weighted=T)
adonis2(bray_dist_genus ~ sample_data(AUE2021_ITS_rarefied_rel_genus)$Site+sample_data(AUE2021_ITS_rarefied_rel_genus)$Time+sample_data(AUE2021_ITS_rarefied_rel_genus)$source)

bray_dist_family = phyloseq::distance(AUE2021_ITS_rarefied_rel_family, method="bray", weighted=T)
adonis2(bray_dist_family ~ sample_data(AUE2021_ITS_rarefied_rel_family)$Site+sample_data(AUE2021_ITS_rarefied_rel_family)$Time+sample_data(AUE2021_ITS_rarefied_rel_family)$source)

bray_dist_order = phyloseq::distance(AUE2021_ITS_rarefied_rel_order, method="bray", weighted=T)
adonis2(bray_dist_order ~ sample_data(AUE2021_ITS_rarefied_rel_order)$Plot+sample_data(AUE2021_ITS_rarefied_rel_order)$Time+sample_data(AUE2021_ITS_rarefied_rel_order)$source)
