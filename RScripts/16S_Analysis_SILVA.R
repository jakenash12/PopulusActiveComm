library(lme4)
library(tidyverse)
library(vegan)
library(phyloseq)
library(magrittr)
library(tibble)
library(car)

#reads in SILVA taxonomy and formats it by splitting the single
#taxonomy column into separate columns for each rank
silva_tax =
  read.delim("./InputFiles/taxonomy_SILVA-138.2.tsv") %>%
  separate(Taxon,
           into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";\\s*",
           fill = "right") %>%
  mutate(across(c(Domain, Phylum, Class, Order, Family, Genus, Species),
                ~ str_remove(., "^[a-z]__"))) %>%
  dplyr::select(-Confidence) %>%
  column_to_rownames("Feature.ID")

#generates list of mitochondrial and chloroplast OTUs
#based on SILVA annotations
mito_chloro_otus_silva=
  silva_tax %>%
  filter((Family=="Mitochondria" | Order == "Chloroplast")) %>% 
  rownames(.)

#reads in base metadata, then joins with comprehensive 
#sample metadata (soils, sensors, etc)
samples_df2_16S <- read.delim("./InputFiles/Metadata_16S.txt") %>% 
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

#reads in OTU table
otu_mat_16S <- read.delim("./InputFiles/OTU_table_16S.txt") 

#creates version of otu table only with samples for which
#there is available soil sensor data
otu_mat_temp_moist_16S <- 
  read.delim("./InputFiles/OTU_table_16S.txt") %>%
  dplyr::select(filter(samples_df2_16S, !is.na(plotavg_soiltemp))$sample,otu)

#second copy of the sample metadata that will be formatted for phyloseq
samples_df_16S <- read.delim("./InputFiles/Metadata_16S.txt") %>% 
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

#filters out mitochondrial and chloroplast ASVs
#from otu tables and taxonomy tables
otu_mat_temp_moist_16S %<>% 
  filter(!(otu %in% mito_chloro_otus_silva)) %>%
  column_to_rownames("otu")
otu_mat_16S %<>% 
  filter(!(otu %in% mito_chloro_otus_silva)) %>%
  column_to_rownames("otu")
silva_tax %<>% 
  filter(!(rownames(.) %in% mito_chloro_otus_silva))
samples_df_16S %<>% column_to_rownames("sample")

#converts otu tables and tax table into matrices for phyloseq
otu_mat_temp_moist_16S <-as.matrix(otu_mat_temp_moist_16S )
otu_mat_16S <- as.matrix(otu_mat_16S)
silva_tax <- as.matrix(silva_tax)

#creates phyloseq objects from otu tables, tax table, and sample meteadata
OTU_moist_16S = otu_table(otu_mat_temp_moist_16S, taxa_are_rows = TRUE )
OTU_16S = otu_table(otu_mat_16S, taxa_are_rows = TRUE)
TAX_SILVA_16S = tax_table(silva_tax)
samples_16S = sample_data(samples_df_16S)

#creates full phyloseq project objects
AUE2021_16S <- phyloseq(OTU_16S, TAX_SILVA_16S, samples_16S)
AUE2021_16S_moist <- phyloseq(OTU_moist_16S, TAX_SILVA_16S, samples_16S)

#calculates sequencing depth per sample
seq_depth_16S <-
  sample_sums(AUE2021_16S) %>%
  as.data.frame %>%
  rename(SeqCount=".") %>%
  mutate(sample=rownames(.)) %>%
  left_join(samples_df2_16S)

#creates otu table in format for rarefaction curve
AUE2021_16S_mat_t=
  AUE2021_16S %>%
  otu_table() %>%
  `@`(., ".Data") %>% 
  t()

#creates rarefaction curve dataframe
rarecurve_16S=
  rarecurve(AUE2021_16S_mat_t, step=100, cex=0.5,label=FALSE,tidy=TRUE)

#plots rarefaction curve with chosen rarefaction
#depth of 13407
ggplot(rarecurve_16S, aes(x=Sample, y=Species, group=Site)) +
  geom_line() +
  geom_vline(xintercept=13407, color="red") +
  xlab("Number of Reads") +
  ylab("Number of OTUs") +
  ggtitle("16S")

#rarefies otu table to 13407 using random seed to ensure
#reproducibility
AUE2021_16S_rarefied = rarefy_even_depth(AUE2021_16S, 13407, rngseed=TRUE)

#creates relativized rarefied dataframe
AUE2021_16S_rarefied_rel = transform_sample_counts(AUE2021_16S_rarefied, function(x) x / sum(x) )

#creates summary otu tables summed at different taxonomic
#levels
AUE2021_16S_rarefied_rel_genus=
  AUE2021_16S_rarefied_rel %>%
  tax_glom("Genus")
AUE2021_16S_rarefied_rel_family=
  AUE2021_16S_rarefied_rel %>%
  tax_glom("Family")
AUE2021_16S_rarefied_rel_order=
  AUE2021_16S_rarefied_rel %>%
  tax_glom("Order")
AUE2021_16S_rarefied_rel_phylum=
  AUE2021_16S_rarefied_rel %>%
  tax_glom("Phylum")

#generates alpha diversity dataframe 
AUE2021_16S_alphadiv =
  estimate_richness(AUE2021_16S_rarefied) %>%
  mutate(sample=rownames(.)) %>%
  left_join(samples_df2_16S, by="sample") %>%
  mutate(Site=factor(Site, levels=c("RA","RB","RC","MA","MB","MC","HA","HB","HC")),
         Time=factor(Time, levels=c("June","August","October")),
         SiteType=factor(.$SiteType, levels=c("High Elevation","Mid Elevation","Riparian")),
         source=dplyr::recode(source, "cDNA"="RNA")) %>%
  mutate(source=factor(.$source, levels=c("DNA", "RNA")))

#calculates mean and SE of Shannon, Inv Simpson, and Richness
AUE2021_16S_alphadiv_summary=
  AUE2021_16S_alphadiv %>%
  group_by(Time, source, SiteType) %>%
  summarise(mean_Shannon=mean(Shannon),
            se_Shannon=sd(Shannon)/sqrt(n()),
            mean_InvSimpson=mean(InvSimpson),
            se_InvSimpson=sd(InvSimpson)/sqrt(n()),
            mean_Richness=mean(Observed),
            se_Richness=sd(Observed)/sqrt(n()))

#sets dodge for plot
Dodge=0.4

#color palette to match the slightly transparent points w/ dark2 palette
#used in the ordination
Alphadivpalette=c("High Elevation" = "#72b599", 
                  "Mid Elevation" = "#e48c5a",
                  "Riparian"="#9791c5")

#creates plot of Shannon Div by season and site type
#plotted without labels so that they can be manually
#added in Adobe Illustrator
AlphaDivPlot_16S=
  ggplot(AUE2021_16S_alphadiv_summary, aes(x=Time, y=mean_Shannon, color=SiteType, group=SiteType)) +
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
  scale_y_continuous(breaks=c(5.4, 5.8, 6.2)) +
  facet_wrap(.~source) +
  scale_color_manual(values=Alphadivpalette) +
  scale_fill_manual(values=Alphadivpalette)

tiff("./OutputFiles/AlphaDivPlot_16S_SILVA.tiff", 
     width = 13, height = 11, units = "cm", res=4000)
AlphaDivPlot_16S
dev.off()

pdf(file="./OutputFiles/AlphaDivPlot_16S_SILVA.pdf", 
    width = 5.11811, height = 4.33071)
AlphaDivPlot_16S
dev.off()

#generates DNA-only alpha div dataframe
AUE2021_16S_alphadiv_DNA=
  AUE2021_16S_alphadiv %>%
  filter(source=="DNA")

#generates RNA-only alpha div dataframe
AUE2021_16S_alphadiv_RNA=
  AUE2021_16S_alphadiv %>%
  filter(source=="cDNA")

#runs global mixed model for effect of site type, season, and source on shannon div
mod=lmer(Shannon~SiteType*Time*source +(1|Site), AUE2021_16S_alphadiv)
Anova(mod)
summary(mod)

#runs mixed model for effect of site type and season on shannon div in DNA dataset
mod=lmer(Shannon~SiteType*Time +(1|Site), AUE2021_16S_alphadiv_DNA)
Anova(mod)
summary(mod)

#runs mixed model for effect of site type and season on shannon div in RNA dataset
mod=lmer(Shannon~SiteType*Time +(1|Site), AUE2021_16S_alphadiv_RNA)
Anova(mod)
summary(mod)

#writes rarefied OTU table to file
AUE2021_16S_rarefied %>%
  otu_table %>%
  as.data.frame %>% 
  write.table("./OutputFiles/16S_table_rarefied.tsv", quote=FALSE, sep='/t', col.names = NA)

