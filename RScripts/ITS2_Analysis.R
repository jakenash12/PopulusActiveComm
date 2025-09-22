#script to perform data import for fungal ITS data
#Includes initial data formatting (conversion to 
#phyloseq format) and alpha diversity analyses

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
  mutate(across(everything(), ~ na_if(str_trim(.), ""))) %>%
  mutate(TaxClass=ifelse(!is.na(Species),Species,ifelse(!is.na(Genus),Genus,ifelse(!is.na(Family),Family,ifelse(!is.na(Order),Order,ifelse(!is.na(Class),Class,ifelse(!is.na(Phylum),Phylum,Domain)))))))

#reads in custom list of dark septate endophyte taxa
DSE_List=read.delim("./InputFiles/DSE_List.txt")

#runs funguild on taxonomy matrix
tax_mat_fg=
  tax_mat %>%
  as.data.frame %>% 
  fungarium::fg_assign(.,tax_cols=c("Domain", "Phylum", "Class", "Order", "Family", "Genus")) %>%
  as.data.frame

#edits guild assignments
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
  
#sets rownames
otu_mat_temp_moist %<>% column_to_rownames("otu")
otu_mat %<>% column_to_rownames("otu")
tax_mat %<>% column_to_rownames("otu")
samples_df %<>% column_to_rownames("sample")

#converts to matrix
otu_mat_temp_moist <-as.matrix(otu_mat_temp_moist )
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

#converts to phyloseq format objects
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
  left_join(read.delim("./InputFiles/Metadata_ITS.txt"))


#creates otu table in format for rarefaction curve
AUE2021_ITS_mat_t=
  AUE2021_ITS %>%
  otu_table() %>%
  `@`(., ".Data") %>% 
  t()

#creates a rarefaction curve
ITS_rarecurve=
  rarecurve(AUE2021_ITS_mat_t, step=100, cex=0.5,label=FALSE,tidy=TRUE)

#plots rarefaction curve with chosen rarefaction
#depth of 13407
ggplot(ITS_rarecurve, aes(x=Sample, y=Species, group=Site)) +
  geom_line() +
  geom_vline(xintercept=17083, color="red") +
  xlab("Number of Reads") +
  ylab("Number of OTUs") +
  ggtitle("ITS")

#conducts rarefaction at 17083 seqs
AUE2021_ITS_rarefied = rarefy_even_depth(AUE2021_ITS,17083, rngseed=TRUE)

#writes rarefied OTU table to file
AUE2021_ITS_rarefied %>%
  otu_table %>%
  as.data.frame %>% 
  write.table("./OutputFiles/ITS_table_rarefied.tsv", quote=FALSE, sep='/t', col.names = NA)

#creates relativized otu table
AUE2021_ITS_rarefied_rel = transform_sample_counts(AUE2021_ITS_rarefied, function(x) x / sum(x) )

AUE2021_ITS_moist_rarefied = rarefy_even_depth(AUE2021_ITS_moist,17083, rngseed=TRUE)
AUE2021_ITS__moistrarefied_rel = transform_sample_counts(AUE2021_ITS_moist_rarefied, function(x) x / sum(x) )

#calculates alpha diversity
AUE2021_alphadiv =
  estimate_richness(AUE2021_ITS_rarefied) %>%
  mutate(sample=rownames(.)) %>%
  left_join(samples_df2, by="sample") %>%
  mutate(Site=factor(Site, levels=c("RA","RB","RC","MA","MB","MC","HA","HB","HC")),
         Time=factor(Time, levels=c("June","August","October")),
         SiteType=factor(.$SiteType, levels=c("High Elevation","Mid Elevation","Riparian")),
         source=dplyr::recode(source, "cDNA"="RNA")) %>%
  mutate(source=factor(.$source, levels=c("DNA", "RNA")))

#subsets alpha div dataframe to DNA only
AUE2021_alphadiv_DNA=
  AUE2021_alphadiv %>%
  filter(source=="DNA")

#subsets alpha div dataframe to RNA only
AUE2021_alphadiv_RNA=
  AUE2021_alphadiv %>%
  filter(source=="RNA")

#models shannon div with lmer - RNA only
mod=lmer(Shannon~SiteType*Time+(1|Plot), AUE2021_alphadiv_RNA)
Anova(mod)
summary(mod)

#models shannon div with lmer - DNA only
mod=lmer(Shannon~SiteType*Time+(1|Plot), AUE2021_alphadiv_DNA)
Anova(mod)
summary(mod)

#calculates mean and se of alpha div metrics
AUE2021_alphadiv_summary=
  AUE2021_alphadiv %>%
  group_by(Time, source, SiteType) %>%
  summarise(mean_Shannon=mean(Shannon),
            se_Shannon=sd(Shannon)/sqrt(n()),
            mean_InvSimpson=mean(InvSimpson),
            se_InvSimpson=sd(InvSimpson)/sqrt(n()),
            mean_Richness=mean(Observed),
            se_Richness=sd(Observed)/sqrt(n()))

#sets dodge for plotting
Dodge=0.4

#color palette to match the slightly transparent points w/ dark2 palette
#used in the ordination
Alphadivpalette=c("High Elevation" = "#72b599", 
                  "Mid Elevation" = "#e48c5a",
                  "Riparian"="#9791c5")

#plots shannon div
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

#saves shannon div plot as pdf
pdf(file="./OutputFiles/AlphaDivPlot_ITS.pdf", 
     width = 5.11811, height = 4.33071)
AlphaDivPlot_ITS
dev.off()

#agglomerates otu table at different taxonomic levels
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

#creates dataframe from phylum level object
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

#creates genus-level df
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
