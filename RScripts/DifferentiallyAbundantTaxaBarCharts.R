#This script makes a stacked barchart showing which season/site type fungal + bacterial
#OTUs reach their maximum abundance in for the DNA and RNA datasets

library(tidyverse)

########create phyloseq objects split into RNA and DNA############
AUE2021_ITS_rarefied_DNA =
  AUE2021_ITS_rarefied %>%
  subset_samples(source=="DNA")

AUE2021_ITS_rarefied_RNA =
  AUE2021_ITS_rarefied %>%
  subset_samples(source=="cDNA")

AUE2021_16S_rarefied_DNA =
  AUE2021_16S_rarefied %>%
  subset_samples(source=="DNA")

AUE2021_16S_rarefied_RNA =
  AUE2021_16S_rarefied %>%
  subset_samples(source=="cDNA")


################First create dataframes indicating which site type/season##############
##################OTUS peaked in abundance in for DNA/RNA and ITS/16S##################

####################Sitetype dfs
SiteType_ITS_Sort_DNA=
  AUE2021_ITS_rarefied_DNA %>%
  otu_table %>%
  as.data.frame %>%
  filter(rownames(.) %in% SiteType_ITS_variable_DNA$otu) %>%
  t %>%
  as.data.frame %>%
  mutate(sample=rownames(.)) %>%
  left_join(samples_df2) %>%
  group_by(Site) %>%
  summarise_at(vars(SiteType_ITS_variable_DNA$otu), funs(mean(., na.rm=TRUE))) %>%
  left_join(select(samples_df2, Site,SiteType)) %>%
  group_by(SiteType) %>%
  summarise_at(vars(SiteType_ITS_variable_DNA$otu), funs(mean(., na.rm=TRUE))) %>%
  column_to_rownames("SiteType") %>%
  t %>%
  as.data.frame %>%
  mutate(otu=rownames(.),
         Peak=ifelse((`High Elevation`>`Mid Elevation`&`High Elevation`>`Riparian`),"High",
                     ifelse((`Mid Elevation`>Riparian&`Mid Elevation`>`High Elevation`),"Mid","Riparian")),
         Peak=factor(Peak, levels=c("High","Mid","Riparian"))) %>%
  arrange(Peak) %>%
  mutate(otu=factor(otu, levels=.$otu)) %>%
  left_join(tax_mat_df)  %>%
  select(otu, Peak, 
         Phylum, Class, Order, Family, 
         Genus, Species, TaxClass) %>%
  mutate(source="DNA", locus="ITS2", Type="SiteType")

SiteType_ITS_Sort_RNA=
  AUE2021_ITS_rarefied_RNA %>%
  otu_table %>%
  as.data.frame %>%
  filter(rownames(.) %in% SiteType_ITS_variable_RNA$otu) %>%
  t %>%
  as.data.frame %>%
  mutate(sample=rownames(.)) %>%
  left_join(samples_df2) %>%
  group_by(Site) %>%
  summarise_at(vars(SiteType_ITS_variable_RNA$otu), funs(mean(., na.rm=TRUE))) %>%
  left_join(select(samples_df2, Site,SiteType)) %>%
  group_by(SiteType) %>%
  summarise_at(vars(SiteType_ITS_variable_RNA$otu), funs(mean(., na.rm=TRUE))) %>%
  column_to_rownames("SiteType") %>%
  t %>%
  as.data.frame %>%
  mutate(otu=rownames(.),
         Peak=ifelse((`High Elevation`>`Mid Elevation`&`High Elevation`>`Riparian`),"High",
                     ifelse((`Mid Elevation`>Riparian&`Mid Elevation`>`High Elevation`),"Mid","Riparian")),
         Peak=factor(Peak, levels=c("High","Mid","Riparian"))) %>%
  arrange(Peak) %>%
  mutate(otu=factor(otu, levels=.$otu)) %>%
  left_join(tax_mat_df)  %>%
  select(otu, Peak, 
         Phylum, Class, Order, Family, 
         Genus, Species, TaxClass) %>%
  mutate(source="RNA", locus="ITS2", Type="SiteType")

SiteType_16S_Sort_DNA=
  AUE2021_16S_rarefied_DNA %>%
  otu_table %>%
  as.data.frame %>%
  filter(rownames(.) %in% SiteType_16S_variable_DNA$otu) %>%
  t %>%
  as.data.frame %>%
  mutate(sample=rownames(.)) %>%
  left_join(samples_df2_16S) %>%
  group_by(Site) %>%
  summarise_at(vars(SiteType_16S_variable_DNA$otu), funs(mean(., na.rm=TRUE))) %>%
  left_join(select(samples_df2_16S, Site,SiteType)) %>%
  group_by(SiteType) %>%
  summarise_at(vars(SiteType_16S_variable_DNA$otu), funs(mean(., na.rm=TRUE))) %>%
  column_to_rownames("SiteType") %>%
  t %>%
  as.data.frame %>%
  mutate(otu=rownames(.),
         Peak=ifelse((`High Elevation`>`Mid Elevation`&`High Elevation`>`Riparian`),"High",
                     ifelse((`Mid Elevation`>Riparian&`Mid Elevation`>`High Elevation`),"Mid","Riparian")),
         Peak=factor(Peak, levels=c("High","Mid","Riparian"))) %>%
  arrange(Peak) %>%
  mutate(otu=factor(otu, levels=.$otu)) %>%
  left_join(silva_tax_df)  %>%
  select(otu, Peak, 
         Phylum, Class, Order, Family, 
         Genus, Species, TaxClass) %>%
  mutate(source="DNA", locus="16S", Type="SiteType")


SiteType_16S_Sort_RNA=
  AUE2021_16S_rarefied_RNA %>%
  otu_table %>%
  as.data.frame %>%
  filter(rownames(.) %in% SiteType_16S_variable_RNA$otu) %>%
  t %>%
  as.data.frame %>%
  mutate(sample=rownames(.)) %>%
  left_join(samples_df2_16S) %>%
  group_by(Site) %>%
  summarise_at(vars(SiteType_16S_variable_RNA$otu), funs(mean(., na.rm=TRUE))) %>%
  left_join(select(samples_df2_16S, Site,SiteType)) %>%
  group_by(SiteType) %>%
  summarise_at(vars(SiteType_16S_variable_RNA$otu), funs(mean(., na.rm=TRUE))) %>%
  column_to_rownames("SiteType") %>%
  t %>%
  as.data.frame %>%
  mutate(otu=rownames(.),
         Peak=ifelse((`High Elevation`>`Mid Elevation`&`High Elevation`>`Riparian`),"High",
                     ifelse((`Mid Elevation`>Riparian&`Mid Elevation`>`High Elevation`),"Mid","Riparian")),
         Peak=factor(Peak, levels=c("High","Mid","Riparian"))) %>%
  arrange(Peak) %>%
  mutate(otu=factor(otu, levels=.$otu)) %>%
  left_join(silva_tax_df)  %>%
  select(otu, Peak, 
         Phylum, Class, Order, Family, 
         Genus, Species, TaxClass) %>%
  mutate(source="RNA", locus="16S", Type="SiteType")

##################season dfs
Time_ITS_Sort_DNA=
  AUE2021_ITS_rarefied_DNA %>%
  otu_table %>%
  as.data.frame %>%
  filter(rownames(.) %in% Time_ITS_variable_DNA$otu) %>%
  t %>%
  as.data.frame %>%
  mutate(sample=rownames(.)) %>%
  left_join(samples_df2) %>%
  mutate(Time=factor(Time,levels=c("June","August","October")),
         TimeSiteType=paste(Time,SiteType,sep="_")) %>%
  group_by(TimeSiteType) %>%
  summarise_at(vars(Time_ITS_variable_DNA$otu), funs(mean(., na.rm=TRUE))) %>%
  left_join(select(samples_df2, TimeSiteType,Time)) %>%
  group_by(Time) %>%
  summarise_at(vars(Time_ITS_variable_DNA$otu), funs(mean(., na.rm=TRUE))) %>%
  column_to_rownames("Time") %>%
  t %>%
  as.data.frame %>%
  mutate(otu=rownames(.),
         Peak=ifelse((`August`>`June`&`August`>`October`),"August",
                     ifelse((`June`>October&`June`>`August`),"June","October")),
         Peak=factor(Peak, levels=c("June","August","October"))) %>%
  arrange(Peak) %>%
  mutate(otu=factor(otu, levels=.$otu)) %>%
  left_join(tax_mat_df) %>%
  select(otu, Peak, 
         Phylum, Class, Order, Family, 
         Genus, Species, TaxClass) %>%
  mutate(source="DNA", locus="ITS2", Type="Time")

Time_ITS_Sort_RNA=
  AUE2021_ITS_rarefied_RNA %>%
  otu_table %>%
  as.data.frame %>%
  filter(rownames(.) %in% Time_ITS_variable_RNA$otu) %>%
  t %>%
  as.data.frame %>%
  mutate(sample=rownames(.)) %>%
  left_join(samples_df2) %>%
  mutate(Time=factor(Time,levels=c("June","August","October")),
         TimeSiteType=paste(Time,SiteType,sep="_")) %>%
  group_by(TimeSiteType) %>%
  summarise_at(vars(Time_ITS_variable_RNA$otu), funs(mean(., na.rm=TRUE))) %>%
  left_join(select(samples_df2, TimeSiteType,Time)) %>%
  group_by(Time) %>%
  summarise_at(vars(Time_ITS_variable_RNA$otu), funs(mean(., na.rm=TRUE))) %>%
  column_to_rownames("Time") %>%
  t %>%
  as.data.frame %>%
  mutate(otu=rownames(.),
         Peak=ifelse((`August`>`June`&`August`>`October`),"August",
                     ifelse((`June`>October&`June`>`August`),"June","October")),
         Peak=factor(Peak, levels=c("June","August","October"))) %>%
  arrange(Peak) %>%
  mutate(otu=factor(otu, levels=.$otu)) %>%
  left_join(tax_mat_df) %>%
  select(otu, Peak, 
         Phylum, Class, Order, Family, 
         Genus, Species, TaxClass) %>%
  mutate(source="RNA", locus="ITS2", Type="Time")

Time_16S_Sort_DNA=
  AUE2021_16S_rarefied_DNA %>%
  otu_table %>%
  as.data.frame %>%
  filter(rownames(.) %in% Time_16S_variable_DNA$otu) %>%
  t %>%
  as.data.frame %>%
  mutate(sample=rownames(.)) %>%
  left_join(samples_df2_16S) %>%
  mutate(Time=factor(Time,levels=c("June","August","October")),
         TimeSiteType=paste(Time,SiteType,sep="_")) %>%
  group_by(TimeSiteType) %>%
  summarise_at(vars(Time_16S_variable_DNA$otu), funs(mean(., na.rm=TRUE))) %>%
  left_join(select(samples_df2_16S, TimeSiteType,Time)) %>%
  group_by(Time) %>%
  summarise_at(vars(Time_16S_variable_DNA$otu), funs(mean(., na.rm=TRUE))) %>%
  column_to_rownames("Time") %>%
  t %>%
  as.data.frame %>%
  mutate(otu=rownames(.),
         Peak=ifelse((`August`>`June`&`August`>`October`),"August",
                     ifelse((`June`>October&`June`>`August`),"June","October")),
         Peak=factor(Peak, levels=c("June","August","October"))) %>%
  arrange(Peak) %>%
  mutate(otu=factor(otu, levels=.$otu)) %>%
  left_join(silva_tax_df)  %>%
  select(otu, Peak, 
         Phylum, Class, Order, Family, 
         Genus, Species, TaxClass) %>%
  mutate(source="DNA", locus="16S", Type="Time")

Time_16S_Sort_RNA=
  AUE2021_16S_rarefied_RNA %>%
  otu_table %>%
  as.data.frame %>%
  filter(rownames(.) %in% Time_16S_variable_RNA$otu) %>%
  t %>%
  as.data.frame %>%
  mutate(sample=rownames(.)) %>%
  left_join(samples_df2_16S) %>%
  mutate(Time=factor(Time,levels=c("June","August","October")),
         TimeSiteType=paste(Time,SiteType,sep="_")) %>%
  group_by(TimeSiteType) %>%
  summarise_at(vars(Time_16S_variable_RNA$otu), funs(mean(., na.rm=TRUE))) %>%
  left_join(select(samples_df2_16S, TimeSiteType,Time)) %>%
  group_by(Time) %>%
  summarise_at(vars(Time_16S_variable_RNA$otu), funs(mean(., na.rm=TRUE))) %>%
  column_to_rownames("Time") %>%
  t %>%
  as.data.frame %>%
  mutate(otu=rownames(.),
         Peak=ifelse((`August`>`June`&`August`>`October`),"August",
                     ifelse((`June`>October&`June`>`August`),"June","October")),
         Peak=factor(Peak, levels=c("June","August","October"))) %>%
  arrange(Peak) %>%
  mutate(otu=factor(otu, levels=.$otu)) %>%
  left_join(silva_tax_df)  %>%
  select(otu, Peak, 
         Phylum, Class, Order, Family, 
         Genus, Species, TaxClass) %>%
  mutate(source="RNA", locus="16S", Type="Time")
########################################################################

#combines all dfs of site type-variable OTUs (16S and ITS, DNA and RNA) 
SiteTypeDiff=
  rbind(SiteType_ITS_Sort_DNA, 
        SiteType_16S_Sort_DNA,
        SiteType_ITS_Sort_RNA,
        SiteType_16S_Sort_RNA) 

#creates a summary of the # of OTUs that differed by site type in each df and 
#by the site type in which it reached its maximum abundance
SiteTypeDiffSum =
  SiteTypeDiff %>%
  group_by(source, locus, Peak, Type) %>%
  summarise(TaxaCount=n())

#combines all dfs of seasonally-variable OTUs (16S and ITS, DNA and RNA) 
TimeDiff=
  rbind(Time_ITS_Sort_DNA, 
        Time_16S_Sort_DNA,
        Time_ITS_Sort_RNA,
        Time_16S_Sort_RNA)

#creates a summary of the # of OTUs that differed by season in each df and 
#by the season in which it reached its maximum abundance
TimeDiffSum=
  TimeDiff %>%
  group_by(source, locus, Peak, Type) %>%
  summarise(TaxaCount=n())

#combines site type and season summary dfs of number of differential abundant
AllDiffSum=
  rbind(SiteTypeDiffSum, TimeDiffSum) %>%
  mutate(Peak=factor(Peak, levels=c("High","Mid", "Riparian", "October","August", "June")),
         locus=factor(locus, levels=c("ITS2", "16S")))

#creates stacked barplot showing # of differentially abundant OTUs in each datasete by 
#site type and season
DiffAbundantPlot=
  ggplot(AllDiff, aes(x=source, y=TaxaCount, fill=Peak)) +
  geom_bar(position="stack", stat="identity", colour = "gray20", size=0.5) +
  facet_grid(Type~locus, scales="free") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks=seq(0,100, by = 10)) +
  theme_test() +
  scale_fill_manual(values=c("#1B9E77", "#D95F02", "#7570B3", "#B2182B", "#BF812D", "#2166AC")) +
  theme(legend.position="none", text=element_blank()) +
  theme(panel.border = element_rect(linetype = "solid",
                                    colour = "gray20", linewidth = 1.5),
        strip.background = element_blank(),
        axis.ticks=element_line(size=1.5, colour = "gray20")) 
  
#saves stacked barplot as PDF to be made into a publication plot
pdf(file="DiffAbunStackedBarchart.pdf", 
    width = 3.5, height = 5.11811)
DiffAbundantPlot
dev.off()

