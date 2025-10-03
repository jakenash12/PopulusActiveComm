#This script makes stacked bar charts of the sequence abundance of
#dominant bacterial phyla and fungal classes (>0.5% avg seq abundance)
#in the total and active communities, split by season and site type

#This script requires that you run the scripts 16S_Analysis_SILVA.R
#and ITS2_Analysis.R and have the resulting data objects loaded in memory

library(tidyverse)
library(magrittr)

###################16S########################
#calculates average relative abundance of each phylum across the dataset
#(for the case of Pseudomonadota, class-level is used)
PhylumMeans_16S=
  AUE2021_16S_rarefied_rel %>%
  otu_table %>%
  as.data.frame %>%
  mutate(otu=rownames(.)) %>%
  left_join(silva_tax_df) %>% 
  mutate(Phylum_edit=case_when(
    Phylum=="Pseudomonadota" ~ Class,
    .default=Phylum
  )) %>%
  filter(!(is.na(Phylum_edit))) %>%
  group_by(Phylum_edit) %>%
  summarise_at(vars(matches('DNA')), sum) %>% 
  column_to_rownames("Phylum_edit") %>%
  rowMeans() %>%
  as.data.frame %>%
  set_colnames("Abundance") %>%
  mutate(Phylum=rownames(.)) %>%
  arrange(desc(Abundance)) %>%
  mutate(Rank=row_number())

#creates list of phyla with >0.5% sequence abundance
DomPhyla_16S=
  PhylumMeans_16S %>%
  filter(Abundance>0.005) %$%
  Phylum

#creates wide-format dataframe of phylum-level abundance in 16S dataset
AUE2021_16S_phylum_summary_wide=
  AUE2021_16S_rarefied_rel %>%
  otu_table %>%
  as.data.frame %>%
  mutate(otu=rownames(.)) %>%
  left_join(silva_tax_df) %>% 
  mutate(Phylum_edit=case_when(
    (Phylum=="Pseudomonadota"& Class %in% DomPhyla_16S) ~ Class,
    Phylum %in% DomPhyla_16S ~ Phylum,
    .default="Other"
  )) %>% 
  group_by(Phylum_edit) %>%
  summarise_at(vars(matches('DNA')), sum) %>%  
  column_to_rownames("Phylum_edit") %>%
  t %>%
  as.data.frame %>%
  mutate(sample=rownames(.)) 

#generates long-format data of phylum-level abundance in 16S dataset
#sets order for plotting of phyla (more highly abundant phyla will be plotted at bottom
#of stacked bar plot)
AUE2021_16S_phylum_summary_long =
  AUE2021_16S_phylum_summary_wide %>% 
  gather("Phylum","Abundance", -sample) %>%
  left_join(samples_df2_16S) %>%
  mutate(Phylum=factor(Phylum, levels=c("Other", rev(DomPhyla_16S)))) %>%
  group_by(Phylum, SiteType, Time, source) %>%
  summarise(meanAbundance=mean(Abundance)) %>%
  mutate(source=factor(source, levels=c("DNA", "cDNA")))

#creates color palette for plotting, manually setting colors for each phylum
BacPhylumPal=
  c("Actinomycetota"="#4c83b9", #formerly Actinobacteria
    "Alphaproteobacteria"="#bb4d53", 
    "Gammaproteobacteria"="#9ebd56",
    "Bacteroidota"="#80649e", #formerly Bacteroidetes
    "Myxococcota"="#ff9d45", #formerly Deltaproteobacteria
    "Acidobacteriota"="#2f5982", #formerly Acidobacteria
    "Planctomycetota"="#78362b", #formerly Planctomycetes
    "Verrucomicrobiota"="#5d7331", #formerly Verrucomicrobia
    "Chloroflexota"="#4b3a5f", #formerly Chloroflexi
    "Gemmatimonadota"="#236c7a", #formerly Gemmatimonadetes
    "Bacillota"="#ba5502",
    "Thermoproteota"="#799bcb", #formerly Crenarchaeota
    "Patescibacteria"="#cf6e75", #formerly TM7
    "Other"="grey65")

#makes stacked bar chart of phylum-level sequence abundance by site type and season
#in DNA and RNA datasets
Phylum_plot_16S=
  ggplot(AUE2021_16S_phylum_summary_long, aes(x=Time, y=meanAbundance, fill=Phylum)) +
  geom_bar(stat="identity") +
  facet_grid(source~SiteType) +
  theme_test() +
  scale_fill_manual(values = BacPhylumPal) +
  scale_y_continuous(limits = c(0,1.05), expand = c(0, 0)) +
  theme(legend.position="none")

###################ITS########################
#calculates average relative abundance of each fungal class across the dataset
ClassMeans_ITS=
  AUE2021_ITS_rarefied_rel %>%
  otu_table %>%
  as.data.frame %>%
  mutate(otu=rownames(.)) %>%
  left_join(tax_mat_df) %>% 
  filter(!(is.na(Class))) %>%
  group_by(Class) %>%
  summarise_at(vars(matches('DNA')), sum) %>% 
  column_to_rownames("Class") %>%
  rowMeans() %>%
  as.data.frame %>%
  set_colnames("Abundance") %>%
  mutate(Class=rownames(.)) %>%
  arrange(desc(Abundance)) %>%
  mutate(Rank=row_number())

#creates list of fungal classes with >0.5% sequence abundance
DomClass_ITS=
  ClassMeans_ITS %>%
  filter(Abundance>0.005) %$%
  Class

#creates a wide format dataframe of class-level sequence abundance
AUE2021_ITS_class_summary_wide=
  AUE2021_ITS_rarefied_rel %>%
  otu_table %>%
  as.data.frame %>%
  mutate(otu=rownames(.)) %>%
  left_join(tax_mat_df) %>% 
  mutate(Class_edit=
           case_when(Class %in% DomClass_ITS ~ Class,
                     .default="Other")) %>%
  group_by(Class_edit) %>%
  summarise_at(vars(matches('DNA')), sum) %>%  
  column_to_rownames("Class_edit") %>%
  t %>%
  as.data.frame %>%
  mutate(sample=rownames(.)) 

#creates a long-format dataframe of class-level sequence abundance
AUE2021_ITS_class_summary_long =
  AUE2021_ITS_class_summary_wide %>% 
  gather("Class_edit","Abundance", -sample) %>%
  left_join(samples_df2) %>%
  mutate(Class_edit=factor(Class_edit, levels=c("Other", rev(DomClass_ITS)))) %>%
  group_by(Class_edit, SiteType, Time, source) %>%
  summarise(meanAbundance=mean(Abundance)) %>%
  mutate(source=factor(source, levels=c("DNA", "cDNA")))

#defines color pallette for plotting stacked barchart of fungal class sequence abundance
ClassPal_ITS=
  c("Agaricomycetes"="#4c83b9",
    "Leotiomycetes"="#bb4d53", 
    "Dothideomycetes"="#9ebd56",
    "Pezizomycetes"="#80649e",
    "Sordariomycetes"="#ff9d45",
    "Eurotiomycetes"="#2f5982",
    "Glomeromycetes"="#78362b",
    "Other"="grey65")

#creates class-level stacked taxonomic barchart for fungi
Class_plot_ITS=
  ggplot(AUE2021_ITS_class_summary_long, aes(x=Time, y=meanAbundance, fill=Class_edit)) +
  geom_bar(stat="identity") +
  facet_grid(source~SiteType) +
  theme_test() +
  scale_fill_manual(values = ClassPal_ITS) +
  scale_y_continuous(limits = c(0,1.05), expand = c(0, 0)) +
  theme(legend.position="none")
  

#saves stacked taxonomic bar charts as pdfs to make into figures
pdf(file="Class_plot_ITS.pdf", 
    width = 6, height = 5)
Class_plot_ITS
dev.off()

pdf(file="Phylum_plot_16S.pdf", 
    width = 6, height = 5)
Phylum_plot_16S
dev.off()
