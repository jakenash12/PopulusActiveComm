#This script defines the core microbiome in both the total and the active community.
#This is done by using an occupancy threshold of 0.95 (i.e. the OTU/ASV must be present
#in >95% of samples)

#This script requires that you run the scripts 16S_Analysis_SILVA.R and
#ITS2_Analysis.R and have the resulting data objects loaded in memory


library(tidyverse)
library(magrittr)
library(ggrepel)
library(ggplot2)

#####################ITS########################
#creates dataframe of OTU occupancy
AUE2021_ITS_occupancy=
  AUE2021_ITS_rarefied %>% 
  otu_table %>%
  as.data.frame %>%
  mutate(across(everything(), ~ ifelse(. > 0, 1, .))) %>%
  mutate(otu=rownames(.)) %>%
  gather("sample", "Presence", -otu) %>%
  left_join(select(samples_df2, sample, source)) %>%
  group_by(otu, source) %>%
  summarise(Occupancy=sum(Presence)/n()) %>% 
  ungroup() %>%
  mutate(source = ifelse(source == "cDNA", "RNA", source)) %>%
  left_join(tax_mat_df) %>% 
  select(-Domain)

#creates df with new column indicating core status of each otu
#"core" means > 0.95 occupancy in both DNA and RNA datasets
#"DNA Cores" means 0.95 occupancy in DNA dataset
#"RNA Cores" means 0.95 occupancy in RNA dataset
#Then filters to only include OTUs that are core in DNA and/or RNA
Core_ITS=
  AUE2021_ITS_occupancy %>%
  spread(source, Occupancy) %>%
  mutate(CoreStatus=case_when(
    (RNA>0.95 & DNA>0.95) ~ "Core",
    (RNA>0.95 & DNA<0.95) ~ "RNA Core",
    (RNA<0.95 & DNA>0.95) ~ "DNA Core",
    .default="Not Core"
  ))  %>%
  filter(CoreStatus != "Not Core") %>%
  mutate(Dataset="ITS")

#####################16S########################
#creates dataframe of OTU occupancy
AUE2021_16S_occupancy=
  AUE2021_16S_rarefied %>% 
  otu_table %>%
  as.data.frame %>%
  mutate(across(everything(), ~ ifelse(. > 0, 1, .))) %>%
  mutate(otu=rownames(.)) %>%
  gather("sample", "Presence", -otu) %>%
  left_join(select(samples_df2_16S, sample, source)) %>%
  group_by(otu, source) %>%
  summarise(Occupancy=sum(Presence)/n()) %>% 
  ungroup() %>%
  mutate(source = ifelse(source == "cDNA", "RNA", source)) %>%
  left_join(silva_tax_df) %>%
  select(-Domain)

#creates df with new column indicating core status of each otu
#"core" means > 0.95 occupancy in both DNA and RNA datasets
#"DNA Cores" means 0.95 occupancy in DNA dataset
#"RNA Cores" means 0.95 occupancy in RNA dataset
#Then filters to only include OTUs that are core in DNA and/or RNA
Core_16S=
  AUE2021_16S_occupancy %>%
  spread(source, Occupancy) %>%
  mutate(CoreStatus=case_when(
    (RNA>0.95 & DNA>0.95) ~ "Core",
    (RNA>0.95 & DNA<0.95) ~ "RNA Core",
    (RNA<0.95 & DNA>0.95) ~ "DNA Core",
    .default="Not Core"
  )) %>%
  filter(CoreStatus != "Not Core") %>%
  mutate(Dataset="16S") %>%
  mutate(across(Phylum:Species, ~ na_if(.x, "Incertae_Sedis"))) %>%
  mutate(TaxClass=ifelse(!is.na(Species),Species,ifelse(!is.na(Genus),Genus,ifelse(!is.na(Family),Family,ifelse(!is.na(Order),Order,ifelse(!is.na(Class),Class,ifelse(!is.na(Phylum),Phylum,Domain)))))))



Core_ITS_16S=
  rbind(Core_16S, Core_ITS) %>%
  mutate(TaxClass=case_when(is.na(TaxClass)~ "Unassigned",
                            .default=TaxClass)) %>%
  pivot_longer(cols = c(DNA, RNA),
               names_to = "source",
               values_to = "Occupancy")
  
# Get ordering of OTUs: first DNA, then RNA
otu_order = Core_ITS_16S %>%
  pivot_wider(names_from = source, values_from = Occupancy) %>%
  arrange(desc(DNA), desc(RNA)) %>%
  pull(otu)

# Use that order in the plot
CoreHeatmap =
  ggplot(Core_ITS_16S, 
         aes(x = source, y = factor(otu, levels = rev(otu_order)), fill = Occupancy)) +
  geom_tile(color="white", size=0.35) +
  facet_grid(Dataset ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(otu) Core_ITS_16S$TaxClass[match(otu, Core_ITS_16S$otu)]) +
  scale_fill_distiller(palette = "Blues", limits = c(0,1), na.value = "#de2d26",
                       direction = 1) +
  theme_minimal() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank())


#makes heatmap of the occupancy of core OTUs in the ITS and 16S dataset.
#For use in supplement
CoreHeatmap=
  ggplot(Core_ITS_16S, 
       aes(x = source, y = factor(otu, levels = unique(otu)), fill = Occupancy)) +
  geom_tile(color="white", size=0.35) +
  facet_grid(Dataset ~ ., scales = "free_y", space = "free_y") +  # Ensure equal tile heights
  scale_y_discrete(labels = function(otu) Core_ITS_16S$TaxClass[match(otu, Core_ITS_16S$otu)]) +
  scale_fill_distiller(palette = "Blues", limits = c(0,1), na.value = "#de2d26",
                       direction = 1) +
  theme_minimal() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank())

#saves core community heatmap as pdf
pdf(file="CoreHeatmap.pdf", 
    width = 4, height = 10)
CoreHeatmap
dev.off()

#saves list of core microbial members 
write_delim(Core_16S, "Core_16S.tsv", delim="\t")
write_delim(Core_ITS, "Core_ITS.tsv", delim="\t")
