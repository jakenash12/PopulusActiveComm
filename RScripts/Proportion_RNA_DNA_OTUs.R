#This script calculates the proportion of OTUs observed in the total community that are also
#observed in the active community. This shows that the active community is a subset of the
#total community.

#This is done by converting OTU tables of active and total communities into long format,
#filtering to only include samples that had a successfully sequenced active and total community
#and then taking number of OTUs that occurred in both active community and total community
#and dividing by number of OTUs that occurred in active community

#This script requires that you run the scripts 16S_Analysis_SILVA.R, 
#ITS2_Analysis.R, and PercentDormantTaxa.R and have the resulting data 
#objects loaded in memory

library(tidyverse)

#############################ITS#############################################

#converts DNA OTU table to long format then codes as binary
#0=absent taxa, 1=present taxa
AUE2021_ITS_rarefied_df_DNA_long_binary <- 
  AUE2021_ITS_rarefied_df_DNA %>%
  rownames_to_column(var = "sample") %>%   # make rownames into a column called "sample"
  pivot_longer(
    cols = -sample,                        # all columns except "sample"
    names_to = "OTU",                      # column names become "OTU"
    values_to = "DNA_presence"                    # cell values become "value"
  ) %>%
  mutate(
    sample = str_replace(sample, "^DNA", "sample"),
    DNA_presence = if_else(DNA_presence > 0, 1, 0)   # recode to binary
  )

#converts DRA OTU table to long format then codes as binary
#0=absent taxa, 1=present taxa
AUE2021_ITS_rarefied_df_RNA_long_binary <- 
  AUE2021_ITS_rarefied_df_RNA %>%
  rownames_to_column(var = "sample") %>%   # make rownames into a column called "sample"
  pivot_longer(
    cols = -sample,                        # all columns except "sample"
    names_to = "OTU",                      # column names become "OTU"
    values_to = "RNA_presence"                    # cell values become "value"
  ) %>%
  mutate(
    sample = str_replace(sample, "^cDNA", "sample"),
    RNA_presence = if_else(RNA_presence > 0, 1, 0)   # recode to binary
  )

#joins binary long dfs of RNA and DNA together,
#then uses that to calculate proportions of active taxa in each sample
AUE2021_ITS_prop_RNA_DNA=
  AUE2021_ITS_rarefied_df_DNA_long_binary %>%
  left_join(AUE2021_ITS_rarefied_df_RNA_long_binary) %>%
  filter(DNA_presence==1, !(is.na(RNA_presence))) %>%
  group_by(sample) %>%
  summarise(prop_RNA=sum(RNA_presence)/n())


#############################16S#############################################

#converts DNA OTU table to long format then codes as binary
#0=absent taxa, 1=present taxa
AUE2021_16S_rarefied_df_DNA_long_binary <- 
  AUE2021_16S_rarefied_df_DNA %>%
  rownames_to_column(var = "sample") %>%   # make rownames into a column called "sample"
  pivot_longer(
    cols = -sample,                        # all columns except "sample"
    names_to = "OTU",                      # column names become "OTU"
    values_to = "DNA_presence"                    # cell values become "value"
  ) %>%
  mutate(
    sample = str_replace(sample, "^DNA", "sample"),
    DNA_presence = if_else(DNA_presence > 0, 1, 0)   # recode to binary
  )

#converts DRA OTU table to long format then codes as binary
#0=absent taxa, 1=present taxa
AUE2021_16S_rarefied_df_RNA_long_binary <- 
  AUE2021_16S_rarefied_df_RNA %>%
  rownames_to_column(var = "sample") %>%   # make rownames into a column called "sample"
  pivot_longer(
    cols = -sample,                        # all columns except "sample"
    names_to = "OTU",                      # column names become "OTU"
    values_to = "RNA_presence"                    # cell values become "value"
  ) %>%
  mutate(
    sample = str_replace(sample, "^cDNA", "sample"),
    RNA_presence = if_else(RNA_presence > 0, 1, 0)   # recode to binary
  )

#joins binary long dfs of RNA and DNA together,
#then uses that to calculate proportions of active taxa in each sample
AUE2021_16S_prop_RNA_DNA=
  AUE2021_16S_rarefied_df_DNA_long_binary %>%
  left_join(AUE2021_16S_rarefied_df_RNA_long_binary) %>%
  filter(DNA_presence==1, !(is.na(RNA_presence))) %>%
  group_by(sample) %>%
  summarise(prop_RNA=sum(RNA_presence)/n())

