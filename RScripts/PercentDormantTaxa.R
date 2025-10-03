#This script calculates the proportion of dormant taxa in each sample
#and then uses random forests to model environmental factors controlling
#this ratio.
#
#The proportion of dormant taxa is calculated as (1-n_taxa_RNA/n_taxa_total)
#n_taxa_total is calculated as the number of taxa found in either the RNA or
#DNA datasets

#This script requires that you run the scripts 16S_Analysis_SILVA.R,
#ITS2_Analysis.R, PlotBacterialPhyla.R, RandomForests_Taxa.R and 
#SoilVarDimensionalReduction.R and have the resulting data objects loaded in memory


library(tidyverse)
library(magrittr)
library(randomForest)

################################Proportion Dormant Fungi (ITS)###########################
#########################################################################################

#first create filtered dataframes that only contain samples present
#in both RNA and DNA
AUE2021_ITS_rarefied_DNA_df =
  AUE2021_ITS_rarefied_DNA %>%
  otu_table %>%
  as.data.frame

AUE2021_ITS_rarefied_RNA_df =
  AUE2021_ITS_rarefied_RNA %>%
  otu_table %>%
  as.data.frame

#renames RNA and DNA OTU tables so col names match
AUE2021_ITS_rarefied_DNA_df_rename=
  AUE2021_ITS_rarefied_DNA_df %>%
  rename_with(~ gsub("^DNA", "", .x))

AUE2021_ITS_rarefied_RNA_df_rename=
  AUE2021_ITS_rarefied_RNA_df %>%
  rename_with(~ gsub("^cDNA", "", .x))

#pick out samples that had successfully sequenced DNA and RNA
Samples_RNA_DNA_ITS=intersect(colnames(AUE2021_ITS_rarefied_DNA_df_rename),
                          colnames(AUE2021_ITS_rarefied_RNA_df_rename))

#filter DNA df to only include these samples
AUE2021_ITS_rarefied_DNA_df_filtered=
  select(AUE2021_ITS_rarefied_DNA_df_rename, all_of(Samples_RNA_DNA_ITS)) %>%
  rename_with(~ paste0("DNA", .x))

#filter RNA df to only include these samples
AUE2021_ITS_rarefied_RNA_df_filtered=
  select(AUE2021_ITS_rarefied_RNA_df_rename, all_of(Samples_RNA_DNA_ITS)) %>%
  rename_with(~ paste0("cDNA", .x))

#makes list of taxa present in either DNA or RNA dataset (or both)
nonzero_either_ITS <- 
  (AUE2021_ITS_rarefied_DNA_df_filtered != 0) |
  (AUE2021_ITS_rarefied_RNA_df_filtered != 0)

# calculates total # of taxa per sample (observed in either DNA or RNA sample) 
TotalTaxa_ITS <- 
  as.data.frame(colSums(nonzero_either_ITS)) %>%
  set_colnames(c("TotalTaxa"))  %>%
  rownames_to_column("sample") %>%
  mutate(sample = sub("^DNA_", "", sample))

#calculates total of active taxa per sample (i.e. observed in RNA)
ActiveTaxa_ITS=
  as.data.frame(colSums(AUE2021_ITS_rarefied_RNA_df_filtered != 0)) %>%
  set_colnames(c("ActiveTaxa"))  %>%
  rownames_to_column("sample") %>%
  mutate(sample = sub("^cDNA_", "", sample))

#calculates proportion dormant fungal taxa per sample
PropDormant_ITS=
  ActiveTaxa_ITS %>%
  left_join(TotalTaxa_ITS) %>%
  mutate(PropDormant=1-(ActiveTaxa/TotalTaxa)) %>%
  left_join((samples_df2_16S%>% mutate(sample = sub("^cDNA_", "", sample)) %>% select(sample:PlotDate))) %>%
  left_join((Env_cluster_summary_sensor_df %>% mutate(sample = sub("^DNA_", "", sample)) %>% select(sample, MACmoist:meansoiltemp))) %>%
  left_join((Env_cluster_summary_df %>% mutate(sample = sub("^DNA_", "", sample))))


#########################################################################################
############################Proportion dormant Bacteria##################################
#first create filtered dataframes that only contain samples present
#in both RNA and DNA
AUE2021_16S_rarefied_DNA_df =
  AUE2021_16S_rarefied_DNA %>%
  otu_table %>%
  as.data.frame

AUE2021_16S_rarefied_RNA_df =
  AUE2021_16S_rarefied_RNA %>%
  otu_table %>%
  as.data.frame

#renames RNA and DNA OTU tables so col names match
AUE2021_16S_rarefied_DNA_df_rename=
  AUE2021_16S_rarefied_DNA_df %>%
  rename_with(~ gsub("^DNA", "", .x))

AUE2021_16S_rarefied_RNA_df_rename=
  AUE2021_16S_rarefied_RNA_df %>%
  rename_with(~ gsub("^cDNA", "", .x))

#pick out samples that had successfully sequenced DNA and RNA
Samples_RNA_DNA_16S=intersect(colnames(AUE2021_16S_rarefied_DNA_df_rename),
                          colnames(AUE2021_16S_rarefied_RNA_df_rename))

#filter DNA df to only include these samples
AUE2021_16S_rarefied_DNA_df_filtered=
  select(AUE2021_16S_rarefied_DNA_df_rename, all_of(Samples_RNA_DNA_16S)) %>%
  rename_with(~ paste0("DNA", .x))

#filter RNA df to only include these samples
AUE2021_16S_rarefied_RNA_df_filtered=
  select(AUE2021_16S_rarefied_RNA_df_rename, all_of(Samples_RNA_DNA_16S)) %>%
  rename_with(~ paste0("cDNA", .x))

#makes list of taxa present in either DNA or RNA dataset (or both)
nonzero_either_16S <- 
  (AUE2021_16S_rarefied_DNA_df_filtered != 0) |
  (AUE2021_16S_rarefied_RNA_df_filtered != 0)

# calculates total # of taxa per sample (observed in either DNA or RNA sample) 
TotalTaxa_16S <- 
  as.data.frame(colSums(nonzero_either_16S)) %>%
  set_colnames(c("TotalTaxa"))  %>%
  rownames_to_column("sample") %>%
  mutate(sample = sub("^DNA_", "", sample))

#calculates total of active taxa per sample (i.e. observed in RNA)
ActiveTaxa_16S=
  as.data.frame(colSums(AUE2021_16S_rarefied_RNA_df_filtered != 0)) %>%
  set_colnames(c("ActiveTaxa"))  %>%
  rownames_to_column("sample") %>%
  mutate(sample = sub("^cDNA_", "", sample))

#calculates proportion dormant fungal taxa per sample
PropDormant_16S=
  ActiveTaxa_16S %>%
  left_join(TotalTaxa_16S) %>%
  mutate(PropDormant=1-(ActiveTaxa/TotalTaxa)) %>%
  left_join((samples_df2_16S%>% mutate(sample = sub("^cDNA_", "", sample)) %>% select(sample:PlotDate))) %>%
  left_join((Env_cluster_summary_sensor_df %>% mutate(sample = sub("^DNA_", "", sample)) %>% select(sample, MACmoist:meansoiltemp))) %>%
  left_join((Env_cluster_summary_df %>% mutate(sample = sub("^DNA_", "", sample))))

########################################################################################
#########################Random Forest Modelling of % dormant fungi#####################
#rf modelling without sensor data
rf_PropDormant_no_sensor_ITS=randomForest(PropDormant_ITS$PropDormant ~ ., 
                             data=select(PropDormant_ITS, NO3Flux:NConc),
                             importance=TRUE)

#formats importance values of random forest as data frame
PropDormant_no_sensor_ITS_df=
  importance(rf_PropDormant_no_sensor_ITS) %>%
  as.data.frame


####rf modelling with sensor data
#create df with dormancy rates filtered to only include samples with sensor data
PropDormant_ITS_sensor=
  PropDormant_ITS %>%
  filter(!(is.na(MACmoist)))

#do random forest modelling
rf_PropDormant_sensor_ITS=randomForest(PropDormant_ITS_sensor$PropDormant ~ ., 
                                          data=select(PropDormant_ITS_sensor, MACmoist:NConc),
                                          importance=TRUE)

#formats importance values of random forest as data frame
PropDormant_sensor_ITS_df=
  importance(rf_PropDormant_sensor_ITS) %>%
  as.data.frame

#########################Random Forest Modelling of % dormant bacteria#####################
#rf modelling without sensor data
rf_PropDormant_no_sensor_16S=randomForest(PropDormant_16S$PropDormant ~ ., 
                                          data=select(PropDormant_16S, NO3Flux:NConc),
                                          importance=TRUE)

#formats importance values of random forest as data frame
PropDormant_no_sensor_16S_df=
  importance(rf_PropDormant_no_sensor_16S) %>%
  as.data.frame


####rf modelling with sensor data
#create df with dormancy rates filtered to only include samples with sensor data
PropDormant_16S_sensor=
  PropDormant_16S %>%
  filter(!(is.na(MACmoist)))

#do random forest modelling
rf_PropDormant_sensor_16S=randomForest(PropDormant_16S_sensor$PropDormant ~ ., 
                                       data=select(PropDormant_16S_sensor, MACmoist:NConc),
                                       importance=TRUE)

#formats importance values of random forest as data frame
PropDormant_sensor_16S_df=
  importance(rf_PropDormant_sensor_16S) %>%
  as.data.frame





########################################################################################

