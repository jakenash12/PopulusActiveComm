#This script uses random forests to identify soil variables that best
#predict fungal and bacterial total and active community composition.
#This is done by using the first two axes of NMDS as response variables.
#A dimensionally-reduced dataframe of soil variables is used as predictors
#in the random forests

#Some plots did not have soil sensor data (providing soil temp and moisture measurements)
#Missing data is not allowed in random forest models. To deal with this, we first did
#random forest modelling without sensor data (and including all plots).
#Then, we did random forest modelling with the sensor data, excluding the plots that did
#not have sensors installed.

#This script requires that you run the scripts 16S_Analysis_SILVA.R,
#ITS2_Analysis.R, PERMANOVA.R, and SoilVarDimensionalReduction.R
#and have the resulting data objects loaded in memory

library(tidyverse)
library(randomForest)
library(corrplot)
library(dplyr)
library(ggfortify)
library(ggplot2)
library(magrittr)


#####################First do random forests without sensor data#############################
#########################ITS DNA - without sensor data##############
#creates df with soil vars, and NMDS axes
rf_df_ITS_DNA=
  ITS_DNA_PCoA_Points %>%
  dplyr::select(sample,MDS1:MDS5) %>%
  left_join(select(AUE2021_alphadiv, sample, Shannon)) %>%
  left_join(Env_cluster_summary_df_full)

#runs random forest with NMDS 1 as response
rf_ITS_DNA_MDS1=randomForest(rf_df_ITS_DNA$MDS1 ~ ., 
                             data=select(rf_df_ITS_DNA, -c(MDS1:MDS5, Shannon, sample)),
                             importance=TRUE)

#runs random forest with NMDS 2 as response
rf_ITS_DNA_MDS2=randomForest(rf_df_ITS_DNA$MDS2 ~ ., 
                             data=select(rf_df_ITS_DNA, -c(MDS1:MDS5, Shannon, sample)),
                             importance=TRUE)

############################ITS RNA - without sensor data#########
#creates df with soil vars, and NMDS axes
rf_df_ITS_RNA=
  ITS_RNA_PCoA_Points %>%
  dplyr::select(sample,MDS1:MDS5) %>%
  left_join(select(AUE2021_alphadiv, sample, Shannon)) %>%
  left_join(Env_cluster_summary_df_full)

#runs random forest with NMDS 1 as response
rf_ITS_RNA_MDS1=randomForest(rf_df_ITS_RNA$MDS1 ~ ., 
                        data=select(rf_df_ITS_RNA, -c(MDS1:MDS5, Shannon, sample)),
                        importance=TRUE)

#runs random forest with NMDS 2 as response
rf_ITS_RNA_MDS2=randomForest(rf_df_ITS_RNA$MDS2 ~ ., 
                             data=select(rf_df_ITS_RNA, -c(MDS1:MDS5, Shannon, sample)),
                             importance=TRUE)



#########################16S DNA - without sensor data##############

#creates df with soil vars, and NMDS axes
rf_df_16S_DNA=
  PCoA_16S_DNA_Points %>%
  dplyr::select(sample,MDS1:MDS5) %>%
  left_join(select(AUE2021_16S_alphadiv, sample, Shannon)) %>%
  left_join(Env_cluster_summary_df_full)

#runs random forest with NMDS 1 as response
rf_16S_DNA_MDS1=randomForest(rf_df_16S_DNA$MDS1 ~ ., 
                             data=select(rf_df_16S_DNA, -c(MDS1:MDS5, Shannon, sample)),
                             importance=TRUE)

#runs random forest with NMDS 2 as response
rf_16S_DNA_MDS2=randomForest(rf_df_16S_DNA$MDS2 ~ ., 
                             data=select(rf_df_16S_DNA, -c(MDS1:MDS5, Shannon, sample)),
                             importance=TRUE)


############################16S RNA - without sensor data#########

#creates df with soil vars, and NMDS axes
rf_df_16S_RNA=
  PCoA_16S_RNA_Points %>%
  dplyr::select(sample,MDS1:MDS5) %>%
  left_join(select(AUE2021_16S_alphadiv, sample, Shannon)) %>%
  left_join(Env_cluster_summary_df_full)

#runs random forest with NMDS 1 as response
rf_16S_RNA_MDS1=randomForest(rf_df_16S_RNA$MDS1 ~ ., 
                             data=select(rf_df_16S_RNA, -c(MDS1:MDS5, Shannon, sample)),
                             importance=TRUE)

#runs random forest with NMDS 2 as response
rf_16S_RNA_MDS2=randomForest(rf_df_16S_RNA$MDS2 ~ ., 
                             data=select(rf_df_16S_RNA, -c(MDS1:MDS5, Shannon, sample)),
                             importance=TRUE)



####################Now do random forests with sensor data#############################
############ITS DNA - with sensor data

#first create df for random forest model
rf_df_ITS_sensor_DNA=
  ITS_DNA_PCoA_Points %>%
  dplyr::select(sample,MDS1:MDS5) %>%
  left_join(Env_cluster_summary_sensor_df_full) %>%
  filter(!is.na(meansoiltemp))

#run random forest
rf_ITS_sensor_DNA_MDS1=randomForest(rf_df_ITS_sensor_DNA$MDS1 ~ ., 
                                    data=select(rf_df_ITS_sensor_DNA, -c(MDS1:MDS5, sample)),
                                    importance=TRUE)

rf_ITS_sensor_DNA_MDS2=randomForest(rf_df_ITS_sensor_DNA$MDS2 ~ ., 
                                    data=select(rf_df_ITS_sensor_DNA, -c(MDS1:MDS5, sample)),
                                    importance=TRUE)

############ITS RNA - with sensor data
#first create df for random forest model
rf_df_ITS_sensor_RNA=
  ITS_RNA_PCoA_Points %>%
  dplyr::select(sample,MDS1:MDS5) %>%
  left_join(Env_cluster_summary_sensor_df_full) %>%
  filter(!is.na(meansoiltemp))

#run random forest
rf_ITS_sensor_RNA_MDS1=randomForest(rf_df_ITS_sensor_RNA$MDS1 ~ ., 
                                    data=select(rf_df_ITS_sensor_RNA, -c(MDS1:MDS5, sample)),
                                    importance=TRUE)

rf_ITS_sensor_RNA_MDS2=randomForest(rf_df_ITS_sensor_RNA$MDS2 ~ ., 
                                    data=select(rf_df_ITS_sensor_RNA, -c(MDS1:MDS5, sample)),
                                    importance=TRUE)

############16S DNA - with sensor data
#first create df for random forest model
rf_df_16S_sensor_DNA=
  PCoA_16S_DNA_Points %>%
  dplyr::select(sample,MDS1:MDS5) %>%
  left_join(Env_cluster_summary_sensor_df_full) %>%
  filter(!is.na(meansoiltemp))

#run random forest
rf_16S_sensor_DNA_MDS1=randomForest(rf_df_16S_sensor_DNA$MDS1 ~ ., 
                                    data=select(rf_df_16S_sensor_DNA, -c(MDS1:MDS5, sample)),
                                    importance=TRUE)

rf_16S_sensor_DNA_MDS2=randomForest(rf_df_16S_sensor_DNA$MDS2 ~ ., 
                                    data=select(rf_df_16S_sensor_DNA, -c(MDS1:MDS5, sample)),
                                    importance=TRUE)

############16S RNA - with sensor data
#first create df for random forest model
rf_df_16S_sensor_RNA=
  PCoA_16S_RNA_Points %>%
  dplyr::select(sample,MDS1:MDS5) %>%
  left_join(Env_cluster_summary_sensor_df_full) %>%
  filter(!is.na(meansoiltemp))

#run random forest
rf_16S_sensor_RNA_MDS1=randomForest(rf_df_16S_sensor_RNA$MDS1 ~ ., 
                                    data=select(rf_df_16S_sensor_RNA, -c(MDS1:MDS5, sample)),
                                    importance=TRUE)

#run random forest
rf_16S_sensor_RNA_MDS2=randomForest(rf_df_16S_sensor_RNA$MDS2 ~ ., 
                                    data=select(rf_df_16S_sensor_RNA, -c(MDS1:MDS5, sample)),
                                    importance=TRUE)

#########################Creating summary dataframes of random forest models##########
#Creates summary dataframes of random forest models to be used in supplement

#####first for MDS1
rf_16S_DNA_full_MDS1 <- 
  bind_rows(
    importance(rf_16S_DNA_MDS1) %>% 
      as.data.frame() %>%
      tibble::rownames_to_column("Variable"),
    
    importance(rf_16S_sensor_DNA_MDS1) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Variable") %>%
      filter(!(Variable %in% rownames(importance(rf_16S_DNA_MDS1))))
  ) %>%
  mutate(MDS1_DNA_16S=`%IncMSE`) %>% 
  select(Variable, MDS1_DNA_16S)

rf_16S_RNA_full_MDS1 <- 
  bind_rows(
    importance(rf_16S_RNA_MDS1) %>% 
      as.data.frame() %>%
      tibble::rownames_to_column("Variable"),
    
    importance(rf_16S_sensor_RNA_MDS1) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Variable") %>%
      filter(!(Variable %in% rownames(importance(rf_16S_RNA_MDS1))))
  ) %>%
  mutate(MDS1_RNA_16S=`%IncMSE`) %>% 
  select(Variable, MDS1_RNA_16S)

rf_ITS_DNA_full_MDS1 <- 
  bind_rows(
    importance(rf_ITS_DNA_MDS1) %>% 
      as.data.frame() %>%
      tibble::rownames_to_column("Variable"),
    
    importance(rf_ITS_sensor_DNA_MDS1) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Variable") %>%
      filter(!(Variable %in% rownames(importance(rf_ITS_DNA_MDS1))))
  ) %>%
  mutate(MDS1_DNA_ITS=`%IncMSE`) %>% 
  select(Variable, MDS1_DNA_ITS)

rf_ITS_RNA_full_MDS1 <- 
  bind_rows(
    importance(rf_ITS_RNA_MDS1) %>% 
      as.data.frame() %>%
      tibble::rownames_to_column("Variable"),
    
    importance(rf_ITS_sensor_RNA_MDS1) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Variable") %>%
      filter(!(Variable %in% rownames(importance(rf_ITS_RNA_MDS1))))
  )%>%
  mutate(MDS1_RNA_ITS=`%IncMSE`) %>% 
  select(Variable, MDS1_RNA_ITS)


rf_MDS1_full <- rf_ITS_DNA_full_MDS1 %>%
  left_join(rf_ITS_RNA_full_MDS1, by = "Variable") %>%
  left_join(rf_16S_DNA_full_MDS1, by = "Variable") %>%
  left_join(rf_16S_RNA_full_MDS1, by = "Variable") %>%
  mutate(across(
    .cols = -Variable,
    .fns = list(rank = ~ rank(-.x, ties.method = "average")),
    .names = "{.col}_rank"
  )) 


####Now do MDS2

rf_16S_DNA_full_MDS2 <- 
  bind_rows(
    importance(rf_16S_DNA_MDS2) %>% 
      as.data.frame() %>%
      tibble::rownames_to_column("Variable"),
    
    importance(rf_16S_sensor_DNA_MDS2) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Variable") %>%
      filter(!(Variable %in% rownames(importance(rf_16S_DNA_MDS2))))
  ) %>%
  mutate(MDS2_DNA_16S=`%IncMSE`) %>% 
  select(Variable, MDS2_DNA_16S)

rf_16S_RNA_full_MDS2 <- 
  bind_rows(
    importance(rf_16S_RNA_MDS2) %>% 
      as.data.frame() %>%
      tibble::rownames_to_column("Variable"),
    
    importance(rf_16S_sensor_RNA_MDS2) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Variable") %>%
      filter(!(Variable %in% rownames(importance(rf_16S_RNA_MDS2))))
  ) %>%
  mutate(MDS2_RNA_16S=`%IncMSE`) %>% 
  select(Variable, MDS2_RNA_16S)

rf_ITS_DNA_full_MDS2 <- 
  bind_rows(
    importance(rf_ITS_DNA_MDS2) %>% 
      as.data.frame() %>%
      tibble::rownames_to_column("Variable"),
    
    importance(rf_ITS_sensor_DNA_MDS2) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Variable") %>%
      filter(!(Variable %in% rownames(importance(rf_ITS_DNA_MDS2))))
  ) %>%
  mutate(MDS2_DNA_ITS=`%IncMSE`) %>% 
  select(Variable, MDS2_DNA_ITS)

rf_ITS_RNA_full_MDS2 <- 
  bind_rows(
    importance(rf_ITS_RNA_MDS2) %>% 
      as.data.frame() %>%
      tibble::rownames_to_column("Variable"),
    
    importance(rf_ITS_sensor_RNA_MDS2) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Variable") %>%
      filter(!(Variable %in% rownames(importance(rf_ITS_RNA_MDS2))))
  )%>%
  mutate(MDS2_RNA_ITS=`%IncMSE`) %>% 
  select(Variable, MDS2_RNA_ITS)


rf_MDS2_full <- rf_ITS_DNA_full_MDS2 %>%
  left_join(rf_ITS_RNA_full_MDS2, by = "Variable") %>%
  left_join(rf_16S_DNA_full_MDS2, by = "Variable") %>%
  left_join(rf_16S_RNA_full_MDS2, by = "Variable") %>%
  mutate(across(
    .cols = -Variable,
    .fns = list(rank = ~ rank(-.x, ties.method = "average")),
    .names = "{.col}_rank"
  )) 

#creates df with both MDS1 and MDS2 models
rf_MDS1_MDS2=left_join(rf_MDS1_full, rf_MDS2_full)

#saves as tsv
write_delim(rf_MDS1_MDS2, "rf_MDS1_MDS2.tsv", delim="\t")


