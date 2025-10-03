#random forest modelling of top fungal and bacterial taxa
#with inclusion of soil sensor data

library(tidyverse)
library(ggrepel)
library(ggplot2)
library(randomForest)
library(dplyr)
library(tidyr)
library(car)
library(lme4)
library(lmerTest)
library(purrr)
library(broom.mixed)

#################Random forest models with sensor data############################
#####################################################################################

#creates data frame w soil data for PCA with sensor data 
EnvVars_sub_sensor=
  samples_df_16S %>% 
  column_to_rownames("sample") %>%
  filter(source=="DNA",) %>%  #filters to only include DNA samples so that metadata is not redundant
  filter(!is.na(mean_soiltemp)) %>% #filters rows that do not have sensor data
  select(NO3_PRS:Cd_PRS,
         GrassCover:ForbCover,
         pH_UGA, BaseSat_UGA:Zn_UGA, N_UGA, 
         TOC_UGA, MAC_moist, mean_moist,
         MAC_soiltemp, mean_soiltemp) %>%
  select(-c(Cr_UGA, Cu_UGA, Mo_UGA, Pb_UGA)) %>% # removes variables with trace concentrations
  select(-c(B_PRS, Pb_PRS, Al_PRS, Zn_PRS, Mn_PRS,
            Cd_PRS, Cd_UGA, Ni_UGA, Zn_UGA, Mn_UGA, 
            MAC_soiltemp)) %>% #removes variables with little biological relevance
  select(where(~ !any(is.na(.)))) %>% #removes columns with NAs
  rename_with(~ str_replace_all(., c("PRS" = "Flux", "UGA" = "Conc"))) %>% #renames column headers as flux and concentration for PRS strips and bulk soil tests respectively
  rename_with(~ str_replace_all(., "_", ""))

###########Correlation based clustering of variables#########
# Compute correlation matrix
cor_mat_sensor <- cor(EnvVars_sub_sensor, use = "pairwise.complete.obs")

# Convert to distance matrix: distance = 1 - abs(correlation)
dist_mat_sensor <- as.dist(1 - abs(cor_mat_sensor))

# Hierarchical clustering
hc_sensor <- hclust(dist_mat_sensor, method = "average")  # or "complete", "ward.D2", etc.

# Plot dendrogram
plot(hc_sensor, main = "Hierarchical Clustering of Variables")

# Cut tree at desired correlation threshold (e.g., 0.5 → height = 1 - 0.8 = 0.2)
clusters_sensor <- cutree(hc_sensor, h = 0.5)

# Assign each variable to a cluster.
#Clustering of variables is forced to follow the clustering
#determined on the full dataset (i.e. containing rows
#that did not have sensor data). For the sensor data, 
#de novo clustering is performed
cluster_sensor_df <- data.frame(Variable = names(clusters_sensor), ClusterNoSens = clusters_sensor) %>%
  left_join(cluster_df) %>%
  mutate(Cluster=
           case_when(is.na(Cluster) ~ clusters_sensor+1000,
                     .default=Cluster))

#############################################################
#Calculate linear combinations of variables in each cluster.
#clusters with 1 variable are left untransformed

# Function to summarize each cluster
summarize_clusters <- function(data, cluster_df, method = c("pca", "mean")) {
  method <- match.arg(method)
  
  out_list <- list()
  
  for (clust in unique(cluster_df$Cluster)) {
    vars <- cluster_df %>% filter(Cluster == clust) %>% pull(Variable)
    
    # Create a name from the variables, joined by underscores
    cluster_name <- paste(vars, collapse = "_")
    
    if (length(vars) == 1) {
      out_list[[cluster_name]] <- data[[vars]]
    } else {
      sub_data <- data[, vars, drop = FALSE]
      if (method == "pca") {
        pc1 <- prcomp(sub_data, scale. = TRUE)$x[, 1]
        out_list[[cluster_name]] <- pc1
      } else if (method == "mean") {
        mean_score <- rowMeans(sub_data, na.rm = TRUE)
        out_list[[cluster_name]] <- mean_score
      }
    }
  }
  
  # Combine into data frame and preserve rownames
  result_df <- as.data.frame(out_list)
  rownames(result_df) <- rownames(data)
  
  return(result_df)
}


# Run the function to calculate reduced dimensionality dataset
Env_cluster_summary_sensor_df <- 
  summarize_clusters(EnvVars_sub_sensor, cluster_sensor_df, method = "pca") %>%
  rownames_to_column(var = "sample")

#formats the clustered env variable dataset so that DNA
#and cDNA are separate rows (with same values)
Env_cluster_summary_sensor_df_full=
  rbind(Env_cluster_summary_sensor_df,
        mutate(Env_cluster_summary_sensor_df, sample = paste0("c", sample)))

EnvVars_sub_sensor_full=
  mutate(EnvVars_sub_sensor, sample=rownames(EnvVars_sub_sensor)) %>%
  rbind(mutate(EnvVars_sub_sensor, sample = paste0("c", rownames(EnvVars_sub_sensor)))) %>%
  { `rownames<-`(., NULL) }

############ITS DNA Dataset modelling MDS axes#####
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

############ITS RNA Dataset modelling MDS axes#####
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

############16S DNA Dataset modelling MDS axes#####
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

############16S RNA Dataset modelling MDS axes#####
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

rf_16S_DNA_full_MDS2 <- 
  bind_rows(
    importance(rf_16S_DNA_MDS2) %>% 
      as.data.frame() %>%
      tibble::rownames_to_column("Variable"),
    
    importance(rf_16S_sensor_DNA_MDS2) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Variable") %>%
      filter(!(Variable %in% rownames(importance(rf_16S_DNA))))
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
      filter(!(Variable %in% rownames(importance(rf_16S_RNA))))
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
      filter(!(Variable %in% rownames(importance(rf_ITS_DNA))))
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
      filter(!(Variable %in% rownames(importance(rf_ITS_RNA))))
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


RankCutoff=5
ImpoCutoff=10
rf_MDS1_MDS2_Full =
  left_join(rf_MDS1_full, rf_MDS2_full) %>%
  rowwise() %>%
  mutate(
    CountImportance = sum(c_across(contains("MDS1_") & !ends_with("rank")) > ImpoCutoff |
                            c_across(contains("MDS2_") & !ends_with("rank")) > ImpoCutoff),
    
    CountRank = sum(c_across(contains("MDS1_") & ends_with("rank")) < RankCutoff |
                      c_across(contains("MDS2_") & ends_with("rank")) < RankCutoff)
  ) %>%
  ungroup()


##################Random forest on specific Taxa########
#generate list of predictor variables
Rf_vars_sensor=
  select(Env_cluster_summary_sensor_df_full, -sample) %>%
  colnames

#Start by formatting taxon summary dfs
AUE2021_ITS_rarefied_class_sensor_t=
  AUE2021_ITS_rarefied_class %>%
  t %>% 
  as.data.frame %>%
  mutate(sample=rownames(.)) %>%
  left_join(Env_cluster_summary_sensor_df_full) %>%
  left_join(select(samples_df2_16S,source, sample, Site, SiteType, Time, Plot)) %>%
  filter(!is.na(MACmoist))

AUE2021_ITS_rarefied_class_sensor_t_DNA=
  filter(AUE2021_ITS_rarefied_class_sensor_t, source=="DNA")

AUE2021_ITS_rarefied_class_sensor_t_RNA=
  filter(AUE2021_ITS_rarefied_class_sensor_t, source=="cDNA")

#Start by formatting taxon summary dfs
AUE2021_16S_rarefied_phylum_sensor_t=
  AUE2021_16S_rarefied_phylum %>%
  t %>% 
  as.data.frame %>%
  mutate(sample=rownames(.)) %>%
  left_join(Env_cluster_summary_sensor_df_full) %>%
  left_join(select(samples_df2_16S,source, sample, Site, SiteType, Time, Plot)) %>%
  filter(!is.na(MACmoist))

AUE2021_16S_rarefied_phylum_sensor_t_DNA=
  filter(AUE2021_16S_rarefied_phylum_sensor_t, source=="DNA")

AUE2021_16S_rarefied_phylum_sensor_t_RNA=
  filter(AUE2021_16S_rarefied_phylum_sensor_t, source=="cDNA")

#####################ITS DNA Tests#######
rf_results_ITS_sensor_DNA <- map(DomClass_ITS, function(response_var) {
  formula <- as.formula(paste(response_var, "~ ."))
  
  rf_model <- randomForest(
    formula,
    data = select(AUE2021_ITS_rarefied_class_sensor_t_DNA, all_of(c(response_var, Rf_vars_sensor))),
    importance = TRUE
  )
  
  importance_df <- as.data.frame(importance(rf_model)) %>%
    rownames_to_column("Predictor") %>%
    select(Predictor, `%IncMSE`) %>%
    mutate(Response = response_var)
})

# Combine to long format
rf_importance_long_ITS_sensor_DNA <- bind_rows(rf_results_ITS_sensor_DNA)

# Pivot to wide format
rf_importance_wide_ITS_sensor_DNA <- rf_importance_long_ITS_sensor_DNA %>%
  pivot_wider(names_from = Response, values_from = `%IncMSE`)

# Now run lmer models
lmer_results_ITS_sensor_DNA <- map_dfr(DomClass_ITS, function(response_var) {
  map_dfr(Rf_vars_sensor, function(predictor_var) {
    formula <- as.formula(paste(response_var, "~", predictor_var, "+ (1|Site)"))
    
    model <- tryCatch(
      lmer(formula, data = AUE2021_ITS_rarefied_class_sensor_t_DNA),
      error = function(e) NULL
    )
    
    if (!is.null(model)) {
      coefs <- summary(model)$coefficients
      if (predictor_var %in% rownames(coefs)) {
        est <- coefs[predictor_var, "Estimate"]
        se <- coefs[predictor_var, "Std. Error"]
        tval <- coefs[predictor_var, "t value"]
        pval <- coefs[predictor_var, "Pr(>|t|)"]
      } else {
        est <- se <- tval <- pval <- NA
      }
      tibble(
        Response = response_var,
        Predictor = predictor_var,
        estimate = est,
        std.error = se,
        statistic = tval,
        p.value = pval
      )
    } else {
      tibble(
        Response = response_var,
        Predictor = predictor_var,
        estimate = NA,
        std.error = NA,
        statistic = NA,
        p.value = NA
      )
    }
  })
})

# Join with random forest results
combined_results_ITS_sensor_DNA <- left_join(
  rf_importance_long_ITS_sensor_DNA,
  lmer_results_ITS_sensor_DNA,
  by = c("Response", "Predictor")
)


#########################ITS RNA Tests#######
rf_results_ITS_sensor_RNA <- map(DomClass_ITS, function(response_var) {
  formula <- as.formula(paste(response_var, "~ ."))
  
  rf_model <- randomForest(
    formula,
    data = select(AUE2021_ITS_rarefied_class_sensor_t_RNA, all_of(c(response_var, Rf_vars_sensor))),
    importance = TRUE
  )
  
  importance_df <- as.data.frame(importance(rf_model)) %>%
    rownames_to_column("Predictor") %>%
    select(Predictor, `%IncMSE`) %>%
    mutate(Response = response_var)
})

# Combine to long format
rf_importance_long_ITS_sensor_RNA <- bind_rows(rf_results_ITS_sensor_RNA)

# Pivot to wide format
rf_importance_wide_ITS_sensor_RNA <- rf_importance_long_ITS_sensor_RNA %>%
  pivot_wider(names_from = Response, values_from = `%IncMSE`)

# Now run lmer models
lmer_results_ITS_sensor_RNA <- map_dfr(DomClass_ITS, function(response_var) {
  map_dfr(Rf_vars_sensor, function(predictor_var) {
    formula <- as.formula(paste(response_var, "~", predictor_var, "+ (1|Site)"))
    
    model <- tryCatch(
      lmer(formula, data = AUE2021_ITS_rarefied_class_sensor_t_RNA),
      error = function(e) NULL
    )
    
    if (!is.null(model)) {
      coefs <- summary(model)$coefficients
      if (predictor_var %in% rownames(coefs)) {
        est <- coefs[predictor_var, "Estimate"]
        se <- coefs[predictor_var, "Std. Error"]
        tval <- coefs[predictor_var, "t value"]
        pval <- coefs[predictor_var, "Pr(>|t|)"]
      } else {
        est <- se <- tval <- pval <- NA
      }
      tibble(
        Response = response_var,
        Predictor = predictor_var,
        estimate = est,
        std.error = se,
        statistic = tval,
        p.value = pval
      )
    } else {
      tibble(
        Response = response_var,
        Predictor = predictor_var,
        estimate = NA,
        std.error = NA,
        statistic = NA,
        p.value = NA
      )
    }
  })
})

# Join with random forest results
combined_results_ITS_sensor_RNA <- left_join(
  rf_importance_long_ITS_sensor_RNA,
  lmer_results_ITS_sensor_RNA,
  by = c("Response", "Predictor")
)


#####################16S DNA Tests#######
rf_results_16S_sensor_DNA <- map(DomPhyla_16S, function(response_var) {
  formula <- as.formula(paste(response_var, "~ ."))
  
  rf_model <- randomForest(
    formula,
    data = select(AUE2021_16S_rarefied_phylum_sensor_t_DNA, all_of(c(response_var, Rf_vars_sensor))),
    importance = TRUE
  )
  
  importance_df <- as.data.frame(importance(rf_model)) %>%
    rownames_to_column("Predictor") %>%
    select(Predictor, `%IncMSE`) %>%
    mutate(Response = response_var)
})

# Combine to long format
rf_importance_long_16S_sensor_DNA <- bind_rows(rf_results_16S_sensor_DNA)

# Pivot to wide format
rf_importance_wide_16S_sensor_DNA <- rf_importance_long_16S_sensor_DNA %>%
  pivot_wider(names_from = Response, values_from = `%IncMSE`)

# Now run lmer models
lmer_results_16S_sensor_DNA <- map_dfr(DomPhyla_16S, function(response_var) {
  map_dfr(Rf_vars_sensor, function(predictor_var) {
    formula <- as.formula(paste(response_var, "~", predictor_var, "+ (1|Site)"))
    
    model <- tryCatch(
      lmer(formula, data = AUE2021_16S_rarefied_phylum_sensor_t_DNA),
      error = function(e) NULL
    )
    
    if (!is.null(model)) {
      coefs <- summary(model)$coefficients
      if (predictor_var %in% rownames(coefs)) {
        est <- coefs[predictor_var, "Estimate"]
        se <- coefs[predictor_var, "Std. Error"]
        tval <- coefs[predictor_var, "t value"]
        pval <- coefs[predictor_var, "Pr(>|t|)"]
      } else {
        est <- se <- tval <- pval <- NA
      }
      tibble(
        Response = response_var,
        Predictor = predictor_var,
        estimate = est,
        std.error = se,
        statistic = tval,
        p.value = pval
      )
    } else {
      tibble(
        Response = response_var,
        Predictor = predictor_var,
        estimate = NA,
        std.error = NA,
        statistic = NA,
        p.value = NA
      )
    }
  })
})

# Join with random forest results
combined_results_16S_sensor_DNA <- left_join(
  rf_importance_long_16S_sensor_DNA,
  lmer_results_16S_sensor_DNA,
  by = c("Response", "Predictor")
)


##################16S RNA Tests###############
rf_results_16S_sensor_RNA <- map(DomPhyla_16S, function(response_var) {
  formula <- as.formula(paste(response_var, "~ ."))
  
  rf_model <- randomForest(
    formula,
    data = select(AUE2021_16S_rarefied_phylum_sensor_t_RNA, all_of(c(response_var, Rf_vars_sensor))),
    importance = TRUE
  )
  
  importance_df <- as.data.frame(importance(rf_model)) %>%
    rownames_to_column("Predictor") %>%
    select(Predictor, `%IncMSE`) %>%
    mutate(Response = response_var)
})

# Combine to long format
rf_importance_long_16S_sensor_RNA <- bind_rows(rf_results_16S_sensor_RNA)

# Pivot to wide format
rf_importance_wide_16S_sensor_RNA <- rf_importance_long_16S_sensor_RNA %>%
  pivot_wider(names_from = Response, values_from = `%IncMSE`)

# Now run lmer models
lmer_results_16S_sensor_RNA <- map_dfr(DomPhyla_16S, function(response_var) {
  map_dfr(Rf_vars_sensor, function(predictor_var) {
    formula <- as.formula(paste(response_var, "~", predictor_var, "+ (1|Site)"))
    
    model <- tryCatch(
      lmer(formula, data = AUE2021_16S_rarefied_phylum_sensor_t_RNA),
      error = function(e) NULL
    )
    
    if (!is.null(model)) {
      coefs <- summary(model)$coefficients
      if (predictor_var %in% rownames(coefs)) {
        est <- coefs[predictor_var, "Estimate"]
        se <- coefs[predictor_var, "Std. Error"]
        tval <- coefs[predictor_var, "t value"]
        pval <- coefs[predictor_var, "Pr(>|t|)"]
      } else {
        est <- se <- tval <- pval <- NA
      }
      tibble(
        Response = response_var,
        Predictor = predictor_var,
        estimate = est,
        std.error = se,
        statistic = tval,
        p.value = pval
      )
    } else {
      tibble(
        Response = response_var,
        Predictor = predictor_var,
        estimate = NA,
        std.error = NA,
        statistic = NA,
        p.value = NA
      )
    }
  })
})

# Join with random forest results
combined_results_16S_sensor_RNA <- left_join(
  rf_importance_long_16S_sensor_RNA,
  lmer_results_16S_sensor_RNA,
  by = c("Response", "Predictor")
)



