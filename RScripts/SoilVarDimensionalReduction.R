#This script uses PCA to perform dimensional reduction on soil variables. This
#reduced dimensionality dataset is used in a few other analyses: 
#Random forests exploring effects of environmental variables on microbial communities
#Spatial autocorellation tests (partial mantel tests accounting for soil factors)

#Some of the field plots did not have soil sensors installed (sensors were installed at 
#3 out of 5 plots per site). So, we generate two types of reduced dimensionality 
#soil variable dataframes: 1) without sensor data, but retaining all plots and 
#2) with sensor data, but discarding plots that did not have sensors.

#This script requires you to have run 16S_Analysis_SILVA.R, SoilSensorDataProcessing.R,
#Import_IonExchangeData.R, MeanAbsoluteChange.R, and Vegsurvey_SoilTest_DataImport.R 
#and have the resulting data objects loaded in memory

library(tidyverse)

######################################################################################
###############Dimension reduction WITHOUT sensor data####################################

#creates data frame w soil data for PCA without sensor data because
#of missing data (i.e. sensors were only installed at 3 out of 5 plots at each site)
#selects relevant variables
EnvVars_sub_nosensor=
  samples_df_16S %>% 
  column_to_rownames("sample") %>%
  filter(source=="DNA") %>%  #filters to only include DNA samples so that metadata is not redundant
  select(NO3_PRS:Cd_PRS,
         GrassCover:ForbCover,
         pH_UGA, BaseSat_UGA:Zn_UGA, N_UGA, 
         TOC_UGA, MAC_moist, mean_moist,
         MAC_soiltemp, mean_soiltemp) %>%
  select(-c(Cr_UGA, Cu_UGA, Mo_UGA, Pb_UGA)) %>% # removes variables with trace concentrations
  select(-c(B_PRS, Pb_PRS, Al_PRS, Zn_PRS, Mn_PRS,
            Cd_PRS, Cd_UGA, Ni_UGA, Zn_UGA, Mn_UGA)) %>% #removes variables with little biological relevance
  select(where(~ !any(is.na(.)))) %>% #removes columns with NAs (i.e. sensor data)
  rename_with(~ str_replace_all(., c("PRS" = "Flux", "UGA" = "Conc"))) %>% #renames column headers as flux and concentration for PRS strips and bulk soil tests respectively
  rename_with(~ str_replace_all(., "_", ""))

###########Correlation based clustering of variables#########
#a correlation matrix is made between all variables. 
#This matrix is then used to do hierarchical clustering.
#A threshold of 0.5 is set to "cut" clusters of correlated variables.
#Variables that are in the same cluster are then collapsed into index variables using PCA
#These clusters are named with the names of all variables in them separated by "_"

# Compute correlation matrix
cor_mat <- cor(EnvVars_sub_nosensor, use = "pairwise.complete.obs")

#Convert to distance matrix: distance = 1 - abs(correlation)
dist_mat <- as.dist(1 - abs(cor_mat))

#Hierarchical clustering
hc <- hclust(dist_mat, method = "average")  # or "complete", "ward.D2", etc.

# Plot dendrogram
plot(hc, main = "Hierarchical Clustering of Variables")

# Cut tree at correlation threshold of 0.5
clusters <- cutree(hc, h = 0.5)

# Assign each variable to a cluster
cluster_df <- data.frame(Variable = names(clusters), Cluster = clusters)

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
Env_cluster_summary_df <- 
  summarize_clusters(EnvVars_sub_nosensor, cluster_df, method = "pca") %>%
  rownames_to_column(var = "sample")

#formats the clustered env variable dataset so that DNA
#and cDNA are separate rows (with same values)
Env_cluster_summary_df_full=
  rbind(Env_cluster_summary_df,
        mutate(Env_cluster_summary_df, sample = paste0("c", sample)))

#formats the raw unclustered env variable dataset so that DNA
#and cDNA are separate rows (with same values)
EnvVars_sub_nosensor_full=
  mutate(EnvVars_sub_nosensor, sample=rownames(EnvVars_sub_nosensor)) %>%
  rbind(mutate(EnvVars_sub_nosensor, sample = paste0("c", rownames(EnvVars_sub_nosensor)))) %>%
  { `rownames<-`(., NULL) }

######################################################################################

######################################################################################
###############Dimension reduction WITH sensor data####################################


#creates data frame w soil data for PCA with sensor data because
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

######################################################################################