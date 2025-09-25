#This script tests for spatial autocorrelation in fungal and bacterial community composition
#while accounting for the effects of nutrients by doing partial mantel tests 
#that first remove variation associated with a dissimilarity matrix of soil properties
#
#The code here is similar to that used in SpatialAutocorrelation.R with the modification
#of making it a partial Mantel test instead of a traditional Mantel test
#to account for soil properties

#This script requires you to have run 16S_Analysis_SILVA.R,
#ITS_Analysis.R, SpatialDistanceMatrix.R and have the resulting data objects in memory

library(vegan) 

#generates distance matrices from environmental factors
#for each timepoint
Env_dm_June=
  Env_cluster_summary_df %>% 
  left_join(select(samples_df2_16S, sample, Time)) %>%
  left_join(SampleToPlot_16S) %>% 
  filter(Time=="June") %>%
  column_to_rownames("Plot") %>%
  dplyr::select(-Time, -sample) %>%
  metaMDSdist(autotransform=FALSE, distance = "manhattan")

Env_dm_August=
  Env_cluster_summary_df %>% 
  left_join(select(samples_df2_16S, sample, Time)) %>%
  left_join(SampleToPlot_16S) %>% 
  filter(Time=="August") %>%
  column_to_rownames("Plot") %>%
  dplyr::select(-Time, -sample) %>%
  metaMDSdist(autotransform=FALSE, distance = "manhattan")

Env_dm_October=
  Env_cluster_summary_df %>% 
  left_join(select(samples_df2_16S, sample, Time)) %>%
  left_join(SampleToPlot_16S) %>% 
  filter(Time=="October") %>%
  column_to_rownames("Plot") %>%
  dplyr::select(-Time, -sample) %>%
  metaMDSdist(autotransform=FALSE, distance = "manhattan")

###########################################
# Define combinations
sources <- c("DNA", "cDNA")
times <- c("June", "August", "October")

rename_dist_matrix <- function(dist_obj, lookup_df, sample_col, plot_col) {
  if (!inherits(dist_obj, "dist")) {
    stop("Input dist_obj must be of class 'dist'")
  }
  
  labels <- attr(dist_obj, "Labels")
  lookup_df <- lookup_df %>%
    dplyr::filter(!!rlang::sym(sample_col) %in% labels) %>%
    dplyr::distinct(!!rlang::sym(sample_col), !!rlang::sym(plot_col))
  
  # Make named vector: sample name -> plot name
  name_map <- lookup_df %>%
    tibble::deframe()
  
  # Rename labels
  renamed_labels <- name_map[labels]
  
  # If any labels are NA (missing match), keep originals
  renamed_labels[is.na(renamed_labels)] <- labels[is.na(renamed_labels)]
  
  attr(dist_obj, "Labels") <- renamed_labels
  return(dist_obj)
}


# Function: Run partial Mantel test
run_partial_mantel_test <- function(comm_dm, spatial_matrix, control_matrix) {
  if (!inherits(comm_dm, "dist") || !inherits(control_matrix, "dist")) {
    stop("comm_dm and control_matrix must be 'dist' objects.")
  }
  if (!inherits(spatial_matrix, "matrix")) {
    stop("spatial_matrix must be a matrix.")
  }
  
  # Match sample/plot names
  common_ids <- attr(comm_dm, "Labels")
  spatial_subset <- spatial_matrix[common_ids, common_ids]
  control_matrix_subset <- as.matrix(control_matrix)[common_ids, common_ids]
  
  # Convert to dist objects
  spatial_dist <- as.dist(spatial_subset)
  control_dist <- as.dist(control_matrix_subset)
  
  # Run partial Mantel test
  mantel.partial(comm_dm, spatial_dist, control_dist, method = "pearson")
}

# ITS: Partial Mantel Tests (spatial | env)
mantel_partial_ITS_results <- list()
for (s in sources) {
  for (t in times) {
    combo_name <- paste0("ITS_", s, "_", t)
    
    sample_df <- as(sample_data(AUE2021_ITS_rarefied), "data.frame")
    samples_to_keep <- rownames(sample_df[sample_df$source == s & sample_df$Time == t, ])
    phy_sub <- prune_samples(samples_to_keep, AUE2021_ITS_rarefied)
    
    if (nsamples(phy_sub) < 2) next
    
    comm_dm <- phy_sub %>%
      otu_table() %>%
      as("matrix") %>%
      t() %>%
      as.data.frame() %>%
      metaMDSdist(distance = "bray")
    
    comm_dm_named <- rename_dist_matrix(
      comm_dm,
      SampleToPlot,
      sample_col = "sample",
      plot_col = "Plot"
    )
    
    control_env_dm <- switch(t,
                             "June" = Env_dm_June,
                             "August" = Env_dm_August,
                             "October" = Env_dm_October)
    
    mantel_partial_ITS_results[[combo_name]] <- run_partial_mantel_test(
      comm_dm = comm_dm_named,
      spatial_matrix = DistanceMatrixMatrix,
      control_matrix = control_env_dm
    )
  }
}

# 16S: Partial Mantel Tests (spatial | env)
mantel_partial_16S_results <- list()
for (s in sources) {
  for (t in times) {
    combo_name <- paste0("16S_", s, "_", t)
    
    sample_df <- as(sample_data(AUE2021_16S_rarefied), "data.frame")
    samples_to_keep <- rownames(sample_df[sample_df$source == s & sample_df$Time == t, ])
    phy_sub <- prune_samples(samples_to_keep, AUE2021_16S_rarefied)
    
    if (nsamples(phy_sub) < 2) next
    
    comm_dm <- phy_sub %>%
      otu_table() %>%
      as("matrix") %>%
      t() %>%
      as.data.frame() %>%
      metaMDSdist(distance = "bray")
    
    comm_dm_named <- rename_dist_matrix(
      comm_dm,
      SampleToPlot_16S,
      sample_col = "sample",
      plot_col = "Plot"
    )
    
    control_env_dm <- switch(t,
                             "June" = Env_dm_June,
                             "August" = Env_dm_August,
                             "October" = Env_dm_October)
    
    mantel_partial_16S_results[[combo_name]] <- run_partial_mantel_test(
      comm_dm = comm_dm_named,
      spatial_matrix = DistanceMatrixMatrix,
      control_matrix = control_env_dm
    )
  }
}

partial_summary_ITS <- bind_rows(lapply(names(mantel_partial_ITS_results), function(name) {
  result <- mantel_partial_ITS_results[[name]]
  if (!is.null(result)) {
    tibble(
      Dataset = name,
      Mantel_r = result$statistic,
      p_value = result$signif
    )
  }
}), .id = "row_id")

partial_summary_16S <- bind_rows(lapply(names(mantel_partial_16S_results), function(name) {
  result <- mantel_partial_16S_results[[name]]
  if (!is.null(result)) {
    tibble(
      Dataset = name,
      Mantel_r = result$statistic,
      p_value = result$signif
    )
  }
}), .id = "row_id")

