#This script conducts tests of spatial autocorrelation with fungal and bacterial
#community composition. Are communities closer together more similar??
#
#This script requires you to have run 16S_Analysis_SILVA.R,
#ITS_Analysis.R, SpatialDistanceMatrix.R and have the resulting data objects in memory

library(vegan)
library(tidyverse)
library(phyloseq)

#creates file linking sample ID to plot ID for ITS dataset
SampleToPlot_ITS=
  read.delim("./InputFiles/Metadata_ITS.txt") %>%
  select(sample, Plot)

#creates file linking sample ID to plot ID for 16S dataset
SampleToPlot_16S=
  read.delim("./InputFiles/Metadata_16S.txt") %>%
  select(sample, Plot)

#function to rename distance object row/columns
rename_dist_matrix <- function(dist_obj, mapping_df, sample_col = "sample", plot_col = "Plot") {
  # Convert the distance object to a matrix
  dist_mat <- as.matrix(dist_obj)
  
  # Create a named vector from sample to plot
  sample_to_plot <- setNames(mapping_df[[plot_col]], mapping_df[[sample_col]])
  
  # Check all sample names exist in the mapping
  if (!all(rownames(dist_mat) %in% names(sample_to_plot))) {
    stop("Not all row names in the distance matrix are present in the mapping.")
  }
  if (!all(colnames(dist_mat) %in% names(sample_to_plot))) {
    stop("Not all column names in the distance matrix are present in the mapping.")
  }
  
  # Rename rows and columns
  rownames(dist_mat) <- sample_to_plot[rownames(dist_mat)]
  colnames(dist_mat) <- sample_to_plot[colnames(dist_mat)]
  
  # Convert back to a dist object (optional)
  return(as.dist(dist_mat))
}

#custom wrapper function for Mantel tests
run_mantel_test <- function(community_dm, spatial_matrix) {
  # Ensure 'community_dm' is a dist object
  if (!inherits(community_dm, "dist")) {
    stop("The community dissimilarity matrix must be a 'dist' object.")
  }
  
  # Extract labels from community dissimilarity matrix
  common_ids <- attr(community_dm, "Labels")
  
  # Subset and reorder the spatial matrix
  spatial_matrix_subset <- spatial_matrix[common_ids, common_ids]
  
  # Check if any NA values were introduced (indicates mismatch)
  if (any(is.na(spatial_matrix_subset))) {
    stop("Mismatch between sample/plot names in community_dm and spatial_matrix.")
  }
  
  # Convert spatial matrix to dist object
  spatial_dist <- as.dist(spatial_matrix_subset)
  
  # Run the Mantel test
  mantel(community_dm, spatial_dist, method = "pearson")
}

#takes community dissimilarity matrix and spatial dissimilarity matrix
#and formats them for plotting
prep_mantel_plot_df <- function(comm_dm_ITS, spatial_dm) {
  # Convert distance objects to matrices if needed
  if (inherits(comm_dm_ITS, "dist")) {
    comm_dm_ITS <- as.matrix(comm_dm_ITS)
  }
  
  if (inherits(spatial_dm, "dist")) {
    spatial_dm <- as.matrix(spatial_dm)
  }
  
  # Check that both matrices are square and have matching sample names
  common_samples <- intersect(rownames(comm_dm_ITS), rownames(spatial_dm))
  
  comm_dm_ITS <- comm_dm_ITS[common_samples, common_samples]
  spatial_dm <- spatial_dm[common_samples, common_samples]
  
  # Convert to long data frame for plotting
  comm_df <- as.data.frame(as.table(comm_dm_ITS)) %>%
    rename(Sample1 = Var1, Sample2 = Var2, CommunityDist = Freq)
  
  spatial_df <- as.data.frame(as.table(spatial_dm)) %>%
    rename(Sample1 = Var1, Sample2 = Var2, SpatialDist = Freq)
  
  # Join and filter out self comparisons
  mantel_df <- left_join(comm_df, spatial_df, by = c("Sample1", "Sample2")) %>%
    filter(Sample1 != Sample2)
  
  return(mantel_df)
}

###########################################
# Define combinations
sources <- c("DNA", "cDNA")
times <- c("June", "August", "October")

# Store results
mantel_ITS_dfs <- list()
spatial_ITS_dfs <- list()

for (s in sources) {
  for (t in times) {
    
    combo_name_ITS <- paste0("ITS_", s, "_", t)
    
    # Extract sample data as data frame
    sample_df_ITS <- as(sample_data(AUE2021_ITS_rarefied), "data.frame")
    
    # Identify samples with matching source and time
    samples_to_keep_ITS <- rownames(sample_df_ITS[sample_df_ITS$source == s & sample_df_ITS$Time == t, ])
    
    # Subset the phyloseq object
    phy_sub_ITS <- prune_samples(samples_to_keep_ITS, AUE2021_ITS_rarefied)
    
    # Skip if not enough samples
    if (nsamples(phy_sub_ITS) < 2) {
      warning(paste("Skipping", combo_name_ITS, "- not enough samples."))
      next
    }
    
    # Compute community distance matrix
    comm_dm_ITS <- phy_sub_ITS %>%
      otu_table() %>%
      as("matrix") %>%
      t() %>%
      as.data.frame() %>%
      metaMDSdist(distance = "bray")
    
    # Rename row/col names
    comm_dm_ITS_named <- rename_dist_matrix(
      comm_dm_ITS,
      SampleToPlot,
      sample_col = "sample",
      plot_col = "Plot"
    )
    
    # Run Mantel test
    mantel_ITS_dfs[[combo_name_ITS]] <- run_mantel_test(
      community_dm = comm_dm_ITS_named,
      spatial_matrix = DistanceMatrixMatrix
    )
    
    # Prepare data frame for plotting
    spatial_ITS_dfs[[combo_name_ITS]] <- prep_mantel_plot_df(
      comm_dm = comm_dm_ITS_named,
      spatial_dm = DistanceMatrixMatrix
    ) %>%
      mutate(DistClass = case_when(
        SpatialDist < 100 ~ "WithinSite",
        .default = "BetweenSite"
      )) %>%
      mutate(DistClass = factor(DistClass, levels = c("WithinSite", "BetweenSite"))) %>%
      mutate(source=s,
             Time=t)
  }
}

# Store results
mantel_16S_dfs <- list()
spatial_16S_dfs <- list()

for (s in sources) {
  for (t in times) {
    
    combo_name_16S <- paste0("16S_", s, "_", t)
    
    # Extract sample data as data frame
    sample_df_16S <- as(sample_data(AUE2021_16S_rarefied), "data.frame")
    
    # Identify samples with matching source and time
    samples_to_keep_16S <- rownames(sample_df_16S[sample_df_16S$source == s & sample_df_16S$Time == t, ])
    
    # Subset the phyloseq object
    phy_sub_16S <- prune_samples(samples_to_keep_16S, AUE2021_16S_rarefied)
    
    # Skip if not enough samples
    if (nsamples(phy_sub_16S) < 2) {
      warning(paste("Skipping", combo_name_16S, "- not enough samples."))
      next
    }
    
    # Compute community distance matrix
    comm_dm_16S <- phy_sub_16S %>%
      otu_table() %>%
      as("matrix") %>%
      t() %>%
      as.data.frame() %>%
      metaMDSdist(distance = "bray")
    
    # Rename row/col names
    comm_dm_16S_named <- rename_dist_matrix(
      comm_dm_16S,
      SampleToPlot_16S,
      sample_col = "sample",
      plot_col = "Plot"
    )
    
    # Run Mantel test
    mantel_16S_dfs[[combo_name_16S]] <- run_mantel_test(
      community_dm = comm_dm_16S_named,
      spatial_matrix = DistanceMatrixMatrix
    )
    
    # Prepare data frame for plotting
    spatial_16S_dfs[[combo_name_16S]] <- prep_mantel_plot_df(
      comm_dm = comm_dm_16S_named,
      spatial_dm = DistanceMatrixMatrix
    ) %>%
      mutate(DistClass = case_when(
        SpatialDist < 100 ~ "WithinSite",
        .default = "BetweenSite"
      )) %>%
      mutate(DistClass = factor(DistClass, levels = c("WithinSite", "BetweenSite"))) %>%
      mutate(source=s,
             Time=t)
  }
}
#######################Plotting spatial autocorrelation relationships####################

#Fungal ITS2 - binds together dfs of community dissimilarity vs spatial distance
#for DNA/RNA and all timepoints so it can be plotted
Full_Spatial_df_ITS=
  rbind(spatial_ITS_dfs$ITS_DNA_June,
        spatial_ITS_dfs$ITS_DNA_August,
        spatial_ITS_dfs$ITS_DNA_October,
        spatial_ITS_dfs$ITS_cDNA_June,
        spatial_ITS_dfs$ITS_cDNA_August,
        spatial_ITS_dfs$ITS_cDNA_October) %>%
  mutate(source=factor(source, levels=c("DNA", "cDNA")),
         Time=factor(Time, levels=c("June", "August", "October")))

#Bacterial 16S - binds together dfs of community dissimilarity vs spatial distance
#for DNA/RNA and all timepoints so it can be plotted
Full_Spatial_df_16S=
  rbind(spatial_16S_dfs$`16S_DNA_June`,
        spatial_16S_dfs$`16S_DNA_August`,
        spatial_16S_dfs$`16S_DNA_October`,
        spatial_16S_dfs$`16S_cDNA_June`,
        spatial_16S_dfs$`16S_cDNA_August`,
        spatial_16S_dfs$`16S_cDNA_October`) %>%
  mutate(source=factor(source, levels=c("DNA", "cDNA")),
         Time=factor(Time, levels=c("June", "August", "October")))

#Fungal ITS2 - creates plot of community dissimilarity vs spatial distance
#facetted by dataset (DNA vs RNA) and season
SpatAutPlot_ITS=
  ggplot(Full_Spatial_df_ITS, 
       aes(x=SpatialDist, y=CommunityDist)) +
  geom_point(size=4, alpha=0.3, aes(color=DistClass)) +
  theme_test(base_size = 16) +
  scale_color_manual(values = c("#f53600", "#4ba7bc")) +
  geom_smooth(method="lm", 
              se=FALSE, 
              color="black", 
              linewidth=2) +
  scale_x_log10(breaks = c(15, 100, 1000, 10000)) +
  facet_grid(source~Time, scales="free_y")

#Bacterial 16S - creates plot of community dissimilarity vs spatial distance
#facetted by dataset (DNA vs RNA) and season
SpatAutPlot_16S=
  ggplot(Full_Spatial_df_16S, 
       aes(x=SpatialDist, y=CommunityDist)) +
  geom_point(size=4, alpha=0.2, aes(color=DistClass)) +
  theme_test(base_size = 16) +
  scale_color_manual(values = c("#f53600", "#4ba7bc")) +
  geom_smooth(method="lm", 
              se=FALSE, 
              color="black", 
              linewidth=2) +
  scale_x_log10(breaks = c(15, 100, 1000, 10000)) +
  facet_grid(source~Time, scales="free_y")


#saves plot of bacterial 16S spatial autocorellation as pdf to be used in figure
pdf(file="SpatAutPlot_16S.pdf", 
    width = 10, height = 6)
SpatAutPlot_16S
dev.off()
