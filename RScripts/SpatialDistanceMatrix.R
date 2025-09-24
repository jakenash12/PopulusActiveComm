#This script uses GPS coordinates taken at each site to calculate a spatial distance matrix
#between each plots. Only 1 GPS coordinate was taken at each site, so to calculate the
#distances between plots at the same site, their measured distances (with a tape in field)
#along the transect were used

library(readxl)
library(dplyr)
library(tidyr)
library(geosphere)
library(tibble)

#Reads in field data with GPS coordinates
PlotSpatialDistMatrix <- read.delim("./InputFiles/AUE_PlotEstablishment.txt")

# Select only relevant columns
PlotLocations <- PlotSpatialDistMatrix %>%
  select(Site, PlotNum, SiteLat, SiteLon)

# Prepare all pairwise combinations
PlotPairs <- PlotLocations %>%
  rename(Site1 = Site, PlotNum1 = PlotNum, Lat1 = SiteLat, Lon1 = SiteLon) %>%
  crossing(
    PlotLocations %>%
      rename(Site2 = Site, PlotNum2 = PlotNum, Lat2 = SiteLat, Lon2 = SiteLon)
  ) %>%
  filter(!(Site1 == Site2 & PlotNum1 == PlotNum2)) %>%
  mutate(
    PlotID1 = paste0(Site1, PlotNum1),
    PlotID2 = paste0(Site2, PlotNum2)
  )

# Predefine transect distances
transect_dists <- c(30, 15, 15, 30)
plot_order <- c(1, 2, 3, 4, 5)

# Compute distances
PlotPairs <- PlotPairs %>%
  rowwise() %>%
  mutate(
    Distance_m = if (Site1 != Site2) {
      distHaversine(c(Lon1, Lat1), c(Lon2, Lat2))
    } else {
      idx1 <- match(PlotNum1, plot_order)
      idx2 <- match(PlotNum2, plot_order)
      if (!is.na(idx1) && !is.na(idx2)) {
        as.numeric(sum(transect_dists[seq(min(idx1, idx2), max(idx1, idx2) - 1)]))
      } else {
        NA_real_
      }
    }
  ) %>%
  ungroup()


###################
# Start from the correct long format
DistanceMatrixWide <- PlotPairs %>%
  mutate(
    PlotID1 = paste0(Site1, PlotNum1),
    PlotID2 = paste0(Site2, PlotNum2)
  ) %>%
  select(PlotID1, PlotID2, Distance_m)

# Add symmetric reverse entries
DistanceMatrixSymmetric <- DistanceMatrixWide %>%
  bind_rows(
    DistanceMatrixWide %>%
      rename(PlotID1 = PlotID2, PlotID2 = PlotID1)
  ) %>%
  distinct()

# Add diagonal (distance to self = 0)
all_ids <- sort(unique(c(DistanceMatrixWide$PlotID1, DistanceMatrixWide$PlotID2)))
diag_df <- tibble(
  PlotID1 = all_ids,
  PlotID2 = all_ids,
  Distance_m = 0
)

# Combine and pivot to wide format
DistanceMatrixFull <- DistanceMatrixSymmetric %>%
  bind_rows(diag_df) %>%
  mutate(
    PlotID1 = factor(PlotID1, levels = all_ids),
    PlotID2 = factor(PlotID2, levels = all_ids)
  ) %>%
  pivot_wider(
    names_from = PlotID2,
    values_from = Distance_m
  ) %>%
  arrange(PlotID1)

# Convert to a numeric matrix with rownames
DistanceMatrixMatrix <- DistanceMatrixFull %>%
  column_to_rownames("PlotID1") %>%
  as.matrix()

#exports spatial distance matrix as tsv
DistanceMatrixMatrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "Plot") %>%
  write_tsv("./OutputFiles/SpatialDistanceMatrix.txt")
