#Goal: Transfer BANKSY cluster labels from subsampled bins to the full dataset
#   using spatial kNN. At lambda=0.8, domains are primarily spatial, so spatial
#   nearest neighbors are appropriate for label transfer.
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5
library(SpatialExperiment)
library(BiocNeighbors)
library(sessioninfo)
library(escheR)
library(here)

################################################################################
#   Setup
################################################################################

message(Sys.time(), " | Loading data")

#   Load the full filtered SPE
spe <- readRDS(here(
  "processed-data", "HD_Full_Analysis", "SPEs", "spe_norm_binQC_filtered.Rds"
))

spe <- spe[, colSums(counts(spe)) > 0]

#   Load the subsampled BANKSY object (has cluster labels)
spe_sub <- readRDS(here(
  "processed-data", "HD_Full_Analysis", "SPEs",
  "spe_bin_banksy_subsampled_clustered.RDS"
))

#   Load subsample indices
subsample_idx <- readRDS(here("processed-data","HD_Full_Analysis","Banksy","subsample_idx_8umbins.Rds"))

#   Identify cluster columns from BANKSY
clust_cols <- grep("^clust_", colnames(colData(spe_sub)), value = TRUE)
message(sprintf(
  "Full SPE: %d bins | Subsampled: %d bins | Cluster columns: %s",
  ncol(spe), ncol(spe_sub), paste(clust_cols, collapse = ", ")
))

################################################################################
#   Label transfer via spatial kNN
################################################################################

message(Sys.time(), " | Transferring cluster labels via spatial kNN")

k_transfer <- 10
samples <- unique(spe$sample_id)

#   Indices of bins NOT in the subsample
remaining_idx <- setdiff(seq_len(ncol(spe)), subsample_idx)
message(sprintf("  Bins to transfer: %d", length(remaining_idx)))

#   Initialize cluster columns in the full SPE
for (clust_col in clust_cols) {
  spe[[clust_col]] <- NA_character_
  #   Fill in subsampled bins directly
  spe[[clust_col]][subsample_idx] <- as.character(spe_sub[[clust_col]])
}

#   Transfer per sample (spatial coords only comparable within a sample)
for (sid in samples) {
  message(sprintf("  Processing sample: %s", sid))
  
  sub_in_sample <- intersect(subsample_idx, which(spe$sample_id == sid))
  rem_in_sample <- intersect(remaining_idx, which(spe$sample_id == sid))
  
  if (length(rem_in_sample) == 0) {
    message("    No remaining bins for this sample, skipping")
    next
  }
  
  #   Spatial coordinates
  coords_ref <- spatialCoords(spe)[sub_in_sample, , drop = FALSE] #Maintain matrix 
  coords_query <- spatialCoords(spe)[rem_in_sample, , drop = FALSE]
  
  #   Find k nearest subsampled neighbors for each remaining bin
  nn <- queryKNN(X = coords_ref, query = coords_query, k = k_transfer)
  
  #   Majority vote for each cluster resolution
  for (clust_col in clust_cols) {
    ref_labels <- spe[[clust_col]][sub_in_sample] #Pull labels
    #nn$index is a dataframe where rows are unlabeled bins and columns are the nearest cells. Values ar ethe indices
    transferred <- apply(nn$index, 1, function(idx) {
      labels <- ref_labels[idx]
      names(sort(table(labels), decreasing = TRUE))[1]
    })
    
    spe[[clust_col]][rem_in_sample] <- transferred
  }
  
  message(sprintf("    Transferred labels for %d bins", length(rem_in_sample)))
}

#   Convert cluster columns to factors
for (clust_col in clust_cols) {
  spe[[clust_col]] <- factor(spe[[clust_col]])
  message(sprintf("  %s: %d clusters", clust_col, nlevels(spe[[clust_col]])))
}

#   Sanity check: no NAs remaining
for (clust_col in clust_cols) {
  n_na <- sum(is.na(spe[[clust_col]]))
  if (n_na > 0) {
    warning(sprintf("  %s has %d NA labels!", clust_col, n_na))
  }
}

################################################################################
#   Set up colors for plotting
################################################################################

message(Sys.time(), " | Generating Colors")

cluster_colors <- vector(mode = "list", length = 3)
names(cluster_colors) <- clust_cols
### Generate cluster colors for each clust_fol
for(clust_col in clust_cols){
  #Make colors for CellType
  cluster_colors[[clust_col]] <- Polychrome::createPalette(length(unique(spe[[clust_col]])),
                                             c("#D81B60", "#1E88E5","#ffcd14","#ff7814","#004D40"))
  names(cluster_colors[[clust_col]]) <- unique(spe[[clust_col]])
  saveRDS(object = cluster_colors[[clust_col]],
          file = here("processed-data","HD_Full_Analysis","Cluster_colors",paste(clust_col,"_bin_level_colors.Rds")))
}



################################################################################
#   Plot clusters on the samples
################################################################################
message(Sys.time(), " | Plotting clusters on samples")
for(sample in samples){
  print(sample)
  spe_sub2 <- spe[,spe$sample_id == sample]
  for(clust_col in clust_cols){
    print(clust_col)
    p <- make_escheR(spe_sub2) |>
      add_fill(clust_col) +
      scale_fill_manual(values = cluster_colors[[clust_col]])
    ggsave(
      here("plots","HD_Full_Analysis","Banksy","Bin_Level",paste0(sample,"_",clust_col,".png")), 
      p,
      width = 10, height = 8, dpi = 200
      )
  }
}

################################################################################
#   Save
################################################################################

message(Sys.time(), " | Saving full SPE with cluster labels")

saveRDS(spe, here(
  "processed-data", "HD_Full_Analysis", "SPEs",
  "spe_bin_level_banksy_lambda_pt8.RDS"
))

message(sprintf("Done. Full SPE has %d bins with cluster labels.", ncol(spe)))

###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
