# Make heatmaps for 02 and 03 scripts
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(MetaNeighbor)
library(gplots)
library(here)

########################################
########### 02 rat human ##############
#######################################

#Load output from human vs rat
output_02 <- readRDS(here("processed-data","MetaNeighbor","Human_Both_rat_snRNA_Modalities_aurocs.Rds"))

#Generate some side colors
labels <- colnames(output_02)
dataset <- sub("\\|.*", "", labels)
species <- sub("\\_.*", "", labels)
modality <- sub(".*_", "", dataset)

# Create color mappings
species_colors <- setNames(
  c("tomato","dodgerblue"),
  unique(species)
)

modality_colors <- setNames(
  c("grey55","black"),
  unique(modality)
)


stopifnot(identical(rownames(output_02),colnames(output_02)))

pdf(file = here("plots","MetaNeighbor","Human_Both_Rat_snRNA_Heatmap.pdf"),width = 12,height = 12)
plotHeatmap(output_02,
            cex = 0.75,
            ColSideColors = species_colors[species],
            RowSideColor = modality_colors[modality],
            key.title = NA,
            key.xlab = "AUROC",
            keysize = 0.6,
            margins = c(12, 12))
legend("topright",
       legend = names(species_colors),
       fill   = species_colors,
       title  = "Species",
       border = FALSE,
       bty    = "n",
       cex    = 0.7
)
legend("bottomleft",
       legend = names(modality_colors),
       fill   = modality_colors,
       title  = "Modality",
       border = FALSE,
       bty    = "n",
       cex    = 0.7
)
dev.off()

########################################
############## 03 MSNs ################
#######################################


#Load output from human vs rat
output_03 <- readRDS(here("processed-data","MetaNeighbor","All_Species_MSNs_Only_aurocs.Rds"))

#Generate some side colors
labels <- colnames(output_03)
dataset <- sub("\\|.*", "", labels)
species <- sub("\\_.*", "", labels)
modality <- sub(".*_", "", dataset)

# Create color mappings
species_colors <- setNames(
  c("tomato","dodgerblue","green4"),
  unique(species)
)

modality_colors <- setNames(
  c("grey55","black"),
  unique(modality)
)


stopifnot(identical(rownames(output_03),colnames(output_03)))

pdf(file = here("plots","MetaNeighbor","All_Species_MSNs_Heatmap.pdf"),width = 12,height = 12)
plotHeatmap(output_03,
            cex = 0.75,
            ColSideColors = species_colors[species],
            RowSideColor = modality_colors[modality],
            key.title = NA,
            key.xlab = "AUROC",
            keysize = 0.6,
            margins = c(12, 12))
legend("topright",
       legend = names(species_colors),
       fill   = species_colors,
       title  = "Species",
       border = FALSE,
       bty    = "n",
       cex    = 0.7
)
legend("bottomleft",
       legend = names(modality_colors),
       fill   = modality_colors,
       title  = "Modality",
       border = FALSE,
       bty    = "n",
       cex    = 0.7
)
dev.off()

sessionInfo()
