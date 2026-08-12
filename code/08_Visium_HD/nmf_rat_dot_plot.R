library(SpatialFeatureExperiment)
library(SingleCellExperiment)
library(HDF5Array)
library(tidyverse)
library(ggplot2)
library(here)


here::i_am("code/08_Visium_HD/nmf_rat_dot_plot.R")

###### Visium-HD (rat NMF)
# Load object
sfe_dir <- here("processed-data", "HD_Full_Analysis", "sfe_with_rat_NMF")
sfe <- loadHDF5SummarizedExperiment(sfe_dir)

# --- Config: adjust these ---
nmf_factors <- c("nmf10_rat", "nmf18_rat")
group_var   <- "Spatial_Domain"
threshold   <- 0

# --- Pull relevant colData into a data.frame ---
cd <- as.data.frame(colData(sfe))[, c(group_var, nmf_factors)]

# --- Center and scale each NMF factor (z-score across all spots/cells) ---
cd_scaled <- cd
cd_scaled[nmf_factors] <- scale(cd[nmf_factors], center = TRUE, scale = TRUE) # z score

# --- Reshape to long format ---
cd_long <- cd_scaled %>%
  pivot_longer(cols = all_of(nmf_factors), names_to = "Factor", values_to = "ScaledValue")

cd_raw_long <- cd %>%
  pivot_longer(cols = all_of(nmf_factors), names_to = "Factor", values_to = "RawValue")

cd_long$RawValue <- cd_raw_long$RawValue

# --- Summarize per group x factor ---
summary_df <- cd_long %>%
  group_by(.data[[group_var]], Factor) %>%
  summarise(
    MeanScaled = mean(ScaledValue, na.rm = TRUE),
    PctAbove   = mean(RawValue > threshold, na.rm = TRUE) * 100,
    .groups = "drop"
  )

summary_df$Factor <- factor(summary_df$Factor, levels = nmf_factors)

# --- Dot plot ---
dotplot <- ggplot(summary_df, aes(x = Factor, y = .data[[group_var]])) +
  geom_point(aes(size = PctAbove, color = MeanScaled)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_size(range = c(0, 8), name = "% spots > threshold") +
  labs(x = "NMF Factor", y = group_var, color = "Mean\n(scaled)") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey90")
  )

ggsave(plot = dotplot,
       filename = here("plots", "HD_Full_Analysis", "rat_nmf_spatial_domain_dotplot.pdf"),
       height = 8, width = 6)

#-------------------------- Reproducibility information -----------------------#
message("\nReproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
