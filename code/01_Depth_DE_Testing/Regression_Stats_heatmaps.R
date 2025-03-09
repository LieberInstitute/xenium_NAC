#module load r_nac
#Goal: Compile stats from MiddlevsOther analysis and make heatmaps for R-squared and -log10(p)

library(ggplot2)
library(here)

list.files(path = "processed-data/Regression_stats/MiddlevsOther")
# [1] "Astrocyte_AMiddlevsOther_Stats.Rds" "Astrocyte_BMiddlevsOther_Stats.Rds"
# [3] "DRD1_MSN_AMiddlevsOther_Stats.Rds"  "DRD1_MSN_BMiddlevsOther_Stats.Rds" 
# [5] "DRD1_MSN_CMiddlevsOther_Stats.Rds"  "DRD1_MSN_DMiddlevsOther_Stats.Rds" 
# [7] "DRD2_MSN_AMiddlevsOther_Stats.Rds"  "DRD2_MSN_BMiddlevsOther_Stats.Rds" 
# [9] "EndothelialMiddlevsOther_Stats.Rds" "ExcitatoryMiddlevsOther_Stats.Rds" 
# [11] "Inh_AMiddlevsOther_Stats.Rds"       "Inh_BMiddlevsOther_Stats.Rds"      
# [13] "Inh_CMiddlevsOther_Stats.Rds"       "Inh_DMiddlevsOther_Stats.Rds"      
# [15] "Inh_EMiddlevsOther_Stats.Rds"       "Inh_FMiddlevsOther_Stats.Rds"      
# [17] "MicrogliaMiddlevsOther_Stats.Rds"   "OligoMiddlevsOther_Stats.Rds"      
# [19] "OPCMiddlevsOther_Stats.Rds" 

# Get a list of all .rds files in the directory
rds_files <- list.files("processed-data/Regression_stats/MiddlevsOther",full.names = TRUE)

# Read each .rds file into a list
rds_list <- lapply(rds_files, readRDS)

# change names
names(rds_list) <- c("Astrocyte_A","Astrocyte_B",
                     "DRD1_MSN_A","DRD1_MSN_B","DRD1_MSN_C","DRD1_MSN_D",
                     "DRD2_MSN_A","DRD2_MSN_B",
                     "Endothelial","Excitatory",
                     "Inh_A","Inh_B","Inh_C", "Inh_D", "Inh_E", "Inh_F",
                     "Microglia","Oligo","OPC")

#Compile to dataframe
rds_df <- do.call(what = rbind,rds_list)

#Make PCA a factor
rds_df$PCA <- factor(rds_df$PCA,levels=rev(c("PC1","PC2","PC3","PC4","PC5",
                                         "PC6","PC7","PC8","PC9","PC10")))


#R^2 heatmap
rds_df$R_2_rounded <- round(rds_df$`R-squared`,3)
R_squared <- ggplot(rds_df,aes(x = CellType,y = PCA,fill = R_2_rounded)) +
  geom_tile(color = "black") +
  geom_text(aes(label = R_2_rounded), color = "black") +
  scale_fill_gradientn(colours = c("lightgrey","red")) + 
  labs(fill = "R-squared",
       x = "Cell Type",
       y = "Principal Component") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave(plot = R_squared,
       filename = "/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/plots/01_Depth_DE_Testing/Regression/Middle_vs_Other/R_squared_heatmap.pdf",
       height = 12,width = 12)



#p-value HEATMAP
rds_df$neg_log10p <- round(-log10(rds_df$p),3)

neg_log10p <- ggplot(rds_df,aes(x = CellType,y = PCA,fill = neg_log10p)) +
  geom_tile(color = "black") +
  geom_text(aes(label = neg_log10p), color = "black") +
  scale_fill_gradientn(colours = c("lightgrey","red")) + 
  labs(fill = "-log10(p-value)",
       x = "Cell Type",
       y = "Principal Component") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave(plot = neg_log10p,
       filename = "/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/plots/01_Depth_DE_Testing/Regression/Middle_vs_Other/neglog10p_heatmap.pdf",
       height = 12,width = 12)

