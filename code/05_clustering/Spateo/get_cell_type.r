#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
#module load conda_R/4.5

library(SpatialExperiment)
library(here)

celltype_cols <- readRDS("/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/processed-data/05_Clustering/CellType_cols.Rds")

# old annotation
# spe_annot <- readRDS("/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/processed-data/06_label_transfer/Objects/Banksy_CellTypes_transfer.Rds")
# write.csv(spe_annot, file = "/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/processed-data/05_Clustering/Banksy_cell_types.csv", row.names = FALSE)

spe_annot <- readRDS("/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/processed-data/02_build_spe/SPEs/spe_celltype_v2.Rds")
cell_types <- spe_annot$CellTypes
cell_ids <- colnames(spe_annot)
out_df <- data.frame(cell_id   = cell_ids, CellType  = cell_types, stringsAsFactors = FALSE)

write.csv(out_df, file = "/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/processed-data/05_Clustering/Banksy_cell_types.csv", row.names = FALSE)
