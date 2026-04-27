#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
#module load conda_R/4.5

library(SpatialExperiment)
library(here)

CellType_cols <- readRDS(here("processed-data","05_Clustering","CellType_cols.Rds"))

#Rename colors 
names(CellType_cols)[18] <- "Astrocyte_Oligo"
names(CellType_cols)[22] <- "Microglia_Oligo"

#Remove WM_A-C colors
CellType_cols <- CellType_cols[-c(3:5)]

#Add a color for WM
CellType_cols <- c(CellType_cols,"orange")
names(CellType_cols)[20] <- "WM"

#     Excitatory     Microglia_A     D1_Island_B         Astro_A    Fibroblast_A 
#       "#DB1C5F"       "#0D87E4"       "#00FF0D"       "#FD00FD"       "#FFA7E2" 
#    Fibroblast_B        DRD2_MSN        DRD1_MSN         Astro_B       MSN_Oligo 
#       "#2AFECA"       "#7AAA16"       "#9400FF"       "#823526"       "#F5DEC0" 
#       Inh_PVALB             OPC     Microglia_B         Inh_SST Astrocyte_Oligo 
#       "#93E5FF"       "#7A1699"       "#FF0DBC"       "#C4B3FB"       "#F80D2A" 
#     D1_Island_A            CHAT       Ependymal Microglia_Oligo              WM 
#       "#224B82"       "#FBA475"       "#73690D"       "#B62A7A"        "orange" 

spe_annot <- readRDS("/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/processed-data/02_build_spe/SPEs/spe_celltype_v2.Rds")
cell_types <- spe_annot$CellTypes
cell_ids <- colnames(spe_annot)
out_df <- data.frame(cell_id   = cell_ids, CellType  = cell_types, stringsAsFactors = FALSE)

write.csv(out_df, file = "/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/processed-data/05_Clustering/Banksy_cell_types.csv", row.names = FALSE)
