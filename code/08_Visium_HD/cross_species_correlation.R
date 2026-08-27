library(SingleCellExperiment)
library(DeconvoBuddies)
library(orthogene)
library(HDF5Array)
library(ggrepel)
library(ggplot2)
library(Seurat)
library(dplyr)
library(here)


setwd("/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/")

###### Visium-HD
# Load object 
sfe_dir <- "/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/processed-data/HD_Full_Analysis/sfe_spatial_annotated"

sfe <- loadHDF5SummarizedExperiment(sfe_dir)


anno_df <- as.data.frame(unique(colData(sfe)[,c("spatial_0.4","Spatial_Domain")]))

vhd_stats <- readRDS("/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/processed-data/HD_Full_Analysis/spatial_0.4_1vALL.Rds")
vhd_stats$Spatial_Domain <- anno_df[match(vhd_stats$cellType.target,anno_df$spatial_0.4),"Spatial_Domain"]

write.csv(vhd_stats,file = "/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/processed-data/HD_Full_Analysis/VHD_Spatial_Domain_1vALL.csv")


###########################################
###########     NHP   #####################
###########################################

message("Loading NHP object |", Sys.time())

#Plan: Identify 1-to-1 orthologs from NHP to human --> Subset the  NHP object --> Change rownames

nhp_obj <- readRDS("/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/processed-data/NHP_Data/Results_MSNs_processed_final.rds")
DefaultAssay(nhp_obj) <- "RNA"
nhp_sce <- as.SingleCellExperiment(x = nhp_obj)
nhp_sce

nhp_stats <- findMarkers_1vAll(sce = nhp_sce,
                               assay_name = "logcounts",
                               cellType_col = "MSN_type",
                               direction = "up",mod = "~monkey")

write.csv(nhp_stats,file = "/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/processed-data/HD_Full_Analysis/NHP_MSN_Type_1vALL.csv")

#Convert to 121 orthologs
nhp_orthologs <- convert_orthologs(gene_df = counts(nhp_sce),
                                   gene_input = "rownames",
                                   gene_output = "dict",
                                   input_species = "macaque",
                                   output_species = "human",
                                   non121_strategy = "drop_both_species",
                                   method = "gprofiler")

nhp_orthologs <- data.frame(macaque = names(nhp_orthologs),
                            human = unname(nhp_orthologs))

nhp_stats$human_gene <- nhp_orthologs[match(nhp_stats$gene,nhp_orthologs$macaque),"human"]
nhp_stats <- nhp_stats[!is.na(nhp_stats$human_gene),]

vhd_stats <- vhd_stats[vhd_stats$gene %in% unique(nhp_stats$human_gene),]


#Subset to celltypes
isl_vhd <- subset(vhd_stats,subset=(Spatial_Domain %in% c("D1_Island_A","D1_Island_B")))
isl_vhd <- isl_vhd[,c("Spatial_Domain","std.logFC","gene")]
colnames(isl_vhd)[1] <- c("cellType.target")
isl_nhp <- subset(nhp_stats,subset=(cellType.target %in% c("D1-ICj","D1-NUDAP")))
isl_nhp <- isl_nhp[,c("cellType.target","std.logFC","human_gene")]
colnames(isl_nhp)[3] <- "gene"

# ---- Reshape long (gene, cellType.target, std.logFC) -> wide -----------------
isl_vhd_wide <- tidyr::pivot_wider(
  isl_vhd[, c("gene", "cellType.target", "std.logFC")],
  id_cols     = gene,
  names_from  = cellType.target,
  values_from = std.logFC
)

isl_nhp_wide <- tidyr::pivot_wider(
  isl_nhp[, c("gene", "cellType.target", "std.logFC")],
  id_cols     = gene,
  names_from  = cellType.target,
  values_from = std.logFC
)

isl_wide <- merge(x = isl_vhd_wide,y = isl_nhp_wide,by = "gene")

isl_wide$quadrant <- with(isl_wide, ifelse(
  D1_Island_A > 0 & `D1-NUDAP` > 0, "Enriched in both",
  ifelse(D1_Island_A < 0 & `D1-NUDAP` < 0, "Depleted in both",
         ifelse(D1_Island_A > 0 & `D1-NUDAP` <= 0, "D1_Island_A only",
                "D1-NUDAP only"))))

d1_isl_cor <- ggplot(data = isl_wide,aes(x = D1_Island_A,y = `D1-NUDAP`,color = quadrant)) +
  geom_point() +
  scale_color_manual(values = c("Enriched in both"  = "tomato",
                                "Depleted in both"  = "dodgerblue",
                                "D1_Island_A only"  = "grey",
                                "D1-NUDAP only"     = "grey")) +
  geom_label_repel(data = subset(isl_wide,
                                 subset=(gene %in% c("OPRM1","TSHZ1","GABRQ",
                                                     "EYA2","SLC35F1","FOXP2",
                                                     "RXFP1","CHST9","VWC2L"))),
                 aes(label = gene),color = "black") +
  xlim(c(-1,1.5)) +
  ylim(c(-2.5,4)) +
  geom_hline(yintercept = 0,lty = 2) +
  geom_vline(xintercept = 0,lty = 2) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("D1_Island_A std.logFC\n(Human)") +
  ylab("D1-NUDAP std.logFC \n(NHP)")
ggsave(d1_isl_cor,filename = "/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/plots/HD_Full_Analysis/D1_Island_A_vs_D1-NUDAP.pdf",
       height = 8, width = 8)


isl_wide$quadrant <- with(isl_wide, ifelse(
  D1_Island_B > 0 & `D1-ICj` > 0, "Enriched in both",
  ifelse(D1_Island_B < 0 & `D1-ICj` < 0, "Depleted in both",
         ifelse(D1_Island_A > 0 & `D1-ICj` <= 0, "D1_Island_B only",
                "D1_ICJ only"))))


Icj <- ggplot(data = isl_wide,aes(x = D1_Island_B,y = `D1-ICj`,color = quadrant)) +
  geom_point() +
  scale_color_manual(values = c("Enriched in both"  = "tomato",
                                "Depleted in both"  = "dodgerblue",
                                "D1_Island_B only"  = "grey",
                                "D1_ICJ only"     = "grey")) +
  geom_label_repel(data = subset(isl_wide,
                                 subset=(gene %in% c("DRD3","NTN1","CPNE4","KCNT2",
                                                     "ISL1","PROK2","MYO16","CXCL14",
                                                     "TCERG1L","MC4R","VIP"))),
                   aes(label = gene),color = "black") +
  ylim(c(-3,4.5)) +
  geom_hline(yintercept = 0,lty = 2) +
  geom_vline(xintercept = 0,lty = 2) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("D1_Island_B std.logFC\n(Human)") +
  ylab("D1-ICj std.logFC \n(NHP)")
ggsave(Icj,filename = "/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/plots/HD_Full_Analysis/D1_Island_B_vs_D1-ICj.pdf",
       height = 8, width = 8)

#-------------------------- Reproducibility information -----------------------#
message("\nReproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
