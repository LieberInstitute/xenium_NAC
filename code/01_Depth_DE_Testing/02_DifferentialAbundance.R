#module load r_nac
#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
#Goal: Perform differential abundance analysis.
#Following workflow at https://bioconductor.org/books/3.13/OSCA.multisample/differential-abundance.html

library(SingleCellExperiment)
library(sessioninfo)
library(ggrepel)
library(ggplot2)
library(edgeR)
library(here)

#load the sce object
sce <- readRDS("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/sce_CellType_noresiduals.Rds")


#Remove the Neuron_Ambig group
sce <- sce[,sce$CellType.Final != "Neuron_Ambig"]

dim(sce)
#[1]  36601 103339

#Remove the genes with 0 counts
sce <- sce[!rowSums(assay(sce, "counts")) == 0, ]

dim(sce)
#[1]  34977 103339

#For this each sample needs an anterior/middle/posterior designation. 
#Combine the Anterior and Posterior samples into a single group called "Anterior_Posterior"
#Make a dataframe. 
Ant_Mid_Post <- data.frame(Brain_ID = unique(sce$Brain_ID))

Ant_Mid_Post <- cbind(Ant_Mid_Post,c("Anterior_Posterior","Middle",
                                     "Middle","Anterior_Posterior",
                                     "Anterior_Posterior","Anterior_Posterior",
                                     "Anterior_Posterior","Middle",
                                     "Middle","Anterior_Posterior"))

colnames(Ant_Mid_Post)[2] <- "Depth"

Ant_Mid_Post
# Brain_ID              Depth
# 1    Br8325 Anterior_Posterior
# 2    Br8492             Middle
# 3    Br2720             Middle
# 4    Br6423 Anterior_Posterior
# 5    Br2743 Anterior_Posterior
# 6    Br3942 Anterior_Posterior
# 7    Br6432 Anterior_Posterior
# 8    Br6471             Middle
# 9    Br6522             Middle
# 10   Br8667 Anterior_Posterior


#Add depth to the sce object
sce$Depth <- Ant_Mid_Post[match(sce$Brain_ID,Ant_Mid_Post$Brain_ID),"Depth"]

as.data.frame(unique(colData(sce)[,c("Brain_ID","Depth")]))
# Brain_ID              Depth
# 1_AAACCCAAGACCAACG-1    Br8325 Anterior_Posterior
# 3_AAACCCAAGGTGAGCT-1    Br8492             Middle
# 5_AAACCCAGTAATTAGG-1    Br2720             Middle
# 7_AAACCCACACCCTTAC-1    Br6423 Anterior_Posterior
# 9_AAACCCACATTGTAGC-1    Br2743 Anterior_Posterior
# 11_AAACCCAAGACTCCGC-1   Br3942 Anterior_Posterior
# 13_AAACCCAAGGACTGGT-1   Br6432 Anterior_Posterior
# 15_AAACCCAAGGGAGATA-1   Br6471             Middle
# 17_AAACCCAAGATTGCGG-1   Br6522             Middle
# 19_AAACCCACAAGGCCTC-1   Br8667 Anterior_Posterior

#Pull table of cell type abundance by cluster and sample
abundances <- table(sce$CellType.Final, sce$Sample) 
abundances <- unclass(abundances) 
colSums(abundances)
# 1c_NAc_SVB  2c_NAc_SVB  3c_NAc_SVB  4c_NAc_SVB  5c_NAc_SVB  6c_NAc_SVB 
# 5188        4797        3749        5198        5241        2375 
# 7c_NAc_SVB  8c_NAc_SVB  9c_NAc_SVB 10c_NAc_SVB 11c_NAc_SVB 12c_NAc_SVB 
# 5513        5431        3928        6183        4746        5098 
# 13c_NAc_SVB 14c_NAc_SVB 15c_Nac_SVB 16c_Nac_SVB 17c_Nac_SVB 18c_Nac_SVB 
# 6040        5421        5425        6766        5061        6098 
# 19c_Nac_SVB 20c_Nac_SVB 
# 5389        5692 

# Attaching some column metadata.
extra.info <- colData(sce)[match(colnames(abundances), sce$Sample),]
y.ab <- DGEList(abundances, samples=extra.info)

#Per OSCA 6.2 "For a DA analysis of cluster abundances, filtering is generally not required as most clusters will not be of low-abundance"
#Will not filter
#Want to model sample and sort as technical variables and Depth as the biological variable
design <- model.matrix(~factor(Sort) + factor(Depth), y.ab$samples)

#Use estimateDisp to "estimate the NB dispersion for each cluster"
#In the workflow, they turn off the trend because too few samples. n=6 and n=4 for two depth levels so also going to turn off trend. 
y.ab <- estimateDisp(y.ab, design, trend="none")
summary(y.ab$common.dispersion)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2321  0.2321  0.2321  0.2321  0.2321  0.2321 

#The edgeR package provides a function to plot BCV, but want to add some additional 
#information to the plots
# Extract data
bcv_data <- data.frame(
  logCPM = y.ab$AveLogCPM,
  BCV = sqrt(y.ab$tagwise.dispersion), #BCV is squareroot of negative binomial dispersion
  CellType = rownames(y.ab)  
)


# Plot with ggplot2
BCV_plot <- ggplot(bcv_data, aes(x=logCPM, y=BCV, label=CellType)) +
  geom_point(aes(color="Tagwise"),size=2) +  # Scatter plot points
  geom_hline(aes(yintercept=sqrt(y.ab$common.dispersion), color="Common"), linetype="solid") +  # Common BCV line
  geom_text_repel(size=3, max.overlaps=15) +  
  scale_color_manual(values = c("Tagwise"="black","Common"="red"),
                     name = "Dispersion Type") +
  guides(color = guide_legend(override.aes = list(linetype = c("solid","blank"),
                                                  shape = c(NA,16)
                                                  )
                              )
         ) +
  labs(x="Average log CPM", y="Biological Coefficient of Variation") +
  theme_bw()

ggsave(plot = BCV_plot,
       filename = here("plots","01_Depth_DE_Testing","Differential_Abundance",paste0("BCV_plot.pdf")),
       height = 8,width = 8)


fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
summary(fit.ab$var.prior)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3858  0.3858  0.3858  0.3858  0.3858  0.3858 

summary(fit.ab$df.prior)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.067   2.067   2.067   2.175   2.377   2.377 

pdf(file = here("plots","01_Depth_DE_Testing","Differential_Abundance","QLDisp_plot.pdf"))
plotQLDisp(fit.ab, cex=1)
dev.off()

res <- glmQLFTest(fit.ab, coef=ncol(design))
summary(decideTests(res))
# factor(Depth)Middle
# Down                     2
# NotSig                  12
# Up                       6

topTags(res,n=20)
# Coefficient:  factor(Depth)Middle 
# logFC   logCPM           F       PValue          FDR
# DRD1_MSN_B   0.84021542 15.95413 44.97311822 1.869977e-06 3.739953e-05
# DRD1_MSN_D   2.53571837 13.71718 34.95727623 1.066378e-05 1.066378e-04
# DRD1_MSN_C   1.66407712 14.71002 23.53123048 1.099748e-04 7.331655e-04
# Inh_E        0.73088502 12.48620 13.38978926 1.630600e-03 8.153002e-03
# Inh_F        1.00399611 13.68900  9.21237698 6.789890e-03 2.715956e-02
# Astrocyte_B -1.34094027 13.61039  7.41321251 1.348143e-02 4.493811e-02
# Inh_D        0.63605147 11.34882  6.66550951 1.822658e-02 4.651643e-02
# Oligo       -0.75678123 17.41619  6.61858578 1.860657e-02 4.651643e-02
# Microglia    0.47461596 15.44249  2.64033080 1.205975e-01 2.679945e-01
# Inh_B        0.41584873 12.74417  2.37630416 1.396220e-01 2.792440e-01
# OPC         -0.49271491 15.16974  1.79954699 1.955277e-01 3.555049e-01
# Astrocyte_A -0.34217274 16.37992  0.96179760 3.390153e-01 5.459321e-01
# Ependymal   -0.81805726 13.27483  0.89921964 3.548559e-01 5.459321e-01
# DRD2_MSN_B   0.10657443 14.12308  0.65380153 4.285732e-01 5.725272e-01
# Inh_A        0.08787834 14.24425  0.58124451 4.550060e-01 5.725272e-01
# DRD2_MSN_A   0.08905092 17.72113  0.56682091 4.605746e-01 5.725272e-01
# Inh_C       -0.11587585 13.78916  0.50296148 4.866481e-01 5.725272e-01
# Excitatory   0.20620528 13.88534  0.27063562 6.088938e-01 6.745831e-01
# Endothelial -0.12492634 13.14258  0.22471472 6.408540e-01 6.745831e-01
# DRD1_MSN_A  -0.03417355 17.74105  0.06411442 8.027743e-01 8.027743e-01

#Pull results dataframe
tag_df <- as.data.frame(topTags(res,n=20))

#Add a celltype column
tag_df$CellType <- rownames(tag_df)

#Write out the results as a csv
write.csv(x = tag_df,file = here("processed-data","Differential_Abundance","DA_results.csv"),quote = FALSE)

#DA volcano plot
DA_volcano <- ggplot(data = tag_df,aes(x = logFC,y = -log10(FDR),size = F,label = CellType)) +
  geom_point() +
  geom_point(data = subset(tag_df,subset=(FDR<=0.05 & logFC >0)),col = "red") +
  geom_point(data = subset(tag_df,subset=(FDR<=0.05 & logFC <0)),col = "dodgerblue") +
  guides(size = guide_legend(override.aes = list(color = "black"))) +
  geom_text_repel(size=3, max.overlaps=15) +
  labs(size = "F-statistic") +
  xlim(c(-3,3)) +
  geom_hline(yintercept = -log10(0.05),lty = 2) +
  geom_vline(xintercept = 0,lty = 2) +
  theme_bw()
ggsave(plot = DA_volcano,
       filename =  here("plots","01_Depth_DE_Testing","Differential_Abundance","DA_Volcano.pdf"),
       height = 8,width = 8)


###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
# [1] "Reproducibility information:"
# [1] "2025-03-17 14:13:17 EDT"
# user   system  elapsed 
# 53.222    4.761 2212.802 
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.2 (2023-10-31)
# os       Rocky Linux 9.4 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2025-03-17
# pandoc   3.1.3 @ /jhpce/shared/libd/core/r_nac/1.0/nac_env/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# abind                  1.4-5     2016-07-21 [1] CRAN (R 4.3.2)
# Biobase              * 2.62.0    2023-10-24 [1] Bioconductor
# BiocGenerics         * 0.48.1    2023-11-01 [1] Bioconductor
# bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.3.2)
# cli                    3.6.2     2023-12-11 [1] CRAN (R 4.3.2)
# colorspace             2.1-0     2023-01-23 [1] CRAN (R 4.3.0)
# crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
# DelayedArray           0.28.0    2023-10-24 [1] Bioconductor
# dplyr                  1.1.4     2023-11-17 [1] CRAN (R 4.3.2)
# edgeR                * 4.0.3     2023-12-10 [1] Bioconductor 3.18 (R 4.3.2)
# fansi                  1.0.6     2023-12-08 [1] CRAN (R 4.3.2)
# farver                 2.1.1     2022-07-06 [1] CRAN (R 4.3.0)
# generics               0.1.3     2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb         * 1.38.1    2023-11-08 [1] Bioconductor
# GenomeInfoDbData       1.2.11    2023-12-12 [1] Bioconductor
# GenomicRanges        * 1.54.1    2023-10-29 [1] Bioconductor
# ggplot2              * 3.5.1     2024-04-23 [1] CRAN (R 4.3.2)
# ggrepel              * 0.9.4     2023-10-13 [1] CRAN (R 4.3.2)
# glue                   1.7.0     2024-01-09 [1] CRAN (R 4.3.2)
# gtable                 0.3.4     2023-08-21 [1] CRAN (R 4.3.1)
# here                 * 1.0.1     2020-12-13 [1] CRAN (R 4.3.2)
# IRanges              * 2.36.0    2023-10-24 [1] Bioconductor
# labeling               0.4.3     2023-08-29 [1] CRAN (R 4.3.1)
# lattice                0.22-5    2023-10-24 [1] CRAN (R 4.3.1)
# lifecycle              1.0.4     2023-11-07 [1] CRAN (R 4.3.2)
# limma                * 3.58.1    2023-10-31 [1] Bioconductor
# locfit                 1.5-9.8   2023-06-11 [1] CRAN (R 4.3.2)
# magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.3.0)
# Matrix                 1.6-4     2023-11-30 [1] CRAN (R 4.3.2)
# MatrixGenerics       * 1.14.0    2023-10-24 [1] Bioconductor
# matrixStats          * 1.2.0     2023-12-11 [1] CRAN (R 4.3.2)
# munsell                0.5.0     2018-06-12 [1] CRAN (R 4.3.0)
# pillar                 1.9.0     2023-03-22 [1] CRAN (R 4.3.0)
# pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
# R6                     2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
# ragg                   1.2.7     2023-12-11 [1] CRAN (R 4.3.2)
# Rcpp                   1.0.12    2024-01-09 [1] CRAN (R 4.3.2)
# RCurl                  1.98-1.13 2023-11-02 [1] CRAN (R 4.3.2)
# rlang                  1.1.3     2024-01-10 [1] CRAN (R 4.3.2)
# rprojroot              2.0.4     2023-11-05 [1] CRAN (R 4.3.2)
# S4Arrays               1.2.0     2023-10-24 [1] Bioconductor
# S4Vectors            * 0.40.2    2023-11-23 [1] Bioconductor 3.18 (R 4.3.2)
# scales                 1.3.0     2023-11-28 [1] CRAN (R 4.3.2)
# sessioninfo          * 1.2.2     2021-12-06 [1] CRAN (R 4.3.2)
# SingleCellExperiment * 1.24.0    2023-10-24 [1] Bioconductor
# SparseArray            1.2.2     2023-11-07 [1] Bioconductor
# statmod                1.5.0     2023-01-06 [1] CRAN (R 4.3.2)
# SummarizedExperiment * 1.32.0    2023-10-24 [1] Bioconductor
# systemfonts            1.0.5     2023-10-09 [1] CRAN (R 4.3.1)
# textshaping            0.3.7     2023-10-09 [1] CRAN (R 4.3.1)
# tibble                 3.2.1     2023-03-20 [1] CRAN (R 4.3.0)
# tidyselect             1.2.0     2022-10-10 [1] CRAN (R 4.3.0)
# utf8                   1.2.4     2023-10-22 [1] CRAN (R 4.3.1)
# vctrs                  0.6.5     2023-12-01 [1] CRAN (R 4.3.2)
# withr                  2.5.2     2023-10-30 [1] CRAN (R 4.3.1)
# XVector                0.42.0    2023-10-24 [1] Bioconductor
# zlibbioc               1.48.0    2023-10-24 [1] Bioconductor
# 
# [1] /jhpce/shared/libd/core/r_nac/1.0/nac_env/lib/R/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
