# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SpatialFeatureExperiment)
library(SingleCellExperiment)
library(SpatialExperiment)
library(sessioninfo)
library(HDF5Array)
library(ggplot2)
library(escheR)
library(tidyr)
library(dplyr)
library(here)



here::i_am("code/08_Visium_HD/03_LabelTransfer/SpaceRanger/compare_RCTD_singleR.R")
# Load object 
sfe_dir <- here("processed-data", "HD_Full_Analysis", "sfe_with_labels")

sfe <- loadHDF5SummarizedExperiment(sfe_dir)

sfe

table(is.na(sfe$snRNA_label))

sfe <- sfe[,!is.na(sfe$snRNA_label)]

sfe

#Make a bargraoh of cell type by spot_class
x <- as.data.frame(table(sfe$snRNA_label,sfe$rctd_spot_class))
colnames(x) <- c("Cell_Type","spot_class","Freq")
#Make a bargraph of number of cells by 
p1 <- ggplot(data = x,aes(x = Cell_Type, y = Freq,fill = spot_class))+
  geom_bar(stat = "identity",position = "dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = p1,filename = here("plots","HD_Full_Analysis","LabelTransfer","SpaceRanger",
                                 "CellType_SpotClass_Bargraph.pdf"),
       height = 8, width = 12)

y <- as.data.frame(prop.table(table(sfe$sample_id,sfe$rctd_spot_class),margin = 1)*100)
colnames(y) <- c("Sample","spot_class","Percentage")
p1 <- ggplot(data = y,aes(x = Sample, y = Percentage,fill = spot_class))+
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = p1,filename = here("plots","HD_Full_Analysis","LabelTransfer","SpaceRanger",
                                 "sampleid_SpotClass_Bargraph.pdf"),
       height = 8, width = 8)

# pull the relevant columns
meta <- as.data.frame(
  colData(sfe)[, c("rctd_first_type", "rctd_second_type", "rctd_spot_class")]
)

# keep only doublets — first/second pairing is only meaningful here
meta_db <- meta %>%
  filter(rctd_spot_class %in% c("doublet_certain", "doublet_uncertain"))

# When there are doublets, what is the predicted second cell type? 
# I anticipate that doublets will be cells in which neurons are surrounded by many astros or oligos
# This will caause a shared transcriptional signature 
# count pairs, fill empty combos with 0, then convert to % within each first_type
pct <- meta_db %>%
  count(rctd_first_type, rctd_second_type, name = "n") %>%
  complete(rctd_first_type, rctd_second_type, fill = list(n = 0)) %>%
  group_by(rctd_first_type) %>%
  mutate(total = sum(n),
         pct   = ifelse(total > 0, 100 * n / total, 0)) %>%
  ungroup()

# plot
p2 <- ggplot(pct, aes(x = rctd_first_type, y = rctd_second_type, fill = pct)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = ifelse(pct > 0, round(pct), "")), size = 2.5) +
  scale_fill_gradientn(colours = c("white","lightgrey","orange","red") )+
  labs(x = "First type", y = "Second type", fill = "% of first type") +
  coord_equal() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plot = p2,filename = here("plots","HD_Full_Analysis","LabelTransfer","SpaceRanger",
                                 "RCTD_First_Second_Type_Comparison_Dou.pdf"),
       height = 8, width = 12)



# Now for the doublets, does the SingleR call agree with the first_type call from RCTD? 
# This will help determine if we just go ahead and move forward with the RCTD calls. 
# pull the relevant columns
meta <- as.data.frame(
  colData(sfe)[, c("rctd_first_type", "rctd_second_type","pruned.labels","rctd_spot_class")]
)

# keep only doublets — first/second pairing is only meaningful here
meta_db <- meta %>%
  filter(rctd_spot_class %in% c("doublet_certain", "doublet_uncertain"))

# put both calls on the SAME factor levels so the diagonal = agreement
lvls <- sort(union(as.character(meta_db$rctd_first_type),
                   as.character(meta_db$pruned.labels)))
meta_db$rctd_first_type   <- factor(meta_db$rctd_first_type,   levels = lvls)
meta_db$pruned.labels <- factor(meta_db$pruned.labels, levels = lvls)

# confusion matrix -> proportion within each first_type (columns sum to 1)
prop <- meta_db %>%
  count(rctd_first_type, pruned.labels, name = "n") %>%
  complete(rctd_first_type, pruned.labels, fill = list(n = 0)) %>%
  group_by(rctd_first_type) %>%
  mutate(total = sum(n),
         prop  = ifelse(total > 0, n / total, 0)) %>%
  ungroup()

p3 <- ggplot(prop, aes(x = rctd_first_type, y = pruned.labels, fill = prop)) +
  geom_tile(color = "black") +
  # uncomment to label tiles:
  geom_text(aes(label = ifelse(prop > 0.01, round(prop, 2), "")), size = 2.5) +
  scale_fill_gradientn(colours = c("white","lightgrey","orange","red")) +
  labs(x = "First type (doublet mode)",
       y = "Singlet call",
       fill = "Proportion\nof first type") +
  coord_equal() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plot = p3,filename = here("plots","HD_Full_Analysis","LabelTransfer","SpaceRanger",
                                 "RCTD_SingleR_Comparison_Doublets_Only.pdf"),
       height = 8, width = 12)


sessionInfo()
