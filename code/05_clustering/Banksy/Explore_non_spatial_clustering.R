# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5
#Goal: Explore resolutions of clustering. 
#All clusters are at lambda of 0 (non-spatial/cell typing mode)

#load libraries
library(SingleCellExperiment)
library(escheR)
library(dplyr)
library(here)

#How many files? 
files <- list.files(path = here("processed-data","05_Clustering","Banksy_Results","NonSpatial"),
           pattern = "Lambda0_*",full.names = TRUE)

#Read in the files
files_list <- lapply(X = files,FUN = read.csv)
names(files_list) <- basename(files)
names_cols <-  c("row","cluster","cell_id")
files_list <- lapply(files_list, setNames, names_cols)

#How many clusters in each set
p <- lapply(files_list,FUN = function(x){
  length(unique(x$cluster)) 
}) %>% as.data.frame() %>% 
  t() %>% 
  as.data.frame() %>% 
  select("Number_of_clusters" = 1) %>% 
  tibble::rownames_to_column("Settings") %>% 
  mutate(Settings = as.character(lapply(strsplit(Settings,split = ".csv"),"[",1))) %>%
  ggplot(aes(x = Settings, y = Number_of_clusters)) +
  geom_bar(stat = "identity") +
  ggtitle("Number of Clusters Per Resolution Value\nNon-spatial clustering") +
  labs(y = "Number of Clusters") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

ggsave(plot = p,
       file = here("plots","05_clustering","Banksy","Non_spatial_Number_clusters.pdf"))
