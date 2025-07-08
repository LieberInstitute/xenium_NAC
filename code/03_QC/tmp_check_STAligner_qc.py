# conda activate STAligner
import os
import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad("/dcs04/hicks/data/jianing/NAc_AP/STAligner_analysis/STAligner_adata/adata_concat_STAligner_default_mclust.h5ad")
cols_to_save = ["slice_name", "cell_id", "mclust_20"]
adata.obs[cols_to_save].to_csv("/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/processed-data/05_Clustering/STAligner_initial_r15_mclust20.csv", index=True)

adata = sc.read_h5ad("/dcs04/hicks/data/jianing/NAc_AP/STAligner_analysis/STAligner_adata/adata_concat_STAligner_smallR_mclust.h5ad")
cols_to_save = ["slice_name", "cell_id", "mclust_20"]
adata.obs[cols_to_save].to_csv("/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/processed-data/05_Clustering/STAligner_initial_r0_mclust20.csv", index=True)
