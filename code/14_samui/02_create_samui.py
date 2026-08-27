#  /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
import os
os.chdir('/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/')
from pathlib import Path
from pyhere import here
import json
import os
import re
import scanpy as sc

import numpy as np
import pandas as pd
from rasterio import Affine
from tifffile import imread, imwrite
from loopy.sample import Sample
from loopy.utils.utils import remove_dupes, Url
import matplotlib.pyplot as plt

sample_info_path = here('processed-data','14_samui', 'sample_list_n22.txt')
spe_path = here("processed-data", "14_samui", "spe.h5ad")
spg = sc.read(spe_path)

#sample_info = pd.read_excel(sample_info_path)
sample_info = pd.read_csv(sample_info_path, sep='\t')
#sample_info["imagepath"] = (sample_info["obj_id"] + "__" +sample_info["sample_id"] + "_domainbounds.png")

SAMP = int(os.environ['SLURM_ARRAY_TASK_ID'])-1
this_sample = sample_info.iloc[SAMP]["ids"]
spgP = spg[spg.obs['Sample'] == this_sample, :]

gene_df = pd.DataFrame(
    spgP.X.toarray(),
    index = spgP.obs.index,
    columns = spgP.var['Symbol']
)

out_dir = Path('/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/processed-data/14_samui')/ sample_info.iloc[SAMP]["ids"]
out_dir.mkdir(exist_ok = True)
this_sample = Sample(name = sample_info.iloc[SAMP]["ids"], path = out_dir)

this_sample.add_coords(
    spgP.obsm["spatial"].rename(
        columns = {"x_final": "x", "y_final": "y"}
    ),
    name = "coords",  size = 2e-5 , mPerPx = 1e-6
)

this_sample.add_chunked_feature(gene_df, name = "Genes", coordName = "coords", dataType = "quantitative")

this_sample.add_csv_feature(
    spgP.obs["CellTypes"], name = "Cell Types", coordName = "coords",
    dataType = "categorical"
)

#img_path = Path('/dcs04/lieber/marmaypag/spatialHYP_LIBD4195/spatial_HYP/xenium_HYP/processed-data/11-XeniumSPE_imageunderlays') / sample_info.iloc[SAMP]["imagepath"]

##arr = plt.imread(str(img_path))    
#arr_rgb = arr[..., :3]           # (Y, X, 3)
#arr_rgb = (arr_rgb * 255).round().astype(np.uint8)
#arr_rgb = np.moveaxis(arr, 0, -1)        # -> (Y, X, 3)
#print("Original shape:", arr.shape)
#print("RGB shape:", arr_rgb.shape)

#new_path = Path(str(img_path)).with_name(Path(str(img_path)).stem + "_YXC.tif")
#imwrite(str(new_path), arr_rgb)

#this_sample.add_image(tiff = new_path, channels = "rgb", scale = 1e-6)
this_sample.set_default_feature(group = "Genes", feature = "DRD1")
this_sample.write()

features_of_interest = [{"feature":"PPP1R1B","group":"Genes"},
                        {"feature":"DRD2","group":"Genes"}]
#                        {"feature":"SOX10","group":"Cell Types"}]
with open(here(out_dir,'sample.json'), 'r') as json_file:
    data = json.load(json_file)

# Replace the "importantFeatures" value
data['overlayParams']['importantFeatures'] = features_of_interest

# Write the modified data back to the JSON file
with open(here(out_dir,'sample.json'), 'w') as json_file:
    json.dump(data, json_file, indent=4)