import matplotlib.pyplot as plt
import scanpy as sc
import os
from pyhere import here
import session_info
import bin2cell as b2c
import datetime
import pandas as pd
import json

#   Grab this sample ID from the full list, using the array task ID
sample_info_path = here('processed-data','visiumHD_sample_info_NAc.csv')
sample_info = pd.read_csv(sample_info_path)
task_id = int(os.getenv('SLURM_ARRAY_TASK_ID')) - 1
sample_id = sample_info.iloc[task_id]['sample_id']
space_ranger_dir = sample_info.iloc[task_id]['spaceranger_dir']


#   Set directories

stardist_dir = here(
  'processed-data', 'HD_Full_Analysis', 'bin2cell', 'stardist'
)

final_out_path = here(
  'processed-data', 'HD_Full_Analysis', 'bin2cell', f'{sample_id}.h5ad'
)

pre_out_path = here(
  'processed-data','HD_Full_Analysis', 'bin2cell',
  f'{sample_id}_pre_bin2cell.h5ad'
)

sr_dir = here(
  'processed-data', '01_spaceranger', sample_id, 'outs', 'binned_outputs',
  'square_002um'
)

sr_spatial_dir = here(
  'processed-data', '01_spaceranger', sample_id, 'outs', 'spatial'
)

plot_dir = here('plots', 'HD_Full_Analysis')

raw_image_path = sample_info.iloc[task_id]['raw_image_path']

os.makedirs(stardist_dir, exist_ok=True)
os.makedirs(plot_dir, exist_ok=True)

################################################################################
#   Dynamically set mpp from scalefactors_json.json
################################################################################

scalefactors_path = here(
    'processed-data', '01_spaceranger', sample_id, 'outs',
    'binned_outputs', 'square_002um', 'spatial', 'scalefactors_json.json'
)
with open(scalefactors_path) as f:
    scalefactors = json.load(f)

mpp = scalefactors['microns_per_pixel']
print(f"{datetime.datetime.now()} | mpp for {sample_id}: {mpp}")

################################################################################
#   Build and preprocess AnnData
################################################################################

print(f"{datetime.datetime.now()} | Building and preprocessing AnnData")

#   In the tutorial at https://nbviewer.org/github/Teichlab/bin2cell/blob/main/notebooks/demo.ipynb,
#   the filtered feature matrix, with additional gene-filtering steps, is used
#   before using specific settings to segment the gene-expression-based image
#   ('secondary segmentation'). While we're interested in retaining all genes
#   (i.e. using the raw feature matrix as in 'adata'), we want similar secondary
#   segmentation behavior as in the tutorial, hence 'adata_filtered'
adata_filtered = b2c.read_visium(
  sr_dir,
  count_file = 'filtered_feature_bc_matrix.h5',
  source_image_path = raw_image_path,
  spaceranger_image_path = sr_spatial_dir
)

#   Require bins with nonzero counts
sc.pp.filter_cells(adata_filtered, min_counts=1)
sc.pp.filter_genes(adata_filtered, min_cells=3)

#   Note that we intentionally skip b2c.destripe(), which I've observed to
#   dramatically worsen the technical effect it's supposed to account for
#   (https://github.com/Teichlab/bin2cell/issues/45)
adata_filtered.obs['sum_umi'] = adata_filtered.X.sum(axis = 1)

#   Read in spaceranger outputs into an AnnData
adata = b2c.read_visium(
  sr_dir,
  count_file = 'raw_feature_bc_matrix.h5',
  source_image_path = raw_image_path,
  spaceranger_image_path = sr_spatial_dir
)

#   Only keep in-tissue bins
assert adata_filtered.obs.index.isin(adata.obs.index).all()
adata = adata[adata_filtered.obs.index, :]

#   Use Ensembl IDs for var_names
adata.var_names = adata.var['gene_ids']
adata.var_names.name = None

#   Create a scaled H&E image attached to the object (and for segmentation with
#   stardist)
b2c.scaled_he_image(
  adata,
  mpp = mpp,
  save_path = os.path.join(stardist_dir, f'he_{sample_id}.tiff')
)

#   Capture the img_key that bin2cell stored in adata.uns['spatial'] so that
#   sc.pl.spatial() calls later can reference it correctly. 
sample_key = list(adata.uns['spatial'].keys())[0]
img_key = list(adata.uns['spatial'][sample_key]['images'].keys())[0]
print(f"{datetime.datetime.now()} | img_key for {sample_id}: {img_key}")

################################################################################
#   Perform nuclear-based ("primary") segmentation
################################################################################

print(f"{datetime.datetime.now()} | Performing nuclear-based ('primary') segmentation")

#   Segment nuclei on H&E image
b2c.stardist(
  image_path=os.path.join(
    stardist_dir, f'he_{sample_id}.tiff'
  ),
  labels_npz_path=os.path.join(
    stardist_dir, f'he_{sample_id}.npz'
  ),
  stardist_model="2D_versatile_he",
  prob_thresh=0.01
)

#   Add segmentations to object
b2c.insert_labels(
  adata,
  labels_npz_path=os.path.join(
    stardist_dir, f'he_{sample_id}.npz'
  ),
  basis="spatial",
  spatial_key="spatial_cropped_150_buffer",
  mpp=mpp,
  labels_key="labels_he"
)

#   Expand labels to attempt to capture cells and not nuclei
b2c.expand_labels(
  adata,
  labels_key='labels_he',
  expanded_labels_key="labels_he_expanded"
)

################################################################################
#   Perform gene-expression-based ('secondary') segmentation
################################################################################

print(f"{datetime.datetime.now()} | Performing gene-expression-based ('secondary') segmentation")

#   Create an image from gene counts
b2c.grid_image(
  adata_filtered,
  "sum_umi",
  mpp=mpp,
  sigma=5,
  save_path=os.path.join(
    stardist_dir, f'gex_{sample_id}.tiff'
  )
)

#   Hack around a different object ('adata_filtered') being used to create the
#   gene-count image as the object ('adata') we want to insert secondary labels
#   into. 'b2c.grid_image' inserted this column silently
adata.uns['bin2cell']['array_check'] = adata_filtered.uns['bin2cell']['array_check']

#   Segment cells on the gene-count image
b2c.stardist(
  image_path=os.path.join(
    stardist_dir, f'gex_{sample_id}.tiff'
  ),
  labels_npz_path = os.path.join(
    stardist_dir, f'gex_{sample_id}.npz'
  ),
  stardist_model="2D_versatile_fluo",
  prob_thresh=0.05,
  nms_thresh=0.5
)

#   Add segmentations to object
b2c.insert_labels(
  adata,
  labels_npz_path = os.path.join(
    stardist_dir, f'gex_{sample_id}.npz'
  ),
  basis="array",
  mpp=mpp,
  labels_key="labels_gex"
)

#   Take the union of cell labels from both segmentation methods
b2c.salvage_secondary_labels(
  adata,
  primary_label="labels_he_expanded",
  secondary_label="labels_gex",
  labels_key="labels_joint"
)

#-------------------------------------------------------------------------------
#   Plot primary and secondary cells
#-------------------------------------------------------------------------------

#   Plot 3 different subregions to get a representative idea
for i in range(3):
  #   Region for plots
  mask = (
    (adata.obs['array_row'] >= 1000 + 500 * i) &
      (adata.obs['array_row'] <= 1050 + 500 * i) &
      (adata.obs['array_col'] >= 1000 + 500 * i) &
      (adata.obs['array_col'] <= 1050 + 500 * i)
  )

  #   If the region has no cells, try to iterate over other regions until
  #   cells are found
  offset = 150
  while not any(mask) and offset < 1000:
    mask = (
      (adata.obs['array_row'] >= 1000 + 500 * i + offset) &
        (adata.obs['array_row'] <= 1050 + 500 * i + offset) &
        (adata.obs['array_col'] >= 1000 + 500 * i + offset) &
        (adata.obs['array_col'] <= 1050 + 500 * i + offset)
    )
    offset += 150
  assert any(mask), "Failed to find a region with cells for plotting"

  #   Plot union of cell labels
  bdata = adata[mask]
  bdata = bdata[bdata.obs['labels_joint'] > 0]
  bdata.obs['labels_joint'] = bdata.obs['labels_joint'].astype(str)
  sc.pl.spatial(
    bdata, color=[None, "labels_joint_source", "labels_joint"],
    img_key=img_key, basis="spatial_cropped_150_buffer"
  )
  plt.savefig(
    os.path.join(plot_dir, f'{sample_id}_cells{i+1}.png')
  )
  plt.close('all')

  #   Plot primary segmentations
  crop = b2c.get_crop(
    adata[mask], basis="spatial", spatial_key="spatial_cropped_150_buffer",
    mpp=mpp
  )
  rendered = b2c.view_labels(
    image_path = os.path.join(
      stardist_dir, f'he_{sample_id}.tiff'
    ),
    labels_npz_path = os.path.join(
      stardist_dir, f'he_{sample_id}.npz'
    ),
    crop = crop
  )
  plt.imshow(rendered)
  plt.savefig(
    os.path.join(plot_dir, f'{sample_id}_primary_segmentation{i+1}.png')
  )
  plt.close('all')

  #   Plot secondary segmentations
  crop = b2c.get_crop(adata[mask], basis="array", mpp=mpp)
  rendered = b2c.view_labels(
    image_path = os.path.join(
      stardist_dir, f'gex_{sample_id}.tiff'
    ),
    labels_npz_path = os.path.join(
      stardist_dir, f'gex_{sample_id}.npz'
    ),
    crop = crop,
    stardist_normalize = True
  )
  plt.imshow(rendered)
  plt.savefig(
    os.path.join(plot_dir, f'{sample_id}_secondary_segmentation{i+1}.png')
  )
  plt.close('all')

#   Keep a copy of the AnnData before aggregation (to enable interactive
#   plotting later, for example)
sc.write(pre_out_path, adata)

################################################################################
#   Aggregate bins into cells
################################################################################

print(f"{datetime.datetime.now()} | Aggregating bins into cells")

adata = b2c.bin_to_cell(
  adata, labels_key="labels_joint",
  spatial_keys=["spatial", "spatial_cropped_150_buffer"]
)

cell_mask = (
  (adata.obs['array_row'] >= 1450) &
    (adata.obs['array_row'] <= 1550) &
    (adata.obs['array_col'] >= 250) &
    (adata.obs['array_col'] <= 450)
)

#   Plot counts within cells after aggregation of bins
bdata = adata[cell_mask]
sc.pl.spatial(
  bdata, color="bin_count", img_key=img_key,
  basis="spatial_cropped_150_buffer"
)
plt.savefig(
  os.path.join(plot_dir, f'{sample_id}_cells_aggregated.png')
)
plt.close('all')

sc.write(final_out_path, adata)
session_info.show()
