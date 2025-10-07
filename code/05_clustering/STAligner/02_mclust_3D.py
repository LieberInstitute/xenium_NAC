# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC
# conda activate STAligner

import sys
DONOR = sys.argv[1]
NUM_CLUST = int(sys.argv[2])

import warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
import subprocess
from pathlib import Path

# import STAligner.ST_utils
# import train_STAligner
import STAligner

# the location of R (used for the mclust clustering)
import os
os.environ['R_HOME'] = "/users/jyao/.conda/envs/proust-310/lib/R"
os.environ['R_USER'] = "/users/jyao/.conda/envs/STAligner/lib/python3.9/site-packages/rpy2"
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri

import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
import scipy.linalg
from scipy.sparse import csr_matrix

git_root = Path(subprocess.check_output(
    ["git", "rev-parse", "--show-toplevel"], text=True).strip()
)

# mclust
adata_concat = sc.read_h5ad(git_root / "processed-data" / "05_Clustering" / f"STAligner_r0_{DONOR}.h5ad")
STAligner.mclust_R(adata_concat, num_cluster=NUM_CLUST, used_obsm='STAligner')
adata_concat.obs['mclust'].to_csv(
    git_root / "processed-data" / "05_Clustering" / f"STAligner_r0_{DONOR}_{NUM_CLUST}_clusters.csv",
    header=True,
    index=True,
    index_label="cell"
)
STAligner_emb = adata_concat.obsm['STAligner']
np.savetxt(git_root / "processed-data" / "05_Clustering" / f"STAligner_r0_{DONOR}_{NUM_CLUST}_embeddings.tsv", STAligner_emb, delimiter='\t')

# Batch effect correction
sc.pp.neighbors(adata_concat, use_rep='STAligner', random_state=666)
sc.tl.umap(adata_concat, random_state=666)
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams["figure.figsize"] = (3, 3)
plt.rcParams['font.size'] = 10
ax = sc.pl.umap(adata_concat, color="Sample", show=False)
plt.savefig(git_root / "plots" / "05_clustering" / "STAligner" / f"STAligner_r0_{DONOR}_{NUM_CLUST}_umap.png", dpi=300, bbox_inches="tight")
plt.close()


# Plot individual slice
sp = adata_concat.obsm.get('spatial')
if isinstance(sp, pd.DataFrame):
    adata_concat.obsm['spatial'] = sp.values

section_ids = adata_concat.obs['Sample'].cat.remove_unused_categories().cat.categories.tolist()
spot_size = 70
for section_id in section_ids:
    print(section_id)
    ad = adata_concat[adata_concat.obs['Sample'] == section_id].copy()
    fig, ax = plt.subplots(figsize=(12, 12), dpi=300)
    sc.pl.spatial(ad, title=section_id, color="mclust", frameon=False, spot_size=spot_size, legend_fontsize=20, show=False, ax=ax)
    ax.set_title(section_id, fontsize=20)
    fig.savefig(git_root / "plots" / "05_clustering" / "STAligner" / "Pre_Annotation" / f"{DONOR}" / f"k{NUM_CLUST}" /f"{section_id}_STAligner_r0_PreAnnotation.png", bbox_inches="tight")
    plt.close(fig)


# Perform 3D reconstruction: use white matter as landmark
if DONOR == "Br6436":
    landmark_domain = 4

if DONOR == "Br6660":
    landmark_domain = 3

if DONOR == "All":
    landmark_domain = 4 #WM in Br6660 but not in Br6436

for section_id in section_ids:
    print(section_id)
    ad = adata_concat[adata_concat.obs['Sample'] == section_id].copy()
    fig, ax = plt.subplots(figsize=(12, 12), dpi=300)
    sc.pl.spatial(ad, title=section_id, color="mclust", frameon=False, spot_size=spot_size, legend_fontsize=20, groups=[landmark_domain], show=False, ax=ax)
    ax.set_title(section_id, fontsize=20)
    fig.savefig(git_root / "plots" / "05_clustering" / "STAligner" / "3D_Reconstruction" / f"{DONOR}" / f"{section_id}_landmark.png", bbox_inches="tight")
    plt.close(fig)

# if DONOR == "All":
#     iter_comb = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 9), (9, 7), (7, 8), (8, 10), 
#                  (11, 0),
#                  (12, 11), (13, 12), (14, 13), (15, 14), (16, 15), (17, 16), (18, 17), (19, 18), (20, 19), (21, 20)]
    
if DONOR == "Br6436":
    iter_comb = [(1, 0), (2, 1), (3, 2), (4, 3), (5, 4), (6, 5),
                 (7, 6), (8, 7), (9, 8), (10, 9)]
    
if DONOR == "Br6660":
    iter_comb = [(1, 0), (2, 1), (3, 2), (4, 3), (5, 4), (6, 5),
                 (9, 6), (7, 9), (8, 7), (10, 8)]
    
adata_concat.obs['louvain'] = adata_concat.obs['mclust'].values.astype('category')
landmark_domain_list = [landmark_domain]
Batch_list = []
base_ref = adata_concat[adata_concat.obs['Sample'] == section_ids[0]].copy()
Batch_list.append(base_ref)
for comb in iter_comb:
    print(comb)
    i, j = comb[0], comb[1]
    adata_target = adata_concat[adata_concat.obs['Sample'] == section_ids[i]].copy()
    adata_ref = adata_concat[adata_concat.obs['Sample'] == section_ids[j]].copy()
    slice_target = section_ids[i]
    slice_ref = section_ids[j]
    aligned_coor = STAligner.ICP_align(adata_concat, adata_target, adata_ref, slice_target, slice_ref, landmark_domain_list)
    adata_target.obsm["spatial"] = aligned_coor
    Batch_list.append(adata_target)


# Plot 3D reconstruction
adata_concat.obs['Z'] = list(Batch_list[0].shape[0] * [0]) + list(Batch_list[1].shape[0] * [100]) \
+ list(Batch_list[2].shape[0] * [200]) + list(Batch_list[3].shape[0] * [300])+ list(Batch_list[4].shape[0] * [400]) \
+ list(Batch_list[5].shape[0] * [500]) + list(Batch_list[6].shape[0] * [600]) + list(Batch_list[7].shape[0] * [700]) \
+ list(Batch_list[8].shape[0] * [800]) + list(Batch_list[9].shape[0] * [900]) + list(Batch_list[10].shape[0] * [1000])
adata_concat.obs['x_centroid'] = adata_concat.obsm['spatial'][:, 0]
adata_concat.obs['y_centroid'] = adata_concat.obsm['spatial'][:, 1]

All_coor = adata_concat.obs[['x_centroid', 'y_centroid', 'Z']].copy()
if DONOR == "Br6436":
    All_coor.loc[adata_concat.obs['batch_name'] == section_ids[0], All_coor.columns[:2]] = Batch_list[0].obsm["spatial"]
    All_coor.loc[adata_concat.obs['batch_name'] == section_ids[1], All_coor.columns[:2]] = Batch_list[1].obsm["spatial"]
    All_coor.loc[adata_concat.obs['batch_name'] == section_ids[2], All_coor.columns[:2]] = Batch_list[2].obsm["spatial"]
    All_coor.loc[adata_concat.obs['batch_name'] == section_ids[3], All_coor.columns[:2]] = Batch_list[3].obsm["spatial"]
    All_coor.loc[adata_concat.obs['batch_name'] == section_ids[4], All_coor.columns[:2]] = Batch_list[4].obsm["spatial"]
    All_coor.loc[adata_concat.obs['batch_name'] == section_ids[5], All_coor.columns[:2]] = Batch_list[5].obsm["spatial"]
    All_coor.loc[adata_concat.obs['batch_name'] == section_ids[6], All_coor.columns[:2]] = Batch_list[6].obsm["spatial"]
    All_coor.loc[adata_concat.obs['batch_name'] == section_ids[7], All_coor.columns[:2]] = Batch_list[7].obsm["spatial"]
    All_coor.loc[adata_concat.obs['batch_name'] == section_ids[8], All_coor.columns[:2]] = Batch_list[8].obsm["spatial"]
    All_coor.loc[adata_concat.obs['batch_name'] == section_ids[9], All_coor.columns[:2]] = Batch_list[9].obsm["spatial"]
    All_coor.loc[adata_concat.obs['batch_name'] == section_ids[10], All_coor.columns[:2]] = Batch_list[10].obsm["spatial"]
if DONOR == "Br6660":
    All_coor.loc[adata_concat.obs['batch_name'] == section_ids[0], All_coor.columns[:2]] = Batch_list[0].obsm["spatial"]
    All_coor.loc[adata_concat.obs['batch_name'] == section_ids[1], All_coor.columns[:2]] = Batch_list[1].obsm["spatial"]
    All_coor.loc[adata_concat.obs['batch_name'] == section_ids[2], All_coor.columns[:2]] = Batch_list[2].obsm["spatial"]
    All_coor.loc[adata_concat.obs['batch_name'] == section_ids[3], All_coor.columns[:2]] = Batch_list[3].obsm["spatial"]
    All_coor.loc[adata_concat.obs['batch_name'] == section_ids[4], All_coor.columns[:2]] = Batch_list[4].obsm["spatial"]
    All_coor.loc[adata_concat.obs['batch_name'] == section_ids[5], All_coor.columns[:2]] = Batch_list[5].obsm["spatial"]
    All_coor.loc[adata_concat.obs['batch_name'] == section_ids[6], All_coor.columns[:2]] = Batch_list[6].obsm["spatial"]
    All_coor.loc[adata_concat.obs['batch_name'] == section_ids[9], All_coor.columns[:2]] = Batch_list[7].obsm["spatial"]
    All_coor.loc[adata_concat.obs['batch_name'] == section_ids[7], All_coor.columns[:2]] = Batch_list[8].obsm["spatial"]
    All_coor.loc[adata_concat.obs['batch_name'] == section_ids[8], All_coor.columns[:2]] = Batch_list[9].obsm["spatial"]
    All_coor.loc[adata_concat.obs['batch_name'] == section_ids[10], All_coor.columns[:2]] = Batch_list[10].obsm["spatial"]

fig = plt.figure(figsize=(20, 20))
ax1 = plt.axes(projection='3d')

for it, label in enumerate(np.unique(adata_concat.obs['mclust'])):
    temp_Coor = All_coor.loc[adata_concat.obs['mclust'] == label, :]
    temp_xd = temp_Coor['x_centroid']
    temp_yd = temp_Coor['y_centroid']
    temp_zd = temp_Coor['Z']
    if label == landmark_domain:
        ax1.scatter3D(temp_xd, temp_yd, temp_zd,
                      s=0.02, marker="o", label=label, alpha=1)
    else:
        ax1.scatter3D(temp_xd, temp_yd, temp_zd,
                      s=0.02, marker="o", label=label, alpha=0.05)
        
plt.legend(bbox_to_anchor=(1.2, 0.8), markerscale=10, frameon=False)
plt.title('3D stacked slices')
ax1.elev = 20
ax1.azim = -40

ax1.set_xlabel('')
ax1.set_ylabel('')
ax1.set_zlabel('')
ax1.set_xticklabels([])
ax1.set_yticklabels([])
ax1.set_zticklabels([])   

out_path = git_root / "plots" / "05_clustering" / "STAligner" / "3D_Reconstruction" / f"{DONOR}" / f"{DONOR}_STAligner_3D.png"
fig.savefig(out_path, dpi=300, bbox_inches="tight")
plt.close(fig)

