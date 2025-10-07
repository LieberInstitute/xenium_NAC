# srun --partition=gpu --mem=500G --gres=gpu:tesh100:1 --pty /bin/bash
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC
# conda activate STAligner

import sys
DONOR = sys.argv[1]
RADIUS = float(sys.argv[2])

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

import threading
import psutil
import pynvml
import time

import torch
used_device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

git_root = Path(subprocess.check_output(
    ["git", "rev-parse", "--show-toplevel"], text=True).strip()
)

# Load the data
adata_all = ad.read_h5ad(git_root / "processed-data" / "02_build_spe" / "h5ad" / "spe_NormCounts_nucleus_normcounts.h5ad")
gene_mask = adata_all.var["Type"] == "Gene Expression"
adata_all = adata_all[:, gene_mask].copy()
# subset to donor level data
if DONOR != "All":
    adata_all = adata_all[adata_all.obs["Donor"] == DONOR].copy()

Batch_list = []
adj_list = []
section_ids = adata_all.obs['Sample'].cat.remove_unused_categories().cat.categories.tolist()

# Prepare the data for STAligner
for section_id in section_ids:
    print(section_id)
    adata = adata_all[adata_all.obs["Sample"] == section_id].copy()
    adata.X = csr_matrix(adata.X)
    adata.var_names_make_unique(join="++")    
    STAligner.Cal_Spatial_Net(adata, rad_cutoff=RADIUS) 
    adj_list.append(adata.uns['adj'])
    Batch_list.append(adata)

adata_concat = ad.concat(Batch_list, label="slice_name", keys=section_ids)
adata_concat.obs["batch_name"] = adata_concat.obs["slice_name"].astype('category')

# Concat the spatial network for multiple slices
# adj_concat = np.asarray(adj_list[0].todense())
# for batch_id in range(1,len(section_ids)):
#     adj_concat = scipy.linalg.block_diag(adj_concat, np.asarray(adj_list[batch_id].todense()))
# adata_concat.uns['edgeList'] = np.nonzero(adj_concat)

# STAligner training
# helper function to trace time and memory usage
pynvml.nvmlInit()
gpu_handle = pynvml.nvmlDeviceGetHandleByIndex(0)  

def train_with_monitoring(*args, log_csv_path="resource_log.csv", interval=5.0, **kwargs):
    metrics = []
    max_proc_ram = 0.0
    max_gpu_mem = 0.0
    monitoring = True
    proc = psutil.Process(os.getpid())
    def _monitor():
        nonlocal max_proc_ram, max_gpu_mem
        t0 = time.time()
        while monitoring:
            elapsed = time.time() - t0
            cpu_pct  = psutil.cpu_percent(interval=None)
            ram_used = psutil.virtual_memory().used    / 1024**3
            rss      = proc.memory_info().rss         / 1024**3
            util     = pynvml.nvmlDeviceGetUtilizationRates(gpu_handle).gpu
            mem_info = pynvml.nvmlDeviceGetMemoryInfo(gpu_handle)
            gpu_used = mem_info.used                   / 1024**3
            metrics.append({
                "elapsed_s":       elapsed,
                "cpu_pct":         cpu_pct,
                "ram_used_GiB":    ram_used,
                "proc_rss_GiB":    rss,
                "gpu_pct":         util,
                "gpu_used_GiB":    gpu_used,
            })
            if rss      > max_proc_ram: max_proc_ram = rss
            if gpu_used > max_gpu_mem:   max_gpu_mem  = gpu_used
            time.sleep(interval)
    mon_thread = threading.Thread(target=_monitor, daemon=True)
    mon_thread.start()
    start_time = time.time()
    try:
        result = STAligner.train_STAligner_subgraph(*args, **kwargs)
    finally:
        monitoring = False
        mon_thread.join()
        total_time = time.time() - start_time
        os.makedirs(os.path.dirname(log_csv_path), exist_ok=True)
        pd.DataFrame(metrics).to_csv(log_csv_path, index=False)
        print(f"Resource log saved to: {log_csv_path}")
        print(f"Total training time: {total_time/60:.2f} minutes")
        print(f"Peak process RSS:     {max_proc_ram:.2f} GiB")
        print(f"Peak GPU VRAM used:   {max_gpu_mem:.2f} GiB")
    return result

if DONOR == "All":
    iter_comb = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 9), (9, 7), (7, 8), (8, 10), 
                 (11, 0),
                 (12, 11), (13, 12), (14, 13), (15, 14), (16, 15), (17, 16), (18, 17), (19, 18), (20, 19), (21, 20)]
if DONOR == "Br6436":
    iter_comb = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6),
                 (6, 7), (7, 8), (8, 9), (9, 10)]
if DONOR == "Br6660":
    iter_comb = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6),
                 (6, 9), (9, 7), (7, 8), (8, 10)]

import gc
gc.collect()
torch.cuda.empty_cache()

adata_concat = train_with_monitoring(
    adata_concat,
    log_csv_path=git_root / "processed-data" / "05_Clustering" / f"STAligner_r0_{DONOR}_resources.csv",
    interval=5.0,
    verbose=True,
    knn_neigh=100,
    n_epochs=600,
    iter_comb=iter_comb,
    Batch_list=Batch_list,
    device=used_device
)

adata_concat.write(git_root / "processed-data" / "05_Clustering" / f"STAligner_r0_{DONOR}.h5ad")

