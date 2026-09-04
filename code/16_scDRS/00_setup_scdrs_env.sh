#!/bin/bash
# One-time setup of the conda environment used by the scDRS and drug2cell
# pipelines (code/16_scDRS and code/17_drug_analysis).
#
# Run interactively on a login/transfer node (needs internet access):
#   bash code/16_scDRS/00_setup_scdrs_env.sh

set -euo pipefail

module load conda/3-24.3.0

conda create -y -n scdrs_nac python=3.11
source activate scdrs_nac

pip install \
    "setuptools<81" \
    "scdrs==1.0.2" \
    "scanpy>=1.9" \
    "anndata>=0.9" \
    drug2cell \
    requests \
    session_info

python -c "import scdrs, scanpy, drug2cell; print('scdrs', scdrs.__version__); print('scanpy', scanpy.__version__)"

echo "Environment 'scdrs_nac' ready."
