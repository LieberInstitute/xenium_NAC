#!/bin/bash
# Submit 09 (risk-gene drivers + expression tables) and chain 10 (figures)
# after it with an afterok dependency.
#
# Usage (from this directory, on a login node):
#   bash run_drivers.sh
#
# Use this when the upstream scoring pipeline (01-07) has already run and only
# the driver analysis and its figures need regenerating -- run_all.sh submits
# the same two steps as part of the full chain.

set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs

j_drivers=$(sbatch --parsable 09_risk_gene_drivers.sh)
echo "09_risk_gene_drivers    -> job ${j_drivers}"

# afterok, not afterany: if 09 fails, the expression table is missing and 10
# would otherwise run and quietly produce only the driver figures.
j_driverplots=$(sbatch --parsable --dependency=afterok:${j_drivers} 10_driver_gene_plots.sh)
echo "10_driver_gene_plots    -> job ${j_driverplots} (after ${j_drivers})"

echo
echo "Submitted. Monitor with: squeue -u ${USER}"
