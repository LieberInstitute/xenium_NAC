#!/usr/bin/env bash
set -euo pipefail

OUT="/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/processed-data/Xenium_QC/outlier_summary.csv"
mkdir -p "$(dirname "$OUT")"
echo "Sample,exclude_low_lib,exclude_any_neg,cell_area_outliers,subsets_any_neg_percent_outliers" > "$OUT"

for GF in /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/processed-data/Xenium_QC/Sample_outliers/*_global_outliers.csv; do
  SAMPLE=$(basename "$GF" _global_outliers.csv)

  # find column indices in the header
  mapfile -t COLS < <(head -1 "$GF" | tr ',' '\n')
  for FLAG in exclude_low_lib exclude_any_neg cell_area_outliers; do
    # find index (1-based)
    declare "IDX_${FLAG}"=$((  $(printf "%s\n" "${COLS[@]}" | grep -n "^${FLAG}$" | cut -d: -f1)  ))
  done

  # count TRUE rows for each flag
  excl_low=$(awk -F, -v c="${IDX_exclude_low_lib}" 'NR>1 && $c=="TRUE"{n++} END{print n+0}' "$GF")
  excl_neg=$(awk -F, -v c="${IDX_exclude_any_neg}"   'NR>1 && $c=="TRUE"{n++} END{print n+0}' "$GF")
  excl_cell=$(awk -F, -v c="${IDX_cell_area_outliers}" 'NR>1 && $c=="TRUE"{n++} END{print n+0}' "$GF")

  # now the corresponding local file
  LF="${GF%_global_outliers.csv}_local_outliers.csv"
  if [[ ! -f "$LF" ]]; then
    echo "WARNING: no local‐outliers file for $SAMPLE" >&2
    local_count=0
  else
    # subtract 1 to drop the header
    local_count=$(( $(wc -l < "$LF") - 1 ))
  fi

  # append to summary
  echo "${SAMPLE},${excl_low},${excl_neg},${excl_cell},${local_count}" >> "$OUT"
done

echo "Done → $OUT"
