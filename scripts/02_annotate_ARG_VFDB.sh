#!/usr/bin/env bash
set -euo pipefail

# Annotate MAGs for ARGs using RGI/CARD and virulence-associated genes using ABRicate/VFDB.
# The script expects a MAG table with columns: mag_id, fasta.

CONFIG=${1:-config/config_template.env}
source "$CONFIG"

ANNOT_OUT="$OUT_DIR/annotations"
mkdir -p "$ANNOT_OUT"/{prodigal,rgi,abricate_vfdb,filtered}

while IFS=$'\t' read -r mag_id fasta; do
  [[ "$mag_id" == "mag_id" ]] && continue
  echo "[INFO] Annotating MAG: $mag_id"

  prodigal \
    -i "$fasta" \
    -a "$ANNOT_OUT/prodigal/${mag_id}.faa" \
    -d "$ANNOT_OUT/prodigal/${mag_id}.ffn" \
    -o "$ANNOT_OUT/prodigal/${mag_id}.gff" \
    -p meta

  rgi main \
    --input_sequence "$ANNOT_OUT/prodigal/${mag_id}.faa" \
    --output_file "$ANNOT_OUT/rgi/${mag_id}" \
    --input_type protein \
    --alignment_tool DIAMOND \
    --clean \
    --num_threads "$THREADS" \
    --include_loose false

  abricate \
    --db vfdb \
    --minid "$VF_MIN_IDENTITY" \
    --mincov "$VF_MIN_COVERAGE" \
    "$fasta" > "$ANNOT_OUT/abricate_vfdb/${mag_id}.vfdb.tsv"

done < "$MAG_TABLE"

# Create ABRicate summary table for VFDB annotations.
abricate --summary "$ANNOT_OUT"/abricate_vfdb/*.vfdb.tsv > "$ANNOT_OUT/abricate_vfdb/vfdb_summary.tsv"
