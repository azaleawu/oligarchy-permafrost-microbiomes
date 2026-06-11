#!/usr/bin/env bash
set -euo pipefail

# Estimate MAG relative abundance/TPM with CoverM.
# MAG FASTA records are renamed with MAG identifiers to preserve contig-to-genome mapping.

CONFIG=${1:-config/config_template.env}
source "$CONFIG"

ABUND_OUT="$OUT_DIR/abundance"
mkdir -p "$ABUND_OUT/renamed_mags" "$ABUND_OUT/mapping"

COMBINED_MAG_FASTA="$ABUND_OUT/all_mags_prefixed.fna"
: > "$COMBINED_MAG_FASTA"

while IFS=$'\t' read -r mag_id fasta; do
  [[ "$mag_id" == "mag_id" ]] && continue
  awk -v prefix="$mag_id" '/^>/ {sub(/^>/, ">" prefix "|"); print; next} {print}' "$fasta" \
    > "$ABUND_OUT/renamed_mags/${mag_id}.prefixed.fna"
  cat "$ABUND_OUT/renamed_mags/${mag_id}.prefixed.fna" >> "$COMBINED_MAG_FASTA"
done < "$MAG_TABLE"

bowtie2-build "$COMBINED_MAG_FASTA" "$ABUND_OUT/mapping/mag_catalogue_index"

while IFS=$'\t' read -r sample_id read1 read2; do
  [[ "$sample_id" == "sample_id" ]] && continue
  echo "[INFO] Mapping reads for MAG abundance: $sample_id"

  bowtie2 \
    -x "$ABUND_OUT/mapping/mag_catalogue_index" \
    -1 "$read1" \
    -2 "$read2" \
    -p "$THREADS" \
    -S "$ABUND_OUT/mapping/${sample_id}.sam"

  samtools view -@ "$THREADS" -bS "$ABUND_OUT/mapping/${sample_id}.sam" | \
    samtools sort -@ "$THREADS" -o "$ABUND_OUT/mapping/${sample_id}.sorted.bam"
  samtools index "$ABUND_OUT/mapping/${sample_id}.sorted.bam"
  rm -f "$ABUND_OUT/mapping/${sample_id}.sam"
done < "$SAMPLE_TABLE"

coverm genome \
  --bam-files "$ABUND_OUT"/mapping/*.sorted.bam \
  --genome-fasta-files "$ABUND_OUT"/renamed_mags/*.prefixed.fna \
  --methods tpm relative_abundance \
  --threads "$THREADS" \
  --output-file "$ABUND_OUT/mag_abundance_coverm.tsv"
