#!/usr/bin/env bash
set -euo pipefail

# Estimate MAG abundance by mapping reads to a combined, genome-prefixed MAG reference and summarizing TPM with CoverM.
# Input: MAG manifest with columns: sample_id, mag_dir, read1, read2
# Usage: bash scripts/03_MAG_abundance_coverm.sh config/mags.tsv config/config.env

MAG_MANIFEST=${1:?"Usage: bash 03_MAG_abundance_coverm.sh config/mags.tsv config/config.env"}
CONFIG=${2:?"Missing config.env"}
source "${CONFIG}"

mkdir -p "${OUT_DIR}/abundance" "${LOG_DIR}"

rename_fasta_with_genome_prefix() {
    local genome_file="$1"
    local genome_id="$2"
    awk -v gid="${genome_id}" '
        /^>/ {sub(/^>/, ""); print ">" gid "|" $0; next}
        {print}
    ' "${genome_file}"
}

while IFS=$'\t' read -r SAMPLE MAG_DIR R1 R2; do
    [[ "${SAMPLE}" == "sample_id" ]] && continue
    [[ -z "${SAMPLE}" ]] && continue

    echo "[INFO] Estimating MAG abundance for ${SAMPLE}"
    SAMPLE_OUT="${OUT_DIR}/abundance/${SAMPLE}"
    mkdir -p "${SAMPLE_OUT}/logs" "${SAMPLE_OUT}/prefixed_mags"

    COMBINED_FASTA="${SAMPLE_OUT}/${SAMPLE}.MAGs.prefixed.fna"
    : > "${COMBINED_FASTA}"

    for GENOME in "${MAG_DIR}"/*.fna "${MAG_DIR}"/*.fa; do
        [[ -e "${GENOME}" ]] || continue
        BASE=$(basename "${GENOME}")
        BASE=${BASE%.fna}
        BASE=${BASE%.fa}
        PREFIXED="${SAMPLE_OUT}/prefixed_mags/${BASE}.prefixed.fna"
        rename_fasta_with_genome_prefix "${GENOME}" "${BASE}" > "${PREFIXED}"
        cat "${PREFIXED}" >> "${COMBINED_FASTA}"
    done

    bowtie2-build "${COMBINED_FASTA}" "${SAMPLE_OUT}/MAG_db" \
      > "${SAMPLE_OUT}/logs/bowtie2_build.log" 2>&1

    bowtie2 \
      -x "${SAMPLE_OUT}/MAG_db" \
      -1 "${R1}" \
      -2 "${R2}" \
      --sensitive \
      -p "${THREADS}" \
      -S "${SAMPLE_OUT}/MAG_alignment.sam" \
      > "${SAMPLE_OUT}/logs/bowtie2_mapping.log" 2>&1

    samtools view -@ "${THREADS}" -F 4 -bS "${SAMPLE_OUT}/MAG_alignment.sam" | \
      samtools sort -@ "${THREADS}" -o "${SAMPLE_OUT}/MAG_alignment.sorted.bam" -
    samtools index "${SAMPLE_OUT}/MAG_alignment.sorted.bam"
    rm -f "${SAMPLE_OUT}/MAG_alignment.sam"

    coverm genome \
      -m tpm \
      -b "${SAMPLE_OUT}/MAG_alignment.sorted.bam" \
      --genome-fasta-directory "${SAMPLE_OUT}/prefixed_mags" \
      > "${SAMPLE_OUT}/${SAMPLE}.MAG_TPM.tsv"

    echo "[INFO] Finished MAG abundance for ${SAMPLE}"
done < "${MAG_MANIFEST}"
