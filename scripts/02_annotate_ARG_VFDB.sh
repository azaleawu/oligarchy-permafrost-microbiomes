#!/usr/bin/env bash
set -euo pipefail

# ORF prediction, ARG annotation with RGI/CARD and virulence-associated gene annotation with ABRicate/VFDB.
# Input: MAG manifest with columns: sample_id, mag_dir, read1, read2
# Usage: bash scripts/02_annotate_ARG_VFDB.sh config/mags.tsv config/config.env

MAG_MANIFEST=${1:?"Usage: bash 02_annotate_ARG_VFDB.sh config/mags.tsv config/config.env"}
CONFIG=${2:?"Missing config.env"}
source "${CONFIG}"

mkdir -p "${OUT_DIR}/annotations" "${LOG_DIR}"

# Load CARD database locally for RGI if needed.
# Run once before annotation if the local RGI database has not already been configured.
if [[ -f "${CARD_JSON}" ]]; then
    rgi load --card_json "${CARD_JSON}" --local > "${LOG_DIR}/rgi_load.log" 2>&1 || true
fi

while IFS=$'\t' read -r SAMPLE MAG_DIR R1 R2; do
    [[ "${SAMPLE}" == "sample_id" ]] && continue
    [[ -z "${SAMPLE}" ]] && continue

    echo "[INFO] Annotating MAGs for ${SAMPLE}"
    SAMPLE_OUT="${OUT_DIR}/annotations/${SAMPLE}"
    mkdir -p "${SAMPLE_OUT}/prodigal" "${SAMPLE_OUT}/rgi" "${SAMPLE_OUT}/vfdb" "${SAMPLE_OUT}/logs"

    for GENOME in "${MAG_DIR}"/*.fna "${MAG_DIR}"/*.fa; do
        [[ -e "${GENOME}" ]] || continue
        BASE=$(basename "${GENOME}")
        BASE=${BASE%.fna}
        BASE=${BASE%.fa}

        # 1. ORF prediction in metagenomic mode
        prodigal \
          -i "${GENOME}" \
          -a "${SAMPLE_OUT}/prodigal/${BASE}.faa" \
          -d "${SAMPLE_OUT}/prodigal/${BASE}.fna" \
          -p meta \
          > "${SAMPLE_OUT}/logs/${BASE}.prodigal.log" 2>&1

        # 2. RGI/CARD annotation on predicted proteins
        rgi main \
          --input_sequence "${SAMPLE_OUT}/prodigal/${BASE}.faa" \
          --output_file "${SAMPLE_OUT}/rgi/${BASE}" \
          --input_type protein \
          --alignment_tool DIAMOND \
          --local \
          --clean \
          > "${SAMPLE_OUT}/logs/${BASE}.rgi.log" 2>&1

        # 3. ABRicate/VFDB annotation on MAG nucleotide sequence
        abricate \
          --db "${VFDB_DB_NAME}" \
          --threads "${THREADS}" \
          --minid "${VFDB_MIN_IDENTITY}" \
          --mincov "${VFDB_MIN_COVERAGE}" \
          "${GENOME}" \
          > "${SAMPLE_OUT}/vfdb/${BASE}.vfdb.tsv"
    done

    # 4. Combine and summarize VFDB outputs for the sample
    if compgen -G "${SAMPLE_OUT}/vfdb/*.vfdb.tsv" > /dev/null; then
        cat "${SAMPLE_OUT}"/vfdb/*.vfdb.tsv > "${SAMPLE_OUT}/vfdb/${SAMPLE}.vfdb.raw.tsv"
        abricate --summary "${SAMPLE_OUT}/vfdb/${SAMPLE}.vfdb.raw.tsv" \
          > "${SAMPLE_OUT}/vfdb/${SAMPLE}.vfdb.summary.tsv"
    fi

    echo "[INFO] Finished annotations for ${SAMPLE}"
done < "${MAG_MANIFEST}"
