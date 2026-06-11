#!/usr/bin/env bash
set -euo pipefail

# Assembly, read mapping, MAG binning, dereplication, quality assessment and taxonomy.
# Input: tab-separated sample manifest with columns: sample_id, read1, read2
# Usage: bash scripts/01_metagenome_to_mags.sh config/samples.tsv config/config.env

SAMPLE_MANIFEST=${1:?"Usage: bash 01_metagenome_to_mags.sh config/samples.tsv config/config.env"}
CONFIG=${2:?"Missing config.env"}
source "${CONFIG}"

mkdir -p "${WORK_DIR}" "${OUT_DIR}" "${LOG_DIR}"

while IFS=$'\t' read -r SAMPLE R1 R2; do
    [[ "${SAMPLE}" == "sample_id" ]] && continue
    [[ -z "${SAMPLE}" ]] && continue

    echo "[INFO] Processing ${SAMPLE}"
    SAMPLE_DIR="${WORK_DIR}/${SAMPLE}"
    mkdir -p "${SAMPLE_DIR}" "${SAMPLE_DIR}/logs"

    ASSEMBLY_DIR="${SAMPLE_DIR}/assembly"
    CONTIGS="${ASSEMBLY_DIR}/assembly.contigs.fa"
    CONTIG_INDEX="${ASSEMBLY_DIR}/assembly_db"
    BAM="${SAMPLE_DIR}/assembly.sorted.bam"

    # 1. Single-sample metagenome assembly
    megahit \
      -1 "${R1}" \
      -2 "${R2}" \
      -o "${ASSEMBLY_DIR}" \
      --out-prefix assembly \
      --min-contig-len "${MIN_CONTIG_LEN}" \
      --presets meta-sensitive \
      -t "${THREADS}" \
      > "${SAMPLE_DIR}/logs/megahit.log" 2>&1

    # 2. Build Bowtie2 index for contigs
    bowtie2-build "${CONTIGS}" "${CONTIG_INDEX}" \
      > "${SAMPLE_DIR}/logs/bowtie2_build.log" 2>&1

    # 3. Map reads back to contigs
    bowtie2 \
      -x "${CONTIG_INDEX}" \
      -1 "${R1}" \
      -2 "${R2}" \
      --sensitive \
      -p "${THREADS}" \
      -S "${SAMPLE_DIR}/assembly.sam" \
      > "${SAMPLE_DIR}/logs/bowtie2_mapping.log" 2>&1

    # 4. Convert, sort and index BAM
    samtools view -@ "${THREADS}" -F 4 -bS "${SAMPLE_DIR}/assembly.sam" | \
      samtools sort -@ "${THREADS}" -o "${BAM}" -
    samtools index "${BAM}"
    rm -f "${SAMPLE_DIR}/assembly.sam"

    # 5. MAG binning with VAMB
    vamb \
      --outdir "${SAMPLE_DIR}/vamb" \
      --fasta "${CONTIGS}" \
      --bamfiles "${BAM}" \
      -o C \
      --minfasta "${VAMB_MINFASTA}" \
      -p "${THREADS}" \
      --model vae-aae \
      > "${SAMPLE_DIR}/logs/vamb.log" 2>&1

    # 6. Dereplicate bins within sample. Quality filtering is applied after CheckM2.
    mkdir -p "${SAMPLE_DIR}/drep"
    dRep dereplicate "${SAMPLE_DIR}/drep" \
      -g "${SAMPLE_DIR}"/vamb/bins/k*/*.fna \
      --S_algorithm ANImf \
      -pa "${DREP_PRIMARY_ANI}" \
      -sa "${DREP_SECONDARY_ANI}" \
      -nc "${DREP_COVERAGE}" \
      --clusterAlg single \
      --ignoreGenomeQuality \
      > "${SAMPLE_DIR}/logs/drep.log" 2>&1

    # 7. MAG quality assessment with CheckM2
    checkm2 predict \
      --threads "${THREADS}" \
      --input "${SAMPLE_DIR}/drep/dereplicated_genomes" \
      --output-directory "${SAMPLE_DIR}/checkm2" \
      --database_path "${CHECKM2_DB}" \
      > "${SAMPLE_DIR}/logs/checkm2.log" 2>&1

    # 8. Taxonomic assignment with GTDB-Tk
    export GTDBTK_DATA_PATH="${GTDBTK_DATA_PATH}"
    gtdbtk classify_wf \
      --genome_dir "${SAMPLE_DIR}/drep/dereplicated_genomes" \
      --out_dir "${SAMPLE_DIR}/gtdbtk" \
      --cpus "${THREADS}" \
      --skip_ani_screen \
      > "${SAMPLE_DIR}/logs/gtdbtk.log" 2>&1

    echo "[INFO] Finished ${SAMPLE}"
done < "${SAMPLE_MANIFEST}"
