#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# User-defined inputs
# -----------------------------
SAMPLE="sample1"
R1="data/${SAMPLE}_1.fastq.gz"
R2="data/${SAMPLE}_2.fastq.gz"
OUTDIR="results/${SAMPLE}"

THREADS=32

mkdir -p "${OUTDIR}"

# -----------------------------
# 1. Assembly with MEGAHIT
# -----------------------------
megahit \
  -1 "${R1}" \
  -2 "${R2}" \
  -o "${OUTDIR}/megahit" \
  --out-prefix assembly \
  --min-contig-len 500 \
  --presets meta-sensitive \
  -t "${THREADS}"

# -----------------------------
# 2. Build Bowtie2 index
# -----------------------------
bowtie2-build \
  "${OUTDIR}/megahit/assembly.contigs.fa" \
  "${OUTDIR}/megahit/assembly_db"

# -----------------------------
# 3. Map reads back to contigs
# -----------------------------
bowtie2 \
  -x "${OUTDIR}/megahit/assembly_db" \
  -1 "${R1}" \
  -2 "${R2}" \
  --sensitive \
  -S "${OUTDIR}/alignment.sam" \
  -p "${THREADS}"

# -----------------------------
# 4. Convert / sort / index BAM
# -----------------------------
samtools view -F 4 -bS "${OUTDIR}/alignment.sam" | \
  samtools sort -o "${OUTDIR}/alignment.sorted.bam"

samtools index "${OUTDIR}/alignment.sorted.bam"
