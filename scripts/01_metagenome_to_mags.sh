#!/usr/bin/env bash
set -euo pipefail

# Generic upstream workflow for single-sample metagenome assembly and MAG reconstruction.
# Edit config/config_template.env and config/sample_table_template.tsv before use.

CONFIG=${1:-config/config_template.env}
source "$CONFIG"

mkdir -p "$OUT_DIR" "$TMP_DIR"

while IFS=$'\t' read -r sample_id read1 read2; do
  [[ "$sample_id" == "sample_id" ]] && continue
  echo "[INFO] Processing sample: $sample_id"

  SAMPLE_OUT="$OUT_DIR/$sample_id"
  mkdir -p "$SAMPLE_OUT"/{assembly,mapping,vamb,drep,checkm2,gtdbtk}

  # Assembly
  megahit \
    -1 "$read1" \
    -2 "$read2" \
    -o "$SAMPLE_OUT/assembly" \
    --out-prefix "$sample_id" \
    --min-contig-len "$MIN_CONTIG_LEN" \
    --presets "$MEGAHIT_PRESET" \
    -t "$THREADS"

  CONTIGS="$SAMPLE_OUT/assembly/${sample_id}.contigs.fa"

  # Read mapping back to assembled contigs
  bowtie2-build "$CONTIGS" "$SAMPLE_OUT/mapping/${sample_id}_assembly_index"
  bowtie2 \
    -x "$SAMPLE_OUT/mapping/${sample_id}_assembly_index" \
    -1 "$read1" \
    -2 "$read2" \
    -p "$THREADS" \
    -S "$SAMPLE_OUT/mapping/${sample_id}.sam"

  samtools view -@ "$THREADS" -bS "$SAMPLE_OUT/mapping/${sample_id}.sam" | \
    samtools sort -@ "$THREADS" -o "$SAMPLE_OUT/mapping/${sample_id}.sorted.bam"
  samtools index "$SAMPLE_OUT/mapping/${sample_id}.sorted.bam"
  rm -f "$SAMPLE_OUT/mapping/${sample_id}.sam"

  # MAG binning with VAMB
  vamb \
    --outdir "$SAMPLE_OUT/vamb" \
    --fasta "$CONTIGS" \
    --bamfiles "$SAMPLE_OUT/mapping/${sample_id}.sorted.bam" \
    --minfasta 200000 \
    --cuda false

  # Dereplication of bins for this sample or batch, if bins are available
  if compgen -G "$SAMPLE_OUT/vamb/bins/*.fna" > /dev/null; then
    dRep dereplicate "$SAMPLE_OUT/drep" \
      -g "$SAMPLE_OUT/vamb/bins/*.fna" \
      -comp 50 \
      -con 10 \
      -p "$THREADS"
  fi

  # MAG quality assessment
  if [ -d "$SAMPLE_OUT/drep/dereplicated_genomes" ]; then
    checkm2 predict \
      --threads "$THREADS" \
      --input "$SAMPLE_OUT/drep/dereplicated_genomes" \
      --output-directory "$SAMPLE_OUT/checkm2"

    # Taxonomic assignment
    gtdbtk classify_wf \
      --genome_dir "$SAMPLE_OUT/drep/dereplicated_genomes" \
      --out_dir "$SAMPLE_OUT/gtdbtk" \
      --extension fna \
      --cpus "$THREADS"
  fi

done < "$SAMPLE_TABLE"
