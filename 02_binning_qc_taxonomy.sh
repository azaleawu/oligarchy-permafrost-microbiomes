# -----------------------------
# 5. MAG binning with VAMB
#    Note: paper used contigs >=50 kb
# -----------------------------
vamb \
  --outdir "${OUTDIR}/vamb" \
  --fasta "${OUTDIR}/megahit/assembly.contigs.fa" \
  --bamfiles "${OUTDIR}/alignment.sorted.bam" \
  -o C \
  --minfasta 50000 \
  -p "${THREADS}" \
  --model vae-aae

# -----------------------------
# 6. Dereplication with dRep
# -----------------------------
dRep dereplicate "${OUTDIR}/drep" \
  -g "${OUTDIR}/vamb/bins/"*.fna \
  --S_algorithm ANImf \
  -pa 0.9 \
  -sa 0.95 \
  -nc 0.1 \
  -comp 50 \
  -con 10 \
  --clusterAlg single \
  --ignoreGenomeQuality

# -----------------------------
# 7. Genome quality with CheckM2
# -----------------------------
checkm2 predict \
  --threads "${THREADS}" \
  --input "${OUTDIR}/drep/dereplicated_genomes" \
  --output-directory "${OUTDIR}/checkm2"

# -----------------------------
# 8. Taxonomy with GTDB-Tk
# -----------------------------
gtdbtk classify_wf \
  --genome_dir "${OUTDIR}/drep/dereplicated_genomes" \
  --out_dir "${OUTDIR}/gtdbtk" \
  --cpus "${THREADS}"




  
