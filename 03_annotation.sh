# -----------------------------
# 9. ORF prediction with Prodigal
# -----------------------------
mkdir -p "${OUTDIR}/prodigal"

for genome in "${OUTDIR}"/drep/dereplicated_genomes/*.fa; do
    base=$(basename "${genome}" .fa)
    prodigal \
      -i "${genome}" \
      -a "${OUTDIR}/prodigal/${base}.faa" \
      -d "${OUTDIR}/prodigal/${base}.fna" \
      -p meta
done

# -----------------------------
# 10. ARG annotation with RGI
# -----------------------------
mkdir -p "${OUTDIR}/rgi"

for faa in "${OUTDIR}"/prodigal/*.faa; do
    base=$(basename "${faa}" .faa)
    rgi main \
      --input_sequence "${faa}" \
      --output_file "${OUTDIR}/rgi/${base}" \
      --input_type protein \
      --alignment_tool DIAMOND \
      --clean
done

# -----------------------------
# 11. Virulence annotation with ABRicate + VFDB
# -----------------------------
mkdir -p "${OUTDIR}/abricate_vfdb"

for genome in "${OUTDIR}"/drep/dereplicated_genomes/*.fa; do
    base=$(basename "${genome}" .fa)
    abricate \
      --db vfdb \
      --minid 70 \
      --mincov 70 \
      "${genome}" \
      > "${OUTDIR}/abricate_vfdb/${base}.vfdb.tsv"
done
