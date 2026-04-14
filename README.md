# oligarchy-permafrost-microbiomes
# Code for: An oligarchy of microbes and genes control the fate of ARGs and pathogens in permafrost soils worldwide

This repository contains analysis code for the genome-resolved study of antibiotic resistance and virulence-associated genes across global permafrost metagenomes.

The repository is currently being cleaned and documented. A public release of the full reproducible workflow will be made available upon acceptance of the manuscript.

Corresponding workflow components include:
1. Assemble reads with MEGAHIT
2. Map reads back to contigs with Bowtie2
3. Bin contigs into MAGs with VAMB
4. Dereplicate MAGs with dRep
5. Assess MAG quality with CheckM2
6. Assign taxonomy with GTDB-Tk
7. Predict ORFs with Prodigal
8. Annotate ARGs with RGI/CARD
9. Annotate virulence genes with ABRicate/VFDB
