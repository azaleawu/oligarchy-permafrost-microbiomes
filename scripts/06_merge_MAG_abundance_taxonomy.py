#!/usr/bin/env python3
"""Merge MAG gene-count, abundance and taxonomy tables.

Usage:
  python scripts/06_merge_MAG_abundance_taxonomy.py \
    --gene-counts results/combined_vfdb_counts.csv \
    --abundance results/combined_mags_tpm.tsv \
    --taxonomy results/merged_MAG_taxonomy.tsv \
    --output results/combined_MAG_table.csv
"""

import argparse
import os
import pandas as pd


def read_table(path):
    ext = os.path.splitext(path)[1].lower()
    if ext in [".xlsx", ".xls"]:
        return pd.read_excel(path)
    if ext in [".tsv", ".txt"]:
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene-counts", required=True)
    parser.add_argument("--abundance", required=True)
    parser.add_argument("--taxonomy", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--genome-column", default="genome")
    parser.add_argument("--abundance-column", default="tpm")
    args = parser.parse_args()

    genes = read_table(args.gene_counts)
    abundance = read_table(args.abundance)
    taxonomy = read_table(args.taxonomy)

    for name, df in [("gene-counts", genes), ("abundance", abundance), ("taxonomy", taxonomy)]:
        if args.genome_column not in df.columns:
            raise ValueError(f"Column '{args.genome_column}' not found in {name} table")

    abundance_keep = [args.genome_column, args.abundance_column]
    abundance = abundance[abundance_keep].rename(columns={args.abundance_column: "abundance_tpm"})

    keep_tax = [c for c in [args.genome_column, "domain", "phylum", "class", "order", "family", "genus", "species", "sample"] if c in taxonomy.columns]
    taxonomy = taxonomy[keep_tax]

    merged = genes.merge(abundance, on=args.genome_column, how="left")
    merged = merged.merge(taxonomy, on=args.genome_column, how="left")

    if "abundance_tpm" in merged.columns:
        merged["abundance_tpm"] = merged["abundance_tpm"].fillna(0)

    tax_cols = [c for c in ["domain", "phylum", "class", "order", "family", "genus", "species"] if c in merged.columns]
    for c in tax_cols:
        merged[c] = merged[c].fillna("unclassified")
    if "sample" in merged.columns:
        merged["sample"] = merged["sample"].fillna("unknown")

    merged.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
