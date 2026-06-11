#!/usr/bin/env python3
"""Merge multiple ABRicate-derived gene-count matrices.

This script is useful when several ABRicate databases were screened and a combined
presence/count matrix is needed. For the manuscript workflow, VFDB is used for
virulence-associated gene annotation; ARG annotation is primarily based on RGI/CARD.

Usage:
  python scripts/05_merge_abricate_matrices.py \
    --input-dir results/abricate_counts \
    --pattern 'all_*_MAG_summary_abundance.csv' \
    --output results/combined_abricate_counts.csv
"""

import argparse
import glob
import os
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--pattern", default="*.csv")
    parser.add_argument("--output", required=True)
    parser.add_argument("--id-column", default=None,
                        help="Identifier column name. If omitted, the first column is used.")
    args = parser.parse_args()

    files = sorted(glob.glob(os.path.join(args.input_dir, args.pattern)))
    if not files:
        raise FileNotFoundError(f"No files matched {args.pattern} in {args.input_dir}")

    frames = []
    for f in files:
        df = pd.read_csv(f)
        id_col = args.id_column or df.columns[0]
        df = df.rename(columns={id_col: "genome"})
        df = df.drop(columns=[c for c in ["NUM_FOUND", "database"] if c in df.columns], errors="ignore")
        frames.append(df)

    combined = pd.concat(frames, axis=0, ignore_index=True)
    numeric_cols = [c for c in combined.columns if c != "genome"]
    combined[numeric_cols] = combined[numeric_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    final = combined.groupby("genome", as_index=False)[numeric_cols].max()
    final[numeric_cols] = final[numeric_cols].astype(int)
    final.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
