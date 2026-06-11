#!/usr/bin/env python3
"""Convert an ABRicate --summary table into a numeric gene-count matrix.

ABRicate summary cells contain either '.' or semicolon-separated hit IDs.
This script converts '.' to 0 and all other cells to the number of entries.

Usage:
  python scripts/04_parse_abricate_summary.py \
    --input results/annotations/sample/vfdb/sample.vfdb.summary.tsv \
    --output results/annotations/sample/vfdb/sample.vfdb.counts.csv
"""

import argparse
import pandas as pd


def convert_count(cell):
    if pd.isna(cell):
        return 0
    text = str(cell).strip()
    if text == "." or text == "":
        return 0
    return len([x for x in text.split(";") if x.strip()])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="ABRicate summary TSV file")
    parser.add_argument("--output", required=True, help="Output CSV count matrix")
    parser.add_argument("--metadata-columns", type=int, default=2,
                        help="Number of leading metadata columns in ABRicate summary; default=2")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")
    if df.empty:
        raise ValueError(f"Input file is empty: {args.input}")

    for col in df.columns[args.metadata_columns:]:
        df[col] = df[col].apply(convert_count)

    df.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
