#!/usr/bin/env python3
"""Convert an ABRicate summary table into a gene-count matrix.

The script expects ABRicate summary output where each database column contains
semicolon-separated gene names or '.' for no hit.
"""

import argparse
from pathlib import Path
import pandas as pd


def count_hits(value: str) -> int:
    if pd.isna(value):
        return 0
    value = str(value).strip()
    if value in {"", ".", "NA", "nan"}:
        return 0
    return len([x for x in value.split(";") if x.strip()])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="ABRicate summary TSV")
    parser.add_argument("--output", required=True, help="Output count matrix TSV")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")
    if df.empty:
        raise ValueError("Input table is empty")

    id_col = df.columns[0]
    out = pd.DataFrame({id_col: df[id_col]})
    for col in df.columns[1:]:
        out[col] = df[col].map(count_hits)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
