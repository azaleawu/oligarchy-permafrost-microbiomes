#!/usr/bin/env python3
"""Merge multiple ABRicate-derived count matrices by genome identifier."""

import argparse
from functools import reduce
from pathlib import Path
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", required=True, help="Input TSV count matrices")
    parser.add_argument("--output", required=True, help="Merged output TSV")
    parser.add_argument("--id-column", default=None, help="Identifier column; defaults to first column")
    args = parser.parse_args()

    tables = []
    id_col = args.id_column
    for path in args.inputs:
        df = pd.read_csv(path, sep="\t")
        if id_col is None:
            id_col = df.columns[0]
        if id_col not in df.columns:
            raise ValueError(f"Identifier column {id_col!r} not found in {path}")
        tables.append(df)

    merged = reduce(lambda left, right: pd.merge(left, right, on=id_col, how="outer"), tables)
    merged = merged.fillna(0)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
