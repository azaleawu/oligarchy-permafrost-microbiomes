#!/usr/bin/env python3
"""Merge MAG annotation, abundance and taxonomy tables."""

import argparse
from pathlib import Path
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mag-metadata", required=True)
    parser.add_argument("--abundance", required=True)
    parser.add_argument("--annotation", required=False)
    parser.add_argument("--output", required=True)
    parser.add_argument("--id-column", default="mag_id")
    args = parser.parse_args()

    meta = pd.read_csv(args.mag_metadata, sep="\t")
    abundance = pd.read_csv(args.abundance, sep="\t")
    if args.id_column not in meta.columns:
        raise ValueError(f"{args.id_column} not found in MAG metadata")
    if args.id_column not in abundance.columns:
        raise ValueError(f"{args.id_column} not found in abundance table")

    merged = meta.merge(abundance, on=args.id_column, how="left")
    if args.annotation:
        annot = pd.read_csv(args.annotation, sep="\t")
        if args.id_column not in annot.columns:
            raise ValueError(f"{args.id_column} not found in annotation table")
        merged = merged.merge(annot, on=args.id_column, how="left")

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
