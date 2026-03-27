#!/usr/bin/env python3
"""
Convert a Parquet file to CSV.

Usage:
    python parquet_to_csv.py input.parquet [output.csv]

If no output file is specified, it will be created with the same name as the input.
"""

import sys
import pandas as pd
from pathlib import Path

def parquet_to_csv(parquet_path, csv_path=None):
    parquet_path = Path(parquet_path)
    if csv_path is None:
        csv_path = parquet_path.with_suffix(".csv")

    print(f"📂 Reading: {parquet_path}")
    df = pd.read_parquet(parquet_path)
    print(f"✅ Loaded {len(df):,} rows and {len(df.columns):,} columns")

    print(f"💾 Writing to: {csv_path}")
    df = df[["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]]
    df.to_csv(csv_path, index=False)
    print("✅ Conversion complete!")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python parquet_to_csv.py input.parquet [output.csv]")
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2] if len(sys.argv) > 2 else None
    parquet_to_csv(input_path, output_path)
