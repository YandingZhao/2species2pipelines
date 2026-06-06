#!/usr/bin/env python3
"""Standalone Python normalization script.

Reads a single .h5ad file, applies normalization, writes output .h5ad.
Preserves raw counts in layers['counts'] so downstream integration scripts
can always recover the original matrix.

Usage
-----
python3 run_normalize_python.py \
    --input Dog.h5ad \
    --output Dog_lognorm.h5ad \
    --method log_norm
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import scanpy as sc
from normalization import NORM_METHODS, apply_normalization

PYTHON_NORM_METHODS = ("log_norm", "pearson_residuals", "raw_counts")


def parse_args():
    p = argparse.ArgumentParser(description="Standalone Python normalization")
    p.add_argument("--input", required=True, help="Input .h5ad file")
    p.add_argument("--output", required=True, help="Output .h5ad file")
    p.add_argument(
        "--method",
        required=True,
        choices=PYTHON_NORM_METHODS,
        help="Normalization method",
    )
    return p.parse_args()


def main():
    args = parse_args()

    adata = sc.read_h5ad(args.input)

    apply_normalization(adata, args.method)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    adata.write_h5ad(args.output)

    print(
        f"status: ok  method={args.method}  "
        f"cells={adata.n_obs}  genes={adata.n_vars}  "
        f"output={args.output}"
    )


if __name__ == "__main__":
    main()
