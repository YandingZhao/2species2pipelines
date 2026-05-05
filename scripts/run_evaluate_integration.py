"""Evaluate integrated AnnData outputs with scib-metrics benchmarker."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scib_metrics.benchmark import BatchCorrection, Benchmarker, BioConservation


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate one integrated .h5ad result")
    parser.add_argument("--input_h5ad", required=True, help="Integrated .h5ad file to evaluate")
    parser.add_argument("--batch_key", default="batch", help="Batch column in adata.obs")
    parser.add_argument("--label_key", default="celltype", help="Cell type column in adata.obs")
    parser.add_argument(
        "--n_jobs",
        type=int,
        default=1,
        help="Number of worker processes used by scib-metrics",
    )
    return parser.parse_args()


def _sanitize_obs(adata, batch_key: str, label_key: str) -> None:
    if batch_key not in adata.obs.columns:
        adata.obs[batch_key] = "unknown"
    if label_key not in adata.obs.columns:
        adata.obs[label_key] = "unknown"

    adata.obs[batch_key] = adata.obs[batch_key].astype(str).astype("category")
    adata.obs[label_key] = adata.obs[label_key].astype(str).astype("category")


def _candidate_embeddings(adata) -> list[str]:
    keys = []
    for key in adata.obsm_keys():
        emb = np.asarray(adata.obsm[key])
        if emb.ndim != 2:
            continue
        if emb.shape[0] != adata.n_obs:
            continue
        if emb.shape[1] < 2:
            continue
        keys.append(key)

    return sorted(keys)


def _write_skip(metrics_path: Path, report_path: Path, reason: str) -> None:
    pd.DataFrame(
        {
            "metric": ["status", "reason"],
            "value": ["skipped", reason],
        }
    ).to_csv(metrics_path, sep="\t", index=False)

    with open(report_path, "w", encoding="utf-8") as handle:
        handle.write("status: skipped\n")
        handle.write(f"reason: {reason}\n")


def main() -> None:
    args = parse_args()

    input_path = Path(args.input_h5ad)
    stem = input_path.stem

    metrics_path = Path(f"{stem}_scib_metrics.tsv")
    scaled_metrics_path = Path(f"{stem}_scib_metrics_scaled.tsv")
    report_path = Path(f"{stem}_scib_report.txt")

    adata = sc.read_h5ad(input_path)
    _sanitize_obs(adata, args.batch_key, args.label_key)

    embeddings = _candidate_embeddings(adata)
    if not embeddings:
        _write_skip(metrics_path, report_path, "No valid 2D embeddings found in adata.obsm")
        pd.DataFrame().to_csv(scaled_metrics_path, sep="\t", index=False)
        return

    batch_count = adata.obs[args.batch_key].nunique()
    label_count = adata.obs[args.label_key].nunique()
    if batch_count < 2:
        _write_skip(metrics_path, report_path, "Need at least two batches for benchmarking")
        pd.DataFrame().to_csv(scaled_metrics_path, sep="\t", index=False)
        return

#    if label_count < 2:
#        _write_skip(metrics_path, report_path, "Need at least two labels for benchmarking")
#        pd.DataFrame().to_csv(scaled_metrics_path, sep="\t", index=False)
#        return

    benchmarker = Benchmarker(
        adata,
        batch_key=args.batch_key,
        label_key=args.label_key,
        bio_conservation_metrics=BioConservation(),
        batch_correction_metrics=BatchCorrection(),
        embedding_obsm_keys=embeddings,
        n_jobs=args.n_jobs,
    )

    benchmarker.benchmark()
    raw_results = benchmarker.get_results(min_max_scale=False)
    scaled_results = benchmarker.get_results(min_max_scale=True)

    raw_results.to_csv(metrics_path, sep="\t")
    scaled_results.to_csv(scaled_metrics_path, sep="\t")

    with open(report_path, "w", encoding="utf-8") as handle:
        handle.write("status: ok\n")
        handle.write(f"input_h5ad: {input_path.name}\n")
        handle.write(f"cells: {adata.n_obs}\n")
        handle.write(f"genes: {adata.n_vars}\n")
        handle.write(f"batch_key: {args.batch_key}\n")
        handle.write(f"label_key: {args.label_key}\n")
        handle.write(f"batch_levels: {batch_count}\n")
        handle.write(f"label_levels: {label_count}\n")
        handle.write(f"embeddings: {', '.join(embeddings)}\n")


if __name__ == "__main__":
    main()
