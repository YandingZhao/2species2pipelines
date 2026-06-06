#!/usr/bin/env python3
"""
Cross-species scRNA-seq integration orchestrator.

Replaces Nextflow as the execution layer for programmatic / agentic use.
Accepts two input files (.rds or .h5ad), runs one or more integration methods,
evaluates with scIB metrics, generates UMAPs, and writes a JSON summary.

Full pipeline for .rds inputs:
  ortholog_convert → seurat_to_anndata (for Python methods)
                   → R integration → seurat_to_anndata (single) → h5ad
                   → Python integration → h5ad
  all h5ads → evaluate → aggregate → umap

Usage:
  python agentic/orchestrate.py \
    --input_a tests/data/data/Dog_pt15_Immune_Lymphoid_diet.rds \
    --input_b tests/data/data/Human_X00004_Immune_Lymphoid_diet.rds \
    --sample_id task4_demo --species_a dog --species_b human \
    --methods harmony scvi bbknn --outdir results/agentic
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

SCRIPTS_DIR = Path(__file__).parent.parent / "scripts"

RSCRIPT = os.environ.get("RSCRIPT", "Rscript")
PYTHON = sys.executable

R_METHODS = {"harmony", "seurat4", "fastmnn"}
PYTHON_METHODS = {"bbknn", "scanorama", "scvi", "scgen", "samap", "saturn"}
ALL_METHODS = R_METHODS | PYTHON_METHODS

# Feature selection methods available in run_feature_selection.py (Python)
PYTHON_FS_METHODS = (
    "seurat_v3", "seurat_v3_batch", "seurat", "seurat_batch",
    "cell_ranger", "cell_ranger_batch", "pearson", "pearson_batch",
    "mean", "variance", "wilcoxon", "triku", "hotspot", "anticor",
    "random", "all",
)
# Feature selection methods available in run_feature_selection_r.R (R only)
R_FS_METHODS = (
    "seurat_vst", "seurat_mvp", "seurat_disp", "seurat_sct",
    "osca", "brennecke", "nbumi", "dubstepr", "scry",
    "scpnmf", "singlecellhaystack", "scsegindex",
)
ALL_FS_METHODS = PYTHON_FS_METHODS + R_FS_METHODS

# Embedding key each method stores in obsm (used for UMAP + neighbour graph)
EMBEDDING_KEYS = {
    "harmony": "X_harmony",
    "seurat4": "X_pca",
    "fastmnn": "X_mnn",
    "bbknn": "X_pca",       # BBKNN also pre-computes the graph; handled specially
    "scanorama": "X_scanorama",
    "scvi": "X_scvi",
    "scgen": "X_scgen",
    "samap": "X_samap",
    "saturn": "X_saturn",
}

# Methods whose h5ad already has a pre-computed neighbour graph (skip sc.pp.neighbors)
PRECOMPUTED_GRAPH_METHODS = {"bbknn"}


def _run(cmd, cwd=None, label=""):
    """Run a subprocess; raise RuntimeError with captured output on failure."""
    result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"{label} failed (exit {result.returncode}):\n"
            f"STDOUT: {result.stdout[-3000:]}\n"
            f"STDERR: {result.stderr[-3000:]}"
        )
    return result


# ---------------------------------------------------------------------------
# Pipeline steps
# ---------------------------------------------------------------------------

def ortholog_convert(input_a, input_b, sample_id, species_a, species_b, workdir):
    """Convert species_a genes to species_b orthologs using orthogene.
    Returns (rds_a, rds_b) paths in workdir."""
    script = SCRIPTS_DIR / "run_ortholog_convert_pair.R"
    _run(
        [RSCRIPT, str(script),
         "--input_a", str(input_a), "--input_b", str(input_b),
         "--sample_id", sample_id,
         "--species_a", species_a, "--species_b", species_b],
        cwd=str(workdir), label="ortholog_convert",
    )
    return (
        Path(workdir) / f"{sample_id}_a_ortholog.rds",
        Path(workdir) / f"{sample_id}_b_ortholog.rds",
    )


def seurat_to_anndata_pair(input_a, input_b, workdir):
    """Convert a pair of Seurat RDS to h5ad for Python methods.
    Returns (h5ad_a, h5ad_b) paths."""
    script = SCRIPTS_DIR / "run_seurat_to_anndata.R"
    _run(
        [RSCRIPT, str(script),
         "--input_a", str(input_a), "--input_b", str(input_b)],
        cwd=str(workdir), label="seurat_to_anndata_pair",
    )
    return (
        Path(workdir) / f"{Path(input_a).stem}.h5ad",
        Path(workdir) / f"{Path(input_b).stem}.h5ad",
    )


def seurat_to_anndata_single(input_rds, workdir):
    """Convert a single integrated Seurat RDS to h5ad for evaluation.
    Returns h5ad path."""
    script = SCRIPTS_DIR / "run_seurat_to_anndata.R"
    _run(
        [RSCRIPT, str(script), "--input_a", str(input_rds)],
        cwd=str(workdir), label="seurat_to_anndata_single",
    )
    return Path(workdir) / f"{Path(input_rds).stem}.h5ad"


def run_feature_selection_py(input_a, input_b, sample_id, species_a, species_b,
                              fs_method, n_features, workdir):
    """Run Python feature selection. Returns path to features.txt."""
    _run(
        [PYTHON, str(SCRIPTS_DIR / "run_feature_selection.py"),
         "--input_a", str(input_a), "--input_b", str(input_b),
         "--sample_id", sample_id,
         "--species_a", species_a, "--species_b", species_b,
         "--method", fs_method, "--n_features", str(n_features)],
        cwd=str(workdir), label=f"feature_selection_{fs_method}",
    )
    return Path(workdir) / f"{sample_id}_{fs_method}_{n_features}_features.txt"


def run_feature_selection_r(input_a, input_b, sample_id, species_a, species_b,
                             fs_method, n_features, workdir):
    """Run R feature selection. Returns path to features.txt."""
    _run(
        [RSCRIPT, str(SCRIPTS_DIR / "run_feature_selection_r.R"),
         "--input_a", str(input_a), "--input_b", str(input_b),
         "--sample_id", sample_id,
         "--species_a", species_a, "--species_b", species_b,
         "--method", fs_method, "--n_features", str(n_features)],
        cwd=str(workdir), label=f"r_feature_selection_{fs_method}",
    )
    return Path(workdir) / f"{sample_id}_{fs_method}_{n_features}_features.txt"


def run_r_normalize(input_a, input_b, sample_id, norm_method, workdir):
    """Run standalone R normalization. Returns (rds_a_norm, rds_b_norm) paths."""
    out_a = Path(workdir) / f"{sample_id}_a_{norm_method}_norm.rds"
    out_b = Path(workdir) / f"{sample_id}_b_{norm_method}_norm.rds"
    _run(
        [RSCRIPT, str(SCRIPTS_DIR / "run_normalize_r.R"),
         "--input_a", str(input_a), "--input_b", str(input_b),
         "--method", norm_method,
         "--output_a", str(out_a), "--output_b", str(out_b)],
        cwd=str(workdir), label=f"normalize_r_{norm_method}",
    )
    return out_a, out_b


def run_r_integration(method, input_a, input_b, sample_id, species_a, species_b,
                      workdir, normalization="log_norm", features_file=None):
    """Run a Harmony / Seurat4 / fastMNN integration. Returns integrated RDS path."""
    script_map = {
        "harmony": "run_harmony_module.R",
        "seurat4": "run_seurat4_module.R",
        "fastmnn": "run_fastmnn_module.R",
    }
    cmd = [RSCRIPT, str(SCRIPTS_DIR / script_map[method]),
           "--input_a", str(input_a), "--input_b", str(input_b),
           "--sample_id", sample_id,
           "--species_a", species_a, "--species_b", species_b,
           "--normalization", normalization]
    if features_file:
        cmd += ["--features_file", str(features_file)]
    _run(cmd, cwd=str(workdir), label=f"{method}_integration")
    return Path(workdir) / f"{sample_id}_{method}_integration.rds"


def run_python_integration(method, input_a, input_b, sample_id, species_a, species_b,
                           workdir, normalization="log_norm", features_file=None):
    """Run a BBKNN / Scanorama / scVI / scGen integration. Returns integrated h5ad path."""
    script_map = {
        "bbknn": "run_bbknn_module.py",
        "scanorama": "run_scanorama_module.py",
        "scvi": "run_scvi_module.py",
        "scgen": "run_scgen_module.py",
        "samap": "run_samap_module.py",
        "saturn": "run_saturn_module.py",
    }
    cmd = [PYTHON, str(SCRIPTS_DIR / script_map[method]),
           "--input_a", str(input_a), "--input_b", str(input_b),
           "--sample_id", sample_id,
           "--species_a", species_a, "--species_b", species_b,
           "--normalization", normalization]
    if features_file:
        cmd += ["--features_file", str(features_file)]
    _run(cmd, cwd=str(workdir), label=f"{method}_integration")
    return Path(workdir) / f"{sample_id}_{method}_integration.h5ad"


def run_evaluation(h5ad_path, workdir):
    """Run scIB metrics on one integrated h5ad. Returns metrics TSV path."""
    script = SCRIPTS_DIR / "run_evaluate_integration.py"
    _run(
        [PYTHON, str(script), "--input_h5ad", str(h5ad_path)],
        cwd=str(workdir), label=f"evaluate_{h5ad_path.stem}",
    )
    return Path(workdir) / f"{h5ad_path.stem}_scib_metrics.tsv"


def run_aggregation(metrics_files, outdir):
    """Aggregate per-method metrics into combined report + figure.
    Returns (report_tsv, long_tsv, figure_png) paths."""
    script = SCRIPTS_DIR / "aggregate_metrics.py"
    report = outdir / "combined_metrics_report.tsv"
    long = outdir / "combined_metrics_long.tsv"
    figure = outdir / "combined_metrics_figure.png"
    _run(
        [PYTHON, str(script),
         "--metrics_files"] + [str(f) for f in metrics_files] +
        ["--output_report", str(report),
         "--output_long", str(long),
         "--output_figure", str(figure)],
        cwd=str(outdir), label="aggregate_metrics",
    )
    return report, long, figure


def generate_umap(h5ad_path, method, outdir):
    """Generate a 2-panel UMAP (batch | celltype) from integration output.
    Returns PNG path, or None if generation fails."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import scanpy as sc

        adata = sc.read_h5ad(h5ad_path)
        emb_key = EMBEDDING_KEYS.get(method)

        # Fall back to any valid 2-D embedding if expected key is absent
        if emb_key not in adata.obsm:
            candidates = [k for k in adata.obsm if adata.obsm[k].ndim == 2 and adata.obsm[k].shape[1] >= 2]
            if not candidates:
                return None
            emb_key = candidates[0]

        if method not in PRECOMPUTED_GRAPH_METHODS or "connectivities" not in adata.obsp:
            sc.pp.neighbors(adata, use_rep=emb_key, n_neighbors=15)

        sc.tl.umap(adata)

        fig, axes = plt.subplots(1, 2, figsize=(13, 5))
        for ax, color, title in zip(
            axes,
            ["batch", "celltype"],
            [f"{method.upper()} — batch", f"{method.upper()} — cell type"],
        ):
            if color not in adata.obs.columns:
                ax.set_visible(False)
                continue
            sc.pl.umap(adata, color=color, ax=ax, show=False, title=title,
                       frameon=False, legend_fontsize=7)

        fig.tight_layout()
        umap_path = outdir / f"{method}_umap.png"
        fig.savefig(umap_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return umap_path

    except Exception as exc:
        print(f"  WARNING: UMAP failed for {method}: {exc}")
        return None


# ---------------------------------------------------------------------------
# Summary helpers
# ---------------------------------------------------------------------------

def _parse_report_txt(path):
    """Parse a key: value report.txt into a dict."""
    report = {}
    try:
        with open(path) as fh:
            for line in fh:
                if ":" in line:
                    k, _, v = line.partition(":")
                    report[k.strip()] = v.strip()
    except OSError:
        pass
    return report


def _method_result(status, **kwargs):
    return {"status": status, **kwargs}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run(
    input_a,
    input_b,
    sample_id,
    species_a,
    species_b,
    methods=None,
    outdir="results/agentic",
    skip_ortholog=False,
    skip_umap=False,
    normalizations=None,
    feature_selections=None,
    n_features=2000,
):
    """Programmatic entry point (also called by CLI). Returns the summary dict."""
    if methods is None:
        methods = sorted(ALL_METHODS)
    if normalizations is None:
        normalizations = ["log_norm"]
    # feature_selections=None means no standalone FS step (use internal HVG in each script)
    if feature_selections is None:
        feature_selections = [None]

    input_a = Path(input_a).resolve()
    input_b = Path(input_b).resolve()
    outdir = Path(outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    is_rds = input_a.suffix.lower() == ".rds"
    r_methods = [m for m in methods if m in R_METHODS]
    py_methods = [m for m in methods if m in PYTHON_METHODS]

    summary = {
        "sample_id": sample_id,
        "species_a": species_a,
        "species_b": species_b,
        "input_format": "rds" if is_rds else "h5ad",
        "methods_requested": methods,
        "methods": {},
    }

    with tempfile.TemporaryDirectory(prefix="integration_") as _tmp:
        tmpdir = Path(_tmp)

        # ── Step 1: Ortholog conversion ──────────────────────────────────────
        if is_rds and not skip_ortholog:
            print(f"\n[1/4] Ortholog conversion  ({species_a} → {species_b})")
            rds_a, rds_b = ortholog_convert(
                input_a, input_b, sample_id, species_a, species_b, tmpdir
            )
        else:
            rds_a, rds_b = input_a, input_b

        # ── Step 2a: Convert RDS pair → h5ad for Python methods ─────────────
        h5ad_a = h5ad_b = None
        if is_rds and py_methods:
            print(f"\n[2/4] Seurat → AnnData (input pair for Python methods)")
            h5ad_a, h5ad_b = seurat_to_anndata_pair(rds_a, rds_b, tmpdir)
        elif not is_rds:
            h5ad_a, h5ad_b = input_a, input_b

        # ── Step 2b: Run integrations ────────────────────────────────────────
        integrated_h5ads = []
        multi_dim = len(normalizations) > 1 or any(f is not None for f in feature_selections)

        for fs in feature_selections:
            # ── R methods ───────────────────────────────────────────────────
            for method in r_methods:
                r_features_file = None
                if fs is not None:
                    fs_label = f"{fs}_{n_features}"
                    fs_sample_id = f"{sample_id}_{fs_label}" if multi_dim else sample_id
                    try:
                        if fs in R_FS_METHODS:
                            r_features_file = run_feature_selection_r(
                                rds_a, rds_b, fs_sample_id,
                                species_a, species_b, fs, n_features, tmpdir,
                            )
                        else:
                            # Python FS methods work on h5ad; use h5ad_a/b if available
                            if h5ad_a is not None:
                                r_features_file = run_feature_selection_py(
                                    h5ad_a, h5ad_b, fs_sample_id,
                                    species_a, species_b, fs, n_features, tmpdir,
                                )
                    except RuntimeError as exc:
                        print(f"  WARNING: feature selection {fs} failed for {method}: {exc}")

                run_id = (f"{sample_id}_{fs}_{n_features}" if fs else sample_id)
                label = f"{fs}/{method}" if fs else method
                print(f"\n[2/4] Running {label} (R)")
                try:
                    rds_out = run_r_integration(
                        method, rds_a, rds_b, run_id,
                        species_a, species_b, tmpdir,
                        features_file=r_features_file,
                    )
                    h5ad_out = seurat_to_anndata_single(rds_out, tmpdir)
                    final = outdir / h5ad_out.name
                    shutil.copy(h5ad_out, final)
                    integrated_h5ads.append(final)
                    report = _parse_report_txt(tmpdir / f"{run_id}_{method}_report.txt")
                    summary["methods"][label] = _method_result(
                        "ok", h5ad=str(final), run_report=report
                    )
                except RuntimeError as exc:
                    print(f"  ERROR: {exc}")
                    summary["methods"][label] = _method_result("failed", error=str(exc)[:800])

            # ── Python methods ───────────────────────────────────────────────
            for norm in normalizations:
                py_features_file = None
                if fs is not None and h5ad_a is not None and fs in PYTHON_FS_METHODS:
                    fs_label = f"{fs}_{n_features}"
                    fs_sample_id = (
                        f"{sample_id}_{norm}_{fs_label}" if multi_dim else sample_id
                    )
                    try:
                        py_features_file = run_feature_selection_py(
                            h5ad_a, h5ad_b, fs_sample_id,
                            species_a, species_b, fs, n_features, tmpdir,
                        )
                    except RuntimeError as exc:
                        print(f"  WARNING: Python feature selection {fs} failed: {exc}")

                # Build run-level sample_id and summary key
                parts = [sample_id]
                if len(normalizations) > 1:
                    parts.append(norm)
                if fs is not None:
                    parts.append(f"{fs}_{n_features}")
                norm_sample_id = "_".join(parts) if len(parts) > 1 else sample_id

                for method in py_methods:
                    label_parts = []
                    if len(normalizations) > 1:
                        label_parts.append(norm)
                    if fs is not None:
                        label_parts.append(fs)
                    label_parts.append(method)
                    run_label = "/".join(label_parts) if len(label_parts) > 1 else method

                    print(f"\n[2/4] Running {run_label} (Python)")
                    try:
                        h5ad_out = run_python_integration(
                            method, h5ad_a, h5ad_b, norm_sample_id,
                            species_a, species_b, tmpdir,
                            normalization=norm,
                            features_file=py_features_file,
                        )
                        final = outdir / h5ad_out.name
                        shutil.copy(h5ad_out, final)
                        integrated_h5ads.append(final)
                        report = _parse_report_txt(
                            tmpdir / f"{norm_sample_id}_{method}_report.txt"
                        )
                        summary["methods"][run_label] = _method_result(
                            "ok", h5ad=str(final), run_report=report
                        )
                    except RuntimeError as exc:
                        print(f"  ERROR: {exc}")
                        summary["methods"][run_label] = _method_result(
                            "failed", error=str(exc)[:800]
                        )

        # ── Step 3: Evaluate ─────────────────────────────────────────────────
        print(f"\n[3/4] scIB evaluation  ({len(integrated_h5ads)} outputs)")
        metrics_files = []
        for h5ad_path in integrated_h5ads:
            try:
                metrics_tsv = run_evaluation(h5ad_path, tmpdir)
                final_metrics = outdir / metrics_tsv.name
                shutil.copy(metrics_tsv, final_metrics)
                metrics_files.append(final_metrics)
                method_key = (
                    h5ad_path.stem
                    .removeprefix(f"{sample_id}_")
                    .removesuffix("_integration")
                )
                if method_key in summary["methods"]:
                    summary["methods"][method_key]["metrics_tsv"] = str(final_metrics)
            except RuntimeError as exc:
                print(f"  WARNING: evaluation failed for {h5ad_path.name}: {exc}")

        if metrics_files:
            print(f"\n[3/4] Aggregating {len(metrics_files)} metrics files")
            try:
                report_p, long_p, fig_p = run_aggregation(metrics_files, outdir)
                summary["aggregated"] = {
                    "report_tsv": str(report_p),
                    "long_tsv": str(long_p),
                    "figure_png": str(fig_p),
                }
            except RuntimeError as exc:
                print(f"  WARNING: aggregation failed: {exc}")

        # ── Step 4: UMAP ─────────────────────────────────────────────────────
        if not skip_umap:
            print(f"\n[4/4] Generating UMAPs")
            for h5ad_path in integrated_h5ads:
                method_key = (
                    h5ad_path.stem
                    .removeprefix(f"{sample_id}_")
                    .removesuffix("_integration")
                )
                umap_path = generate_umap(h5ad_path, method_key, outdir)
                if umap_path and method_key in summary["methods"]:
                    summary["methods"][method_key]["umap_png"] = str(umap_path)

    # ── Write JSON summary ───────────────────────────────────────────────────
    summary_path = outdir / "integration_summary.json"
    with open(summary_path, "w") as fh:
        json.dump(summary, fh, indent=2)

    n_ok = sum(1 for v in summary["methods"].values() if v["status"] == "ok")
    n_total = len(summary["methods"])
    print(f"\nDone — {n_ok}/{n_total} methods succeeded.")
    print(f"Results : {outdir}")
    print(f"Summary : {summary_path}")
    return summary


def main():
    parser = argparse.ArgumentParser(
        description="Cross-species scRNA-seq integration orchestrator",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input_a", required=True, help="Species A data (.rds or .h5ad)")
    parser.add_argument("--input_b", required=True, help="Species B data (.rds or .h5ad)")
    parser.add_argument("--sample_id", required=True, help="Sample identifier")
    parser.add_argument("--species_a", required=True, help="Species A name (e.g. dog)")
    parser.add_argument("--species_b", required=True, help="Species B name (e.g. human)")
    parser.add_argument(
        "--methods", nargs="+", default=sorted(ALL_METHODS),
        choices=sorted(ALL_METHODS),
        help="Integration methods to run",
    )
    parser.add_argument("--outdir", default="results/agentic", help="Output directory")
    parser.add_argument(
        "--skip_ortholog", action="store_true",
        help="Skip ortholog conversion (inputs already share gene space)",
    )
    parser.add_argument("--skip_umap", action="store_true", help="Skip UMAP generation")
    parser.add_argument(
        "--normalizations", nargs="+",
        default=["log_norm"],
        choices=["log_norm", "pearson_residuals", "scran", "sctransform", "raw_counts"],
        help="Normalization methods to sweep. Multiple values produce one run each.",
    )
    parser.add_argument(
        "--feature_selections", nargs="+",
        default=None,
        choices=list(ALL_FS_METHODS),
        help="Feature selection methods to sweep via standalone scripts. "
             "When omitted each integration script uses its own internal HVG selection. "
             "Multiple values produce one run per method.",
    )
    parser.add_argument(
        "--n_features", type=int, default=2000,
        help="Number of features to select (default: 2000).",
    )
    args = parser.parse_args()

    run(
        input_a=args.input_a,
        input_b=args.input_b,
        sample_id=args.sample_id,
        species_a=args.species_a,
        species_b=args.species_b,
        methods=args.methods,
        outdir=args.outdir,
        skip_ortholog=args.skip_ortholog,
        skip_umap=args.skip_umap,
        normalizations=args.normalizations,
        feature_selections=args.feature_selections,
        n_features=args.n_features,
    )


if __name__ == "__main__":
    main()
