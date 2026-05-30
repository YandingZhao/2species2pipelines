"""Aggregate per-integration unscaled scIB metrics into combined report tables."""

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Combine individual *_scib_metrics.tsv files into report tables"
    )
    parser.add_argument(
        "--metrics_files",
        nargs="+",
        required=True,
        help="Input unscaled metrics TSV files",
    )
    parser.add_argument(
        "--output_report",
        default="combined_unscaled_metrics_report.tsv",
        help="Output wide report table path",
    )
    parser.add_argument(
        "--output_long",
        default="combined_unscaled_metrics_long.tsv",
        help="Output long-format metrics table path",
    )
    parser.add_argument(
        "--output_figure",
        default="combined_unscaled_metrics_report.png",
        help="Output Nature-style heatmap figure path",
    )
    return parser.parse_args()


def integration_name_from_path(path: Path) -> str:
    suffix = "_scib_metrics.tsv"
    if path.name.endswith(suffix):
        return path.name[: -len(suffix)]
    return path.stem


def _display_integration_name(name: str) -> str:
    cleaned = re.sub(r"_integration$", "", name)
    cleaned = re.sub(r"^task\d+_", "", cleaned)
    cleaned = re.sub(r"^demo_", "", cleaned)

    token_map = {
        "bbknn": "BBKNN",
        "fastmnn": "fastMNN",
        "harmony": "Harmony",
        "scanorama": "Scanorama",
        "scgen": "scGen",
        "scvi": "scVI",
        "seurat4": "Seurat4",
    }
    return token_map.get(cleaned, cleaned.replace("_", " "))


def _coerce_numeric(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for col in out.columns:
        try:
            out[col] = pd.to_numeric(out[col])
        except (ValueError, TypeError):
            continue
    return out


def _read_status_metrics(raw: pd.DataFrame, integration: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    kv = dict(zip(raw["metric"].astype(str), raw["value"].astype(str)))

    report = pd.DataFrame(
        {"NA": [kv.get("status", "unknown"), kv.get("reason", "")]},
        index=pd.MultiIndex.from_tuples(
            [("status", "Run status"), ("reason", "Run status")],
            names=["Metric", "Metric Type"],
        ),
    )
    report.columns = pd.MultiIndex.from_tuples(
        [(integration, "NA")], names=["Integration", "Embedding"]
    )

    long_rows = raw[["metric", "value"]].copy()
    long_rows.columns = ["Metric", "value"]
    long_rows.insert(1, "Metric Type", "Run status")
    long_rows.insert(2, "Embedding", "NA")
    long_rows.insert(0, "Integration", integration)
    return report, long_rows


def _read_matrix_metrics(raw: pd.DataFrame, integration: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    metric_col = raw.columns[0]
    embedding_cols = [col for col in raw.columns if col not in {metric_col, "Metric Type"}]

    matrix = raw.rename(columns={metric_col: "Metric"}).set_index(["Metric", "Metric Type"])
    matrix = matrix[embedding_cols]
    matrix = _coerce_numeric(matrix)
    matrix.index.names = ["Metric", "Metric Type"]

    long_rows = matrix.reset_index().melt(
        id_vars=["Metric", "Metric Type"],
        var_name="Embedding",
        value_name="value",
    )
    long_rows.insert(0, "Integration", integration)
    long_rows = long_rows[["Integration", "Embedding", "Metric", "Metric Type", "value"]]

    matrix.columns = pd.MultiIndex.from_tuples(
        [(integration, embedding) for embedding in matrix.columns],
        names=["Integration", "Embedding"],
    )

    return matrix, long_rows


def read_one_metrics(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    integration = integration_name_from_path(path)

    raw = pd.read_csv(path, sep="\t")
    if set(raw.columns) >= {"metric", "value"}:
        return _read_status_metrics(raw, integration)

    return _read_matrix_metrics(raw, integration)


def _format_annotation(value: float) -> str:
    if pd.isna(value):
        return ""
    if abs(value) >= 100:
        return f"{value:.0f}"
    if abs(value) >= 10:
        return f"{value:.1f}"
    return f"{value:.2f}"


def _palette(metric_types: list[str]) -> dict[str, str]:
    base = {
        "Bio conservation": "#4C78A8",
        "Batch correction": "#D17C3F",
        "Aggregate score": "#5A9E6F",
        "Run status": "#8A8F98",
    }
    extra = ["#7E6AA2", "#9B6A4D", "#4F8C8A", "#A5647A"]
    missing = [name for name in metric_types if name not in base]
    for idx, name in enumerate(missing):
        base[name] = extra[idx % len(extra)]
    return {name: base[name] for name in metric_types}


def render_report_figure(combined_report: pd.DataFrame, output_figure: str) -> None:
    import matplotlib.pyplot as plt
    from matplotlib import colors
    from matplotlib.patches import Rectangle
    numeric_report = combined_report.apply(pd.to_numeric, errors="coerce")

    row_metric_types = numeric_report.index.get_level_values("Metric Type").astype(str)
    valid_rows = ~(
        numeric_report.isna().all(axis=1)
        & row_metric_types.isin(["Run status"])
    )
    numeric_report = numeric_report.loc[valid_rows]

    if numeric_report.empty:
        fig, ax = plt.subplots(figsize=(8, 3), dpi=300)
        ax.axis("off")
        ax.text(
            0.5,
            0.5,
            "No numeric metrics available for visualization",
            ha="center",
            va="center",
            fontsize=12,
            color="#374151",
        )
        fig.savefig(output_figure, dpi=300, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        return

    metric_types = list(dict.fromkeys(numeric_report.index.get_level_values("Metric Type").astype(str)))
    type_colors = _palette(metric_types)
    values = numeric_report.to_numpy(dtype=float)
    n_rows, n_cols = values.shape

    width = max(9.5, 2.8 + 0.75 * n_cols)
    height = max(5.5, 2.4 + 0.42 * n_rows)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 14,
            "axes.labelsize": 10,
        }
    )

    fig = plt.figure(figsize=(width, height), dpi=300, constrained_layout=True)
    gs = fig.add_gridspec(
        nrows=2,
        ncols=4,
        height_ratios=[0.12, 0.88],
        width_ratios=[0.12, 0.06, 0.76, 0.06],
        wspace=0.03,
        hspace=0.02,
    )

    ax_type_label = fig.add_subplot(gs[1, 0])
    ax_type = fig.add_subplot(gs[1, 1])
    ax = fig.add_subplot(gs[1, 2])
    ax_top = fig.add_subplot(gs[0, 2], sharex=ax)
    cax = fig.add_subplot(gs[1, 3])

    cmap = colors.LinearSegmentedColormap.from_list(
        "nature_scores",
        ["#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"],
    ).copy()
    cmap.set_bad("#eceff3")
    norm = colors.Normalize(vmin=0, vmax=1)

    im = ax.imshow(values, aspect="auto", cmap=cmap, norm=norm)

    ax.set_xticks(np.arange(n_cols))
    embedding_labels = [str(embedding) for _, embedding in numeric_report.columns]
    ax.set_xticklabels(embedding_labels, rotation=65, ha="right", rotation_mode="anchor")
    row_labels = [str(metric) for metric, _ in numeric_report.index]
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(row_labels)
    ax.tick_params(length=0, pad=4)

    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.7)
    ax.tick_params(which="minor", bottom=False, left=False)
    for spine in ax.spines.values():
        spine.set_visible(False)

    aggregate_rows = numeric_report.index.get_level_values("Metric Type") == "Aggregate score"
    for row_idx in range(n_rows):
        for col_idx in range(n_cols):
            value = values[row_idx, col_idx]
            if np.isnan(value):
                continue
            if n_rows <= 16 or aggregate_rows[row_idx]:
                text_color = "white" if value >= 0.55 else "#102a43"
                ax.text(
                    col_idx,
                    row_idx,
                    _format_annotation(value),
                    ha="center",
                    va="center",
                    fontsize=7.2,
                    color=text_color,
                )

    integration_labels = [str(integration) for integration, _ in numeric_report.columns]
    integration_spans = []
    start = 0
    current = integration_labels[0]
    for idx, label in enumerate(integration_labels[1:], start=1):
        if label != current:
            integration_spans.append((current, start, idx - 1))
            current = label
            start = idx
    integration_spans.append((current, start, n_cols - 1))

    ax_top.set_xlim(-0.5, n_cols - 0.5)
    ax_top.set_ylim(0, 1)
    ax_top.axis("off")
    for label, start_idx, end_idx in integration_spans:
        center = (start_idx + end_idx) / 2
        ax_top.text(
            center,
            0.55,
            _display_integration_name(label),
            ha="center",
            va="center",
            fontsize=10,
            fontweight="semibold",
            color="#1f2937",
        )
        ax_top.plot([start_idx - 0.45, end_idx + 0.45], [0.18, 0.18], color="#64748b", linewidth=1.0)

    ax_type.set_xlim(0, 1)
    ax_type.set_ylim(n_rows - 0.5, -0.5)
    ax_type.axis("off")
    metric_type_values = numeric_report.index.get_level_values("Metric Type").astype(str)
    for row_idx, metric_type in enumerate(metric_type_values):
        ax_type.add_patch(
            Rectangle((0, row_idx - 0.5), 1, 1, facecolor=type_colors[metric_type], edgecolor="white", linewidth=0.7)
        )

    ax_type_label.set_xlim(0, 1)
    ax_type_label.set_ylim(n_rows - 0.5, -0.5)
    ax_type_label.axis("off")
    type_spans = []
    start = 0
    current = metric_type_values[0]
    for idx, metric_type in enumerate(metric_type_values[1:], start=1):
        if metric_type != current:
            type_spans.append((current, start, idx - 1))
            current = metric_type
            start = idx
    type_spans.append((current, start, n_rows - 1))

    for metric_type, start_idx, end_idx in type_spans:
        center = (start_idx + end_idx) / 2
        ax_type_label.text(
            0.96,
            center,
            metric_type,
            rotation=90,
            ha="right",
            va="center",
            fontsize=9,
            fontweight="semibold",
            color=type_colors[metric_type],
        )

    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("Score", rotation=90)
    cbar.outline.set_visible(False)
    fig.savefig(output_figure, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _compute_aggregate_score(combined_report: pd.DataFrame) -> pd.DataFrame:
    """Append one Overall score row using the Zhong et al. 2025 weighting: 60% bio + 40% batch."""
    numeric = combined_report.apply(pd.to_numeric, errors="coerce")
    metric_types = numeric.index.get_level_values("Metric Type").astype(str)

    # Idempotent — skip if already computed
    if "Aggregate score" in metric_types.values:
        return combined_report

    bio_rows   = numeric.loc[metric_types == "Bio conservation"]
    batch_rows = numeric.loc[metric_types == "Batch correction"]

    if bio_rows.empty or batch_rows.empty:
        return combined_report

    agg_score = 0.6 * bio_rows.mean(axis=0) + 0.4 * batch_rows.mean(axis=0)

    # Use "Total" to match scib_metrics naming convention
    agg_index = pd.MultiIndex.from_tuples(
        [("Total", "Aggregate score")],
        names=["Metric", "Metric Type"],
    )
    agg_row = pd.DataFrame([agg_score.values], index=agg_index, columns=combined_report.columns)
    return pd.concat([combined_report, agg_row])


def _aggregate_score_long_rows(combined_report: pd.DataFrame) -> pd.DataFrame:
    """Return the aggregate score row in long format to append to combined_long."""
    numeric = combined_report.apply(pd.to_numeric, errors="coerce")
    metric_types = numeric.index.get_level_values("Metric Type").astype(str)
    agg_rows = numeric.loc[metric_types == "Aggregate score"]
    if agg_rows.empty:
        return pd.DataFrame()

    # Iterate over MultiIndex columns explicitly — melt doesn't handle them reliably
    rows = []
    for (metric, metric_type), row in agg_rows.iterrows():
        for (integration, embedding), value in row.items():
            rows.append({
                "Integration": integration,
                "Embedding": embedding,
                "Metric": metric,
                "Metric Type": metric_type,
                "value": value,
            })
    return pd.DataFrame(rows, columns=["Integration", "Embedding", "Metric", "Metric Type", "value"])


def main() -> None:
    args = parse_args()

    report_tables = []
    long_tables = []

    for file_str in args.metrics_files:
        path = Path(file_str)
        report_part, long_part = read_one_metrics(path)
        report_tables.append(report_part)
        long_tables.append(long_part)

    combined_report = pd.concat(report_tables, axis=1, sort=False)
    combined_long = pd.concat(long_tables, ignore_index=True, sort=False)

    combined_report = _compute_aggregate_score(combined_report)
    agg_long = _aggregate_score_long_rows(combined_report)
    if not agg_long.empty:
        combined_long = pd.concat([combined_long, agg_long], ignore_index=True, sort=False)

    combined_report.to_csv(args.output_report, sep="\t")
    combined_long.to_csv(args.output_long, sep="\t", index=False)
    render_report_figure(combined_report, args.output_figure)


if __name__ == "__main__":
    main()
