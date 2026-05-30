"""Unit tests for the weighted aggregate score in aggregate_metrics.py."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))
from aggregate_metrics import _compute_aggregate_score, _aggregate_score_long_rows


def _make_report(bio_vals: list[float], batch_vals: list[float]) -> pd.DataFrame:
    """Build a minimal combined_report DataFrame with Bio and Batch rows."""
    n_bio   = len(bio_vals)
    n_batch = len(batch_vals)

    bio_index = pd.MultiIndex.from_arrays(
        [[f"bio_metric_{i}" for i in range(n_bio)], ["Bio conservation"] * n_bio],
        names=["Metric", "Metric Type"],
    )
    batch_index = pd.MultiIndex.from_arrays(
        [[f"batch_metric_{i}" for i in range(n_batch)], ["Batch correction"] * n_batch],
        names=["Metric", "Metric Type"],
    )
    col = pd.MultiIndex.from_tuples([("method_a", "X_emb")], names=["Integration", "Embedding"])

    bio_df   = pd.DataFrame(np.array(bio_vals).reshape(-1, 1),   index=bio_index,   columns=col)
    batch_df = pd.DataFrame(np.array(batch_vals).reshape(-1, 1), index=batch_index, columns=col)
    return pd.concat([bio_df, batch_df])


class TestComputeAggregateScore:

    def test_weighted_formula(self):
        """0.6 * mean(bio) + 0.4 * mean(batch) is appended as Overall score."""
        bio_vals   = [0.8, 0.6]          # mean = 0.7
        batch_vals = [0.5, 0.5]          # mean = 0.5
        expected   = 0.6 * 0.7 + 0.4 * 0.5   # = 0.62

        report = _make_report(bio_vals, batch_vals)
        result = _compute_aggregate_score(report)

        # Row added
        metric_types = result.index.get_level_values("Metric Type").astype(str)
        assert "Aggregate score" in metric_types.values

        # Value is correct
        agg_row = result.loc[metric_types == "Aggregate score"]
        agg_val = float(agg_row.iloc[0, 0])
        assert abs(agg_val - expected) < 1e-9

    def test_aggregate_row_label(self):
        """The appended row has Metric='Overall score' and Metric Type='Aggregate score'."""
        report = _make_report([0.5], [0.5])
        result = _compute_aggregate_score(report)

        metrics      = result.index.get_level_values("Metric").astype(str)
        metric_types = result.index.get_level_values("Metric Type").astype(str)

        assert "Total"   in metrics.values
        assert "Aggregate score" in metric_types.values

    def test_no_duplicate_aggregate_rows(self):
        """Calling the function twice does not add a second aggregate row."""
        report = _make_report([0.7, 0.8], [0.4, 0.6])
        result = _compute_aggregate_score(_compute_aggregate_score(report))

        metric_types = result.index.get_level_values("Metric Type").astype(str)
        agg_count = (metric_types == "Aggregate score").sum()
        assert agg_count == 1

    def test_missing_bio_rows_returns_unchanged(self):
        """If Bio conservation rows are absent the report is returned unchanged."""
        batch_index = pd.MultiIndex.from_arrays(
            [["batch_metric_0"], ["Batch correction"]],
            names=["Metric", "Metric Type"],
        )
        col = pd.MultiIndex.from_tuples([("method_a", "X_emb")], names=["Integration", "Embedding"])
        report = pd.DataFrame([[0.5]], index=batch_index, columns=col)

        result = _compute_aggregate_score(report)
        assert len(result) == len(report)

    def test_missing_batch_rows_returns_unchanged(self):
        """If Batch correction rows are absent the report is returned unchanged."""
        bio_index = pd.MultiIndex.from_arrays(
            [["bio_metric_0"], ["Bio conservation"]],
            names=["Metric", "Metric Type"],
        )
        col = pd.MultiIndex.from_tuples([("method_a", "X_emb")], names=["Integration", "Embedding"])
        report = pd.DataFrame([[0.8]], index=bio_index, columns=col)

        result = _compute_aggregate_score(report)
        assert len(result) == len(report)

    def test_multiple_methods(self):
        """Aggregate score is computed per column (per integration × embedding)."""
        col = pd.MultiIndex.from_tuples(
            [("method_a", "X_emb"), ("method_b", "X_emb")],
            names=["Integration", "Embedding"],
        )
        bio_index = pd.MultiIndex.from_arrays(
            [["bio_metric_0"], ["Bio conservation"]], names=["Metric", "Metric Type"]
        )
        batch_index = pd.MultiIndex.from_arrays(
            [["batch_metric_0"], ["Batch correction"]], names=["Metric", "Metric Type"]
        )
        bio_df   = pd.DataFrame([[0.8, 0.4]], index=bio_index,   columns=col)
        batch_df = pd.DataFrame([[0.6, 0.2]], index=batch_index, columns=col)
        report   = pd.concat([bio_df, batch_df])

        result = _compute_aggregate_score(report)
        metric_types = result.index.get_level_values("Metric Type").astype(str)
        agg_row = result.loc[metric_types == "Aggregate score"]

        # method_a: 0.6*0.8 + 0.4*0.6 = 0.72
        # method_b: 0.6*0.4 + 0.4*0.2 = 0.32
        assert abs(float(agg_row.iloc[0, 0]) - 0.72) < 1e-9
        assert abs(float(agg_row.iloc[0, 1]) - 0.32) < 1e-9

    def test_nan_handling(self):
        """NaN values are ignored when computing per-column means."""
        report = _make_report([0.8, float("nan")], [0.6, float("nan")])
        result = _compute_aggregate_score(report)

        metric_types = result.index.get_level_values("Metric Type").astype(str)
        agg_row = result.loc[metric_types == "Aggregate score"]
        agg_val = float(agg_row.iloc[0, 0])

        expected = 0.6 * 0.8 + 0.4 * 0.6   # NaN rows excluded from mean
        assert abs(agg_val - expected) < 1e-9


class TestAggregateLongRows:

    def test_long_rows_produced(self):
        """After _compute_aggregate_score, _aggregate_score_long_rows returns non-empty DF."""
        report = _make_report([0.7], [0.5])
        report_with_agg = _compute_aggregate_score(report)
        long = _aggregate_score_long_rows(report_with_agg)

        assert not long.empty
        assert set(long.columns) >= {"Integration", "Embedding", "Metric", "Metric Type", "value"}

    def test_long_rows_metric_type(self):
        """All rows in the long output have Metric Type = 'Aggregate score'."""
        report = _make_report([0.7, 0.8], [0.4, 0.5])
        report_with_agg = _compute_aggregate_score(report)
        long = _aggregate_score_long_rows(report_with_agg)

        assert (long["Metric Type"] == "Aggregate score").all()

    def test_no_long_rows_without_agg(self):
        """Without calling _compute_aggregate_score first, long rows are empty."""
        report = _make_report([0.7], [0.5])
        long = _aggregate_score_long_rows(report)
        assert long.empty
