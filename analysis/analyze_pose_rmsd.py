#!/usr/bin/env python3
"""Analyze and plot ligand pose RMSD results."""

import argparse
import re
from pathlib import Path
from typing import List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


RANK_RE = re.compile(r"(?:rank|idx_)(\d+)", re.IGNORECASE)


def extract_rank(pose_name: str) -> int:
    """Extract rank index from pose_name, fallback to large value."""
    if not isinstance(pose_name, str):
        return 9999
    match = RANK_RE.search(pose_name)
    if not match:
        return 9999
    return int(match.group(1))


def load_rmsd_table(csv_path: Path, include_low_coverage: bool, diffdock_rank1_only: bool) -> pd.DataFrame:
    """Load RMSD CSV and keep only usable rows."""
    df = pd.read_csv(csv_path)
    if "rmsd_angstrom" not in df.columns:
        raise ValueError("Input CSV must contain 'rmsd_angstrom'.")

    allowed_status = ["ok", "ok_low_coverage"] if include_low_coverage else ["ok"]
    df = df[df["status"].isin(allowed_status)].copy()

    if diffdock_rank1_only:
        keep = (df["source"] != "diffdock_sdf") | (df["pose_name"].apply(extract_rank) == 1)
        df = df[keep].copy()

    df["rmsd_angstrom"] = pd.to_numeric(df["rmsd_angstrom"], errors="coerce")
    df = df.dropna(subset=["rmsd_angstrom"])
    return df


def exclude_sources(df: pd.DataFrame, sources_to_exclude: List[str]) -> pd.DataFrame:
    """Exclude rows by source labels."""
    if not sources_to_exclude:
        return df
    excluded = {s.strip() for s in sources_to_exclude if s.strip()}
    if not excluded:
        return df
    return df[~df["source"].astype(str).isin(excluded)].copy()


def enforce_common_pdb_across_groups(df: pd.DataFrame) -> pd.DataFrame:
    """Keep only pdb_id values present in every source|scenario group."""
    if df.empty:
        return df

    grouped_pdb_sets = (
        df.groupby(["source", "scenario"], dropna=False)["pdb_id"]
        .apply(lambda s: set(s.astype(str)))
        .tolist()
    )
    if not grouped_pdb_sets:
        return df

    common_pdb = set.intersection(*grouped_pdb_sets)
    if not common_pdb:
        return df.iloc[0:0].copy()

    return df[df["pdb_id"].astype(str).isin(common_pdb)].copy()


def summarize_groups(df: pd.DataFrame) -> pd.DataFrame:
    """Create summary metrics grouped by source/scenario."""
    summary = (
        df.groupby(["source", "scenario"], dropna=False)["rmsd_angstrom"]
        .agg(
            n="count",
            mean="mean",
            median="median",
            std="std",
            p25=lambda x: np.percentile(x, 25),
            p75=lambda x: np.percentile(x, 75),
            p90=lambda x: np.percentile(x, 90),
            min="min",
            max="max",
        )
        .reset_index()
        .sort_values(["source", "scenario"])
    )
    return summary


def summarize_best_per_pdb(df: pd.DataFrame) -> pd.DataFrame:
    """For each pdb_id and group, keep the best (lowest) RMSD."""
    best = (
        df.groupby(["pdb_id", "source", "scenario"], dropna=False)["rmsd_angstrom"]
        .min()
        .reset_index()
    )
    return best


def plot_box_by_group(df: pd.DataFrame, out_path: Path) -> None:
    """Boxplot of RMSD by source+scenario."""
    grouped = df.copy()
    grouped["group"] = grouped["source"].astype(str) + " | " + grouped["scenario"].astype(str)
    groups = sorted(grouped["group"].unique())
    values: List[np.ndarray] = [grouped[grouped["group"] == g]["rmsd_angstrom"].to_numpy() for g in groups]

    fig, ax = plt.subplots(figsize=(max(8, len(groups) * 1.4), 5))
    ax.boxplot(values, tick_labels=groups, showfliers=False)
    ax.set_title("Pose RMSD Distribution by Group")
    ax.set_ylabel("RMSD (Angstrom)")
    ax.grid(True, axis="y", alpha=0.3)
    plt.xticks(rotation=35, ha="right")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_best_per_pdb(best_df: pd.DataFrame, out_path: Path) -> None:
    """Plot best-per-PDB RMSD as violin+box per group."""
    best = best_df.copy()
    best["group"] = best["source"].astype(str) + " | " + best["scenario"].astype(str)
    groups = sorted(best["group"].unique())
    values: List[np.ndarray] = [best[best["group"] == g]["rmsd_angstrom"].to_numpy() for g in groups]

    fig, ax = plt.subplots(figsize=(max(8, len(groups) * 1.4), 5))
    parts = ax.violinplot(values, showmeans=False, showmedians=True, showextrema=False)
    for body in parts["bodies"]:
        body.set_alpha(0.35)
    ax.boxplot(values, tick_labels=groups, widths=0.2, showfliers=False)

    ax.set_title("Best Pose RMSD per PDB (Lower is Better)")
    ax.set_ylabel("Best RMSD (Angstrom)")
    ax.grid(True, axis="y", alpha=0.3)
    plt.xticks(rotation=35, ha="right")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_cdf(best_df: pd.DataFrame, out_path: Path) -> None:
    """Empirical CDF of best RMSD per group."""
    best = best_df.copy()
    best["group"] = best["source"].astype(str) + " | " + best["scenario"].astype(str)
    groups = sorted(best["group"].unique())

    fig, ax = plt.subplots(figsize=(8, 5))
    for group in groups:
        vals = np.sort(best[best["group"] == group]["rmsd_angstrom"].to_numpy())
        if len(vals) == 0:
            continue
        y = np.arange(1, len(vals) + 1) / len(vals)
        ax.plot(vals, y, label=f"{group} (n={len(vals)})")

    ax.set_title("CDF of Best Pose RMSD per PDB")
    ax.set_xlabel("Best RMSD (Angstrom)")
    ax.set_ylabel("Fraction <= RMSD")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Analyze pose RMSD results and generate plots")
    parser.add_argument(
        "--input",
        default="./results/analysis/pose_rmsd_summary.csv",
        help="Path to pose RMSD summary CSV",
    )
    parser.add_argument(
        "--output-dir",
        default="./results/analysis",
        help="Output directory for plots and summaries",
    )
    parser.add_argument(
        "--include-low-coverage",
        action="store_true",
        help="Include rows with status=ok_low_coverage",
    )
    parser.add_argument(
        "--diffdock-rank1-only",
        action="store_true",
        help="For source=diffdock_sdf, keep only rank1 poses",
    )
    parser.add_argument(
        "--exclude-sources",
        default="",
        help="Comma-separated source names to exclude (e.g. diffdock_sdf)",
    )
    parser.add_argument(
        "--common-pdb-only",
        action="store_true",
        help="Use only pdb_id values present in all displayed groups (consistent n)",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"[ERROR] Input CSV not found: {input_path}")
        return

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df = load_rmsd_table(
        input_path,
        include_low_coverage=args.include_low_coverage,
        diffdock_rank1_only=args.diffdock_rank1_only,
    )
    sources_to_exclude = [s.strip() for s in args.exclude_sources.split(",")] if args.exclude_sources else []
    df = exclude_sources(df, sources_to_exclude)
    if args.common_pdb_only:
        before = len(df)
        df = enforce_common_pdb_across_groups(df)
        print(f"[INFO] common-pdb filter active: {before} -> {len(df)} rows")
    if df.empty:
        print("[ERROR] No usable rows after filtering status/NaN values.")
        return

    summary_df = summarize_groups(df)
    best_df = summarize_best_per_pdb(df)

    summary_path = output_dir / "pose_rmsd_group_summary.csv"
    best_path = output_dir / "pose_rmsd_best_per_pdb.csv"
    summary_df.to_csv(summary_path, index=False)
    best_df.to_csv(best_path, index=False)
    print(f"✓ Saved group summary: {summary_path}")
    print(f"✓ Saved best-per-pdb table: {best_path}")

    plot_box_by_group(df, output_dir / "pose_rmsd_boxplot_by_group.png")
    print(f"✓ Saved plot: {output_dir / 'pose_rmsd_boxplot_by_group.png'}")
    plot_best_per_pdb(best_df, output_dir / "pose_rmsd_best_violin_by_group.png")
    print(f"✓ Saved plot: {output_dir / 'pose_rmsd_best_violin_by_group.png'}")
    plot_cdf(best_df, output_dir / "pose_rmsd_best_cdf.png")
    print(f"✓ Saved plot: {output_dir / 'pose_rmsd_best_cdf.png'}")

    print("\nTopline metrics:")
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
