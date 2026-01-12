#!/usr/bin/env python3
"""Affinity prediction analysis for PLIP Constraints Pipeline.

Compares Boltz-2 affinity predictions across three scenarios:
  - default: no constraints
  - crystal_pocket: cocrystal structure constraints
  - diffdock_pocket: DiffDock pose constraints

Generates regression plots and statistical metrics.
"""
from pathlib import Path
from typing import List, Dict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCENARIOS: List[Dict[str, str]] = [
    {"key": "default", "label": "Default (No Constraints)"},
    {"key": "crystal_pocket", "label": "Cocrystal Constraints"},
    {"key": "diffdock_pocket", "label": "DiffDock Constraints"},
]


def to_log10_um(ki_nanomolar: float) -> float:
    """Convert Ki in nM to log10(IC50) in µM."""
    return float(np.log10(ki_nanomolar / 1000.0))


def to_pic50_kcal(log10_um: float) -> float:
    """Convert log10(IC50) to the pIC50 kcal/mol scale used by Boltz."""
    return (6.0 - log10_um) * 1.364


def load_predictions(predictions_csv: Path) -> pd.DataFrame:
    """Load and convert predictions summary."""
    df = pd.read_csv(predictions_csv)
    
    # Convert experimental Ki to log10(µM) and pIC50 kcal/mol
    df["exp_log10_um"] = df["exp_Ki_nM"].apply(to_log10_um)
    df["exp_pic50_kcal"] = df["exp_log10_um"].apply(to_pic50_kcal)
    
    # Boltz outputs are already in log10(µM) scale
    df["pred_log10_um"] = df["boltz_pred_value_mean"]
    df["pred_pic50_kcal"] = df["pred_log10_um"].apply(to_pic50_kcal)
    
    # Calculate deltas
    df["delta_log10_um"] = df["pred_log10_um"] - df["exp_log10_um"]
    df["delta_pic50_kcal"] = df["pred_pic50_kcal"] - df["exp_pic50_kcal"]
    
    return df


def calculate_regression_metrics(df: pd.DataFrame) -> pd.DataFrame:
    """Calculate regression metrics for each scenario."""
    metrics = []
    
    for scenario in SCENARIOS:
        subset = df[df["scenario"] == scenario["key"]].copy()
        
        # Drop NaN values
        subset = subset.dropna(subset=["pred_log10_um", "exp_log10_um"])
        
        if len(subset) < 2:
            print(f"[WARNING] Not enough data for scenario '{scenario['key']}' ({len(subset)} samples)")
            continue
        
        x = subset["pred_log10_um"].to_numpy()
        y = subset["exp_log10_um"].to_numpy()
        
        # Linear regression
        slope, intercept = np.polyfit(x, y, 1)
        y_fit = slope * x + intercept
        
        # R²
        ss_res = np.sum((y - y_fit) ** 2)
        ss_tot = np.sum((y - np.mean(y)) ** 2)
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
        
        # Error metrics (pred - exp)
        residuals = x - y
        rmse = float(np.sqrt(np.mean(residuals ** 2)))
        mae = float(np.mean(np.abs(residuals)))
        bias = float(np.mean(residuals))
        
        # Correlation coefficients
        pearson = float(np.corrcoef(x, y)[0, 1])
        
        # Spearman (rank correlation)
        rank_pred = pd.Series(x).rank(method="average").to_numpy()
        rank_exp = pd.Series(y).rank(method="average").to_numpy()
        spearman = float(np.corrcoef(rank_pred, rank_exp)[0, 1])
        
        metrics.append({
            "scenario": scenario["key"],
            "scenario_label": scenario["label"],
            "samples": len(subset),
            "slope": slope,
            "intercept": intercept,
            "r2": r2,
            "rmse_log10_um": rmse,
            "mae_log10_um": mae,
            "bias_log10_um": bias,
            "pearson_r": pearson,
            "spearman_rho": spearman,
        })
    
    return pd.DataFrame(metrics)


def plot_regressions(summary_df: pd.DataFrame, metrics_df: pd.DataFrame, out_path: Path) -> None:
    """Create 3-panel regression plot comparing scenarios."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5), sharex=True, sharey=True)
    
    # Calculate plot range
    exp_min = summary_df["exp_log10_um"].min() - 0.5
    exp_max = summary_df["exp_log10_um"].max() + 0.5
    grid = np.linspace(exp_min, exp_max, 100)
    
    for ax, scenario in zip(axes, SCENARIOS):
        subset = summary_df[summary_df["scenario"] == scenario["key"]].dropna(subset=["pred_log10_um", "exp_log10_um"])
        
        if len(subset) == 0:
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(scenario["label"])
            continue
        
        # Scatter plot
        ax.scatter(
            subset["pred_log10_um"],
            subset["exp_log10_um"],
            alpha=0.7,
            s=80,
            label="Predictions",
            color="#1f77b4",
            edgecolors="white",
            linewidth=0.5,
        )
        
        # 1:1 reference line
        ax.plot(grid, grid, linestyle="--", color="#555555", linewidth=1.5, label="y = x", alpha=0.6)
        
        # Regression line
        m_row = metrics_df[metrics_df["scenario"] == scenario["key"]]
        if not m_row.empty:
            m_row = m_row.iloc[0]
            y_reg = m_row["slope"] * grid + m_row["intercept"]
            ax.plot(grid, y_reg, color="#d62728", linewidth=2, label=f"Fit (R²={m_row['r2']:.2f})")
            
            # Title with metrics
            ax.set_title(
                f"{scenario['label']}\n"
                f"R² = {m_row['r2']:.3f} | RMSE = {m_row['rmse_log10_um']:.2f} | n = {m_row['samples']}",
                fontsize=11,
                pad=10
            )
        else:
            ax.set_title(scenario["label"])
        
        # Labels and styling
        ax.set_xlabel("Predicted log₁₀(IC50) [µM]", fontsize=10)
        ax.set_xlim(exp_min, exp_max)
        ax.set_ylim(exp_min, exp_max)
        ax.grid(True, alpha=0.3, linestyle=":", linewidth=0.5)
        ax.set_aspect("equal", adjustable="box")
        
        if ax is axes[0]:
            ax.set_ylabel("Experimental log₁₀(Ki) [µM]", fontsize=10)
            ax.legend(loc="upper left", fontsize=9, framealpha=0.9)
    
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"✓ Saved regression plot: {out_path}")


def plot_error_distribution(summary_df: pd.DataFrame, out_path: Path) -> None:
    """Plot error distribution for each scenario."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 4), sharex=True, sharey=True)
    
    for ax, scenario in zip(axes, SCENARIOS):
        subset = summary_df[summary_df["scenario"] == scenario["key"]].dropna(subset=["delta_log10_um"])
        
        if len(subset) == 0:
            continue
        
        errors = subset["delta_log10_um"].to_numpy()
        
        # Histogram
        ax.hist(errors, bins=15, alpha=0.7, color="#1f77b4", edgecolor="black")
        
        # Mean line
        mean_err = np.mean(errors)
        ax.axvline(mean_err, color="#d62728", linestyle="--", linewidth=2, label=f"Mean = {mean_err:.2f}")
        ax.axvline(0, color="#2ca02c", linestyle="-", linewidth=1.5, label="Zero error", alpha=0.6)
        
        ax.set_title(f"{scenario['label']}\nBias = {mean_err:.2f}, MAE = {np.mean(np.abs(errors)):.2f}")
        ax.set_xlabel("Prediction Error (Pred - Exp) [log₁₀(µM)]")
        ax.legend()
        ax.grid(True, alpha=0.3, axis="y")
        
        if ax is axes[0]:
            ax.set_ylabel("Frequency")
    
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"✓ Saved error distribution plot: {out_path}")


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="Analyze Boltz affinity predictions")
    parser.add_argument(
        "--predictions",
        default="./results/predictions_summary.csv",
        help="Path to predictions summary CSV"
    )
    parser.add_argument(
        "--output-dir",
        default="./results/analysis",
        help="Output directory for analysis results"
    )
    args = parser.parse_args()
    
    predictions_path = Path(args.predictions)
    if not predictions_path.exists():
        print(f"[ERROR] Predictions file not found: {predictions_path}")
        print("Run parse_boltz_predictions.py first!")
        return
    
    print(f"Loading predictions from: {predictions_path}")
    summary_df = load_predictions(predictions_path)
    
    print(f"Loaded {len(summary_df)} predictions")
    print(f"  PDB IDs: {summary_df['pdb_id'].nunique()}")
    print(f"  Scenarios: {', '.join(summary_df['scenario'].unique())}")
    
    # Calculate metrics
    print("\nCalculating regression metrics...")
    metrics_df = calculate_regression_metrics(summary_df)
    
    # Prepare output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save enhanced summary with conversions
    summary_path = output_dir / "affinity_predictions_summary.csv"
    summary_df.to_csv(summary_path, index=False)
    print(f"✓ Saved summary: {summary_path}")
    
    # Save metrics
    metrics_path = output_dir / "affinity_regression_metrics.csv"
    metrics_df.to_csv(metrics_path, index=False)
    print(f"✓ Saved metrics: {metrics_path}")
    
    # Print metrics table
    print("\n" + "=" * 80)
    print("REGRESSION METRICS")
    print("=" * 80)
    print(metrics_df.to_string(index=False))
    print("=" * 80)
    
    # Generate plots
    print("\nGenerating plots...")
    plot_regressions(summary_df, metrics_df, output_dir / "affinity_regression.png")
    plot_error_distribution(summary_df, output_dir / "error_distribution.png")
    
    print(f"\n✓ Analysis complete! Results saved to: {output_dir}")


if __name__ == "__main__":
    main()
