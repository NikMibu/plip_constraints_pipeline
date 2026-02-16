#!/usr/bin/env python3
"""Parse Boltz-2 affinity predictions and merge with experimental data.

Scans Boltz prediction outputs for three scenarios:
  - default: no constraints
  - crystal_pocket: constraints from cocrystal structure
  - diffdock_pocket: constraints from DiffDock pose

Outputs a summary CSV with predictions and experimental Ki values.
"""
import argparse
import json
import os
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd


def safe_mean_from_prefix(d: dict, prefix: str) -> float:
    """Extract mean from JSON fields matching prefix."""
    vals = [float(v) for k, v in d.items() if k.startswith(prefix) and isinstance(v, (int, float))]
    return float(np.mean(vals)) if vals else float("nan")


def find_predictions_dir(root: str) -> Path:
    """Find the predictions directory (may be nested)."""
    root = Path(root)
    
    # Direct predictions folder
    direct = root / "predictions"
    if direct.is_dir():
        return direct
    
    # Nested: boltz_output/boltz_results_*/predictions
    for entry in root.iterdir():
        if entry.is_dir():
            nested = entry / "predictions"
            if nested.is_dir():
                return nested
    
    return root


def parse_predictions(pred_dir: Path, scenarios: List[str]) -> Dict[str, Dict[str, Dict[str, float]]]:
    """Parse Boltz predictions for all PDB IDs and scenarios.
    
    Returns:
        {pdb_id: {scenario: {metric: value}}}
    """
    predictions = {}
    
    for entry in pred_dir.iterdir():
        if not entry.is_dir():
            continue
        
        name = entry.name
        
        # Parse folder name: {pdb_id}_{scenario}
        matched_scenario = None
        pdb_id = None
        
        for scenario in scenarios:
            if name.endswith(f"_{scenario}"):
                matched_scenario = scenario
                pdb_id = name[: -len(f"_{scenario}")]
                break
        
        if not matched_scenario or not pdb_id:
            continue
        
        # Find affinity JSON
        affinity_files = [f for f in entry.iterdir() if f.name.startswith("affinity_") and f.name.endswith(".json")]
        if not affinity_files:
            continue
        
        affinity_path = affinity_files[0]
        try:
            with open(affinity_path) as f:
                data = json.load(f)
        except Exception as e:
            print(f"[WARNING] Failed to parse {affinity_path}: {e}")
            continue
        
        # Extract metrics
        if pdb_id not in predictions:
            predictions[pdb_id] = {}
        
        predictions[pdb_id][matched_scenario] = {
            "prob_mean": safe_mean_from_prefix(data, "affinity_probability_binary"),
            "pred_value_mean": safe_mean_from_prefix(data, "affinity_pred_value"),
        }
    
    return predictions


def main():
    parser = argparse.ArgumentParser(description="Parse Boltz predictions and merge with experimental data")
    parser.add_argument("--predictions-dir", required=True, help="Path to Boltz predictions directory")
    parser.add_argument("--metadata", default="./output/metadata.csv", help="Path to metadata.csv with experimental Ki values")
    parser.add_argument("--output", default="./results/predictions_summary.csv", help="Output CSV path")
    args = parser.parse_args()
    
    # Load metadata
    metadata_path = Path(args.metadata)
    if not metadata_path.exists():
        print(f"[ERROR] Metadata file not found: {metadata_path}")
        return
    
    print(f"Loading metadata from: {metadata_path}")
    metadata_df = pd.read_csv(metadata_path)
    
    if "pdb_id" not in metadata_df.columns:
        print("[ERROR] metadata.csv must contain 'pdb_id' column")
        return
    
    # Parse predictions
    pred_dir = find_predictions_dir(args.predictions_dir)
    print(f"Scanning predictions in: {pred_dir}")
    
    scenarios = ["default", "crystal_pocket", "diffdock_pocket"]
    predictions = parse_predictions(pred_dir, scenarios)
    
    print(f"Found predictions for {len(predictions)} PDB IDs")
    
    # Build output dataframe
    records = []
    for _, row in metadata_df.iterrows():
        pdb_id = row["pdb_id"]
        
        if pdb_id not in predictions:
            print(f"[WARNING] No predictions found for {pdb_id}")
            continue
        
        for scenario in scenarios:
            if scenario not in predictions[pdb_id]:
                print(f"[WARNING] Missing scenario '{scenario}' for {pdb_id}")
                continue
            
            pred = predictions[pdb_id][scenario]
            
            record = {
                "pdb_id": pdb_id,
                "ligand": row.get("ligand", ""),
                "scenario": scenario,
                "exp_Ki_nM": float(row.get("value", np.nan)),
                "exp_type": row.get("type", ""),
                "exp_unit": row.get("unit", ""),
                "boltz_prob_mean": pred["prob_mean"],
                "boltz_pred_value_mean": pred["pred_value_mean"],
            }
            records.append(record)
    
    if not records:
        print("[ERROR] No predictions could be matched with metadata")
        return
    
    output_df = pd.DataFrame(records)
    
    # Save output
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_df.to_csv(output_path, index=False)
    
    print(f"\n✓ Created predictions summary: {output_path}")
    print(f"  Total records: {len(output_df)}")
    print(f"  PDB IDs: {output_df['pdb_id'].nunique()}")
    print(f"  Scenarios: {', '.join(output_df['scenario'].unique())}")
    
    # Statistics
    for scenario in scenarios:
        count = (output_df["scenario"] == scenario).sum()
        print(f"    - {scenario}: {count} predictions")


if __name__ == "__main__":
    main()
