#!/usr/bin/env python3
"""
PLIP Constraints Pipeline - Main CLI

Automated pipeline for generating Boltz-2 YAML inputs with PLIP-based pocket constraints.
"""
import argparse
import os
import sys
import yaml
import pandas as pd
from typing import Dict, List

from plip_pipeline import (
    PDBFetcher, 
    PDBCleaner, 
    PLIPAnalyzer, 
    BoltzGenerator,
    LigandPreparator,
    DiffDockRunner
)


def load_config(config_path: str) -> Dict:
    """Load configuration from YAML file."""
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        return config
    except FileNotFoundError:
        print(f"[ERROR] Config file not found: {config_path}")
        sys.exit(1)
    except yaml.YAMLError as e:
        print(f"[ERROR] Failed to parse config: {e}")
        sys.exit(1)


def run_pipeline(config: Dict, steps: List[str]) -> None:
    """
    Run the complete pipeline or specific steps.
    
    Args:
        config: Configuration dictionary
        steps: List of steps to run (fetch, clean, plip, generate, all)
    """
    output_dir = config['output']['base_dir']
    metadata_path = os.path.join(output_dir, "metadata.csv")
    
    print("=" * 70)
    print("  PLIP CONSTRAINTS PIPELINE")
    print("=" * 70)
    print(f"\nTarget: {config['target']['uniprot_id']}")
    print(f"Output: {output_dir}")
    print(f"Steps: {', '.join(steps)}")
    print()
    
    df = None
    constraints = None
    
    # Step 1: Fetch
    if "fetch" in steps or "all" in steps:
        fetcher = PDBFetcher(config)
        df = fetcher.fetch_all()
        
        if df.empty:
            print("\n[ERROR] No structures found. Aborting.")
            sys.exit(1)
        
        print(f"\n✓ Fetched {len(df)} structures")
    
    # Load existing data if not fetching
    if df is None and os.path.exists(metadata_path):
        df = pd.read_csv(metadata_path)
        print(f"\n[INFO] Loaded {len(df)} structures from metadata")
    elif df is None:
        print("\n[ERROR] No data available. Run with --steps fetch first.")
        sys.exit(1)
    
    # Step 2: Clean
    if "clean" in steps or "all" in steps:
        cleaner = PDBCleaner(config)
        df = cleaner.clean_all(df)
        
        # Update metadata
        df.to_csv(metadata_path, index=False)
        print(f"\n✓ Cleaned {len(df)} structures")
    
    # CRYSTAL WORKFLOW
    # Step 3: PLIP Analysis (Crystal)
    if "plip" in steps or "all" in steps or "crystal" in steps:
        analyzer = PLIPAnalyzer(config, source='crystal')
        constraints = analyzer.analyze_all(df)
        print(f"\n✓ Analyzed {len(constraints)} structures (Crystal)")
    
    # Step 4: Generate YAMLs (Crystal)
    if "generate" in steps or "all" in steps or "crystal" in steps:
        if constraints is None:
            analyzer = PLIPAnalyzer(config, source='crystal')
            constraints = analyzer.analyze_all(df)
        
        generator = BoltzGenerator(config, source='crystal')
        generator.generate_all(df, constraints, include_default=True)
        print(f"\n✓ Generated Crystal YAML files (incl. default)")
    
    # DIFFDOCK WORKFLOW
    # Step 5: Prepare Ligands (SMILES → SDF)
    if "prep_ligands" in steps or "all" in steps or "diffdock_full" in steps:
        prep = LigandPreparator(config)
        df = prep.prepare_all(df)
        df.to_csv(metadata_path, index=False)
        print(f"\n✓ Prepared ligands (SMILES → SDF)")
    
    # Step 6: Run DiffDock
    if "diffdock" in steps or "all" in steps or "diffdock_full" in steps:
        runner = DiffDockRunner(config)
        df = runner.run_all(df)
        df.to_csv(metadata_path, index=False)
        print(f"\n✓ DiffDock docking completed")
    
    # Step 7: PLIP Analysis (DiffDock poses)
    if "diffdock_plip" in steps or "all" in steps or "diffdock_full" in steps:
        analyzer_dd = PLIPAnalyzer(config, source='diffdock')
        constraints_dd = analyzer_dd.analyze_all(df)
        print(f"\n✓ Analyzed {len(constraints_dd)} structures (DiffDock)")
    else:
        constraints_dd = None
    
    # Step 8: Generate YAMLs (DiffDock)
    if "diffdock_yamls" in steps or "all" in steps or "diffdock_full" in steps:
        if constraints_dd is None:
            analyzer_dd = PLIPAnalyzer(config, source='diffdock')
            constraints_dd = analyzer_dd.analyze_all(df)
        
        generator_dd = BoltzGenerator(config, source='diffdock')
        generator_dd.generate_all(df, constraints_dd, include_default=False)
        print(f"\n✓ Generated DiffDock YAML files (no default)")
    
    # Final summary
    print("\n" + "=" * 70)
    print("  PIPELINE COMPLETED")
    print("=" * 70)
    
    if constraints:
        with_constraints = sum(1 for c in constraints if len(c.get('contacts', [])) > 0)
        without_constraints = len(constraints) - with_constraints
        
        print(f"\nStructures with constraints: {with_constraints}")
        print(f"Structures without constraints: {without_constraints}")
    
    print(f"\nOutput directory: {output_dir}")
    print(f"  - Raw PDBs: {os.path.join(output_dir, 'raw_pdb')}")
    print(f"  - Clean PDBs: {os.path.join(output_dir, 'clean_pdb')}")
    print(f"  - Crystal Workflow: {os.path.join(output_dir, 'crystal_workflow')}")
    print(f"  - DiffDock Workflow: {os.path.join(output_dir, 'diffdock_workflow')}")
    print(f"  - Metadata: {metadata_path}")
    print()


def show_stats(config: Dict) -> None:
    """Show statistics from existing data."""
    output_dir = config['output']['base_dir']
    metadata_path = os.path.join(output_dir, "metadata.csv")
    
    if not os.path.exists(metadata_path):
        print("[ERROR] No metadata found. Run pipeline first.")
        sys.exit(1)
    
    df = pd.read_csv(metadata_path)
    
    print("=" * 70)
    print("  PIPELINE STATISTICS")
    print("=" * 70)
    print()
    print(f"Total structures: {len(df)}")
    print()
    
    # Affinity statistics
    print("Affinity Data:")
    print(f"  Types: {', '.join(df['type'].unique())}")
    print()
    
    for aff_type in df['type'].unique():
        subset = df[df['type'] == aff_type]
        print(f"  {aff_type}:")
        print(f"    Count: {len(subset)}")
        print(f"    Min: {subset['value'].min():.2f} {subset['unit'].iloc[0]}")
        print(f"    Max: {subset['value'].max():.2f} {subset['unit'].iloc[0]}")
        print(f"    Mean: {subset['value'].mean():.2f} {subset['unit'].iloc[0]}")
        print()
    
    # File counts
    boltz_dir = os.path.join(output_dir, "boltz_inputs")
    if os.path.exists(boltz_dir):
        yaml_files = [f for f in os.listdir(boltz_dir) if f.endswith('.yaml')]
        default_files = [f for f in yaml_files if 'default' in f]
        constrained_files = [f for f in yaml_files if 'constrained' in f]
        
        print(f"Generated Files:")
        print(f"  Default YAMLs: {len(default_files)}")
        print(f"  Constrained YAMLs: {len(constrained_files)}")
        print()


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="PLIP Constraints Pipeline - Generate Boltz-2 inputs with pocket constraints",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run complete pipeline (both workflows)
  python pipeline.py --config config.yaml
  
  # Crystal workflow only
  python pipeline.py --steps fetch,clean,plip,generate
  
  # DiffDock workflow only
  python pipeline.py --steps prep_ligands,diffdock,diffdock_plip,diffdock_yamls
  
  # Or use shortcuts
  python pipeline.py --steps crystal
  python pipeline.py --steps diffdock_full
  
  # Show statistics
  python pipeline.py --stats
  
Available steps:
  fetch, clean, plip, generate          - Crystal workflow
  prep_ligands, diffdock,               - DiffDock workflow
  diffdock_plip, diffdock_yamls
  crystal, diffdock_full, all           - Shortcuts
        """
    )
    
    parser.add_argument(
        "--config",
        default="config.yaml",
        help="Path to configuration file (default: config.yaml)"
    )
    
    parser.add_argument(
        "--steps",
        default="all",
        help="Comma-separated list of steps (default: all)"
    )
    
    parser.add_argument(
        "--output",
        help="Override output directory from config"
    )
    
    parser.add_argument(
        "--stats",
        action="store_true",
        help="Show statistics only (no pipeline run)"
    )
    
    args = parser.parse_args()
    
    # Load config
    config = load_config(args.config)
    
    # Override output dir if specified
    if args.output:
        config['output']['base_dir'] = args.output
    
    # Show stats or run pipeline
    if args.stats:
        show_stats(config)
    else:
        steps = [s.strip() for s in args.steps.split(',')]
        run_pipeline(config, steps)


if __name__ == "__main__":
    main()

