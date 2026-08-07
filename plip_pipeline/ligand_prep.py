"""Ligand Preparator - Convert SMILES to 3D SDF files"""
import os
from typing import Dict
import pandas as pd
from .utils import ensure_dir


class LigandPreparator:
    """Convert SMILES strings to 3D SDF files using RDKit."""
    
    def __init__(self, config: Dict):
        """
        Initialize Ligand Preparator.
        
        Args:
            config: Configuration dictionary
        """
        base_dir = config['output']['base_dir']
        self.output_dir = os.path.join(base_dir, "diffdock_workflow", "ligands_sdf")
        ensure_dir(self.output_dir)
    
    def prepare_all(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Convert all SMILES to SDF files.
        
        Args:
            df: DataFrame with pdb_id, ligand, smiles columns
            
        Returns:
            DataFrame with added sdf_path column
        """
        print("\n[LIGAND PREP] Converting SMILES to 3D SDF files...")
        
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
        except ImportError:
            print("[ERROR] RDKit not installed. Install with: pip install rdkit")
            return df
        
        successful = []
        sdf_paths = []
        
        for idx, row in df.iterrows():
            pdb_id = row['pdb_id']
            ligand_id = row['ligand']
            smiles = row['smiles']
            
            # The guard has to survive pandas: when every SMILES in the CSV is
            # "N/A", the column is read as float NaN. `nan == "N/A"` is False
            # and `not nan` is False too, so a bare check lets NaN through and
            # RDKit fails deep inside Boost with an unreadable converter error.
            if not isinstance(smiles, str) or not smiles.strip() or smiles.strip() == "N/A":
                print(f"  [SKIP] {pdb_id}: no usable SMILES in the metadata "
                      f"(got {smiles!r})")
                sdf_paths.append(None)
                continue
            smiles = smiles.strip()
            
            try:
                # Parse SMILES
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    print(f"  [ERROR] {pdb_id}: Invalid SMILES")
                    sdf_paths.append(None)
                    continue
                
                # Add hydrogens
                mol = Chem.AddHs(mol)
                
                # Generate 3D coordinates
                result = AllChem.EmbedMolecule(mol, randomSeed=42)
                if result != 0:
                    print(f"  [ERROR] {pdb_id}: Could not embed molecule")
                    sdf_paths.append(None)
                    continue
                
                # Optimize geometry (try both old and new RDKit API)
                try:
                    AllChem.UFFOptimizeMolecule(mol)
                except AttributeError:
                    AllChem.UFFOptimize(mol)
                
                # Save as SDF
                sdf_path = os.path.join(self.output_dir, f"{pdb_id}_{ligand_id}.sdf")
                writer = Chem.SDWriter(sdf_path)
                writer.write(mol)
                writer.close()
                
                print(f"  ✓ {pdb_id}: {ligand_id} → SDF")
                successful.append(idx)
                sdf_paths.append(sdf_path)
                
            except Exception as e:
                print(f"  [ERROR] {pdb_id}: {e}")
                sdf_paths.append(None)
        
        # Add SDF paths to dataframe
        df['sdf_path'] = sdf_paths
        
        print(f"\n[LIGAND PREP] Successfully prepared {len(successful)}/{len(df)} ligands")
        
        # Return full dataframe (with sdf_path column added, even if None for failures)
        return df

