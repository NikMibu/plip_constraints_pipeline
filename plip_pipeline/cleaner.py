"""PDB Cleaner - Removes unwanted ligands and prepares structures for PLIP"""
import os
from typing import List, Dict
from Bio import PDB
from Bio.PDB import Select
import pandas as pd
from .utils import ensure_dir


class CleanSelect(Select):
    """BioPython Select class for filtering residues."""
    
    def __init__(self, target_ligand: str, whitelist: List[str], blacklist: List[str]):
        """
        Initialize selector.
        
        Args:
            target_ligand: 3-letter code of ligand to keep
            whitelist: Residues to always keep (e.g., ["ZN"])
            blacklist: Residues to remove
        """
        self.target_ligand = target_ligand
        self.whitelist = whitelist
        self.blacklist = blacklist
    
    def accept_residue(self, residue):
        """
        Decide whether to keep a residue.
        
        Returns:
            1 to keep, 0 to remove
        """
        resname = residue.get_resname().strip()
        
        # Keep target ligand and whitelist
        if resname == self.target_ligand or resname in self.whitelist:
            return 1
        
        # Remove blacklist and water
        if resname in self.blacklist or resname == "HOH":
            return 0
        
        # Keep protein (standard amino acids)
        if residue.id[0] == " ":
            return 1
        
        # Remove other heteroatoms
        return 0


class PDBCleaner:
    """Clean PDB structures for PLIP analysis."""
    
    def __init__(self, config: Dict):
        """
        Initialize PDB Cleaner.
        
        Args:
            config: Configuration dictionary
        """
        base_dir = config['output']['base_dir']
        self.raw_dir = os.path.join(base_dir, "raw_pdb")
        self.clean_dir = os.path.join(base_dir, "clean_pdb")
        ensure_dir(self.clean_dir)
        
        self.whitelist = config['target']['ligand_whitelist']
        self.blacklist = config['target']['forbidden_ligands']
    
    def clean_all(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Clean all PDB structures from DataFrame.
        
        Args:
            df: DataFrame with 'pdb_id' and 'ligand' columns
            
        Returns:
            DataFrame with successfully cleaned structures
        """
        print("\n[CLEAN] Cleaning PDB structures...")
        
        parser = PDB.PDBParser(QUIET=True)
        io = PDB.PDBIO()
        
        successful = []
        
        for idx, row in df.iterrows():
            pdb_id = row['pdb_id']
            ligand = row['ligand']
            
            try:
                # Load structure
                input_path = os.path.join(self.raw_dir, f"{pdb_id}.pdb")
                structure = parser.get_structure(pdb_id, input_path)
                
                # Save cleaned structure
                output_path = os.path.join(self.clean_dir, f"{pdb_id}_clean.pdb")
                io.set_structure(structure)
                io.save(output_path, CleanSelect(ligand, self.whitelist, self.blacklist))
                
                print(f"  ✓ {pdb_id}: Kept {ligand} + {self.whitelist}")
                successful.append(idx)
                
            except Exception as e:
                print(f"  [ERROR] Failed to clean {pdb_id}: {e}")
        
        # Return only successfully cleaned structures
        cleaned_df = df.loc[successful].copy()
        print(f"\n[CLEAN] Successfully cleaned {len(cleaned_df)}/{len(df)} structures")
        
        return cleaned_df

