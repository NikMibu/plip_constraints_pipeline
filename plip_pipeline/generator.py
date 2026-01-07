"""Boltz Generator - Creates Boltz-2 YAML input files"""
import os
from typing import List, Dict, Optional
import pandas as pd
from .utils import ensure_dir, yaml_quote, get_uniprot_sequence


class BoltzGenerator:
    """Generate Boltz-2 YAML input files with optional pocket or contact constraints."""
    
    def __init__(self, config: Dict, source: str = 'crystal'):
        """
        Initialize Boltz Generator.
        
        Args:
            config: Configuration dictionary
            source: 'crystal' or 'diffdock'
        """
        base_dir = config['output']['base_dir']
        
        if source == 'crystal':
            self.output_dir = os.path.join(base_dir, "crystal_workflow", "boltz_yamls")
        elif source == 'diffdock':
            self.output_dir = os.path.join(base_dir, "diffdock_workflow", "boltz_yamls")
        
        ensure_dir(self.output_dir)
        self.source = source
        
        self.uniprot_id = config['target']['uniprot_id']
        self.max_distance = config['constraints']['max_distance']
        
        # Fetch protein sequence once
        self.sequence = get_uniprot_sequence(self.uniprot_id)
        if not self.sequence:
            print("[WARNING] Could not fetch protein sequence - using placeholder!")
            self.sequence = "SEQUENCE_PLACEHOLDER"
    
    def generate_all(self, df: pd.DataFrame, constraints_list: List[Dict]) -> None:
        """
        Generate YAML files for all structures.
        
        Creates three files per structure:
        - {pdb_id}_default.yaml: Without constraints
        - {pdb_id}_pocket.yaml: With pocket constraints (if available)
        - {pdb_id}_contact.yaml: With contact constraints (if available)
        
        Args:
            df: DataFrame with pdb_id, ligand, smiles columns
            constraints_list: List of constraint dictionaries
        """
        print("\n[GENERATE] Creating Boltz-2 YAML inputs...")
        
        # Create constraint lookup
        constraints_map = {c['pdb_id']: c for c in constraints_list}
        
        total_files = 0
        
        for idx, row in df.iterrows():
            pdb_id = row['pdb_id']
            ligand = row['ligand']
            smiles = row['smiles']
            
            # Get constraints for this structure
            constraint_data = constraints_map.get(pdb_id, {})
            contacts = constraint_data.get('contacts', [])
            
            # Prefix for output files
            prefix = f"{pdb_id}_{self.source}"
            
            # 1. Generate default YAML (no constraints)
            yaml_default = self._make_yaml(ligand, smiles, constraint_type=None)
            path_default = os.path.join(self.output_dir, f"{prefix}_default.yaml")
            with open(path_default, 'w') as f:
                f.write(yaml_default)
            total_files += 1
            
            # 2. Generate pocket YAML (if contacts exist)
            if contacts:
                yaml_pocket = self._make_yaml(ligand, smiles, constraint_type='pocket', contacts=contacts)
                path_pocket = os.path.join(self.output_dir, f"{prefix}_pocket.yaml")
                with open(path_pocket, 'w') as f:
                    f.write(yaml_pocket)
                total_files += 1
                
                # 3. Generate contact YAML (if contacts exist)
                yaml_contact = self._make_yaml(ligand, smiles, constraint_type='contact', contacts=contacts)
                path_contact = os.path.join(self.output_dir, f"{prefix}_contact.yaml")
                with open(path_contact, 'w') as f:
                    f.write(yaml_contact)
                total_files += 1
                
                print(f"  ✓ {pdb_id}: Generated 3 YAMLs ({len(contacts)} constraints)")
            else:
                print(f"  ✓ {pdb_id}: Generated default only (no constraints)")
        
        print(f"\n[GENERATE] Created {total_files} YAML files in {self.output_dir}")
    
    def _make_yaml(self, ligand_id: str, smiles: str, 
                   constraint_type: Optional[str] = None,
                   contacts: Optional[List[List[str]]] = None) -> str:
        """
        Create YAML text in Boltz-2 format.
        
        Args:
            ligand_id: Ligand identifier (3-letter code)
            smiles: SMILES string
            constraint_type: 'pocket', 'contact', or None
            contacts: Optional list of [chain, residue] contacts
            
        Returns:
            YAML string
        """
        lines = []
        
        # Header
        lines.append("version: 1")
        lines.append("")
        
        # Sequences
        lines.append("sequences:")
        lines.append("  - protein:")
        lines.append('      id: ["A"]')
        lines.append(f"      sequence: {yaml_quote(self.sequence)}")
        lines.append("      msa: empty")
        lines.append("")
        lines.append("  - ligand:")
        lines.append(f"      id: [{yaml_quote(ligand_id)}]")
        lines.append(f"      smiles: {yaml_quote(smiles)}")
        lines.append("")
        
        # Properties
        lines.append("properties:")
        lines.append("  - affinity:")
        lines.append(f"      binder: {yaml_quote(ligand_id)}")
        
        # Constraints (if provided)
        if constraint_type and contacts:
            lines.append("")
            lines.append("constraints:")
            
            if constraint_type == 'pocket':
                # Pocket constraint: binder + list of contact residues
                lines.append("  - pocket:")
                lines.append(f"      binder: {yaml_quote(ligand_id)}")
                
                # Format contacts
                contacts_str = "["
                for i, contact in enumerate(contacts):
                    if i > 0:
                        contacts_str += ", "
                    chain = contact[0]
                    residue = contact[1]
                    # Ensure residue is integer format
                    try:
                        residue_int = int(str(residue))
                        contacts_str += f'["{chain}", {residue_int}]'
                    except ValueError:
                        contacts_str += f'["{chain}", "{residue}"]'
                contacts_str += "]"
                
                lines.append(f"      contacts: {contacts_str}")
                lines.append(f"      max_distance: {self.max_distance}")
                lines.append("      force: false")
                
            elif constraint_type == 'contact':
                # Contact constraints: individual pairwise contacts between ligand and each residue
                for contact in contacts:
                    chain = contact[0]
                    residue = contact[1]
                    
                    lines.append("  - contact:")
                    lines.append(f"      token1: [{yaml_quote(ligand_id)}, 1]")  # Ligand, residue 1
                    
                    # Format token2
                    try:
                        residue_int = int(str(residue))
                        lines.append(f'      token2: ["{chain}", {residue_int}]')
                    except ValueError:
                        lines.append(f'      token2: ["{chain}", "{residue}"]')
                    
                    lines.append(f"      max_distance: {self.max_distance}")
                    lines.append("      force: false")
        
        lines.append("")
        return "\n".join(lines)

