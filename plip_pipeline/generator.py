"""Boltz Generator - Creates Boltz-2 YAML input files"""
import os
from typing import List, Dict, Optional
import pandas as pd
from .utils import ensure_dir, yaml_quote, get_uniprot_sequence, read_msa_query_sequence


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
        
        # All YAMLs go to a central boltz_inputs directory
        self.output_dir = os.path.join(base_dir, "boltz_inputs")
        ensure_dir(self.output_dir)
        self.source = source
        
        self.uniprot_id = config['target']['uniprot_id']
        self.max_distance = config['constraints']['max_distance']
        
        # MSA path (optional)
        self.msa_path = config['target'].get('msa_path', None)
        if self.msa_path and not os.path.isabs(self.msa_path):
            # Convert relative path to absolute
            base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            self.msa_path = os.path.join(base_path, self.msa_path)
        
        # Fetch protein sequence once. Order: UniProt, then the query sequence of
        # the configured MSA, then abort.
        #
        # This used to fall back to the literal string "SEQUENCE_PLACEHOLDER" and
        # carry on. A transient network error therefore produced a full set of
        # YAMLs that looked complete, passed every later step, and encoded no
        # protein at all. Failing loudly is the only safe option here.
        self.sequence = get_uniprot_sequence(self.uniprot_id)

        if not self.sequence and self.msa_path:
            self.sequence = read_msa_query_sequence(self.msa_path)
            if self.sequence:
                print(f"[INFO] UniProt unreachable - using the query sequence of "
                      f"{os.path.basename(self.msa_path)} ({len(self.sequence)} residues)")

        if not self.sequence:
            raise RuntimeError(
                f"Could not determine the protein sequence for {self.uniprot_id}.\n"
                f"UniProt was unreachable and no usable MSA is configured "
                f"(target.msa_path = {self.msa_path!r}).\n"
                f"Provide the sequence via an MSA file or restore network access - "
                f"generating YAMLs without it would silently produce unusable inputs."
            )
    
    def generate_all(self, df: pd.DataFrame, constraints_list: List[Dict], 
                     include_contact: bool = False, include_default: bool = True) -> None:
        """
        Generate YAML files for all structures.
        
        Creates files per structure:
        - {pdb_id}_default.yaml: Without constraints (if include_default=True)
        - {pdb_id}_{source}_pocket.yaml: With pocket constraints (if available)
        - {pdb_id}_{source}_contact.yaml: With contact constraints (if include_contact=True)
        
        Args:
            df: DataFrame with pdb_id, ligand, smiles columns
            constraints_list: List of constraint dictionaries
            include_contact: Generate contact YAMLs (requires atom-level data)
            include_default: Generate default YAMLs (set False to avoid duplicates)
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
            
            # 1. Generate default YAML (no constraints) - optional
            if include_default:
                yaml_default = self._make_yaml(ligand, smiles, constraint_type=None)
                path_default = os.path.join(self.output_dir, f"{pdb_id}_default.yaml")
                with open(path_default, 'w') as f:
                    f.write(yaml_default)
                total_files += 1
            
            # Prefix for constraint files
            prefix = f"{pdb_id}_{self.source}"
            
            # 2. Generate pocket YAML (if contacts exist)
            if contacts:
                yaml_pocket = self._make_yaml(ligand, smiles, constraint_type='pocket', contacts=contacts)
                path_pocket = os.path.join(self.output_dir, f"{prefix}_pocket.yaml")
                with open(path_pocket, 'w') as f:
                    f.write(yaml_pocket)
                total_files += 1
                
                # 3. Generate contact YAML (optional, disabled by default)
                if include_contact:
                    yaml_contact = self._make_yaml(ligand, smiles, constraint_type='contact', contacts=contacts)
                    path_contact = os.path.join(self.output_dir, f"{prefix}_contact.yaml")
                    with open(path_contact, 'w') as f:
                        f.write(yaml_contact)
                    total_files += 1
                
                num_yamls = 3 if include_contact else 2
                print(f"  ✓ {pdb_id}: Generated {num_yamls} YAMLs ({len(contacts)} constraints)")
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
        
        # MSA: use path if provided, otherwise empty
        if self.msa_path and os.path.exists(self.msa_path):
            lines.append(f"      msa: {yaml_quote(self.msa_path)}")
        else:
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
                # For ligands: use [ligand_id] without atom specification (entire molecule)
                for contact in contacts:
                    chain = contact[0]
                    residue = contact[1]
                    
                    lines.append("  - contact:")
                    lines.append(f'      token1: [{yaml_quote(ligand_id)}]')  # Entire ligand molecule
                    
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

