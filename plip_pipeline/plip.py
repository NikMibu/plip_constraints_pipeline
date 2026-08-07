"""PLIP Analyzer - Runs PLIP and extracts pocket constraints"""
import os
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET
from typing import List, Dict, Tuple, Optional
import pandas as pd
from .utils import ensure_dir, resolve_micromamba


def resolve_plip_command(micromamba_exe: str, micromamba_env: str) -> List[str]:
    """Return the command prefix that runs PLIP on this machine.

    Prefers the micromamba environment, because that is where PLIP and its
    OpenBabel are pinned. Falls back to a plain `plip` on PATH - but loudly:
    OpenBabel places the hydrogens PLIP derives hydrogen bonds from, so a
    different one there silently changes the generated constraints.
    """
    mamba = resolve_micromamba(micromamba_exe)
    if mamba:
        return [mamba, "run", "-n", micromamba_env, "plip"]

    on_path = shutil.which("plip")
    if on_path:
        print(
            "\n" + "!" * 72 +
            f"\n[PLIP] micromamba not found - falling back to {on_path}"
            f"\n[PLIP] The '{micromamba_env}' environment pins PLIP and OpenBabel."
            "\n[PLIP] OpenBabel places the hydrogens the hydrogen bonds are derived"
            "\n[PLIP] from, so a different one here yields different constraints."
            "\n[PLIP] Check the versions before trusting these results.\n"
            + "!" * 72 + "\n"
        )
        return [on_path]

    raise RuntimeError(
        "PLIP not found.\n"
        "  - micromamba was not located. A standard install defines it as a shell\n"
        "    function, so it is invisible to Python; its binary path is exported as\n"
        "    $MAMBA_EXE. Either export that variable or set micromamba.executable\n"
        "    in the config to the absolute path.\n"
        "  - no 'plip' on PATH either.\n"
        "See SETUP.md section 2."
    )


class PLIPAnalyzer:
    """Run PLIP analysis and extract ligand interaction constraints."""
    
    def __init__(self, config: Dict, source: str = 'crystal'):
        """
        Initialize PLIP Analyzer.
        
        Args:
            config: Configuration dictionary
            source: 'crystal' or 'diffdock'
        """
        base_dir = config['output']['base_dir']
        
        if source == 'crystal':
            self.pdb_dir = os.path.join(base_dir, "clean_pdb")
            self.output_dir = os.path.join(base_dir, "crystal_workflow", "plip_reports")
        elif source == 'diffdock':
            self.pdb_dir = os.path.join(base_dir, "diffdock_workflow", "complexes")
            self.output_dir = os.path.join(base_dir, "diffdock_workflow", "plip_reports")
        
        ensure_dir(self.output_dir)
        
        self.source = source
        self.max_distance = config['constraints']['max_distance']
        self.ignore_metal = config['constraints']['ignore_metal_interactions']
        self.timeout = config['docker']['plip_timeout']
        self.micromamba_exe = config.get('micromamba', {}).get('executable', 'micromamba')
        self.micromamba_env = config.get('micromamba', {}).get('plip_env', 'plip')
        self.plip_cmd = resolve_plip_command(self.micromamba_exe, self.micromamba_env)
        print(f"[PLIP] Using: {' '.join(self.plip_cmd)}")
    
    def analyze_all(self, df: pd.DataFrame) -> List[Dict]:
        """
        Run PLIP analysis for all structures.
        
        Args:
            df: DataFrame with 'pdb_id' and 'ligand' columns
            
        Returns:
            List of constraint dictionaries
        """
        print("\n[PLIP] Running interaction analysis...")
        
        constraints_list = []
        
        for idx, row in df.iterrows():
            pdb_id = row['pdb_id']
            ligand = row['ligand']
            
            # Run PLIP
            xml_path = self._run_plip(pdb_id)
            
            if xml_path and os.path.exists(xml_path):
                # Extract constraints
                constraints = self._extract_constraints(xml_path, pdb_id, ligand)
                constraints_list.append(constraints)
                
                # Report
                num_contacts = len(constraints.get('contacts', []))
                if num_contacts > 0:
                    print(f"  ✓ {pdb_id}: {num_contacts} pocket contacts found")
                else:
                    print(f"  ⚠ {pdb_id}: No ligand interactions (possibly metal-only binding)")
            else:
                print(f"  [ERROR] {pdb_id}: PLIP failed")
                constraints_list.append({"pdb_id": pdb_id, "ligand": ligand, "contacts": []})
        
        return constraints_list
    
    def _run_plip(self, pdb_id: str) -> Optional[str]:
        """
        Run PLIP analysis using micromamba.
        
        Args:
            pdb_id: PDB ID
            
        Returns:
            Path to XML report or None on error
        """
        if self.source == 'crystal':
            pdb_file = os.path.join(self.pdb_dir, f"{pdb_id}_clean.pdb")
        elif self.source == 'diffdock':
            pdb_file = os.path.join(self.pdb_dir, f"{pdb_id}_diffdock_complex.pdb")
        
        if not os.path.exists(pdb_file):
            print(f"  [ERROR] {pdb_id}: Clean PDB not found")
            return None
        
        # Output directory for this structure
        structure_output = os.path.join(self.output_dir, pdb_id)
        ensure_dir(structure_output)
        
        # The PDB file now contains both protein and ligand
        cmd = self.plip_cmd + [
            "-f", os.path.abspath(pdb_file),
            "-x",
            "-o", os.path.abspath(structure_output)
        ]
        
        try:
            result = subprocess.run(
                cmd,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=self.timeout
            )
            
            # PLIP creates either report.xml or {filename}_report.xml
            xml_path1 = os.path.join(structure_output, "report.xml")
            xml_path2 = os.path.join(structure_output, f"{os.path.basename(pdb_file).replace('.pdb', '')}_report.xml")
            
            if os.path.exists(xml_path1):
                return xml_path1
            elif os.path.exists(xml_path2):
                return xml_path2
            else:
                print(f"  [ERROR] {pdb_id}: XML report not found")
                return None
            
        except subprocess.TimeoutExpired:
            print(f"  [ERROR] {pdb_id}: PLIP timeout")
            return None
        except subprocess.CalledProcessError as e:
            stderr = e.stderr.decode() if e.stderr else "No error output"
            print(f"  [ERROR] {pdb_id}: PLIP failed - {stderr[:200]}")
            return None
        except Exception as e:
            print(f"  [ERROR] {pdb_id}: PLIP failed - {e}")
            return None
    
    def _extract_constraints(self, xml_path: str, pdb_id: str, target_ligand: str) -> Dict:
        """
        Extract pocket constraints from PLIP XML.
        
        Args:
            xml_path: Path to PLIP XML report
            pdb_id: PDB ID
            target_ligand: Target ligand ID
            
        Returns:
            Dictionary with pdb_id, ligand, and contacts
        """
        try:
            tree = ET.parse(xml_path)
            root = tree.getroot()
        except Exception as e:
            print(f"  [ERROR] {pdb_id}: Failed to parse XML - {e}")
            return {"pdb_id": pdb_id, "ligand": target_ligand, "contacts": []}
        
        all_contacts = []
        
        # Iterate through binding sites
        for site in root.findall(".//bindingsite"):
            # Check if this is our target ligand
            hetid_elem = site.find(".//hetid")
            if hetid_elem is None:
                continue
            
            hetid = hetid_elem.text.strip() if hetid_elem.text else ""
            
            # For DiffDock source, RDKit generates ligands with hetid "UNL"
            # For crystal source, check exact match
            if self.source == 'crystal':
                if hetid != target_ligand:
                    continue
            elif self.source == 'diffdock':
                # Accept UNL or matching hetid
                if hetid not in ['UNL', target_ligand]:
                    continue
            
            # Check ligand type - skip metals if configured
            if self.ignore_metal:
                ligtype_elem = site.find(".//ligtype")
                ligtype = ligtype_elem.text.strip().upper() if ligtype_elem is not None and ligtype_elem.text else ""
                
                if ligtype in {"ION", "METAL", "COFACTOR"}:
                    continue
            
            # Extract contacts from binding site residues
            contacts = self._parse_bs_residues(site)
            
            # Fallback: Parse individual interactions
            if not contacts:
                contacts = self._parse_interactions(site)
            
            all_contacts.extend(contacts)
        
        # Remove duplicates
        unique_contacts = []
        seen = set()
        for contact in all_contacts:
            key = tuple(contact)
            if key not in seen:
                seen.add(key)
                unique_contacts.append(contact)
        
        return {
            "pdb_id": pdb_id,
            "ligand": target_ligand,
            "contacts": unique_contacts
        }
    
    def _parse_bs_residues(self, site) -> List[List[str]]:
        """Parse binding site residues from bs_residues section."""
        contacts = []
        
        bs_residues = site.find("bs_residues")
        if bs_residues is None:
            return contacts
        
        for bs_res in bs_residues.findall("bs_residue"):
            # Only include residues marked as contact
            contact_flag = bs_res.get("contact", "False").lower() == "true"
            if not contact_flag:
                continue
            
            chain, resnr = self._parse_residue_id(bs_res.text)
            if chain and resnr:
                contacts.append([chain, str(resnr)])
        
        return contacts
    
    def _parse_interactions(self, site) -> List[List[str]]:
        """Parse individual interactions (H-bonds, hydrophobic, etc.)."""
        contacts = []
        
        interaction_types = [
            ".//hydrogen_bond",
            ".//hydrophobic_interaction",
            ".//pi_stacking",
            ".//salt_bridge",
            ".//halogen_bond"
        ]
        
        for interaction_type in interaction_types:
            for interaction in site.findall(interaction_type):
                chain_elem = interaction.find("reschain")
                resnr_elem = interaction.find("resnr")
                
                if chain_elem is not None and resnr_elem is not None:
                    chain = chain_elem.text.strip()
                    resnr = resnr_elem.text.strip()
                    contacts.append([chain, resnr])
        
        return contacts
    
    def _parse_residue_id(self, residue_text: str) -> Tuple[Optional[str], Optional[int]]:
        """
        Parse residue identifier like '199A' or '199:A'.
        
        Args:
            residue_text: Residue string from PLIP
            
        Returns:
            (chain, residue_number) or (None, None)
        """
        if not residue_text:
            return None, None
        
        residue_text = residue_text.strip()
        
        # Handle format '199:A'
        if ":" in residue_text:
            parts = residue_text.split(":")
            try:
                return parts[1], int(parts[0])
            except (ValueError, IndexError):
                return None, None
        
        # Handle format '199A'
        match = re.match(r"(-?\d+)([A-Za-z]+)$", residue_text)
        if match:
            resnr = int(match.group(1))
            chain = match.group(2)
            return chain, resnr
        
        return None, None

