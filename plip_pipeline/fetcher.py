"""PDB Fetcher - Downloads structures and metadata from RCSB"""
import os
import time
import requests
import pandas as pd
from typing import List, Dict, Optional
from .utils import ensure_dir


class PDBFetcher:
    """Fetch PDB structures and affinity data from RCSB API."""
    
    def __init__(self, config: Dict):
        """
        Initialize PDB Fetcher.
        
        Args:
            config: Configuration dictionary
        """
        self.uniprot_id = config['target']['uniprot_id']
        self.limit = config['fetch']['limit']
        self.allowed_types = config['fetch']['allowed_types']
        self.pdb_ids_file = config['fetch'].get('pdb_ids_file') or None
        self.output_dir = os.path.join(config['output']['base_dir'], "raw_pdb")
        ensure_dir(self.output_dir)

    def _read_pdb_ids_file(self) -> List[str]:
        """Read a fixed list of PDB IDs, one per line, '#' for comments."""
        path = os.path.expanduser(os.path.expandvars(self.pdb_ids_file))
        if not os.path.isabs(path):
            root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            path = os.path.join(root, path)
        if not os.path.exists(path):
            raise FileNotFoundError(
                f"fetch.pdb_ids_file points to {path}, which does not exist."
            )
        with open(path) as f:
            ids = [line.split("#")[0].strip().upper() for line in f]
        return [i for i in ids if i]

    def fetch_all(self) -> pd.DataFrame:
        """
        Fetch PDB IDs, download structures, and collect metadata.

        Returns:
            DataFrame with columns: pdb_id, ligand, type, value, unit, smiles
        """
        # A fixed list reproduces an earlier run exactly. The RCSB search is
        # run against live annotations, so the same query returns a different
        # set as entries and affinity annotations are added over time.
        if self.pdb_ids_file:
            pdb_ids = self._read_pdb_ids_file()
            print(f"\n[FETCH] Using the fixed list {self.pdb_ids_file}: "
                  f"{len(pdb_ids)} structures (no search)")
        else:
            print(f"\n[FETCH] Searching structures for UniProt ID: {self.uniprot_id}")
            pdb_ids = self._search_pdb_ids()
            print(f"[FETCH] Found {len(pdb_ids)} structures with affinity data")
        
        # 2. Fetch metadata for each structure
        data = self._fetch_metadata(pdb_ids)
        
        if not data:
            print("[FETCH] No valid structures found!")
            return pd.DataFrame()
        
        df = pd.DataFrame(data)

        # Say it here, not three steps later. Without SMILES the DiffDock
        # workflow silently produces nothing, and the cause is this call.
        usable = df['smiles'].apply(lambda s: isinstance(s, str) and s.strip() not in ("", "N/A"))
        if not usable.any():
            print(f"\n[FETCH] WARNING: no SMILES for any of the {len(df)} ligands. "
                  f"The crystal workflow still works; the DiffDock workflow needs "
                  f"them and will find nothing to dock.")
        elif (~usable).any():
            missing = df.loc[~usable, 'ligand'].tolist()
            print(f"\n[FETCH] WARNING: no SMILES for {(~usable).sum()} of {len(df)} "
                  f"ligands: {', '.join(map(str, missing[:10]))}"
                  f"{' ...' if len(missing) > 10 else ''}")

        # 3. Download PDB files
        print(f"\n[FETCH] Downloading PDB files...")
        for pdb_id in df['pdb_id'].unique():
            self._download_pdb(pdb_id)
        
        # Save metadata
        csv_path = os.path.join(os.path.dirname(self.output_dir), "metadata.csv")
        df.to_csv(csv_path, index=False)
        print(f"[FETCH] Metadata saved to: {csv_path}")
        
        return df
    
    def _search_pdb_ids(self) -> List[str]:
        """Search RCSB for PDB IDs matching criteria."""
        url = "https://search.rcsb.org/rcsbsearch/v2/query"
        
        query = {
            "query": {
                "type": "group",
                "logical_operator": "and",
                "nodes": [
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                            "operator": "exact_match",
                            "value": self.uniprot_id
                        }
                    },
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_binding_affinity.comp_id",
                            "operator": "exists"
                        }
                    }
                ]
            },
            "request_options": {"return_all_hits": True},
            "return_type": "entry"
        }
        
        try:
            response = requests.post(url, json=query, timeout=30)
            response.raise_for_status()
            results = response.json().get("result_set", [])
            return [item["identifier"] for item in results]
        except Exception as e:
            print(f"[ERROR] PDB search failed: {e}")
            return []
    
    def _fetch_metadata(self, pdb_ids: List[str]) -> List[Dict]:
        """Fetch affinity metadata for PDB IDs."""
        data = []
        
        for i, pdb_id in enumerate(pdb_ids):
            if len(data) >= self.limit:
                break
            
            try:
                url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
                response = requests.get(url, timeout=10)
                response.raise_for_status()
                entry_data = response.json()
                
                if "rcsb_binding_affinity" in entry_data:
                    for aff in entry_data["rcsb_binding_affinity"]:
                        aff_type = aff.get("type")
                        
                        if aff_type in self.allowed_types:
                            ligand = aff.get("comp_id")
                            smiles = self._fetch_smiles(ligand)
                            
                            data.append({
                                "pdb_id": pdb_id,
                                "ligand": ligand,
                                "type": aff_type,
                                "value": aff.get("value"),
                                "unit": aff.get("unit", "nM"),
                                "smiles": smiles
                            })
                            
                            print(f"  [{i+1}/{len(pdb_ids)}] {pdb_id}: {aff_type}={aff.get('value')} {aff.get('unit', 'nM')}")
                            break  # Only take first matching affinity per structure
                            
            except Exception as e:
                print(f"  [WARNING] Failed to fetch {pdb_id}: {e}")
            
            time.sleep(0.1)  # Rate limiting
        
        return data
    
    # Field names in rcsb_chem_comp_descriptor, most preferred first.
    # RCSB renamed "smiles" to "SMILES"; the lower-case name is kept so older
    # deployments still work. SMILES_stereo is a last resort - it carries
    # stereochemistry the plain field does not, so it can yield a different 3D
    # ligand and different docking poses. It is not the field the published
    # results were produced with.
    _SMILES_KEYS = ("SMILES", "smiles", "SMILES_stereo", "smiles_stereo")

    def _fetch_smiles(self, ligand_code: str) -> str:
        """
        Fetch the SMILES string for a ligand from RCSB.

        Args:
            ligand_code: 3-letter ligand code

        Returns:
            SMILES string, or "N/A" if none could be retrieved
        """
        url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand_code}"
        try:
            response = requests.get(url, timeout=15)
            response.raise_for_status()
            descriptor = response.json().get("rcsb_chem_comp_descriptor") or {}
        except Exception as e:
            print(f"  [WARN] {ligand_code}: could not query RCSB - {e}")
            return "N/A"

        for key in self._SMILES_KEYS:
            value = descriptor.get(key)
            if isinstance(value, str) and value.strip():
                if key not in ("SMILES", "smiles"):
                    print(f"  [WARN] {ligand_code}: falling back to '{key}'")
                return value.strip()

        # Silence here is expensive: without SMILES the whole DiffDock workflow
        # is skipped later, and nothing upstream says why.
        print(f"  [WARN] {ligand_code}: no SMILES in the RCSB response "
              f"(fields present: {sorted(descriptor) or 'none'})")
        return "N/A"
    
    def _download_pdb(self, pdb_id: str) -> None:
        """
        Download PDB file from RCSB.
        
        Args:
            pdb_id: 4-letter PDB ID
        """
        filepath = os.path.join(self.output_dir, f"{pdb_id}.pdb")
        
        if os.path.exists(filepath):
            return
        
        try:
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            
            with open(filepath, 'wb') as f:
                f.write(response.content)
            
            print(f"  ✓ Downloaded {pdb_id}.pdb")
            
        except Exception as e:
            print(f"  [ERROR] Failed to download {pdb_id}: {e}")

