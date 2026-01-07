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
        self.output_dir = os.path.join(config['output']['base_dir'], "raw_pdb")
        ensure_dir(self.output_dir)
        
    def fetch_all(self) -> pd.DataFrame:
        """
        Fetch PDB IDs, download structures, and collect metadata.
        
        Returns:
            DataFrame with columns: pdb_id, ligand, type, value, unit, smiles
        """
        print(f"\n[FETCH] Searching structures for UniProt ID: {self.uniprot_id}")
        
        # 1. Search for PDB IDs
        pdb_ids = self._search_pdb_ids()
        print(f"[FETCH] Found {len(pdb_ids)} structures with affinity data")
        
        # 2. Fetch metadata for each structure
        data = self._fetch_metadata(pdb_ids)
        
        if not data:
            print("[FETCH] No valid structures found!")
            return pd.DataFrame()
        
        df = pd.DataFrame(data)
        
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
    
    def _fetch_smiles(self, ligand_code: str) -> Optional[str]:
        """
        Fetch SMILES string from RCSB.
        
        Args:
            ligand_code: 3-letter ligand code
            
        Returns:
            SMILES string or "N/A" on error
        """
        try:
            url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand_code}"
            response = requests.get(url, timeout=5)
            response.raise_for_status()
            data = response.json()
            return data['rcsb_chem_comp_descriptor']['smiles']
        except Exception:
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

