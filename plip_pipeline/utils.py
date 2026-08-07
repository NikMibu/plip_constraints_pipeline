"""Utility Functions"""
import os
import requests
from typing import Optional


def ensure_dir(path: str) -> None:
    """Create directory if it doesn't exist."""
    os.makedirs(path, exist_ok=True)


def yaml_quote(s: str) -> str:
    """YAML-safe double-quote escaping."""
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def get_uniprot_sequence(uniprot_id: str) -> Optional[str]:
    """
    Fetch protein sequence from UniProt.
    
    Args:
        uniprot_id: UniProt accession (e.g., "P00918")
        
    Returns:
        Protein sequence or None on error
    """
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
        response = requests.get(url, timeout=10)
        
        if response.status_code == 200:
            lines = response.text.strip().split('\n')
            # Skip header line, join sequence
            sequence = ''.join(line.strip() for line in lines[1:] if line.strip())
            return sequence if sequence else None
        else:
            print(f"[WARNING] Could not fetch sequence for {uniprot_id}: HTTP {response.status_code}")
            return None
            
    except Exception as e:
        print(f"[ERROR] Failed to fetch UniProt sequence: {e}")
        return None


def read_msa_query_sequence(msa_path: str) -> Optional[str]:
    """Read the query sequence from an a3m/FASTA multiple sequence alignment.

    The first record of an a3m file is the query itself, so a configured MSA
    already contains the target sequence. Used as an offline fallback when
    UniProt cannot be reached.

    Lowercase letters mark insertions relative to the query and dashes mark
    gaps; both are stripped so the result is the plain sequence.
    """
    if not msa_path or not os.path.exists(msa_path):
        return None

    try:
        with open(msa_path) as f:
            sequence_lines = []
            seen_header = False
            for line in f:
                line = line.strip()
                if not line:
                    continue
                # hhsuite writes a "#<ncol>\t<nseq>" line before the first
                # record. It is not a record header and must not be counted.
                if line.startswith("#"):
                    continue
                if line.startswith(">"):
                    if seen_header:
                        break  # second record reached, query is complete
                    seen_header = True
                    continue
                sequence_lines.append(line)

        sequence = "".join(sequence_lines)
        sequence = "".join(c for c in sequence if c.isupper())
        return sequence or None

    except Exception as e:
        print(f"[WARNING] Could not read MSA {msa_path}: {e}")
        return None


def format_affinity(value: float, unit: str = "nM") -> str:
    """
    Format affinity value with unit.
    
    Args:
        value: Affinity value
        unit: Unit (nM, uM, mM)
        
    Returns:
        Formatted string
    """
    if value < 1:
        return f"{value:.3f} {unit}"
    elif value < 100:
        return f"{value:.1f} {unit}"
    else:
        return f"{value:.0f} {unit}"

