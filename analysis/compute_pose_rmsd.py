#!/usr/bin/env python3
"""Compute ligand pose RMSD for DiffDock and Boltz outputs."""

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser
from rdkit import Chem
from rdkit.Chem import rdMolAlign


AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "ASX", "GLX", "SEC", "PYL",
}
IGNORE_RESNAMES = {"HOH", "WAT", "DOD", "SO4", "PO4", "CL", "NA", "K", "MG", "CA", "ZN"}
POSE_INDEX_RE = re.compile(r"(?:pose|model|rank|sample)[_-]?(\d+)", re.IGNORECASE)


@dataclass
class AtomData:
    atom_name: str
    element: str
    coord: np.ndarray


def normalize_element(atom_name: str, element: str) -> str:
    if element and element.strip():
        return element.strip().upper()
    letters = "".join(ch for ch in atom_name if ch.isalpha())
    return letters[:1].upper() if letters else "X"


def is_heavy(atom: AtomData) -> bool:
    return atom.element != "H"


def kabsch_rmsd(ref: np.ndarray, target: np.ndarray) -> float:
    ref_centered = ref - ref.mean(axis=0)
    target_centered = target - target.mean(axis=0)
    cov = target_centered.T @ ref_centered
    u, _, vt = np.linalg.svd(cov)
    d = np.linalg.det(vt.T @ u.T)
    corr = np.eye(3)
    corr[2, 2] = 1.0 if d >= 0 else -1.0
    rotation = vt.T @ corr @ u.T
    aligned = target_centered @ rotation
    diff = aligned - ref_centered
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


def atom_name_mapping(ref_atoms: Sequence[AtomData], target_atoms: Sequence[AtomData]) -> Tuple[np.ndarray, np.ndarray, str]:
    ref_by_name = {a.atom_name.upper(): a for a in ref_atoms}
    target_by_name = {a.atom_name.upper(): a for a in target_atoms}
    common = sorted(set(ref_by_name).intersection(target_by_name))
    if len(common) >= 3:
        ref_coords = np.array([ref_by_name[n].coord for n in common], dtype=float)
        target_coords = np.array([target_by_name[n].coord for n in common], dtype=float)
        return ref_coords, target_coords, "atom_name"

    if len(ref_atoms) != len(target_atoms):
        raise ValueError("No reliable mapping: atom names differ and atom count differs.")

    ref_by_element: Dict[str, List[AtomData]] = {}
    target_by_element: Dict[str, List[AtomData]] = {}
    for atom in ref_atoms:
        ref_by_element.setdefault(atom.element, []).append(atom)
    for atom in target_atoms:
        target_by_element.setdefault(atom.element, []).append(atom)

    if set(ref_by_element) != set(target_by_element):
        raise ValueError("No reliable mapping: element composition differs.")
    for element, ref_list in ref_by_element.items():
        if len(ref_list) != len(target_by_element[element]):
            raise ValueError(f"No reliable mapping: element count differs for {element}.")

    ref_coords_list: List[np.ndarray] = []
    target_coords_list: List[np.ndarray] = []
    for element in sorted(ref_by_element):
        ref_sorted = sorted(ref_by_element[element], key=lambda a: (a.coord[0], a.coord[1], a.coord[2]))
        target_sorted = sorted(target_by_element[element], key=lambda a: (a.coord[0], a.coord[1], a.coord[2]))
        ref_coords_list.extend(a.coord for a in ref_sorted)
        target_coords_list.extend(a.coord for a in target_sorted)

    return np.array(ref_coords_list, dtype=float), np.array(target_coords_list, dtype=float), "element_sorted_fallback"


def pick_best_residue(residue_atoms: Dict[str, List[AtomData]], preferred_resname: Optional[str]) -> Tuple[str, List[AtomData]]:
    if not residue_atoms:
        raise ValueError("No ligand-like residues found.")

    if preferred_resname:
        preferred = preferred_resname.strip().upper()
        candidates = [(res_id, atoms) for res_id, atoms in residue_atoms.items() if res_id.startswith(f"{preferred}|")]
        if candidates:
            return max(candidates, key=lambda item: len(item[1]))

    return max(residue_atoms.items(), key=lambda item: len(item[1]))


def parse_pdb_ligand_atoms(path: Path, ligand_resname: Optional[str]) -> Tuple[str, List[AtomData]]:
    residues: Dict[str, List[AtomData]] = {}
    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            record = line[:6].strip()
            if record not in {"HETATM", "ATOM"}:
                continue
            resname = line[17:20].strip().upper()
            chain = line[21:22].strip()
            resseq = line[22:26].strip()
            if record == "ATOM" and resname in AMINO_ACIDS:
                continue
            if resname in IGNORE_RESNAMES:
                continue
            atom_name = line[12:16].strip()
            element = normalize_element(atom_name, line[76:78].strip())
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            key = f"{resname}|{chain}|{resseq}"
            residues.setdefault(key, []).append(AtomData(atom_name=atom_name, element=element, coord=np.array([x, y, z])))
    return pick_best_residue(residues, ligand_resname)


def parse_pdb_ligand_mol(path: Path, ligand_resname: Optional[str]) -> Tuple[str, Chem.Mol]:
    residues: Dict[str, Dict[str, object]] = {}
    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            record = line[:6].strip()
            if record not in {"HETATM", "ATOM"}:
                continue
            resname = line[17:20].strip().upper()
            chain = line[21:22].strip()
            resseq = line[22:26].strip()
            if record == "ATOM" and resname in AMINO_ACIDS:
                continue
            if resname in IGNORE_RESNAMES:
                continue
            key = f"{resname}|{chain}|{resseq}"
            serial = int(line[6:11])
            residues.setdefault(key, {"atoms": [], "serials": set()})
            residues[key]["atoms"].append(line.rstrip("\n"))
            residues[key]["serials"].add(serial)

    if not residues:
        raise ValueError("No ligand-like residues found in PDB.")

    preferred = ligand_resname.strip().upper() if ligand_resname else None
    candidate_keys = list(residues.keys())
    if preferred:
        preferred_keys = [k for k in candidate_keys if k.startswith(f"{preferred}|")]
        if preferred_keys:
            candidate_keys = preferred_keys
    selected_key = max(candidate_keys, key=lambda k: len(residues[k]["atoms"]))  # type: ignore[index]
    atom_lines = residues[selected_key]["atoms"]  # type: ignore[index]
    serials = residues[selected_key]["serials"]  # type: ignore[index]

    conect_lines: List[str] = []
    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("CONECT"):
                continue
            tokens = line.split()
            if len(tokens) < 3:
                continue
            try:
                source = int(tokens[1])
                targets = [int(t) for t in tokens[2:]]
            except ValueError:
                continue
            if source in serials and any(t in serials for t in targets):
                conect_lines.append(line.rstrip("\n"))

    pdb_block = "\n".join(atom_lines + conect_lines + ["END", ""])
    mol = Chem.MolFromPDBBlock(pdb_block, removeHs=False, sanitize=True, proximityBonding=True)
    if mol is None:
        mol = Chem.MolFromPDBBlock(pdb_block, removeHs=False, sanitize=False, proximityBonding=True)
    if mol is None:
        raise ValueError("Could not build RDKit molecule from PDB ligand.")
    return selected_key, mol


def parse_cif_ligand_atoms(path: Path, ligand_resname: Optional[str]) -> Tuple[str, List[AtomData]]:
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(path.stem, str(path))
    residues: Dict[str, List[AtomData]] = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname().strip().upper()
                if resname in AMINO_ACIDS or resname in IGNORE_RESNAMES:
                    continue
                hetflag = residue.id[0].strip()
                if not hetflag and resname in AMINO_ACIDS:
                    continue
                resseq = str(residue.id[1]).strip()
                key = f"{resname}|{chain.id}|{resseq}"
                atom_list: List[AtomData] = []
                for atom in residue:
                    atom_name = atom.get_name().strip()
                    element = normalize_element(atom_name, getattr(atom, "element", ""))
                    coord = atom.coord.astype(float)
                    atom_list.append(AtomData(atom_name=atom_name, element=element, coord=coord))
                if atom_list:
                    residues[key] = atom_list
        break
    return pick_best_residue(residues, ligand_resname)


def parse_cif_ligand_mol(path: Path, ligand_resname: Optional[str]) -> Tuple[str, Chem.Mol]:
    """Build an RDKit ligand molecule from a CIF residue without modifying input files."""
    residue_id, atoms = parse_cif_ligand_atoms(path, ligand_resname)
    resname = residue_id.split("|")[0][:3] if residue_id else "LIG"

    pdb_lines: List[str] = []
    for idx, atom in enumerate(atoms, start=1):
        atom_name = atom.atom_name.strip()[:4]
        element = atom.element.strip()[:2].upper() or "C"
        if len(atom_name) < 4 and len(element) == 1:
            atom_field = f" {atom_name:<3}"
        else:
            atom_field = f"{atom_name:>4}"
        x, y, z = atom.coord
        line = (
            f"HETATM{idx:5d} {atom_field} {resname:>3s} L{1:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{20.00:6.2f}          {element:>2s}"
        )
        pdb_lines.append(line)

    pdb_block = "\n".join(pdb_lines + ["END", ""])

    mol = Chem.MolFromPDBBlock(pdb_block, removeHs=False, sanitize=True, proximityBonding=True)
    if mol is None:
        mol = Chem.MolFromPDBBlock(pdb_block, removeHs=False, sanitize=False, proximityBonding=True)
    if mol is None:
        raise ValueError("Could not build RDKit molecule from CIF ligand.")
    return str(residue_id), mol


def parse_sdf_ligand_atoms(path: Path) -> Tuple[str, List[AtomData]]:
    supplier = Chem.SDMolSupplier(str(path), removeHs=False)
    mol = next((m for m in supplier if m is not None), None)
    if mol is None:
        raise ValueError("Could not parse SDF molecule.")
    conf = mol.GetConformer()
    atoms: List[AtomData] = []
    for idx, atom in enumerate(mol.GetAtoms()):
        pos = conf.GetAtomPosition(idx)
        element = atom.GetSymbol().upper()
        atom_name = f"{element}{idx + 1}"
        atoms.append(AtomData(atom_name=atom_name, element=element, coord=np.array([pos.x, pos.y, pos.z], dtype=float)))
    return "SDF|LIG|1", atoms


def parse_sdf_ligand_mol(path: Path) -> Chem.Mol:
    supplier = Chem.SDMolSupplier(str(path), removeHs=False)
    mol = next((m for m in supplier if m is not None), None)
    if mol is None:
        raise ValueError("Could not parse SDF molecule.")
    return mol


def load_ligand_atoms(path: Path, ligand_resname: Optional[str]) -> Tuple[str, List[AtomData]]:
    suffix = path.suffix.lower()
    if suffix == ".pdb":
        return parse_pdb_ligand_atoms(path, ligand_resname)
    if suffix == ".cif":
        return parse_cif_ligand_atoms(path, ligand_resname)
    if suffix == ".sdf":
        return parse_sdf_ligand_atoms(path)
    raise ValueError(f"Unsupported file format: {path.suffix}")


def best_rmsd_rdkit(ref_mol: Chem.Mol, pose_mol: Chem.Mol, include_hydrogens: bool,
                    max_matches: int = 1000000) -> float:
    """Symmetry-corrected RMSD.

    GetBestRMS enumerates the substructure matches between the two molecules,
    which for symmetric ligands grows fast enough to look like a hang. RDKit's
    own limit is 1e6 and is kept as the default so results do not shift;
    --max-matches lowers it when a structure stalls.
    """
    ref = Chem.Mol(ref_mol)
    pose = Chem.Mol(pose_mol)
    if not include_hydrogens:
        ref = Chem.RemoveHs(ref)
        pose = Chem.RemoveHs(pose)
    try:
        return float(rdMolAlign.GetBestRMS(pose, ref, maxMatches=max_matches))
    except TypeError:
        # Older RDKit without the maxMatches keyword.
        return float(rdMolAlign.GetBestRMS(pose, ref))


def build_diffdock_candidates(diffdock_poses_dir: Path, diffdock_complex_dir: Path, pdb_id: str, max_poses: int) -> List[Tuple[str, str, Path]]:
    candidates: List[Tuple[str, str, Path]] = []

    pose_dir = diffdock_poses_dir / pdb_id / "complex_0"
    if pose_dir.is_dir():
        pose_files = sorted(pose_dir.glob("rank*.sdf"), key=lambda p: int(POSE_INDEX_RE.search(p.stem).group(1)) if POSE_INDEX_RE.search(p.stem) else 9999)
        for pose_file in pose_files[:max_poses]:
            pose_name = pose_file.stem
            candidates.append(("diffdock_sdf", "diffdock", pose_file))

    complex_file = diffdock_complex_dir / f"{pdb_id}_diffdock_complex.pdb"
    if complex_file.exists():
        candidates.append(("diffdock_complex_pdb", "diffdock", complex_file))

    return candidates


def infer_scenario(path: Path, pdb_id: str, scenarios: Sequence[str]) -> Optional[str]:
    lower = str(path).lower()
    pdb_lower = pdb_id.lower()
    for scenario in scenarios:
        token = f"{pdb_lower}_{scenario.lower()}"
        if token in lower or scenario.lower() in lower:
            return scenario
    return None


def build_boltz_candidates(boltz_dir: Path, pdb_id: str, scenarios: Sequence[str]) -> List[Tuple[str, str, Path]]:
    if not boltz_dir.exists():
        return []
    candidates: List[Tuple[str, str, Path]] = []
    for cif_path in sorted(boltz_dir.rglob("*.cif")):
        lower = str(cif_path).lower()
        if pdb_id.lower() not in lower:
            continue
        scenario = infer_scenario(cif_path, pdb_id, scenarios)
        if not scenario:
            continue
        candidates.append(("boltz_cif", scenario, cif_path))
    return candidates


def pose_label(path: Path) -> str:
    match = POSE_INDEX_RE.search(path.stem)
    if match:
        return f"{path.stem}|idx_{match.group(1)}"
    return path.stem


def resolve_reference_path(reference_dir: Path, pdb_id: str) -> Path:
    """Resolve reference file for a PDB ID in raw or clean folders."""
    candidates = [
        reference_dir / f"{pdb_id}.pdb",
        reference_dir / f"{pdb_id}_clean.pdb",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute ligand pose RMSD for DiffDock and Boltz outputs")
    parser.add_argument("--metadata", default="./output/metadata.csv", help="Path to metadata.csv")
    parser.add_argument("--raw-pdb-dir", default="./output/raw_pdb", help="Reference crystal PDB directory")
    parser.add_argument("--diffdock-poses-dir", default="./output/diffdock_workflow/poses", help="DiffDock poses directory")
    parser.add_argument("--diffdock-complex-dir", default="./output/diffdock_workflow/complexes", help="DiffDock complex PDB directory")
    parser.add_argument("--boltz-dir", default="", help="Boltz results root directory with CIF files (optional)")
    parser.add_argument("--max-matches", type=int, default=1000000,
                        help="Cap on substructure matches in the symmetry-corrected "
                             "RMSD. RDKit's own default; lower it (e.g. 10000) if a "
                             "symmetric ligand makes a structure stall")
    parser.add_argument("--max-diffdock-poses", type=int, default=5, help="Max rank*.sdf poses to evaluate per complex")
    parser.add_argument("--scenarios", default="default,crystal_pocket,diffdock_pocket", help="Comma-separated Boltz scenarios")
    parser.add_argument("--include-hydrogens", action="store_true", help="Include hydrogens in RMSD")
    parser.add_argument("--output", default="./results/analysis/pose_rmsd_summary.csv", help="Output CSV path")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    metadata_path = Path(args.metadata)
    raw_pdb_dir = Path(args.raw_pdb_dir)
    diffdock_poses_dir = Path(args.diffdock_poses_dir)
    diffdock_complex_dir = Path(args.diffdock_complex_dir)
    boltz_dir = Path(args.boltz_dir) if args.boltz_dir else None
    scenarios = [s.strip() for s in args.scenarios.split(",") if s.strip()]

    if not metadata_path.exists():
        raise FileNotFoundError(f"Metadata not found: {metadata_path}")

    df = pd.read_csv(metadata_path)
    required_cols = {"pdb_id", "ligand"}
    missing = required_cols.difference(df.columns)
    if missing:
        raise ValueError(f"metadata.csv missing required columns: {sorted(missing)}")

    records: List[Dict[str, object]] = []

    total = len(df)
    for position, (_, row) in enumerate(df.iterrows(), start=1):
        pdb_id = str(row["pdb_id"]).strip()
        ligand = str(row["ligand"]).strip()
        if not pdb_id or not ligand:
            continue

        # Symmetry-corrected RMSD can take a long time on a single ligand, so
        # say which one is being worked on rather than going quiet.
        print(f"  [{position}/{total}] {pdb_id} ({ligand})", flush=True)

        ref_path = resolve_reference_path(raw_pdb_dir, pdb_id)
        if not ref_path.exists():
            records.append(
                {
                    "pdb_id": pdb_id,
                    "ligand": ligand,
                    "source": "reference",
                    "scenario": "reference",
                    "pose_name": "",
                    "pose_path": str(ref_path),
                    "status": "missing_reference",
                    "method": "",
                    "rmsd_angstrom": np.nan,
                    "n_atoms_ref": np.nan,
                    "n_atoms_pose": np.nan,
                    "n_atoms_mapped": np.nan,
                    "ref_residue": "",
                    "pose_residue": "",
                    "error": "reference pdb missing",
                }
            )
            continue

        try:
            ref_residue, ref_atoms = load_ligand_atoms(ref_path, ligand)
            _, ref_mol = parse_pdb_ligand_mol(ref_path, ligand)
        except Exception as exc:
            records.append(
                {
                    "pdb_id": pdb_id,
                    "ligand": ligand,
                    "source": "reference",
                    "scenario": "reference",
                    "pose_name": "",
                    "pose_path": str(ref_path),
                    "status": "reference_parse_failed",
                    "method": "",
                    "rmsd_angstrom": np.nan,
                    "n_atoms_ref": np.nan,
                    "n_atoms_pose": np.nan,
                    "n_atoms_mapped": np.nan,
                    "ref_residue": "",
                    "pose_residue": "",
                    "error": str(exc),
                }
            )
            continue

        if not args.include_hydrogens:
            ref_atoms = [a for a in ref_atoms if is_heavy(a)]

        candidates = build_diffdock_candidates(diffdock_poses_dir, diffdock_complex_dir, pdb_id, args.max_diffdock_poses)
        if boltz_dir is not None:
            candidates.extend(build_boltz_candidates(boltz_dir, pdb_id, scenarios))

        if not candidates:
            records.append(
                {
                    "pdb_id": pdb_id,
                    "ligand": ligand,
                    "source": "none",
                    "scenario": "none",
                    "pose_name": "",
                    "pose_path": "",
                    "status": "no_candidate_poses",
                    "method": "",
                    "rmsd_angstrom": np.nan,
                    "n_atoms_ref": len(ref_atoms),
                    "n_atoms_pose": np.nan,
                    "n_atoms_mapped": np.nan,
                    "ref_residue": ref_residue,
                    "pose_residue": "",
                    "error": "",
                }
            )
            continue

        for source, scenario, pose_path in candidates:
            pose_name = pose_label(pose_path)
            try:
                pose_residue, pose_atoms = load_ligand_atoms(pose_path, ligand if pose_path.suffix.lower() != ".sdf" else None)
                if not args.include_hydrogens:
                    pose_atoms = [a for a in pose_atoms if is_heavy(a)]
                rmsd = np.nan
                method = ""
                n_mapped = np.nan
                result_status = "ok"
                mapping_coverage = np.nan

                rdkit_error = None
                if pose_path.suffix.lower() in {".sdf", ".pdb", ".cif"}:
                    try:
                        ref_mol_for_rdkit = ref_mol
                        if pose_path.suffix.lower() == ".sdf":
                            pose_mol = parse_sdf_ligand_mol(pose_path)
                        elif pose_path.suffix.lower() == ".cif":
                            _, pose_mol = parse_cif_ligand_mol(pose_path, ligand)
                        else:
                            _, pose_mol = parse_pdb_ligand_mol(pose_path, ligand)
                        rmsd = best_rmsd_rdkit(ref_mol_for_rdkit, pose_mol,
                                               args.include_hydrogens, args.max_matches)
                        method = "rdkit_best_rms"
                        n_mapped = float(
                            Chem.Mol(ref_mol_for_rdkit).GetNumAtoms()
                            if args.include_hydrogens
                            else Chem.RemoveHs(Chem.Mol(ref_mol_for_rdkit)).GetNumAtoms()
                        )
                    except Exception as exc:
                        rdkit_error = str(exc)

                if not method:
                    ref_coords, pose_coords, method = atom_name_mapping(ref_atoms, pose_atoms)
                    rmsd = kabsch_rmsd(ref_coords, pose_coords)
                    n_mapped = float(len(ref_coords))
                    min_atoms = min(len(ref_atoms), len(pose_atoms))
                    mapping_coverage = (len(ref_coords) / min_atoms) if min_atoms > 0 else np.nan
                    if min_atoms > 0 and mapping_coverage < 0.7:
                        result_status = "ok_low_coverage"
                else:
                    mapping_coverage = 1.0

                records.append(
                    {
                        "pdb_id": pdb_id,
                        "ligand": ligand,
                        "source": source,
                        "scenario": scenario,
                        "pose_name": pose_name,
                        "pose_path": str(pose_path),
                        "status": result_status,
                        "method": method,
                        "rmsd_angstrom": rmsd,
                        "n_atoms_ref": len(ref_atoms),
                        "n_atoms_pose": len(pose_atoms),
                        "n_atoms_mapped": n_mapped,
                        "mapping_coverage": mapping_coverage,
                        "ref_residue": ref_residue,
                        "pose_residue": pose_residue,
                        "error": rdkit_error if rdkit_error and method != "rdkit_best_rms" else "",
                    }
                )
            except Exception as exc:
                records.append(
                    {
                        "pdb_id": pdb_id,
                        "ligand": ligand,
                        "source": source,
                        "scenario": scenario,
                        "pose_name": pose_name,
                        "pose_path": str(pose_path),
                        "status": "failed",
                        "method": "",
                        "rmsd_angstrom": np.nan,
                        "n_atoms_ref": len(ref_atoms),
                        "n_atoms_pose": np.nan,
                        "n_atoms_mapped": np.nan,
                        "mapping_coverage": np.nan,
                        "ref_residue": ref_residue,
                        "pose_residue": "",
                        "error": str(exc),
                    }
                )

    out_df = pd.DataFrame(records)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(output_path, index=False)

    ok_count = int((out_df["status"] == "ok").sum()) if not out_df.empty else 0
    low_cov_count = int((out_df["status"] == "ok_low_coverage").sum()) if not out_df.empty else 0
    fail_count = int((out_df["status"] == "failed").sum()) if not out_df.empty else 0
    print(f"Saved RMSD results: {output_path}")
    print(f"Rows: {len(out_df)} | OK: {ok_count} | OK(low coverage): {low_cov_count} | Failed: {fail_count}")


if __name__ == "__main__":
    main()
