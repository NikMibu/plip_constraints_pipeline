"""DiffDock Runner - Molecular docking via micromamba"""
import os
import shutil
import subprocess
from typing import Dict, Optional, Tuple
import pandas as pd
from Bio import PDB
from .utils import ensure_dir


class DiffDockRunner:
    """Run DiffDock docking via micromamba and manage outputs."""
    
    def __init__(self, config: Dict):
        """
        Initialize DiffDock Runner.
        
        Args:
            config: Configuration dictionary
        """
        base_dir = config['output']['base_dir']
        self.diffdock_dir = os.path.join(base_dir, "diffdock_workflow")
        
        self.protein_dir = os.path.join(self.diffdock_dir, "proteins_apo")
        self.poses_dir = os.path.join(self.diffdock_dir, "poses")
        self.complex_dir = os.path.join(self.diffdock_dir, "complexes")
        
        ensure_dir(self.protein_dir)
        ensure_dir(self.poses_dir)
        ensure_dir(self.complex_dir)
        
        self.clean_pdb_dir = os.path.join(base_dir, "clean_pdb")
        
        # DiffDock settings
        self.diffdock_repo = config.get('diffdock', {}).get('repo_path', '/mnt/e/master_thesis/DiffDock')
        self.timeout = config['docker']['diffdock_timeout']
        self.samples = config['diffdock']['samples_per_complex']
        self.top_n = config['diffdock']['top_n_poses']
        self.micromamba_exe = config.get('micromamba', {}).get('executable', 'micromamba')
        self.micromamba_env = config.get('micromamba', {}).get('diffdock_env', 'diffdock')
        
    def run_all(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Run DiffDock for all protein-ligand pairs.
        
        Args:
            df: DataFrame with pdb_id, ligand, sdf_path columns
            
        Returns:
            DataFrame with added pose_path and complex_path columns
        """
        print("\n[DIFFDOCK] Running molecular docking...")
        
        successful = []
        pose_paths = []
        complex_paths = []
        
        for idx, row in df.iterrows():
            pdb_id = row['pdb_id']
            ligand_id = row['ligand']
            sdf_path = row.get('sdf_path')
            
            if not sdf_path or not os.path.exists(sdf_path):
                print(f"  [SKIP] {pdb_id}: No SDF file")
                pose_paths.append(None)
                complex_paths.append(None)
                continue
            
            try:
                # 1. Prepare apo protein (remove ligand)
                apo_path = self._prepare_apo_protein(pdb_id, ligand_id)
                if not apo_path:
                    pose_paths.append(None)
                    complex_paths.append(None)
                    continue
                
                # 2. Run DiffDock
                pose_path = self._run_docking(pdb_id, apo_path, sdf_path)
                if not pose_path:
                    pose_paths.append(None)
                    complex_paths.append(None)
                    continue
                
                # 3. Create complex (protein + pose)
                complex_path = self._create_complex(pdb_id, apo_path, pose_path)
                
                if complex_path:
                    print(f"  ✓ {pdb_id}: Docking completed")
                    successful.append(idx)
                    pose_paths.append(pose_path)
                    complex_paths.append(complex_path)
                else:
                    pose_paths.append(None)
                    complex_paths.append(None)
                    
            except Exception as e:
                print(f"  [ERROR] {pdb_id}: {e}")
                pose_paths.append(None)
                complex_paths.append(None)
        
        df['pose_path'] = pose_paths
        df['complex_path'] = complex_paths
        
        print(f"\n[DIFFDOCK] Successfully docked {len(successful)}/{len(df)} ligands")
        
        # Return full dataframe (with new columns added, even if None for failures)
        return df
    
    def _prepare_apo_protein(self, pdb_id: str, ligand_id: str) -> Optional[str]:
        """
        Remove ligand from clean PDB to create apo structure.
        
        Args:
            pdb_id: PDB ID
            ligand_id: Ligand to remove
            
        Returns:
            Path to apo PDB or None on error
        """
        clean_pdb = os.path.join(self.clean_pdb_dir, f"{pdb_id}_clean.pdb")
        apo_pdb = os.path.join(self.protein_dir, f"{pdb_id}_apo.pdb")
        
        if not os.path.exists(clean_pdb):
            print(f"  [ERROR] {pdb_id}: Clean PDB not found")
            return None
        
        try:
            parser = PDB.PDBParser(QUIET=True)
            structure = parser.get_structure(pdb_id, clean_pdb)
            
            # Remove ligand residues
            for model in structure:
                for chain in model:
                    residues_to_remove = []
                    for residue in chain:
                        if residue.get_resname().strip() == ligand_id:
                            residues_to_remove.append(residue.id)
                    
                    for res_id in residues_to_remove:
                        chain.detach_child(res_id)
            
            # Save apo structure
            io = PDB.PDBIO()
            io.set_structure(structure)
            io.save(apo_pdb)
            
            return apo_pdb
            
        except Exception as e:
            print(f"  [ERROR] {pdb_id}: Failed to create apo - {e}")
            return None
    
    def _run_docking(self, pdb_id: str, protein_path: str, ligand_sdf: str) -> Optional[str]:
        """
        Run DiffDock via micromamba.
        
        Args:
            pdb_id: PDB ID
            protein_path: Path to apo protein PDB
            ligand_sdf: Path to ligand SDF
            
        Returns:
            Path to best pose or None on error
        """
        output_dir = os.path.join(self.poses_dir, pdb_id)
        ensure_dir(output_dir)
        
        # Absolute paths
        protein_abs = os.path.abspath(protein_path)
        ligand_abs = os.path.abspath(ligand_sdf)
        output_abs = os.path.abspath(output_dir)
        
        # DiffDock command via micromamba
        cmd = [
            self.micromamba_exe, "run", "-n", self.micromamba_env,
            "python", "-m", "inference",
            "--protein_path", protein_abs,
            "--ligand", ligand_abs,
            "--out_dir", output_abs,
            "--inference_steps", "20",
            "--samples_per_complex", str(self.samples),
            "--batch_size", "10",
            "--actual_steps", "18"
        ]
        
        try:
            subprocess.run(
                cmd,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=self.timeout,
                cwd=self.diffdock_repo  # Run from DiffDock repo directory
            )
            
            # Find best pose (rank1.sdf or similar)
            best_pose = self._find_best_pose(output_dir)
            return best_pose
            
        except subprocess.TimeoutExpired:
            print(f"  [ERROR] {pdb_id}: DiffDock timeout")
            return None
        except subprocess.CalledProcessError as e:
            stderr = e.stderr.decode() if e.stderr else "No error output"
            print(f"  [ERROR] {pdb_id}: DiffDock failed - {stderr[:100]}")
            return None
        except Exception as e:
            print(f"  [ERROR] {pdb_id}: {e}")
            return None
    
    def _find_best_pose(self, output_dir: str) -> Optional[str]:
        """
        Find the best pose from DiffDock output.
        
        DiffDock creates output in complex_0/ subdirectory with rank1.sdf, rank2.sdf, etc.
        
        Args:
            output_dir: DiffDock output directory
            
        Returns:
            Path to best pose SDF
        """
        # DiffDock creates complex_0/ subdirectory
        complex_dir = os.path.join(output_dir, "complex_0")
        
        if os.path.exists(complex_dir):
            # Look for rank1.sdf (best pose)
            rank1_path = os.path.join(complex_dir, "rank1.sdf")
            if os.path.exists(rank1_path):
                return rank1_path
            
            # Fallback: find rank1_confidence*.sdf
            try:
                sdf_files = [f for f in os.listdir(complex_dir) 
                           if f.startswith('rank1') and f.endswith('.sdf')]
                if sdf_files:
                    return os.path.join(complex_dir, sdf_files[0])
            except Exception:
                pass
        
        return None
    
    def _create_complex(self, pdb_id: str, protein_path: str, pose_sdf: str) -> Optional[str]:
        """
        Combine protein and docked pose into complex PDB.
        
        Args:
            pdb_id: PDB ID
            protein_path: Path to apo protein
            pose_sdf: Path to docked pose SDF
            
        Returns:
            Path to complex PDB
        """
        complex_path = os.path.join(self.complex_dir, f"{pdb_id}_diffdock_complex.pdb")
        
        try:
            # Use RDKit to convert SDF to PDB format
            from rdkit import Chem
            
            # Read protein PDB
            with open(protein_path, 'r') as f:
                protein_lines = f.readlines()
            
            # Read ligand from SDF and convert to PDB format
            suppl = Chem.SDMolSupplier(pose_sdf, removeHs=False)
            mol = next(suppl)
            
            if mol is None:
                print(f"  [ERROR] {pdb_id}: Could not read ligand from SDF")
                return None
            
            # Convert ligand to PDB block
            ligand_pdb = Chem.MolToPDBBlock(mol)
            
            # Combine protein and ligand
            with open(complex_path, 'w') as f:
                # Write protein (remove END if present)
                for line in protein_lines:
                    if not line.startswith('END'):
                        f.write(line)
                
                # Write ligand with chain L
                for line in ligand_pdb.split('\n'):
                    if line.startswith(('HETATM', 'ATOM')):
                        # Replace chain ID with 'L' for ligand
                        line = line[:21] + 'L' + line[22:]
                        f.write(line + '\n')
                
                # Write END
                f.write('END\n')
            
            # Also save the SDF pose separately for reference
            pose_copy = os.path.join(self.complex_dir, f"{pdb_id}_diffdock_complex_pose.sdf")
            shutil.copy(pose_sdf, pose_copy)
            
            return complex_path
            
        except Exception as e:
            print(f"  [ERROR] {pdb_id}: Failed to create complex - {e}")
            return None

