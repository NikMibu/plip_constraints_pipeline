# Setup

The pipeline drives two external tools that cannot share one Python
environment: **PLIP** needs OpenBabel and Python 3.9, **DiffDock** needs
PyTorch 1.13.1 and its own dependency set. Each lives in its own micromamba
environment; the pipeline calls into them with `micromamba run -n <env>`.

The pipeline itself runs in whatever Python you invoke it with — it only needs
the packages in [`requirements.txt`](requirements.txt).

Verify a setup at any point with:

```bash
python validate_setup.py            # add --skip-network if you are offline
```

It reads `config.yaml`, so it checks the environments the pipeline will
actually use.

---

## 1. Pipeline environment

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## 2. PLIP

PLIP's output determines the constraints, so the version matters. The results
in the thesis were produced with **PLIP 3.0.0**; version 3.0.1 was observed to
report one interaction fewer for at least one structure.

```bash
micromamba create -n plip -c conda-forge python=3.9 openbabel=3.1.1
micromamba run -n plip pip install plip==3.0.0
```

Check it:

```bash
micromamba run -n plip plip -f <structure>.pdb -x -o /tmp/plip_check
```

**Without micromamba.** PLIP and OpenBabel both publish wheels, so a plain
`pip install plip==3.0.0 openbabel==3.1.1` also works. The pipeline falls back
to a `plip` on `PATH` when micromamba is absent and prints which one it uses.
Micromamba remains the tested route.

## 3. DiffDock

Only needed for the DiffDock workflow. The crystal workflow runs without it.

```bash
git clone https://github.com/gcorso/DiffDock.git
cd DiffDock
git checkout v1.1                  # the version used for the thesis

micromamba create -n diffdock -c conda-forge python=3.9.18
micromamba run -n diffdock pip install torch==1.13.1
micromamba run -n diffdock pip install -r requirements.txt
```

Then tell the pipeline where it is:

```bash
export DIFFDOCK_HOME=/path/to/DiffDock
```

`DIFFDOCK_HOME` takes precedence over `diffdock.repo_path` in `config.yaml`,
which defaults to `${DIFFDOCK_HOME}`. Nothing in the repository hard-codes a
machine-specific path.

> **First run takes longer.** DiffDock builds its SO(3) lookup tables on first
> use, which takes roughly 5–10 minutes and looks like a hang. Let it finish;
> subsequent runs reuse the tables.

Check it:

```bash
cd "$DIFFDOCK_HOME"
micromamba run -n diffdock python -m inference \
    --protein_path protein.pdb --ligand ligand.sdf --out_dir /tmp/dd_check
```

## 4. Configuration

The relevant section of `config.yaml`:

```yaml
micromamba:
  executable: "micromamba"   # or an absolute path
  plip_env: "plip"
  diffdock_env: "diffdock"

diffdock:
  repo_path: "${DIFFDOCK_HOME}"
  samples_per_complex: 40
```

## 5. Running

```bash
python3 pipeline.py --steps crystal          # fetch, clean, PLIP, YAMLs   (CPU)
python3 pipeline.py --steps diffdock_full    # ligand prep, DiffDock, PLIP, YAMLs (GPU)
python3 pipeline.py --steps all              # both
python3 pipeline.py --stats                  # summary of an existing run
```

Start with `--steps crystal`. It needs no GPU and no DiffDock, produces
visible output, and exercises everything except the docking itself.

## Troubleshooting

| Symptom | Cause and fix |
|---|---|
| `PLIP not found` | Neither micromamba nor a `plip` on `PATH`. See section 2. |
| `DiffDock repository not found` | `DIFFDOCK_HOME` unset and `diffdock.repo_path` unusable. See section 3. |
| `Could not determine the protein sequence` | UniProt unreachable and no usable MSA configured. Set `target.msa_path` to an a3m file — its first record is the target sequence — or restore network access. The pipeline refuses to emit YAMLs without a real sequence. |
| DiffDock appears to hang on first use | SO(3) lookup tables are being built. Wait 5–10 minutes. |
| No constraints found for a structure | Expected for pure metal chelators. Metal coordination is excluded on purpose (`constraints.ignore_metal_interactions`), because the coordinating residues are identical for every ligand and carry no ligand-specific information. |
