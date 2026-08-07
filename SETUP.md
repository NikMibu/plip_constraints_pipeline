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

**Pin OpenBabel, not just PLIP.** PLIP does not place hydrogens — OpenBabel
does, and PLIP derives the hydrogen bonds from them. A different OpenBabel
therefore changes individual constraints. Measured on PDB 3KWA: residue His94
is contacted under OpenBabel 3.1.1 and not under 3.2.1, with PLIP 2.3.1 and
3.0.1 agreeing once the protonation is fixed. The reports behind the published
results were produced with PLIP 2.3.1.

```bash
micromamba create -n plip -c conda-forge python=3.9 openbabel=3.1.1
micromamba run -n plip pip install plip==2.3.1
```

Check it:

```bash
micromamba run -n plip plip -f <structure>.pdb -x -o /tmp/plip_check
```

**Without micromamba.** PLIP and OpenBabel both publish wheels, so a plain
`pip install plip==2.3.1 openbabel==3.1.1` also works. The pipeline falls back
to a `plip` on `PATH` when micromamba is absent, and says so loudly, because
that fallback bypasses the pinned environment. Micromamba remains the tested
route.

> **`which micromamba` finds nothing but `micromamba run` works?**
> That is the normal install. The shell hook defines `micromamba` as a shell
> *function* and leaves the binary off `PATH`, so it works when you type it but
> is invisible to Python's `subprocess`. The hook exports the real path as
> `$MAMBA_EXE`, which the pipeline uses. If your shell does not set it:
>
> ```bash
> export MAMBA_EXE=/path/to/micromamba          # find it with: type micromamba
> ```
>
> or put the absolute path into `micromamba.executable` in `config.yaml`.
> `python validate_setup.py` tells you which one it resolved to.

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

## 4. Boltz-2

Only needed to run the generated YAMLs — the pipeline does not call Boltz
itself. A plain venv is the tested route:

```bash
python3 -m venv ~/boltz-venv
~/boltz-venv/bin/pip install boltz==2.2.1
```

Then point the config at it, so `validate_setup.py` can find it:

```yaml
boltz:
  executable: "~/boltz-venv/bin/boltz"
```

`$BOLTZ_EXE` overrides the config; a micromamba environment named in
`micromamba.boltz_env` and a `boltz` on `PATH` are also tried, in that order.

Model weights and the CCD (~2–3 GB) download into `~/.boltz` on first use.

> **On WSL2, pass `--no_kernels` to `boltz predict`.** The cuEquivariance
> Triton kernels crash there because WSL's NVML does not implement
> `nvmlDeviceGetNumGpuCores`. The flag disables them at some cost in speed.
> Native Linux and cloud GPUs do not need it.

## 5. Configuration

The relevant section of `config.yaml`:

```yaml
micromamba:
  executable: "micromamba"   # or an absolute path
  plip_env: "plip"
  diffdock_env: "diffdock"

diffdock:
  repo_path: "${DIFFDOCK_HOME}"
  samples_per_complex: 40

boltz:
  executable: null           # or a path, e.g. "~/boltz-venv/bin/boltz"

fetch:
  pdb_ids_file: null         # or "data/pdb_ids.txt" to reproduce the benchmark
```

## 6. Running

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
| `micromamba not found` although it works in your shell | It is a shell function; the binary is off `PATH`. Export `$MAMBA_EXE` or set `micromamba.executable`. See section 2. |
| Loud `falling back to .../plip` banner | micromamba was not resolved, so PLIP is taken from `PATH` and bypasses the pinned environment. Constraints can differ. Fix the micromamba resolution before using the results. |
| `boltz not found` in the validator | Set `boltz.executable` to the binary — a venv is fine — or `$BOLTZ_EXE`. See section 4. |
| `boltz predict` crashes on WSL2 with an NVML error | Add `--no_kernels`. See section 4. |
| `DiffDock repository not found` | `DIFFDOCK_HOME` unset and `diffdock.repo_path` unusable. See section 3. |
| `Could not determine the protein sequence` | UniProt unreachable and no usable MSA configured. Set `target.msa_path` to an a3m file — its first record is the target sequence — or restore network access. The pipeline refuses to emit YAMLs without a real sequence. |
| DiffDock appears to hang on first use | SO(3) lookup tables are being built. Wait 5–10 minutes. |
| No constraints found for a structure | Expected for pure metal chelators. Metal coordination is excluded on purpose (`constraints.ignore_metal_interactions`), because the coordinating residues are identical for every ligand and carry no ligand-specific information. |
