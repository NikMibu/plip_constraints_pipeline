# Example outputs

Outputs of the full benchmark run behind the numbers in the top-level README,
not of a smoke test. 200 human carbonic anhydrase II co-crystal structures,
three scenarios each, Boltz-2 2.2.1 at default inference parameters
(3 / 200 / 5 / 200).

## Affinity

| File | |
|---|---|
| `affinity_regression_metrics_common_200.csv` | R², slope and intercept per scenario |
| `affinity_regression_common_200.png` | predicted vs. experimental, one panel per scenario |
| `error_distribution_common_200.png` | residual distributions |
| `affinity_common_200_ids.txt` | the 200 PDB IDs |

```
default          R² 0.586
diffdock_pocket  R² 0.563
crystal_pocket   R² 0.531
```

## Pose quality

| File | |
|---|---|
| `pose_rmsd_group_summary.csv` | n, mean, median, quartiles per scenario |
| `pose_rmsd_best_per_pdb.csv` | best RMSD per structure |
| `pose_rmsd_best_violin_by_group.png` | distribution per scenario |
| `pose_rmsd_best_cdf.png` | cumulative distribution |
| `rmsd_common_193_ids.txt` | the 193 PDB IDs |

```
default          median 0.477 Å
crystal_pocket   median 0.491 Å
diffdock_pocket  median 0.538 Å
```

**Why 193 and not 200.** RMSD needs a comparable pose in every scenario. Seven
structures are missing one, so they are excluded to keep the comparison
paired — otherwise the scenarios would be scored on different structures.
Affinity has no such requirement and uses all 200.

## `yaml_example/`

The three Boltz-2 inputs generated for PDB 3DC3 (ligand acetazolamide), so the
output format is visible without running anything:

| File | |
|---|---|
| `3DC3_default.yaml` | baseline, no constraints |
| `3DC3_crystal_pocket.yaml` | pocket constraint from the crystal interactions |
| `3DC3_diffdock_pocket.yaml` | pocket constraint from the top DiffDock pose |

The three differ only in the `constraints:` block — same sequence, same ligand,
same MSA — which is what makes the scenarios comparable.

> These are verbatim from the original run, including the absolute `msa:` path
> of the compute server that produced them. The pipeline writes whatever
> `target.msa_path` resolves to; on your machine that will be your own path.
