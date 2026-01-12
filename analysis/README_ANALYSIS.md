# Boltz Affinity Prediction Analysis

Analyse-Pipeline für Boltz-2 Affinitätsvorhersagen mit drei Constraint-Szenarien.

## Szenarien

1. **default**: Keine Constraints
2. **crystal_pocket**: Constraints aus Cocrystal-Struktur (PLIP-analysiert)
3. **diffdock_pocket**: Constraints aus DiffDock-Pose (PLIP-analysiert)

## Workflow

### 1. Predictions parsen

Parse Boltz-Outputs und merge mit experimentellen Ki-Werten aus `metadata.csv`:

```bash
python parse_boltz_predictions.py \
    --predictions-dir /path/to/boltz/results \
    --metadata ./output/metadata.csv \
    --output ./results/predictions_summary.csv
```

**Output**: `results/predictions_summary.csv` mit Spalten:
- `pdb_id`, `ligand`, `scenario`
- `exp_Ki_nM`, `exp_type`, `exp_unit`
- `boltz_prob_mean`, `boltz_pred_value_mean`

### 2. Analyse durchführen

Berechne Regression-Metriken und erstelle Plots:

```bash
python analyze_affinity_predictions.py \
    --predictions ./results/predictions_summary.csv \
    --output-dir ./results/analysis
```

**Outputs**:
- `results/analysis/affinity_predictions_summary.csv`: Erweiterte Summary mit Conversions
- `results/analysis/affinity_regression_metrics.csv`: R², RMSE, MAE, Pearson, Spearman
- `results/analysis/affinity_regression.png`: 3-Panel Regression-Plot
- `results/analysis/error_distribution.png`: Fehlerverteilung

## Metriken

### Conversions
- Experimentelle Ki (nM) → log₁₀(µM)
- Boltz Predictions: bereits in log₁₀(µM) Scale
- pIC50 kcal/mol: `(6.0 - log10_um) * 1.364`

### Regression Metrics
- **R²**: Bestimmtheitsmaß der linearen Regression
- **RMSE**: Root Mean Squared Error (Pred - Exp)
- **MAE**: Mean Absolute Error
- **Bias**: Systematischer Fehler (Mean Error)
- **Pearson r**: Lineare Korrelation
- **Spearman ρ**: Rank-Korrelation

## Beispiel: Vollständiger Run

```bash
# 1. Parse Predictions vom Server
python parse_boltz_predictions.py \
    --predictions-dir /mnt/server/boltz_results \
    --metadata ./output/metadata.csv \
    --output ./results/predictions_summary.csv

# 2. Analyse durchführen
python analyze_affinity_predictions.py \
    --predictions ./results/predictions_summary.csv \
    --output-dir ./results/analysis
```

## Erwartete Boltz-Output-Struktur

```
boltz_results/
├── predictions/
│   ├── 3DC3_default/
│   │   └── affinity_3DC3_default.json
│   ├── 3DC3_crystal_pocket/
│   │   └── affinity_3DC3_crystal_pocket.json
│   ├── 3DC3_diffdock_pocket/
│   │   └── affinity_3DC3_diffdock_pocket.json
│   └── ...
```

Oder nested:
```
boltz_results/
├── boltz_results_timestamp/
│   └── predictions/
│       └── ...
```

## Dependencies

- pandas
- numpy
- matplotlib

(bereits in `requirements.txt` enthalten)
