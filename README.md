# PLIP Constraints Pipeline

Automatisierte Pipeline zur Generierung von Boltz-2 YAML Inputs mit PLIP-basierten Pocket-Constraints.

## 🎯 Übersicht

Die Pipeline führt 4 Schritte aus:
1. **Fetch**: PDB-Strukturen mit Affinitätsdaten (Ki/Kd/IC50) von RCSB laden
2. **Clean**: Strukturen bereinigen (Kristallisationshilfen entfernen)
3. **PLIP**: Liganden-Interaktionen analysieren
4. **Generate**: Boltz-2 YAML-Dateien erstellen (mit/ohne Constraints)

## 🚀 Quick Start

### Installation

```bash
# Python Dependencies
pip install -r requirements.txt

# Docker für PLIP
docker pull pharmai/plip

# Setup validieren
python validate_setup.py
```

### Erste Verwendung

```bash
# 1. Config anpassen
vim config.yaml  # uniprot_id und limit einstellen

# 2. Pipeline ausführen
python pipeline.py

# Oder schrittweise
python pipeline.py --steps fetch,clean
python pipeline.py --steps plip,generate

# Statistiken
python pipeline.py --stats
```

## 📁 Struktur

```
plip_constraints_pipeline/
├── pipeline.py           # Main CLI
├── validate_setup.py     # Setup-Validator
├── example_run.sh        # Beispiel-Script
├── config.yaml           # Konfiguration
├── requirements.txt
├── plip_pipeline/        # Core Module
│   ├── fetcher.py       # PDB Download
│   ├── cleaner.py       # PDB Cleaning
│   ├── plip.py          # PLIP Analysis
│   ├── generator.py     # YAML Generation
│   └── utils.py
└── output/              # Erstellt bei Run
    ├── raw_pdb/
    ├── clean_pdb/
    ├── plip_reports/
    ├── boltz_inputs/    # ← Hier sind die YAMLs
    └── metadata.csv
```

## ⚙️ Konfiguration (`config.yaml`)

```yaml
target:
  uniprot_id: "P00918"                    # Dein Protein (UniProt ID)
  ligand_whitelist: ["ZN"]                # Residuen behalten
  forbidden_ligands: ["GOL", "DMS", "SO4", "HOH"]  # Entfernen

fetch:
  limit: 50                               # Max. Strukturen
  allowed_types: ["Ki", "IC50"]

output:
  base_dir: "./output"

constraints:
  max_distance: 6.0                       # Pocket-Distanz (Å)
  ignore_metal_interactions: true         # Metall-Komplexe ignorieren

docker:
  plip_image: "pharmai/plip"
  timeout: 300
```

## 📊 Output

Pro Struktur werden **bis zu 3 YAML-Dateien** generiert:

**`{pdb_id}_default.yaml`** - Baseline ohne Constraints  
**`{pdb_id}_pocket.yaml`** - Mit Pocket-Constraints (wenn Kontakte vorhanden)  
**`{pdb_id}_contact.yaml`** - Mit Contact-Constraints (wenn Kontakte vorhanden)

### Constraint-Typen

**Pocket:** Binder zu Liste von Residuen
```yaml
constraints:
  - pocket:
      binder: "AZM"
      contacts: [["A", 199], ["A", 200], ["A", 131]]
      max_distance: 6.0
      force: false
```

**Contact:** Paarweise Kontakte zwischen Ligand und jedem Residue
```yaml
constraints:
  - contact:
      token1: ["AZM", 1]
      token2: ["A", 199]
      max_distance: 6.0
      force: false
  - contact:
      token1: ["AZM", 1]
      token2: ["A", 200]
      max_distance: 6.0
      force: false
```

`output/metadata.csv` enthält PDB ID, Ligand, SMILES, Affinitätswerte

## 🔧 CLI Optionen

```bash
python pipeline.py [OPTIONS]

  --config PATH    Config-Datei (default: config.yaml)
  --steps STEPS    Schritte: fetch,clean,plip,generate,all (default: all)
  --output DIR     Output-Verzeichnis überschreiben
  --stats          Nur Statistiken anzeigen
  --help           Hilfe
```

## 🔬 Verwendung mit Boltz-2

```bash
# YAMLs kopieren
cp output/boltz_inputs/*.yaml /path/to/boltz/inputs/

# Boltz-2 ausführen
boltz predict /path/to/boltz/inputs/ --out_dir results/
```

## 📝 Wichtiges

### Constraint-Typen: Pocket vs. Contact

**Pocket-Constraint:** Definiert eine Bindungstasche als Ganzes
- Ligand muss innerhalb von `max_distance` zu **mindestens einem** der Kontakt-Residuen sein
- Weniger restriktiv, flexibler

**Contact-Constraint:** Definiert einzelne paarweise Kontakte
- Ligand muss innerhalb von `max_distance` zu **jedem** Kontakt-Residue sein
- Restriktiver, spezifischer

→ Vergleiche beide Ansätze um zu sehen, welcher für dein System besser funktioniert!

### Metall-Komplexe werden ignoriert!

Die Pipeline verwendet **nur Liganden-Interaktionen** (H-Bonds, Hydrophobic, Pi-Stacking, Salt Bridges), **nicht** Metall-Komplexe.

**Warum?** Bei vielen Enzymen (z.B. CA2) binden Inhibitoren über Metall-Koordination. Die Metall-Residuen sind immer gleich → nicht informativ für liganden-spezifische Constraints.

### Troubleshooting

| Problem | Lösung |
|---------|--------|
| Docker Permission Denied | `sudo usermod -aG docker $USER` |
| Keine UniProt-Sequenz | Pipeline nutzt Placeholder, manuell anpassen |
| Keine Constraints gefunden | Normal bei reinen Metall-Chelatoren |
| PLIP Timeout | `docker.timeout` in config.yaml erhöhen |

## 📚 Referenzen

- **PLIP**: https://plip-tool.biotec.tu-dresden.de/
- **Boltz-2**: https://github.com/jwohlwend/boltz
- **RCSB PDB**: https://www.rcsb.org/

