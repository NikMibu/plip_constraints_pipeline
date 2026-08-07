#!/usr/bin/env bash
# Example run: the crystal workflow end to end.
#
# Fetches co-crystal structures for the target in config.yaml, cleans them,
# profiles the ligand interactions with PLIP and writes Boltz-2 YAML inputs
# with and without pocket constraints.
#
# CPU only. Running the resulting YAMLs through Boltz-2 is a separate step,
# see the end of this script.
set -euo pipefail

CONFIG="${1:-config.yaml}"

usage() {
    cat <<EOF
Usage: bash example_run.sh [CONFIG]

  CONFIG   pipeline config (default: config.yaml)

Runs fetch -> clean -> plip -> generate and prints statistics.
EOF
}

if [ "${1:-}" = "-h" ] || [ "${1:-}" = "--help" ]; then
    usage
    exit 0
fi

if [ ! -f "$CONFIG" ]; then
    echo "ERROR: config not found: $CONFIG" >&2
    exit 1
fi

echo "=================================="
echo "PLIP Constraints Pipeline"
echo "=================================="

# Validate first. This checks micromamba, the plip environment and the RCSB
# and UniProt endpoints - everything the steps below depend on.
echo
echo "[1/6] Validating setup..."
python3 validate_setup.py --config "$CONFIG"

echo
echo "[2/6] Fetching PDB structures..."
python3 pipeline.py --config "$CONFIG" --steps fetch

echo
echo "[3/6] Cleaning structures..."
python3 pipeline.py --config "$CONFIG" --steps clean

echo
echo "[4/6] Running PLIP interaction analysis..."
python3 pipeline.py --config "$CONFIG" --steps plip

echo
echo "[5/6] Generating Boltz-2 YAML inputs..."
python3 pipeline.py --config "$CONFIG" --steps generate

echo
echo "[6/6] Statistics"
python3 pipeline.py --config "$CONFIG" --stats

OUTPUT_DIR=$(python3 -c "import yaml,sys; print(yaml.safe_load(open('$CONFIG'))['output']['base_dir'])")

cat <<EOF

==================================
Done.
==================================

YAML inputs: ${OUTPUT_DIR}/boltz_inputs/

Three variants per structure, so the effect of the constraints is measurable
against a baseline from the same run:
  <pdb>_default.yaml         no constraints
  <pdb>_crystal_pocket.yaml  pocket constraint from the crystal interactions
  <pdb>_crystal_contact.yaml pairwise contact constraints

To predict, in the Boltz environment (see SETUP.md):

  micromamba run -n boltz boltz predict ${OUTPUT_DIR}/boltz_inputs \\
      --out_dir ${OUTPUT_DIR}/boltz_results

For the DiffDock workflow, set \$DIFFDOCK_HOME and run:

  python3 pipeline.py --config $CONFIG --steps diffdock_full
EOF
