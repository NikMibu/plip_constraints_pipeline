#!/bin/bash
# Example Pipeline Run Script
# This demonstrates a typical pipeline workflow

set -e  # Exit on error

echo "=================================="
echo "PLIP Constraints Pipeline Example"
echo "=================================="
echo ""

# Configuration
CONFIG_FILE="config.yaml"
OUTPUT_DIR="./output"

# Check if config exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "ERROR: Config file not found: $CONFIG_FILE"
    exit 1
fi

# Check Docker
echo "[1/5] Checking Docker..."
if ! docker run --rm pharmai/plip --help > /dev/null 2>&1; then
    echo "ERROR: Docker or PLIP image not available"
    echo "Run: docker pull pharmai/plip"
    exit 1
fi
echo "✓ Docker ready"
echo ""

# Step 1: Fetch PDB data
echo "[2/5] Fetching PDB structures..."
python pipeline.py --config "$CONFIG_FILE" --steps fetch
echo ""

# Step 2: Clean structures
echo "[3/5] Cleaning structures..."
python pipeline.py --config "$CONFIG_FILE" --steps clean
echo ""

# Step 3: PLIP analysis
echo "[4/5] Running PLIP analysis..."
python pipeline.py --config "$CONFIG_FILE" --steps plip
echo ""

# Step 4: Generate Boltz-2 YAMLs
echo "[5/5] Generating Boltz-2 YAML files..."
python pipeline.py --config "$CONFIG_FILE" --steps generate
echo ""

# Show statistics
echo "=================================="
echo "Pipeline Statistics"
echo "=================================="
python pipeline.py --config "$CONFIG_FILE" --stats
echo ""

# Success message
echo "=================================="
echo "✓ Pipeline completed successfully!"
echo "=================================="
echo ""
echo "Output directory: $OUTPUT_DIR"
echo ""
echo "Next steps:"
echo "  1. Review generated YAMLs in: $OUTPUT_DIR/boltz_inputs/"
echo "  2. Copy to Boltz-2 input directory"
echo "  3. Run Boltz-2 predictions"
echo ""

