# scripts/02_detect_pockets.sh
#!/bin/bash

KRAS_DIR="data/structures/kras"

# ── Automatically find the best ranked relaxed PDB ──────────────
# ColabFold always puts "relaxed_rank_001" in the filename of the best model
PDB_FILE=$(find "$KRAS_DIR" -name "*relaxed_rank_001*.pdb" -type f | head -1)

# If no relaxed file, fall back to any rank_001 pdb
if [ -z "$PDB_FILE" ]; then
    PDB_FILE=$(find "$KRAS_DIR" -name "*rank_001*.pdb" -type f | head -1)
fi

# If still nothing found, list what IS there and exit clearly
if [ -z "$PDB_FILE" ]; then
    echo "ERROR: Could not find a rank_001 PDB file in $KRAS_DIR"
    echo ""
    echo "Files present in $KRAS_DIR:"
    ls "$KRAS_DIR"
    echo ""
    echo "Fix: Pass the correct PDB path manually:"
    echo "  bash scripts/02_detect_pockets.sh <path_to_your.pdb>"
    exit 1
fi

# If a path was passed manually as argument, use that instead
if [ ! -z "$1" ]; then
    PDB_FILE=$1
fi

echo "Using PDB file: $PDB_FILE"

# ── Run fpocket ──────────────────────────────────────────────────
fpocket -f "$PDB_FILE"

# ── Derive output directory name ────────────────────────────────
# fpocket creates: <same_directory>/<basename_without_extension>_out/
PDB_DIR=$(dirname "$PDB_FILE")
BASENAME=$(basename "$PDB_FILE" .pdb)
OUTPUT_DIR="${PDB_DIR}/${BASENAME}_out"

echo ""
echo "fpocket output folder: $OUTPUT_DIR"
echo ""

# ── Verify output was created ───────────────────────────────────
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "ERROR: fpocket did not create the output folder."
    echo "This usually means fpocket couldn't parse the PDB."
    echo "Try: fpocket -f \"$PDB_FILE\" and check for errors above."
    exit 1
fi

echo "Files generated:"
ls "$OUTPUT_DIR/"

echo ""
echo "Quick summary of top 5 pockets:"
INFO_FILE="$OUTPUT_DIR/${BASENAME}_info.txt"
grep -A 10 "Pocket [1-5] :" "$INFO_FILE" | grep -E "Score|Druggability|Volume|Center"

# ── Save the paths for the next script ──────────────────────────
# Write the output dir path to a temp file so 02_parse_pockets.py can find it
echo "$OUTPUT_DIR" > data/pockets/.last_fpocket_output
echo ""
echo "Path saved to data/pockets/.last_fpocket_output"
echo "Next: python scripts/02_parse_pockets.py \"$OUTPUT_DIR\""