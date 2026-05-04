#!/bin/bash
# ============================================================
# Stage 1: Protein Structure Prediction using ColabFold
# Usage: bash scripts/01_fold_protein.sh data/inputs/sequences/egfr.fasta
# ============================================================

INPUT_FASTA=$1

# Derive output dir name from the FASTA filename (without extension)
PROTEIN_NAME=$(basename "$INPUT_FASTA" .fasta)
OUTPUT_DIR="data/structures/${PROTEIN_NAME}"

mkdir -p "$OUTPUT_DIR"

echo "=================================================="
echo "  Folding protein: $PROTEIN_NAME"
echo "  Input: $INPUT_FASTA"
echo "  Output: $OUTPUT_DIR"
echo "=================================================="

colabfold_batch \
  --num-recycle 3 \
  --model-type alphafold2_ptm \
  --use-gpu-relax \
  "$INPUT_FASTA" \
  "$OUTPUT_DIR"

echo ""
echo "Folding complete. Files in: $OUTPUT_DIR"
echo "Look for: *_relaxed_rank_001_*.pdb  ← this is your best structure"