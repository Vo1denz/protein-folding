# scripts/04a_prepare_receptor.py

import os
import subprocess

# ── Paths ────────────────────────────────────────────────────────────────────
PDB_IN  = "data/structures/kras/KRAS_G12D_Human_KRAS_proto-oncogene_GTPase_G12D_mutant_unrelaxed_rank_001_alphafold2_model_1_seed_000.pdb"
PDB_CLEAN   = "data/structures/kras/receptor_clean.pdb"
PDBQT_OUT   = "data/structures/kras/receptor.pdbqt"

os.makedirs("data/structures/kras", exist_ok=True)

# ── Step 1: Clean the PDB ─────────────────────────────────────────────────────
# Why: AlphaFold PDBs sometimes have OXT atoms and non-standard records.
# We keep only ATOM lines (backbone + sidechains) and TER (chain terminators).
# HETATM = heteroatoms (waters, ligands, ions) — AutoDock doesn't want these
# in the receptor file; they confuse the charge assignment.

print("Step 1: Cleaning PDB...")
kept = 0
removed = 0

with open(PDB_IN) as f_in, open(PDB_CLEAN, 'w') as f_out:
    for line in f_in:
        record = line[:6].strip()
        if record in ("ATOM", "TER", "END"):
            f_out.write(line)
            kept += 1
        else:
            removed += 1

print(f"  Kept    : {kept} lines")
print(f"  Removed : {removed} lines (HETATM, REMARK, etc.)")
print(f"  Saved   : {PDB_CLEAN}")

# ── Step 2: Convert to PDBQT with Open Babel ─────────────────────────────────
# Why obabel: it adds Gasteiger partial charges (needed by AutoDock's scoring
# function) and assigns AutoDock4 atom types (A, C, NA, OA, N, SA, HD, etc.)
# -xr flag = receptor mode (rigid, no torsions added)

print("\nStep 2: Converting to PDBQT...")
cmd = [
    "obabel",
    PDB_CLEAN,
    "-O", PDBQT_OUT,
    "-xr",          # receptor mode — no rotatable bonds
    "-h",           # add polar hydrogens (needed for H-bond scoring)
]

result = subprocess.run(cmd, capture_output=True, text=True)

if result.returncode != 0:
    print("ERROR from obabel:")
    print(result.stderr)
else:
    print(f"  Saved   : {PDBQT_OUT}")
    # Quick sanity check — count ATOM lines in PDBQT
    with open(PDBQT_OUT) as f:
        atom_lines = sum(1 for l in f if l.startswith("ATOM"))
    print(f"  ATOM records in PDBQT: {atom_lines}")
    if atom_lines > 100:
        print("  ✓ Receptor PDBQT looks good")
    else:
        print("  ✗ WARNING: Very few atoms — check the PDB input")

print("\nStderr from obabel (informational):")
print(result.stderr[:300] if result.stderr else "  (none)")