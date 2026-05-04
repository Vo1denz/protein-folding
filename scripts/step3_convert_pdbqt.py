#!/usr/bin/env python3
"""
Phase 3 — Step 3: Convert SDF to PDBQT for AutoDock-GPU
=========================================================
AutoDock-GPU does NOT read SDF or standard PDB files.
It requires PDBQT format — PDB + partial charges + AutoDock atom types.

What is PDBQT?
  - PDB: standard format with ATOM/HETATM lines, (x, y, z) per atom
  - Q: partial charges added to each atom line (column 71-76)
  - T: AutoDock atom types (e.g., A=aromatic C, OA=H-bond acceptor O, HD=donor H)

Why can't AutoDock read SDF directly?
  The scoring function (force field) uses its own atom typing system.
  It needs to know "this oxygen is an H-bond acceptor" explicitly.
  Open Babel knows how to translate standard atom types → AutoDock types.

Open Babel:
  The Swiss Army knife of chemical file conversion.
  Supports 100+ formats. We use it specifically for the -x flags that
  control partial charge calculation (Gasteiger method) and H-atom handling.

Usage:
  python scripts/step3_convert_pdbqt.py

Input:  data/inputs/ligands/prepared/filtered_ligands.sdf
Output: data/inputs/ligands/prepared/pdbqt/ligand_XXXXX.pdbqt (one per molecule)
"""

import os
import subprocess
import sys
import glob

# ── Configuration ───────────────────────────────────────────────────────────────
INPUT_SDF  = "data/inputs/ligands/prepared/filtered_ligands.sdf"
OUTPUT_DIR = "data/inputs/ligands/prepared/pdbqt"
LOG_FILE   = "data/inputs/ligands/prepared/conversion_log.txt"

os.makedirs(OUTPUT_DIR, exist_ok=True)


# ════════════════════════════════════════════════════════════════════════════════
# CONCEPT: Gasteiger Partial Charges
# ════════════════════════════════════════════════════════════════════════════════
#
# Every atom in a molecule has a partial charge — a fractional electronic charge
# caused by electronegativity differences between bonded atoms.
# Example in water (H-O-H): O is electronegative → pulls charge from H
#   O gets charge ≈ -0.8, each H gets ≈ +0.4
#
# In AutoDock's scoring function, electrostatic interaction is:
#   E_elec = q_ligand × q_receptor / distance
#
# So we need partial charges for every atom.
# Gasteiger charges are computed iteratively based on electronegativity:
#   Simple, fast, good enough for docking (AM1-BCC charges are better but slower)
#
# Open Babel computes these automatically when you use the -xh flag.


def check_openbabel() -> bool:
    import shutil
    path = shutil.which("obabel")   
    if path:
        print(f"[Check] Open Babel found: {path}")
        return True
    return False


def split_sdf_and_convert(input_sdf: str, output_dir: str) -> dict:
    """
    Strategy: Use obabel's built-in split (-m flag) to generate one PDBQT per molecule.
    
    The -m flag tells obabel: "split multi-molecule input into individual files"
    Each output file is named: {output_prefix}NNN.pdbqt where NNN is the index.
    
    Key obabel flags we use:
      -isdf          input format is SDF
      -opdbqt        output format is PDBQT
      -m             split into multiple files (one per molecule)
      -h             add hydrogens (ensure explicit H for docking)
      --partialcharge gasteiger   compute Gasteiger partial charges
      -xr            rigid (no torsion tree) — use for receptor; for ligands omit this
    
    Wait — why NOT use -xr for ligands?
      -xr makes the molecule rigid (one fixed conformation).
      For docking, ligands need to be FLEXIBLE — AutoDock will rotate bonds.
      Without -xr, obabel defines which bonds are rotatable (the torsion tree).
    """
    
    output_prefix = os.path.join(output_dir, "ligand_")
    
    cmd = [
        "obabel",
        input_sdf,               # Input SDF (multi-molecule)
        "-opdbqt",               # Output format: PDBQT
        f"-O{output_prefix}.pdbqt", # Output filename prefix (obabel adds index)
        "-m",                    # Split into individual files
        "-h",                    # Add missing hydrogens
        "--partialcharge", "gasteiger",   # Compute Gasteiger charges
        "--gen3d",               # Ensure 3D coords (our SDF already has them, but safety)
    ]
    
    print(f"[Convert] Running Open Babel to split + convert SDF → PDBQT...")
    print(f"[Convert] Command: {' '.join(cmd)}\n")
    
    result = subprocess.run(
        cmd,
        capture_output=True, text=True, timeout=300  # 5 min max
    )
    
    stats = {
        "returncode": result.returncode,
        "stdout": result.stdout,
        "stderr": result.stderr,
    }
    
    return stats


def validate_pdbqt_files(output_dir: str) -> dict:
    """
    After conversion, check each PDBQT file is valid.
    
    A valid ligand PDBQT must have:
      1. At least one HETATM or ATOM line (atom coordinates)
      2. At least one ROOT line (defines the rigid core for docking)
      3. No ATOM lines with charge = 999 (conversion error)
    
    Returns dict with validation stats.
    """
    pdbqt_files = sorted(glob.glob(os.path.join(output_dir, "*.pdbqt")))
    
    valid = []
    invalid = []
    
    for filepath in pdbqt_files:
        issues = []
        
        with open(filepath) as f:
            content = f.read()
            lines = content.strip().split("\n")
        
        # Check 1: Has atom coordinates
        atom_lines = [l for l in lines if l.startswith("HETATM") or l.startswith("ATOM")]
        if len(atom_lines) == 0:
            issues.append("no ATOM/HETATM lines")
        
        # Check 2: Has ROOT (torsion tree defined)
        if "ROOT" not in content:
            issues.append("missing ROOT (torsion tree not defined)")
        
        # Check 3: Has at least some charge info in the file
        # PDBQT charge is in column 71-76. Check it's not all zeros.
        charges = []
        for line in atom_lines:
            if len(line) >= 76:
                try:
                    charge = float(line[70:76].strip())
                    charges.append(charge)
                except ValueError:
                    pass
        
        if charges and all(c == 0.0 for c in charges):
            issues.append("all partial charges are zero (charge calculation may have failed)")
        
        fname = os.path.basename(filepath)
        if issues:
            invalid.append((fname, issues))
        else:
            valid.append(fname)
    
    return {
        "total_files":   len(pdbqt_files),
        "valid":         len(valid),
        "invalid":       len(invalid),
        "invalid_files": invalid,
        "valid_files":   valid,
    }


def fallback_individual_convert(input_sdf: str, output_dir: str) -> int:
    """
    Alternative approach: Use RDKit to split SDF, then convert each individually.
    Use this if the batch obabel -m approach fails.
    
    Why this might be needed:
    - Some obabel versions handle -m differently
    - Large molecules can cause obabel to hang on batch mode
    """
    from rdkit import Chem
    
    print("[Fallback] Using molecule-by-molecule conversion...")
    
    suppl = Chem.SDMolSupplier(input_sdf, removeHs=False)
    success_count = 0
    
    # Temp SDF for individual molecules
    tmp_sdf = os.path.join(output_dir, "_tmp_single.sdf")
    
    for i, mol in enumerate(suppl):
        if mol is None:
            continue
        
        mol_name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"mol_{i:05d}"
        out_pdbqt = os.path.join(output_dir, f"ligand_{i:05d}.pdbqt")
        
        # Write this single molecule to temp SDF
        writer = Chem.SDWriter(tmp_sdf)
        writer.write(mol)
        writer.close()
        
        # Convert this one molecule
        cmd = [
            "obabel", tmp_sdf,
            "-opdbqt", f"-O{out_pdbqt}",
            "-h",
            "--partialcharge", "gasteiger",
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        
        if result.returncode == 0 and os.path.exists(out_pdbqt):
            success_count += 1
        
        if (i + 1) % 50 == 0:
            print(f"  Converted {i+1} molecules, {success_count} succeeded...")
    
    # Clean up temp file
    if os.path.exists(tmp_sdf):
        os.remove(tmp_sdf)
    
    return success_count


def main():
    print("=" * 65)
    print("Phase 3 — Step 3: SDF → PDBQT Conversion")
    print("=" * 65)
    print(f"Input:  {INPUT_SDF}")
    print(f"Output: {OUTPUT_DIR}/\n")
    
    # ── Check prerequisites ──
    if not os.path.exists(INPUT_SDF):
        print(f"[Error] Input SDF not found: {INPUT_SDF}")
        print("[Error] Run step2_filter_ligands.py first.")
        sys.exit(1)
    
    babel_available = check_openbabel()
    
    if not babel_available:
        print("\n[Error] Open Babel (obabel) not found!")
        print("[Fix]   Install with: conda install -c conda-forge openbabel")
        print("[Fix]   Or: sudo apt-get install openbabel")
        sys.exit(1)
    
    # ── Count input molecules ──
    from rdkit import Chem
    suppl = Chem.SDMolSupplier(INPUT_SDF, removeHs=False)
    mol_count = sum(1 for m in suppl if m is not None)
    print(f"\n[Info] Input SDF contains {mol_count} molecules\n")
    
    # ── Convert ──
    stats = split_sdf_and_convert(INPUT_SDF, OUTPUT_DIR)
    
    # Print obabel output
    if stats["stdout"]:
        print("[obabel stdout]")
        print(stats["stdout"][:1000])
    if stats["stderr"]:
        print("[obabel stderr]")
        print(stats["stderr"][:1000])
    
    # ── Check if conversion worked ──
    pdbqt_files = glob.glob(os.path.join(OUTPUT_DIR, "*.pdbqt"))
    
    if len(pdbqt_files) == 0:
        print("\n[Warning] Batch conversion produced no files. Trying individual mode...")
        count = fallback_individual_convert(INPUT_SDF, OUTPUT_DIR)
        print(f"[Fallback] Converted {count} molecules individually.")
        pdbqt_files = glob.glob(os.path.join(OUTPUT_DIR, "*.pdbqt"))
    
    print(f"\n[Result] PDBQT files created: {len(pdbqt_files)}")
    
    # ── Validate ──
    print("\n[Validate] Checking PDBQT file integrity...")
    val = validate_pdbqt_files(OUTPUT_DIR)
    
    print(f"  Total files:  {val['total_files']}")
    print(f"  Valid:        {val['valid']}")
    print(f"  Invalid:      {val['invalid']}")
    
    if val["invalid_files"]:
        print(f"\n  Invalid files (first 5):")
        for fname, issues in val["invalid_files"][:5]:
            print(f"    {fname}: {', '.join(issues)}")
    
    # ── Example: Print first PDBQT file to understand the format ──
    if val["valid_files"]:
        first_file = os.path.join(OUTPUT_DIR, val["valid_files"][0])
        print(f"\n[Example] First 20 lines of {val['valid_files'][0]}:")
        print("-" * 65)
        with open(first_file) as f:
            for i, line in enumerate(f):
                if i >= 20:
                    break
                print(f"  {line}", end="")
        print("-" * 65)
        print("""
[Reading a PDBQT file:]
  Column 1-6:   Record type (ATOM/HETATM/ROOT/BRANCH/etc.)
  Column 13-16: Atom name
  Column 18-20: Residue name (UNL = unknown ligand)
  Column 31-38: X coordinate
  Column 39-46: Y coordinate  
  Column 47-54: Z coordinate
  Column 71-76: Partial charge (Gasteiger)
  Column 78-79: AutoDock atom type (C, A, N, OA, HD, etc.)
""")
    
    # ── Save log ──
    log_content = f"""Phase 3 Step 3 — Conversion Log
Input: {INPUT_SDF}
Output: {OUTPUT_DIR}
Total PDBQT files: {len(pdbqt_files)}
Valid: {val['valid']}
Invalid: {val['invalid']}
Invalid files: {val['invalid_files']}
"""
    with open(LOG_FILE, "w") as f:
        f.write(log_content)
    
    print(f"\n{'=' * 65}")
    print("PHASE 3 — STEP 3 COMPLETE ✅")
    print(f"{'=' * 65}")
    print(f"  {val['valid']} ready-to-dock PDBQT files in {OUTPUT_DIR}/")
    print(f"\nNext steps:")
    print(f"  1. Run step4_verify_library.py to get a visual summary")
    print(f"  2. Then move to Phase 4: AutoDock-GPU molecular docking")
    print(f"{'=' * 65}")


if __name__ == "__main__":
    main()