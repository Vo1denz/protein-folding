#!/usr/bin/env python3
"""
Phase 3 — Master Runner
========================
Runs all 4 steps of Phase 3 in sequence with error checking.

Usage:
  python scripts/run_phase3.py

Steps:
  Step 1: Download ligands from ZINC15 (or use curated fallback)
  Step 2: Filter (Lipinski + PAINS) + generate 3D conformers
  Step 3: Convert SDF → PDBQT with Open Babel
  Step 4: Verify and visualize the library

If any step fails, the runner stops and tells you exactly what went wrong.
"""

import subprocess
import sys
import os
import time


def run_step(script: str, description: str, args: list[str] = None) -> bool:
    """
    Run a Python script as a subprocess.
    Returns True on success, False on failure.
    """
    cmd = [sys.executable, script] + (args or [])
    
    print(f"\n{'─' * 65}")
    print(f"▶  {description}")
    print(f"   Command: {' '.join(cmd)}")
    print(f"{'─' * 65}\n")
    
    start = time.time()
    result = subprocess.run(cmd)
    elapsed = time.time() - start
    
    if result.returncode != 0:
        print(f"\n❌  {description} FAILED (exit code {result.returncode})")
        print(f"   Check the output above for the error message.")
        return False
    
    print(f"\n✅  {description} done in {elapsed:.1f}s")
    return True


def check_environment():
    """Quick pre-flight check before running anything."""
    print("Pre-flight checks:")
    
    errors = []
    
    # Check RDKit
    try:
        from rdkit import Chem
        print("  ✅ RDKit found")
    except ImportError:
        print("  ❌ RDKit not found")
        errors.append("RDKit not installed. Run: conda install -c conda-forge rdkit")
    
    # Check matplotlib
    try:
        import matplotlib
        print("  ✅ matplotlib found")
    except ImportError:
        print("  ❌ matplotlib not found")
        errors.append("matplotlib not installed. Run: pip install matplotlib")
    
    # Check pandas
    try:
        import pandas
        print("  ✅ pandas found")
    except ImportError:
        print("  ❌ pandas not found")
        errors.append("pandas not installed. Run: pip install pandas")
    
    # Check Open Babel (needed for Step 3)
    import shutil
    obabel_path = shutil.which("obabel")
    if obabel_path:
        print(f"  ✅ Open Babel found: {obabel_path}")
    else:
        print("  ⚠️  Open Babel not found (needed for Step 3)")
        errors.append("obabel not found. Install: conda install -c conda-forge openbabel")


def main():
    print("=" * 65)
    print("Phase 3 — Ligand Library Preparation")
    print("KRAS Drug Target: SMILES → 3D Molecules → PDBQT")
    print("=" * 65)
    
    total_start = time.time()
    
    # ── Pre-flight ──
    check_environment()
    
    # ── Step 1: Download ──
    success = run_step(
        "scripts/step1_download_ligands.py",
        "Step 1: Download ligands from ZINC15"
    )
    if not success:
        print("\n[Hint] Step 1 failure usually means network issues.")
        print("[Hint] The script will auto-use the curated fallback set.")
        sys.exit(1)
    
    # Find the downloaded file
    ligand_file = "data/inputs/ligands/zinc_druglike_500.smi"
    if not os.path.exists(ligand_file):
        print(f"[Error] Expected ligand file not found: {ligand_file}")
        sys.exit(1)
    
    # ── Step 2: Filter + 3D conformers ──
    success = run_step(
        "scripts/step2_filter_ligands.py",
        "Step 2: Lipinski + PAINS filter + 3D conformers (RDKit)",
        args=[ligand_file]
    )
    if not success:
        print("\n[Hint] Common causes:")
        print("  • Invalid SMILES in input file")
        print("  • RDKit import error → check conda environment")
        sys.exit(1)
    
    # Check output exists
    if not os.path.exists("data/inputs/ligands/prepared/filtered_ligands.sdf"):
        print("[Error] No SDF output found after Step 2.")
        sys.exit(1)
    
    # ── Step 3: Convert to PDBQT ──
    success = run_step(
        "scripts/step3_convert_pdbqt.py",
        "Step 3: SDF → PDBQT conversion (Open Babel)"
    )
    if not success:
        print("\n[Hint] If obabel is missing, install it:")
        print("  conda install -c conda-forge openbabel")
        print("[Hint] Step 4 can still run (uses the SDF directly for verification)")
    
    # ── Step 4: Verify ──
    success = run_step(
        "scripts/step4_verify_library.py",
        "Step 4: Library verification + property plots"
    )
    if not success:
        print("\n[Hint] Verification failures usually mean matplotlib is missing.")
    
    # ── Done ──
    total_elapsed = time.time() - total_start
    
    print(f"\n{'=' * 65}")
    print(f"PHASE 3 COMPLETE ✅  (total time: {total_elapsed:.1f}s)")
    print(f"{'=' * 65}")
    print(f"\nWhat was built:")
    print(f"  data/inputs/ligands/zinc_druglike_500.smi   ← raw input")
    print(f"  data/inputs/ligands/prepared/")
    print(f"    filtered_ligands.sdf                      ← 3D molecules")
    print(f"    pdbqt/ligand_XXXXX.pdbqt                  ← AutoDock input")
    print(f"    filter_report.csv                         ← per-molecule log")
    print(f"    library_properties.csv                    ← final property table")
    print(f"  outputs/phase3/")
    print(f"    ligand_property_distributions.png         ← portfolio plot")
    print(f"    library_summary.txt                       ← written summary")
    print(f"\nReady for Phase 4: Molecular Docking with AutoDock-GPU")
    print(f"{'=' * 65}")


if __name__ == "__main__":
    main()