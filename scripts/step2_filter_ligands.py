#!/usr/bin/env python3
"""
Phase 3 — Step 2: Filter, Prepare & Minimize Ligands with RDKit
================================================================
This is the core of Phase 3. We take raw SMILES strings and turn them
into physically realistic 3D molecular structures ready for docking.

Pipeline inside this script:
  SMILES string
      │
      ▼ Parse with RDKit
      │
      ▼ Lipinski Rule of Five filter  ← removes non-drug-like molecules
      │
      ▼ PAINS filter                  ← removes assay-interference compounds  
      │
      ▼ Add hydrogens                 ← explicit H atoms needed for 3D geometry
      │
      ▼ ETKDG conformer generation    ← creates a 3D shape from 2D SMILES
      │
      ▼ MMFF energy minimization      ← physically realistic geometry
      │
      ▼ Save to .sdf file             ← standard 3D molecule format

Usage:
  python scripts/step2_filter_ligands.py data/inputs/ligands/zinc_druglike_500.smi

Output:
  data/inputs/ligands/prepared/filtered_ligands.sdf   ← 3D molecules
  data/inputs/ligands/prepared/filter_report.csv      ← stats per molecule
  data/inputs/ligands/prepared/filter_summary.txt     ← overall stats
"""

import os
import sys
import csv
import time

from fastapi import params

# ─── RDKit imports ──────────────────────────────────────────────────────────────
# RDKit is structured as several sub-modules. Here's what each does:
from rdkit import Chem                           # Core molecule parsing
from rdkit.Chem import AllChem                   # 3D conformer generation, force fields
from rdkit.Chem import Descriptors               # MW, LogP, TPSA, etc.
from rdkit.Chem import rdMolDescriptors          # H-bond donors/acceptors, rings
from rdkit.Chem.FilterCatalog import (           # PAINS filter database
    FilterCatalog,
    FilterCatalogParams
)

# ── Configuration ───────────────────────────────────────────────────────────────
INPUT_SMILES = sys.argv[1] if len(sys.argv) > 1 else "data/inputs/ligands/zinc_druglike_500.smi"
OUTPUT_DIR   = "data/inputs/ligands/prepared"
OUTPUT_SDF   = os.path.join(OUTPUT_DIR, "filtered_ligands.sdf")
OUTPUT_CSV   = os.path.join(OUTPUT_DIR, "filter_report.csv")
OUTPUT_TXT   = os.path.join(OUTPUT_DIR, "filter_summary.txt")

os.makedirs(OUTPUT_DIR, exist_ok=True)


# ════════════════════════════════════════════════════════════════════════════════
# CONCEPT: SMILES Notation
# ════════════════════════════════════════════════════════════════════════════════
#
# SMILES (Simplified Molecular-Input Line-Entry System) encodes a 3D molecule
# as a flat string. Think of it as the "text format" for chemistry.
#
# How to read SMILES:
#   - Atoms are their element symbol:  C=carbon, N=nitrogen, O=oxygen, c=aromatic C
#   - Single bond:  implicit (CC = ethane)
#   - Double bond:  = (C=O = ketone)
#   - Triple bond:  # (C#N = nitrile)
#   - Branches:     () → CC(=O)O = acetic acid (C with branch =O and branch O)
#   - Rings:        numbers mark ring closures → C1CCCCC1 = cyclohexane
#   - Aromatic:     lowercase → c1ccccc1 = benzene
#
# Examples:
#   CC(=O)Oc1ccccc1C(=O)O     → aspirin
#   CC(C)Cc1ccc(cc1)C(C)C(=O)O → ibuprofen
#   c1ccccc1                  → benzene
#
# SMILES is what ZINC15 gives us. RDKit reads it and builds an internal
# molecular graph, from which we can compute anything.


def parse_smiles_file(filepath: str) -> list[tuple[str, str]]:
    """
    Read a .smi file and return (smiles, name) pairs.
    
    .smi format: one molecule per line, "SMILES name"
    Lines starting with # are comments.
    """
    molecules = []
    
    with open(filepath) as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            
            # Handle both space-separated and tab-separated
            parts = line.split("\t") if "\t" in line else line.split()
            
            smiles = parts[0]
            # Name is optional — use index if missing
            name = parts[1] if len(parts) > 1 else f"mol_{i:05d}"
            
            molecules.append((smiles, name))
    
    return molecules


# ════════════════════════════════════════════════════════════════════════════════
# FILTER 1: Lipinski Rule of Five
# ════════════════════════════════════════════════════════════════════════════════
#
# In 1997, Christopher Lipinski analyzed all FDA-approved oral drugs and found
# they share these properties (with a few exceptions like antibiotics):
#
#   MW   < 500 Da      — small enough to be absorbed through the gut wall
#   HBD  ≤ 5           — H-bond donors (NH, OH groups) — too many = can't cross membranes
#   HBA  ≤ 10          — H-bond acceptors (N, O atoms)
#   LogP ≤ 5           — lipophilicity. Too high = insoluble in blood. Too low = can't enter cells
#
# "Rule of Five" because all limits are multiples of 5.
# A molecule violating 2+ rules is unlikely to be orally bioavailable.
#
# For KRAS:
# - KRAS inhibitors tend to be slightly larger (MW 400-500) due to the pocket shape
# - We'll use strict Lipinski but note near-miss compounds too

def lipinski_filter(mol) -> tuple[bool, dict]:
    """
    Apply Lipinski Rule of Five.
    Returns (passed: bool, properties: dict)
    
    We return properties even on failure so we can log WHY it failed.
    """
    mw   = Descriptors.MolWt(mol)
    hbd  = rdMolDescriptors.CalcNumHBD(mol)   # H-bond Donors: NH, OH
    hba  = rdMolDescriptors.CalcNumHBA(mol)   # H-bond Acceptors: N, O
    logp = Descriptors.MolLogP(mol)
    
    # Extra useful properties (not Lipinski, but good to log)
    tpsa = Descriptors.TPSA(mol)              # Topological Polar Surface Area
    rot  = rdMolDescriptors.CalcNumRotatableBonds(mol)
    rings = rdMolDescriptors.CalcNumRings(mol)
    
    props = {
        "mol_weight": round(mw, 2),
        "hbd": hbd,
        "hba": hba,
        "logp": round(logp, 2),
        "tpsa": round(tpsa, 2),
        "rot_bonds": rot,
        "num_rings": rings,
    }
    
    # Count violations (not just pass/fail — lets us track near-misses)
    violations = sum([
        mw   >= 500,
        hbd  >  5,
        hba  >  10,
        logp >  5,
    ])
    
    props["lipinski_violations"] = violations
    passed = violations <= 1   # Allow 1 violation (Lipinski himself allowed this)
    
    return passed, props


# ════════════════════════════════════════════════════════════════════════════════
# FILTER 2: PAINS (Pan-Assay Interference Compounds)
# ════════════════════════════════════════════════════════════════════════════════
#
# Some molecules cause false positives in biological assays — not because they
# actually bind the target, but because they:
#   - Form aggregates (like tiny soap bubbles) that trap proteins
#   - React covalently and non-specifically with any protein
#   - Fluoresce (confusing fluorescence-based assays)
#   - Chelate metal ions needed for enzyme activity
#
# PAINS filters identify these structural motifs.
# The catalog is split into PAINS_A, PAINS_B, PAINS_C (increasing strictness).
# We use all three.
#
# Example PAINS patterns to avoid:
#   - Rhodanine scaffold (C1=CC(=O)SC(=S)N1) — classic assay interference
#   - α,β-unsaturated carbonyl — Michael acceptors that react non-specifically
#   - Frequent hitters: quinones, catechols, oxidizable compounds

def build_pains_catalog() -> FilterCatalog:
    """
    Build RDKit's PAINS filter catalog.
    This only needs to be created once (it's slow to initialize).
    """
    params = FilterCatalogParams()
    # Add all three PAINS sub-catalogs
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
    return FilterCatalog(params)


def pains_filter(mol, catalog: FilterCatalog) -> tuple[bool, str]:
    """
    Returns (passed: bool, reason: str)
    If failed, reason describes which PAINS pattern matched.
    """
    if catalog.HasMatch(mol):
        # Get which pattern(s) matched — useful for debugging
        matches = catalog.GetMatches(mol)
        reason = "; ".join([m.GetDescription() for m in matches])
        return False, reason
    return True, "clean"


# ════════════════════════════════════════════════════════════════════════════════
# STEP 3: 3D Conformer Generation
# ════════════════════════════════════════════════════════════════════════════════
#
# A SMILES string is 2D — it encodes connectivity but not 3D coordinates.
# To dock a molecule, we need (x, y, z) coordinates for every atom.
#
# WHAT IS A CONFORMER?
# The same molecule can adopt different 3D shapes depending on bond rotations.
# These are called "conformers" or "conformations".
# Example: butane (CCCC) can be in gauche or anti conformation.
#
# For docking, we generate ONE starting conformer. AutoDock will then explore
# many conformations during the docking search — so we just need a reasonable start.
#
# ETKDG (Experimental-Torsion Knowledge Distance Geometry):
# - Distance Geometry: place atoms in 3D by first computing pairwise distance
#   constraints (from bond lengths, angles, ring geometry)
# - Then apply torsional corrections from an experimental database of crystal
#   structures (so torsion angles are realistic, not random)
# - ETKDGv3 = the 2022 version, most accurate
#
# WHY NOT JUST USE ANY 3D GEOMETRY?
# Bad starting geometry = docking failure or wrong result.
# MMFF minimization then corrects small geometry errors.

def generate_conformer(mol, mol_name: str) -> tuple[object | None, str]:
    """
    Add explicit H atoms, generate a 3D conformer with ETKDG, minimize with MMFF.
    
    Returns (mol_3d, status_message)
    Returns (None, error_msg) if conformer generation fails.
    """
    # ── Step A: Add explicit hydrogens ──
    # RDKit molecules initially have implicit H (just counted, not placed).
    # For 3D geometry, we need explicit H atoms with coordinates.
    mol_h = Chem.AddHs(mol)
    
    # ── Step B: ETKDG conformer generation ──
    # EmbedMolecule places the heavy atoms + H in 3D space.
    # randomSeed=42 makes results reproducible.
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.numThreads = 0     # Try up to 5 times before giving up
    
    result = AllChem.EmbedMolecule(mol_h, params)
    
    # Return code: 0 = success, -1 = failure
    if result == -1:
        # WHEN DOES THIS FAIL?
        # - Very strained ring systems (e.g., 4-membered rings fused to others)
        # - Extremely large molecules
        # - Molecules with unusual valence
        # Try once more with random coordinates as fallback
        params.useRandomCoords = True
        result = AllChem.EmbedMolecule(mol_h, params)
        
        if result == -1:
            return None, f"ETKDG failed even with random coords"
    
    # ── Step C: MMFF Force Field Minimization ──
    # ETKDG gives good geometry, but MMFF makes it physically correct.
    #
    # WHAT IS MMFF (Merck Molecular Force Field)?
    # A mathematical model of atomic interactions:
    #   E_total = E_bond + E_angle + E_torsion + E_vdW + E_electrostatic
    #
    # Minimization = adjust coordinates to minimize total energy.
    # This corrects: clashing atoms, weird bond angles, unrealistic torsions.
    #
    # maxIters=2000: stop after 2000 steps even if not fully converged
    # (full convergence not needed — just a good starting point for docking)
    
    ff_result = AllChem.MMFFOptimizeMolecule(mol_h, maxIters=2000, mmffVariant="MMFF94")
    
    # ff_result: 0 = converged, 1 = not converged but improved, -1 = MMFF not applicable
    if ff_result == -1:
        # MMFF can't handle some unusual atom types — fall back to UFF force field
        uff_result = AllChem.UFFOptimizeMolecule(mol_h, maxIters=2000)
        if uff_result == -1:
            return None, "Neither MMFF nor UFF could minimize this molecule"
        status = "minimized_with_UFF"
    elif ff_result == 1:
        status = "minimized_MMFF_not_converged_but_ok"
    else:
        status = "minimized_MMFF_converged"
    
    mol_h.SetProp("_Name", mol_name)  # Tag with molecule ID
    return mol_h, status


# ════════════════════════════════════════════════════════════════════════════════
# MAIN PIPELINE
# ════════════════════════════════════════════════════════════════════════════════

def main():
    start_time = time.time()
    
    print("=" * 65)
    print("Phase 3 — Step 2: RDKit Ligand Preparation")
    print("=" * 65)
    print(f"Input:  {INPUT_SMILES}")
    print(f"Output: {OUTPUT_SDF}\n")
    
    # ── Load SMILES ──
    raw_molecules = parse_smiles_file(INPUT_SMILES)
    print(f"[Parse] Loaded {len(raw_molecules)} SMILES entries from file\n")
    
    # ── Initialize PAINS catalog (once) ──
    print("[Setup] Building PAINS filter catalog...")
    pains_catalog = build_pains_catalog()
    print("[Setup] PAINS catalog ready\n")
    
    # ── Stats tracking ──
    stats = {
        "total_input":       len(raw_molecules),
        "invalid_smiles":    0,
        "lipinski_fail":     0,
        "pains_fail":        0,
        "conformer_fail":    0,
        "passed":            0,
    }
    
    passed_molecules = []
    report_rows = []   # Per-molecule report for CSV
    
    # ── Process each molecule ──
    print(f"[Process] Starting molecule pipeline...\n")
    
    for idx, (smiles, name) in enumerate(raw_molecules):
        row = {
            "index": idx,
            "name": name,
            "smiles": smiles,
            "status": "",
            "lipinski_pass": "",
            "pains_pass": "",
            "mol_weight": "",
            "logp": "",
            "hbd": "",
            "hba": "",
            "tpsa": "",
            "lipinski_violations": "",
            "pains_reason": "",
            "conformer_status": "",
        }
        
        # ── Parse SMILES ──
        # MolFromSmiles returns None if SMILES is invalid
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            stats["invalid_smiles"] += 1
            row["status"] = "FAIL_invalid_smiles"
            report_rows.append(row)
            continue
        
        # ── Lipinski filter ──
        lip_pass, props = lipinski_filter(mol)
        row.update(props)
        row["lipinski_pass"] = lip_pass
        
        if not lip_pass:
            stats["lipinski_fail"] += 1
            row["status"] = f"FAIL_lipinski ({props['lipinski_violations']} violations)"
            report_rows.append(row)
            continue
        
        # ── PAINS filter ──
        pains_pass, pains_reason = pains_filter(mol, pains_catalog)
        row["pains_pass"] = pains_pass
        row["pains_reason"] = pains_reason
        
        if not pains_pass:
            stats["pains_fail"] += 1
            row["status"] = f"FAIL_pains ({pains_reason[:60]})"
            report_rows.append(row)
            continue
        
        # ── 3D Conformer Generation ──
        mol_3d, conformer_status = generate_conformer(mol, name)
        row["conformer_status"] = conformer_status
        
        if mol_3d is None:
            stats["conformer_fail"] += 1
            row["status"] = f"FAIL_conformer ({conformer_status})"
            report_rows.append(row)
            continue
        
        # ── All filters passed! ──
        stats["passed"] += 1
        row["status"] = "PASS"
        report_rows.append(row)
        passed_molecules.append(mol_3d)
        
        # Progress update every 50 molecules
        if (idx + 1) % 50 == 0 or (idx + 1) == len(raw_molecules):
            elapsed = time.time() - start_time
            print(f"  [{idx+1:4d}/{len(raw_molecules)}] "
                  f"Passed: {stats['passed']:3d} | "
                  f"Lip fail: {stats['lipinski_fail']:3d} | "
                  f"PAINS fail: {stats['pains_fail']:3d} | "
                  f"Conf fail: {stats['conformer_fail']:3d} | "
                  f"Time: {elapsed:.1f}s")
    
    # ── Save filtered molecules as SDF ──
    print(f"\n[Save] Writing {len(passed_molecules)} molecules to SDF...")
    writer = Chem.SDWriter(OUTPUT_SDF)
    for mol in passed_molecules:
        writer.write(mol)
    writer.close()
    print(f"[Save] ✅ {OUTPUT_SDF}")
    
    # ── Save per-molecule report as CSV ──
    print(f"[Save] Writing per-molecule report...")
    if report_rows:
        fieldnames = report_rows[0].keys()
        with open(OUTPUT_CSV, "w", newline="") as f:
            writer_csv = csv.DictWriter(f, fieldnames=fieldnames)
            writer_csv.writeheader()
            writer_csv.writerows(report_rows)
    print(f"[Save] ✅ {OUTPUT_CSV}")
    
    # ── Final summary ──
    total_time = time.time() - start_time
    pass_rate = stats["passed"] / max(stats["total_input"], 1) * 100
    
    summary = f"""
{'=' * 65}
PHASE 3 — STEP 2 COMPLETE
{'=' * 65}

Input library:          {stats['total_input']:>6d} molecules
Invalid SMILES:         {stats['invalid_smiles']:>6d} ({stats['invalid_smiles']/max(stats['total_input'],1)*100:.1f}%)
Lipinski failures:      {stats['lipinski_fail']:>6d} ({stats['lipinski_fail']/max(stats['total_input'],1)*100:.1f}%)
PAINS failures:         {stats['pains_fail']:>6d} ({stats['pains_fail']/max(stats['total_input'],1)*100:.1f}%)
Conformer failures:     {stats['conformer_fail']:>6d} ({stats['conformer_fail']/max(stats['total_input'],1)*100:.1f}%)
                        ──────
FINAL PASSED:           {stats['passed']:>6d} ({pass_rate:.1f}%)

Time elapsed:           {total_time:.1f}s

Output files:
  3D molecules (SDF):   {OUTPUT_SDF}
  Per-mol report (CSV): {OUTPUT_CSV}

Next step:
  Convert SDF → PDBQT:  python scripts/step3_convert_pdbqt.py
{'=' * 65}
"""
    print(summary)
    
    # Save summary to text file too
    with open(OUTPUT_TXT, "w") as f:
        f.write(summary)


if __name__ == "__main__":
    main()