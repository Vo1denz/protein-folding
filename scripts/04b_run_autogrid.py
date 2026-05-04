# scripts/04b_run_autogrid.py

import os
import subprocess

# ── Paths ─────────────────────────────────────────────────────────────────────
RECEPTOR_PDBQT = "data/structures/kras/receptor.pdbqt"
DOCKING_DIR    = "data/docking"
GPF_FILE       = os.path.join(DOCKING_DIR, "grid.gpf")
GLG_FILE       = os.path.join(DOCKING_DIR, "grid.glg")

os.makedirs(DOCKING_DIR, exist_ok=True)

# ── Pocket coordinates from Phase 2 ───────────────────────────────────────────
CENTER_X = -9.341
CENTER_Y =  2.976
CENTER_Z = -11.641

# ── Box size ──────────────────────────────────────────────────────────────────
# Why 22x22x22 grid points at 0.375 Å spacing?
#   22 × 0.375 = 8.25 Å half-width → total box = 16.5 Å per side
#   Pocket volume is 544 Å³ ≈ cube of side ~8 Å — our box comfortably covers it
#   Too large = slow grid generation + ligand finds false pockets
#   Too small = ligand gets clipped during docking
NPTS    = 22
SPACING = 0.375

# ── Atom types ────────────────────────────────────────────────────────────────
# These must cover every atom type present in both receptor AND ligands.
# A=aromatic C, C=aliphatic C, NA=N-acceptor, OA=O-acceptor,
# N=nitrogen, SA=S-acceptor, HD=H-bond donor hydrogen
ATOM_TYPES = ["A", "C", "NA", "OA", "N", "SA", "HD"]

# ── Build the receptor base name (AutoGrid uses this for output file naming) ──
receptor_base = os.path.abspath(RECEPTOR_PDBQT).replace(".pdbqt", "")

# ── Write the Grid Parameter File (.gpf) ─────────────────────────────────────
# The GPF is AutoGrid's config file. Each 'map' line tells AutoGrid to compute
# one energy grid for one atom type. The ligand will probe all these grids.

map_lines = "\n".join(f"map {receptor_base}.{t}.map" for t in ATOM_TYPES)

gpf_content = f"""npts {NPTS} {NPTS} {NPTS}                  # number of grid points x y z
parameter_file /mnt/d/projects/Protein_Folding/protein-drug-finder/AutoDock-GPU/input/1ac8/derived/AD4.1_bound.dat  # force field params
gridfld {receptor_base}.maps.fld              # output field file
spacing {SPACING}                             # grid spacing in Angstroms
receptor_types {" ".join(ATOM_TYPES)}         # atom types in receptor
ligand_types {" ".join(ATOM_TYPES)}           # atom types to map for ligand
receptor {os.path.abspath(RECEPTOR_PDBQT)}    # receptor PDBQT
gridcenter {CENTER_X} {CENTER_Y} {CENTER_Z}   # pocket center from fpocket
smooth 0.5                                    # potential smoothing
{map_lines}
elecmap {receptor_base}.e.map                 # electrostatic map
dsolvmap {receptor_base}.d.map                # desolvation map
dielectric -0.1465                            # distance-dependent dielectric
"""

with open(GPF_FILE, 'w') as f:
    f.write(gpf_content)

print(f"Grid Parameter File written: {GPF_FILE}")
print(f"  Pocket center : ({CENTER_X}, {CENTER_Y}, {CENTER_Z})")
print(f"  Grid size     : {NPTS}×{NPTS}×{NPTS} points")
print(f"  Spacing       : {SPACING} Å  →  box = {NPTS*SPACING:.1f} Å per side")
print(f"  Atom types    : {ATOM_TYPES}")

# ── Run AutoGrid4 ─────────────────────────────────────────────────────────────
print(f"\nRunning AutoGrid4...")
print("(This generates one energy map per atom type — takes ~30 seconds)\n")

result = subprocess.run(
    ["autogrid4", "-p", GPF_FILE, "-l", GLG_FILE],
    capture_output=True, text=True
)

# ── Parse the log for success/failure ─────────────────────────────────────────
if result.returncode != 0:
    print("ERROR: AutoGrid4 failed")
    print(result.stderr)
else:
    print("AutoGrid4 finished. Checking output...\n")

    # Check that all map files were actually created
    expected_maps = [f"{receptor_base}.{t}.map" for t in ATOM_TYPES]
    expected_maps += [f"{receptor_base}.e.map", f"{receptor_base}.d.map"]
    expected_maps += [f"{receptor_base}.maps.fld"]

    all_ok = True
    for mf in expected_maps:
        exists = os.path.exists(mf)
        status = "✓" if exists else "✗ MISSING"
        print(f"  {status}  {os.path.basename(mf)}")
        if not exists:
            all_ok = False

    print()
    if all_ok:
        print("✓ All grid maps generated successfully")
        print("\nReady for Step 4.3 — running AutoDock-GPU")
    else:
        print("✗ Some maps missing — check grid.glg for errors:")
        print(f"  cat {GLG_FILE} | tail -30")