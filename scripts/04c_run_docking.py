# scripts/04c_run_docking.py

import os
import subprocess
import glob
import time

# ── Paths ─────────────────────────────────────────────────────────────────────
AUTODOCK_BIN  = "/mnt/d/projects/Protein_Folding/protein-drug-finder/AutoDock-GPU/bin/autodock_gpu_128wi"
MAPS_FLD      = "data/structures/kras/receptor.maps.fld"
LIGAND_DIR    = "data/inputs/ligands/prepared/pdbqt"
OUTPUT_DIR    = "data/docking/results"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── What each parameter means ─────────────────────────────────────────────────
# --ffile   : the .maps.fld file AutoGrid generated (defines the search space)
# --lfile   : the ligand PDBQT to dock
# --nrun    : number of independent genetic algorithm runs per ligand
#             More runs = better sampling of pose space = more reliable results
#             20 is standard for screening; 50+ for lead optimization
# --resnam  : output file base name (produces NAME.dlg)
# --gbest   : write only the single best pose per run (keeps files small)
# --lsmet   : local search method — "sw" = Solis-Wets (standard for AutoDock4)
# --seed    : random seed for reproducibility
# ─────────────────────────────────────────────────────────────────────────────

ligand_files = sorted(glob.glob(os.path.join(LIGAND_DIR, "*.pdbqt")))
total = len(ligand_files)

print(f"AutoDock-GPU Docking Run")
print(f"{'='*50}")
print(f"Binary   : {AUTODOCK_BIN}")
print(f"Maps     : {MAPS_FLD}")
print(f"Ligands  : {total}")
print(f"Output   : {OUTPUT_DIR}")
print(f"{'='*50}\n")

# Verify maps file exists
if not os.path.exists(MAPS_FLD):
    print(f"ERROR: Maps file not found: {MAPS_FLD}")
    print("Did AutoGrid4 complete successfully?")
    exit(1)

failed  = []
success = []
start_total = time.time()

for i, ligand_path in enumerate(ligand_files, 1):
    lig_name = os.path.splitext(os.path.basename(ligand_path))[0]
    output_base = os.path.join(OUTPUT_DIR, lig_name)
    dlg_file    = output_base + ".dlg"

    print(f"[{i:2d}/{total}] Docking: {lig_name} ...", end=" ", flush=True)
    t0 = time.time()

    result = subprocess.run(
        [
            AUTODOCK_BIN,
            "--ffile", MAPS_FLD,
            "--lfile", ligand_path,
            "--nrun",  "20",
            "--resnam", output_base,
            "--gbest", "1",
            "--lsmet", "sw",
            "--seed",  "1234",
        ],
        capture_output=True,
        text=True
    )

    elapsed = time.time() - t0

    # ── Validate output ───────────────────────────────────────────────────────
    # A successful run always produces a non-empty .dlg file.
    # Empty or missing DLG = atom type mismatch or grid box too small.
    if os.path.exists(dlg_file) and os.path.getsize(dlg_file) > 0:
        # Quick check: does it contain an energy value?
        with open(dlg_file) as f:
            content = f.read()
        if "Free Energy of Binding" in content:
            print(f"✓  ({elapsed:.1f}s)")
            success.append(lig_name)
        else:
            print(f"✗  DLG produced but no energy found — possible grid miss")
            failed.append(lig_name)
    else:
        print(f"✗  No DLG produced")
        print(f"   stderr: {result.stderr[-200:] if result.stderr else 'none'}")
        failed.append(lig_name)

total_time = time.time() - start_total

print(f"\n{'='*50}")
print(f"Docking complete in {total_time:.1f}s  ({total_time/60:.1f} min)")
print(f"  Successful : {len(success)}/{total}")
print(f"  Failed     : {len(failed)}/{total}")
if failed:
    print(f"\nFailed ligands:")
    for f in failed:
        print(f"  - {f}")
print(f"\nDLG files in: {OUTPUT_DIR}")
print(f"Next step   : python scripts/05_parse_results.py")