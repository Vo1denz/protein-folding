# scripts/extract_pocket_centers.py

import os
import glob
import numpy as np
import pandas as pd

POCKET_DIR = "data/structures/kras/KRAS_G12D_Human_KRAS_proto-oncogene_GTPase_G12D_mutant_unrelaxed_rank_001_alphafold2_model_1_seed_000_out/pockets"

# ── What this does ────────────────────────────────────────────────────────────
# Each pocket_N_atm.pdb contains the protein atoms that line that pocket.
# The geometric center (centroid) of those atoms = the pocket center for docking.
# We average all X, Y, Z coordinates of ATOM records to get that centroid.
# ─────────────────────────────────────────────────────────────────────────────

def get_centroid(pdb_file):
    """Read all ATOM/HETATM lines, return mean x,y,z and atom count."""
    coords = []
    with open(pdb_file) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append((x, y, z))
                except ValueError:
                    continue
    if not coords:
        return None, 0
    arr = np.array(coords)
    return arr.mean(axis=0), len(coords)


pocket_files = sorted(
    glob.glob(os.path.join(POCKET_DIR, "pocket*_atm.pdb")),
    key=lambda f: int(os.path.basename(f).replace("pocket","").replace("_atm.pdb",""))
)

print(f"Found {len(pocket_files)} pocket files\n")

# Load druggability scores from CSV we already have
scores_df = pd.read_csv("data/pockets/all_pockets.csv")

results = []
for pf in pocket_files:
    pid = int(os.path.basename(pf).replace("pocket","").replace("_atm.pdb",""))
    centroid, n_atoms = get_centroid(pf)
    if centroid is None:
        continue

    row = scores_df[scores_df['pocket_id'] == pid]
    drugg  = row['druggability_score'].values[0] if len(row) else 0
    volume = row['volume_A3'].values[0]          if len(row) else 0
    score  = row['score'].values[0]              if len(row) else 0

    results.append({
        'pocket_id':    pid,
        'center_x':     round(centroid[0], 3),
        'center_y':     round(centroid[1], 3),
        'center_z':     round(centroid[2], 3),
        'n_atoms':      n_atoms,
        'druggability': drugg,
        'volume_A3':    volume,
        'score':        score,
    })

df = pd.DataFrame(results).sort_values('druggability', ascending=False).reset_index(drop=True)

print(f"{'='*70}")
print(f"{'ID':>4} {'Drugg':>7} {'Score':>7} {'Vol(Å³)':>9} {'Atoms':>6}   Center (x, y, z)")
print(f"{'='*70}")
for _, r in df.iterrows():
    marker = " ← BEST" if _ == 0 else ""
    print(f"{int(r.pocket_id):>4} {r.druggability:>7.3f} {r.score:>7.3f} "
          f"{r.volume_A3:>9.1f} {int(r.n_atoms):>6}   "
          f"({r.center_x}, {r.center_y}, {r.center_z}){marker}")

best = df.iloc[0]
print(f"\n{'='*70}")
print(f"BEST POCKET  →  Pocket {int(best.pocket_id)}")
print(f"  Druggability : {best.druggability:.3f}")
print(f"  Volume       : {best.volume_A3:.1f} Å³")
print(f"  Center X     : {best.center_x}")
print(f"  Center Y     : {best.center_y}")
print(f"  Center Z     : {best.center_z}")
print(f"{'='*70}")
print(f"\nSave these — you'll paste them into the AutoGrid config in Phase 4.")

# Save updated CSV with coordinates filled in
df.to_csv("data/pockets/pockets_with_centers.csv", index=False)
print(f"\nFull table saved → data/pockets/pockets_with_centers.csv")