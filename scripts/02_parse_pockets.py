# scripts/02_parse_pockets.py
#!/usr/bin/env python3

import re
import sys
import json
import os
import glob
import csv

# ─────────────────────────────────────────────
# SECTION 1: Parse the _info.txt file  (FIXED)
# ─────────────────────────────────────────────

def get_val(block, label):
    """
    Extract a float value after a label in a pocket block.
    
    THE BUG WE FIXED:
    The old version used re.search(f'{label}\\s*:') which meant
    searching for "Score" would match "Druggability Score" first,
    returning the wrong number.
    
    THE FIX:
    We anchor the match to the START OF A LINE using ^ + re.MULTILINE.
    This means "Score :" only matches a line that BEGINS with Score,
    never a line that has "Druggability Score".
    
    Pattern breakdown:
      ^\s+     — start of line, then leading whitespace (tabs/spaces)
      {label}  — the exact label text
      \s*:\s+  — colon with optional space before, required space after
      ([-\d.]+) — capture the number (handles negatives like -4)
    """
    pattern = rf'^\s+{re.escape(label)}\s*:\s+([-\d.]+)'
    match = re.search(pattern, block, re.MULTILINE)
    return float(match.group(1)) if match else None


def parse_fpocket_info(info_file_path):
    """
    Parse the fpocket _info.txt file into a list of pocket dicts.
    
    The file format (confirmed from your output) is:
    
        Pocket 1 :
                Score :         0.110
                Druggability Score :    0.019
                Volume :        1011.014
                ...
        Pocket 2 :
                Score :         0.089
                ...
    
    We split on 'Pocket N :' lines to isolate each block.
    """
    with open(info_file_path) as f:
        content = f.read()

    # Split into blocks at each "Pocket N :" line
    # re.MULTILINE makes ^ match start of each line
    pocket_blocks = re.split(r'(?=^Pocket \d+ :)', content, flags=re.MULTILINE)
    pocket_blocks = [b for b in pocket_blocks if b.strip().startswith("Pocket")]

    print(f"  Raw pocket blocks found in file: {len(pocket_blocks)}")

    pockets = []

    for block in pocket_blocks:
        # Extract pocket number
        pocket_id_match = re.match(r'Pocket (\d+) :', block)
        if not pocket_id_match:
            continue
        pocket_id = int(pocket_id_match.group(1))

        # Extract each field using the fixed get_val()
        # NOTE: label strings must match EXACTLY what's in the file
        score           = get_val(block, "Score")
        druggability    = get_val(block, "Druggability Score")
        volume          = get_val(block, "Volume")
        center_x        = get_val(block, "Center x")
        center_y        = get_val(block, "Center y")
        center_z        = get_val(block, "Center z")
        alpha_spheres   = get_val(block, "Number of Alpha Spheres")
        hydrophobic     = get_val(block, "Mean local hydrophobic density")

        # Debug: print what we extracted per pocket so you can verify
        print(f"  Pocket {pocket_id:>2}: score={score}, "
              f"druggability={druggability}, volume={volume}, "
              f"center=({center_x},{center_y},{center_z})")

        # Skip if the core values are missing
        if None in [score, druggability, volume]:
            print(f"    ↳ SKIPPED — missing core values")
            continue

        pockets.append({
            "pocket_id":          pocket_id,
            "score":              round(score, 4),
            "druggability_score": round(druggability, 4),
            "volume_A3":          round(volume, 2),
            "center_x":           round(center_x, 3) if center_x is not None else None,
            "center_y":           round(center_y, 3) if center_y is not None else None,
            "center_z":           round(center_z, 3) if center_z is not None else None,
            "alpha_spheres":      int(alpha_spheres) if alpha_spheres else None,
            "hydrophobic_density": round(hydrophobic, 4) if hydrophobic is not None else None,
        })

    return pockets


# ─────────────────────────────────────────────
# SECTION 2: Filter pockets
# ─────────────────────────────────────────────

def filter_pockets(pockets, min_druggability=0.5, min_volume=300, max_volume=2000):
    """
    Standard drug-like filters:
    - druggability > 0.5  : fpocket's ML estimate of binding likelihood
    - volume 300–2000 Å³  : drug-sized molecule needs room to fit
    """
    filtered = [
        p for p in pockets
        if p["druggability_score"] >= min_druggability
        and min_volume <= p["volume_A3"] <= max_volume
    ]
    filtered.sort(key=lambda x: x["druggability_score"], reverse=True)
    return filtered


# ─────────────────────────────────────────────
# SECTION 3: Handle KRAS-specific situation
# ─────────────────────────────────────────────

def select_best_pocket(filtered_pockets, all_pockets):
    """
    KRAS G12D SPECIAL CASE — important context:
    
    KRAS G12D is different from KRAS G12C (sotorasib target).
    - G12C has a well-defined Switch II Pocket (S-IIP) that opens reliably
    - G12D does NOT have the same covalent-binding site
    - G12D is harder — the pocket is shallower and more transient
    
    If no pocket passes druggability > 0.5, this is EXPECTED for G12D.
    We lower the threshold and pick the best available pocket,
    which is still scientifically valid for virtual screening.
    
    This is why G12D was called "undruggable" for decades — the pocket
    genuinely scores lower than G12C in computational tools.
    """
    all_sorted = sorted(all_pockets, key=lambda x: x["druggability_score"], reverse=True)

    if not filtered_pockets:
        print("\n" + "="*60)
        print("NOTE: No pocket scored > 0.5 druggability.")
        print("This is EXPECTED for KRAS G12D — the Switch II pocket")
        print("is shallower than in G12C and scores lower in fpocket.")
        print("Selecting best available pocket (highest druggability).")
        print("="*60)
        best = all_sorted[0]
    else:
        best = filtered_pockets[0]

    print("\n" + "="*60)
    print("POCKET SELECTION REPORT")
    print("="*60)
    print(f"\n  Selected: Pocket {best['pocket_id']}")
    print(f"  Druggability score : {best['druggability_score']:.4f}")
    print(f"  Geometric score    : {best['score']:.4f}")
    print(f"  Volume             : {best['volume_A3']:.2f} Å³")
    if best['center_x'] is not None:
        print(f"  Center             : ({best['center_x']}, {best['center_y']}, {best['center_z']})")
    else:
        print(f"  Center             : not found in info file")
        print(f"  → Will extract from pocket PDB file instead (see below)")

    print(f"\nAll pockets ranked by druggability:")
    print(f"  {'ID':>3}  {'Druggability':>12}  {'Score':>7}  {'Volume Å³':>10}")
    print(f"  " + "-"*38)
    for p in all_sorted:
        marker = " ← SELECTED" if p['pocket_id'] == best['pocket_id'] else ""
        print(f"  {p['pocket_id']:>3}  {p['druggability_score']:>12.4f}  "
              f"{p['score']:>7.4f}  {p['volume_A3']:>10.2f}{marker}")

    return best


# ─────────────────────────────────────────────
# SECTION 4: Extract center from pocket PDB
# (fallback if _info.txt has no Center x/y/z)
# ─────────────────────────────────────────────

def extract_center_from_pocket_pdb(fpocket_dir, pocket_id):
    """
    Some fpocket versions don't write Center x/y/z to the info file.
    In that case we compute the center ourselves from the pocket atom file.
    
    The pocket<N>_atm.pdb file contains the protein atoms that LINE the pocket.
    The geometric center of these atoms = the pocket center for docking.
    """
    pocket_pdb = os.path.join(fpocket_dir, "pockets", f"pocket{pocket_id}_atm.pdb")

    if not os.path.exists(pocket_pdb):
        print(f"  WARNING: {pocket_pdb} not found")
        return None, None, None

    xs, ys, zs = [], [], []
    with open(pocket_pdb) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    # PDB format: columns 31-38 = X, 39-46 = Y, 47-54 = Z
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    xs.append(x)
                    ys.append(y)
                    zs.append(z)
                except ValueError:
                    continue

    if not xs:
        return None, None, None

    cx = round(sum(xs) / len(xs), 3)
    cy = round(sum(ys) / len(ys), 3)
    cz = round(sum(zs) / len(zs), 3)

    print(f"  Center computed from {len(xs)} pocket atoms: ({cx}, {cy}, {cz})")
    return cx, cy, cz


# ─────────────────────────────────────────────
# SECTION 5: Save docking config
# ─────────────────────────────────────────────

def save_docking_config(best_pocket, fpocket_dir, output_path):
    cx = best_pocket.get("center_x")
    cy = best_pocket.get("center_y")
    cz = best_pocket.get("center_z")

    # If center wasn't in info file, compute it from the pocket PDB
    if cx is None:
        print("\n  Center not in info file — computing from pocket atom coordinates...")
        cx, cy, cz = extract_center_from_pocket_pdb(fpocket_dir, best_pocket["pocket_id"])

    if cx is None:
        print("  ERROR: Could not determine pocket center by any method.")
        return None

    docking_config = {
        "center_x":          cx,
        "center_y":          cy,
        "center_z":          cz,
        "size_x":            22,
        "size_y":            22,
        "size_z":            22,
        "pocket_id":         best_pocket["pocket_id"],
        "druggability_score": best_pocket["druggability_score"],
        "volume_A3":         best_pocket["volume_A3"],
    }

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(docking_config, f, indent=2)

    print(f"\n  ✓ Docking config saved: {output_path}")
    print(f"    Center : ({cx}, {cy}, {cz})")
    print(f"    Box    : 22 × 22 × 22 Å")
    return docking_config


# ─────────────────────────────────────────────
# SECTION 6: Main
# ─────────────────────────────────────────────

def main():
    if len(sys.argv) < 2:
        print("Usage: python 02_parse_pockets.py <fpocket_output_dir>")
        sys.exit(1)

    fpocket_dir = sys.argv[1].rstrip('/')

    # Find info file
    candidates = glob.glob(os.path.join(fpocket_dir, "*_info.txt"))
    if not candidates:
        print(f"ERROR: No *_info.txt in {fpocket_dir}")
        sys.exit(1)
    info_file = sorted(candidates, key=os.path.getsize, reverse=True)[0]
    print(f"Info file : {os.path.basename(info_file)}")

    # Parse
    print("\nParsing pockets...")
    all_pockets = parse_fpocket_info(info_file)
    print(f"\nSuccessfully parsed: {len(all_pockets)} pockets")

    if not all_pockets:
        print("ERROR: Parsed 0 pockets. Check debug output above.")
        sys.exit(1)

    # Filter
    filtered = filter_pockets(all_pockets)

    # Select
    best = select_best_pocket(filtered, all_pockets)

    # Save config for Stage 4
    config_path = "data/docking/pocket_config.json"
    config = save_docking_config(best, fpocket_dir, config_path)

    # Save full CSV for notebook
    csv_path = "data/pockets/all_pockets.csv"
    os.makedirs("data/pockets", exist_ok=True)
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=all_pockets[0].keys())
        writer.writeheader()
        writer.writerows(all_pockets)
    print(f"\n  ✓ Full pocket table: {csv_path}")

    print("\n" + "="*60)
    print("Phase 2 complete ✓")
    print("Next: Phase 3 — Ligand Library Prep (RDKit)")
    print("="*60)


if __name__ == "__main__":
    main()