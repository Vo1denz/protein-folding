# scripts/05_parse_results.py
# Replace Step 2 and Step 3 with this fixed version
# Add this right after Step 1 (after the docking_results list is built)

import os, re, glob
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

DOCKING_DIR = "data/docking/results"
LIGAND_SDF  = "data/inputs/ligands/prepared/filtered_ligands.sdf"
PDBQT_DIR   = "data/inputs/ligands/prepared/pdbqt"
OUTPUT_CSV  = "data/results/ranked_candidates.csv"
OUTPUT_DIR  = "outputs"

os.makedirs("data/results", exist_ok=True)
os.makedirs(OUTPUT_DIR,     exist_ok=True)

def parse_dlg(dlg_path):
    best = None
    runs_found = 0
    with open(dlg_path) as f:
        for line in f:
            if "Estimated Free Energy of Binding" in line:
                runs_found += 1
                m = re.search(r'=\s*([-\d.]+)', line)
                if m:
                    e = float(m.group(1))
                    if best is None or e < best:
                        best = e
    return best, runs_found

# ── Step 1: Parse DLG files ───────────────────────────────────────────────────
print("Step 1: Parsing DLG files...")
dlg_files = sorted(glob.glob(os.path.join(DOCKING_DIR, "*.dlg")))
dlg_files = [f for f in dlg_files if "test_" not in os.path.basename(f)]

docking_results = []
for dlg in dlg_files:
    lid = os.path.splitext(os.path.basename(dlg))[0]
    energy, n_runs = parse_dlg(dlg)
    if energy is not None:
        docking_results.append({
            'ligand_id': lid,
            'binding_energy_kcal_mol': energy,
            'n_runs': n_runs,
        })
print(f"  Parsed {len(docking_results)}/{len(dlg_files)} ligands")

# ── Step 2: Build ligand_id → real_name mapping from PDBQT REMARK lines ──────
# Why: Phase 3 saved files as ligand_1.pdbqt but the original name is inside
# the file as "REMARK  Name = AMG_510_analog". We extract that here.
print("\nStep 2: Building name mapping from PDBQT files...")

id_to_name = {}
for pdbqt in glob.glob(os.path.join(PDBQT_DIR, "*.pdbqt")):
    lid = os.path.splitext(os.path.basename(pdbqt))[0]
    with open(pdbqt) as f:
        for line in f:
            if line.startswith("REMARK  Name"):
                m = re.search(r'Name\s*=\s*(.+)', line)
                if m:
                    id_to_name[lid] = m.group(1).strip()
                break

print(f"  Mapped {len(id_to_name)} ligands:")
for lid, name in sorted(id_to_name.items(), key=lambda x: int(x[0].split('_')[1])):
    print(f"    {lid:12s} → {name}")

# ── Step 3: Load molecular properties from SDF using real names ───────────────
print("\nStep 3: Loading molecular properties from SDF...")

def get_props(mol):
    return {
        'mol_name':        mol.GetProp("_Name"),
        'smiles':          Chem.MolToSmiles(Chem.RemoveHs(mol)),
        'mol_weight':      round(Descriptors.MolWt(mol), 1),
        'logp':            round(Descriptors.MolLogP(mol), 2),
        'hbd':             rdMolDescriptors.CalcNumHBD(mol),
        'hba':             rdMolDescriptors.CalcNumHBA(mol),
        'tpsa':            round(Descriptors.TPSA(mol), 1),
        'rotatable_bonds': rdMolDescriptors.CalcNumRotatableBonds(mol),
        'heavy_atoms':     rdMolDescriptors.CalcNumHeavyAtoms(mol),
        'rings':           rdMolDescriptors.CalcNumRings(mol),
    }

mol_props = {}
supplier = Chem.SDMolSupplier(LIGAND_SDF, removeHs=False)
for mol in supplier:
    if mol:
        mol_props[mol.GetProp("_Name")] = get_props(mol)
print(f"  Loaded {len(mol_props)} molecules from SDF")

# ── Step 4: Merge everything ──────────────────────────────────────────────────
print("\nStep 4: Merging and ranking...")

df = pd.DataFrame(docking_results)

# Add real name column using the mapping
df['mol_name'] = df['ligand_id'].map(id_to_name)

# Merge molecular properties on real name
props_df = pd.DataFrame.from_dict(mol_props, orient='index').reset_index(drop=True)
df = df.merge(props_df, on='mol_name', how='left')

# Sort and rank
df = df.sort_values('binding_energy_kcal_mol').reset_index(drop=True)
df['rank'] = df.index + 1

# Categories
df['category'] = 'weak'
df.loc[df['binding_energy_kcal_mol'] < -7.0,  'category'] = 'hit'
df.loc[df['binding_energy_kcal_mol'] < -9.0,  'category'] = 'strong hit'
df.loc[df['binding_energy_kcal_mol'] < -10.0, 'category'] = 'excellent'

# Ligand efficiency = binding energy / heavy atoms (good LE > 0.3)
df['ligand_efficiency'] = (df['binding_energy_kcal_mol'] / df['heavy_atoms']).round(3)

# Outlier flag
df['outlier'] = (df['binding_energy_kcal_mol'] > 0) | (df['binding_energy_kcal_mol'] < -15)
outliers = df[df['outlier']]
if len(outliers) > 0:
    print(f"  ⚠ {len(outliers)} outlier(s) detected — check these manually")

df.to_csv(OUTPUT_CSV, index=False)

# ── Step 5: Print summary ─────────────────────────────────────────────────────
print(f"\n{'='*80}")
print(f"RANKED RESULTS — Top 10")
print(f"{'='*80}")
print(f"{'Rk':>3} {'Name':>20} {'Energy':>8} {'MW':>7} {'LogP':>6} "
      f"{'TPSA':>6} {'LE':>6}  Category")
print(f"{'-'*80}")
for _, row in df.head(10).iterrows():
    print(f"{int(row['rank']):>3} {row['mol_name']:>20} "
          f"{row['binding_energy_kcal_mol']:>8.2f} "
          f"{row['mol_weight']:>7.1f} {row['logp']:>6.2f} "
          f"{row['tpsa']:>6.1f} {row['ligand_efficiency']:>6.3f}  "
          f"{row['category']}")

print(f"\nCategory breakdown:")
for cat, count in df['category'].value_counts().items():
    print(f"  {cat:15s}: {count}")

# ── Step 6: Plots ─────────────────────────────────────────────────────────────
colors = {'weak':'#95a5a6','hit':'#3498db','strong hit':'#e67e22','excellent':'#e74c3c'}

fig, axes = plt.subplots(1, 3, figsize=(16, 5))
fig.suptitle('Phase 4 Docking Results — KRAS G12D', fontsize=13, fontweight='bold')

# Plot 1: Binding energy bar chart
ax = axes[0]
bar_colors = [colors[c] for c in df['category']]
ax.bar(range(len(df)), df['binding_energy_kcal_mol'], color=bar_colors, edgecolor='white', lw=0.5)
ax.axhline(-7.0,  color='#3498db', linestyle='--', lw=1)
ax.axhline(-9.0,  color='#e67e22', linestyle='--', lw=1)
ax.axhline(-10.0, color='#e74c3c', linestyle='--', lw=1)
# Label top 3
for _, row in df.head(3).iterrows():
    ax.text(int(row['rank'])-1, row['binding_energy_kcal_mol']-0.15,
            row['mol_name'], fontsize=6, ha='center', rotation=45)
patches = [mpatches.Patch(color=v, label=k) for k, v in colors.items()]
ax.legend(handles=patches, fontsize=8)
ax.set_xlabel('Ligand rank')
ax.set_ylabel('Binding energy (kcal/mol)')
ax.set_title('Binding energy by rank')

# Plot 2: Energy vs MW
ax = axes[1]
for cat in df['category'].unique():
    sub = df[df['category']==cat]
    ax.scatter(sub['mol_weight'], sub['binding_energy_kcal_mol'],
               c=colors[cat], label=cat, s=60, edgecolors='white', lw=0.5)
    for _, r in sub.head(2).iterrows():
        ax.annotate(r['mol_name'], (r['mol_weight'], r['binding_energy_kcal_mol']),
                    fontsize=6, xytext=(3,3), textcoords='offset points')
ax.axhline(-7.0, color='gray', linestyle='--', lw=0.8)
ax.set_xlabel('Molecular weight (Da)')
ax.set_ylabel('Binding energy (kcal/mol)')
ax.set_title('Energy vs Molecular Weight')
ax.legend(fontsize=8)

# Plot 3: Ligand Efficiency vs Energy
ax = axes[2]
for cat in df['category'].unique():
    sub = df[df['category']==cat]
    ax.scatter(sub['ligand_efficiency'], sub['binding_energy_kcal_mol'],
               c=colors[cat], label=cat, s=60, edgecolors='white', lw=0.5)
    for _, r in sub.head(2).iterrows():
        ax.annotate(r['mol_name'], (r['ligand_efficiency'], r['binding_energy_kcal_mol']),
                    fontsize=6, xytext=(3,3), textcoords='offset points')
ax.axhline(-7.0, color='gray', linestyle='--', lw=0.8)
ax.axvline(-0.3, color='gray', linestyle=':', lw=0.8, label='LE=0.3 threshold')
ax.set_xlabel('Ligand Efficiency (kcal/mol/heavy atom)')
ax.set_ylabel('Binding energy (kcal/mol)')
ax.set_title('Energy vs Ligand Efficiency')
ax.legend(fontsize=8)

plt.tight_layout()
plot_path = os.path.join(OUTPUT_DIR, "docking_results.png")
plt.savefig(plot_path, dpi=150, bbox_inches='tight')
plt.show()

print(f"\n{'='*80}")
print(f"Phase 5 complete.")
print(f"  Results CSV : {OUTPUT_CSV}")
print(f"  Best hit    : {df.iloc[0]['mol_name']}  {df.iloc[0]['binding_energy_kcal_mol']:.2f} kcal/mol")
print(f"  Best LE     : {df.nsmallest(1,'ligand_efficiency').iloc[0]['mol_name']}")
print(f"{'='*80}")