import json, glob
import matplotlib.pyplot as plt

# Filter specifically for the scores file
json_files = glob.glob("data/structures/kras/*scores_rank_001*.json")
print("JSON file found:", json_files)

with open(json_files[0]) as f:
    data = json.load(f)

print("Keys in file:", data.keys())  # let's see what's inside

plddt = data['plddt']
residues = list(range(1, len(plddt) + 1))

plt.figure(figsize=(12, 4))
plt.plot(residues, plddt, color='steelblue', linewidth=1.5)
plt.axhline(70, color='orange', linestyle='--', label='Threshold (70)')
plt.axhline(90, color='green',  linestyle='--', label='High confidence (90)')
plt.fill_between(residues, plddt, 70,
                 where=[p > 70 for p in plddt],
                 alpha=0.3, color='green', label='Reliable')
plt.fill_between(residues, plddt, 70,
                 where=[p < 70 for p in plddt],
                 alpha=0.3, color='red', label='Low confidence')
plt.xlabel('Residue position')
plt.ylabel('pLDDT score')
plt.title('KRAS G12D — AlphaFold2 per-residue confidence')
plt.legend()
plt.tight_layout()
plt.savefig('outputs/kras_plddt.png', dpi=150)
plt.show()

high = sum(1 for p in plddt if p > 70)
low  = len(plddt) - high
print(f"Mean pLDDT    : {sum(plddt)/len(plddt):.1f}")
print(f"Residues > 70 : {high}/{len(plddt)}  ← dockable region")
print(f"Residues < 70 : {low}/{len(plddt)}   ← flexible/disordered") n