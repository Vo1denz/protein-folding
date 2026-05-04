# notebooks/02_pocket_analysis.py

import json
import csv
import os
import matplotlib.pyplot as plt
import numpy as np

# ── Make sure outputs folder exists ──────────────────────────
os.makedirs("outputs", exist_ok=True)

# ── Load the pocket data ──────────────────────────────────────
with open("data/pockets/all_pockets.csv") as f:
    reader = csv.DictReader(f)
    pockets = []
    for row in reader:
        pocket = {}
        for k, v in row.items():
            if v == "None" or v == "":
                continue
            if k in ["pocket_id", "alpha_spheres"]:
                pocket[k] = int(float(v))
            else:
                pocket[k] = float(v)
        pockets.append(pocket)

ids          = [p["pocket_id"]         for p in pockets]
druggability = [p["druggability_score"] for p in pockets]
volumes      = [p["volume_A3"]          for p in pockets]
scores       = [p["score"]              for p in pockets]

print(f"Loaded {len(pockets)} pockets from CSV")

# ── Plot 1: Druggability scores per pocket ────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Green = passes threshold, grey = fails
# NOTE: for KRAS G12D none will be green — that's expected and correct
colors = ["#2ecc71" if d >= 0.5 else "#bdc3c7" for d in druggability]

axes[0].bar(ids, druggability, color=colors, edgecolor="white", linewidth=0.8)
axes[0].axhline(0.5, color="red", linestyle="--", linewidth=1.5,
                label="Druggability threshold (0.5)")
axes[0].set_xlabel("Pocket ID")
axes[0].set_ylabel("Druggability Score")
axes[0].set_title("fpocket: Druggability Score per Pocket\n(KRAS G12D)")  # fixed: G12D
axes[0].legend()

# Annotate best pocket
best_idx = druggability.index(max(druggability))
axes[0].annotate(
    f"Best: Pocket {ids[best_idx]}\n({max(druggability):.3f})",
    xy=(ids[best_idx], max(druggability)),
    xytext=(ids[best_idx] + 1.5, max(druggability) + 0.02),
    arrowprops=dict(arrowstyle="->", color="black"),
    fontsize=9
)

# Add a note explaining low scores
axes[0].text(0.98, 0.95,
    "No pocket exceeds 0.5\n(expected for KRAS G12D —\nshallow binding site)",
    transform=axes[0].transAxes,
    fontsize=7.5, color="#555555",
    ha="right", va="top",
    bbox=dict(boxstyle="round,pad=0.3", facecolor="#fff9c4", edgecolor="#cccccc"))

# ── Plot 2: Druggability vs Volume scatter ────────────────────
sc = axes[1].scatter(volumes, druggability,
                     c=scores, cmap="viridis",
                     s=80, edgecolors="white", linewidth=0.5, zorder=3)

axes[1].axhline(0.5, color="red",    linestyle="--", linewidth=1,
                alpha=0.7, label="Druggability ≥ 0.5")
axes[1].axvline(300, color="orange", linestyle="--", linewidth=1,
                alpha=0.7, label="Volume ≥ 300 Å³")

# Shade the ideal zone (top-right quadrant)
axes[1].fill_betweenx([0.5, 1.0], 300, max(volumes) * 1.05,
                       alpha=0.08, color="green", label="Ideal druggable zone")

# Label each point with pocket number
for p in pockets:
    if "volume_A3" in p and "druggability_score" in p:
        axes[1].annotate(
            f"P{p['pocket_id']}",
            (p["volume_A3"], p["druggability_score"]),
            textcoords="offset points", xytext=(5, 5), fontsize=7
        )

# Highlight the selected pocket with a red ring
selected_pocket_id = json.load(open("data/docking/pocket_config.json"))["pocket_id"]
for p in pockets:
    if p["pocket_id"] == selected_pocket_id:
        axes[1].scatter(p["volume_A3"], p["druggability_score"],
                        s=160, facecolors="none", edgecolors="red",
                        linewidth=2, zorder=4, label=f"Selected (P{selected_pocket_id})")
        break

plt.colorbar(sc, ax=axes[1], label="Geometric score")
axes[1].set_xlabel("Volume (Å³)")
axes[1].set_ylabel("Druggability Score")
axes[1].set_title("Druggability vs Volume — KRAS G12D\n(selected pocket circled in red)")
axes[1].legend(fontsize=8)

plt.tight_layout()
plt.savefig("outputs/pocket_analysis.png", dpi=150, bbox_inches="tight")
print("Plot saved: outputs/pocket_analysis.png")
plt.show()

# ── Print final selection ─────────────────────────────────────
with open("data/docking/pocket_config.json") as f:
    config = json.load(f)

print("\n" + "="*50)
print("SELECTED DOCKING BOX FOR STAGE 4")
print("="*50)
print(f"Pocket ID    : {config['pocket_id']}")
print(f"Druggability : {config['druggability_score']}")
print(f"Volume       : {config['volume_A3']} Å³")
print(f"Center       : ({config['center_x']}, {config['center_y']}, {config['center_z']})")
print(f"Box size     : {config['size_x']} × {config['size_y']} × {config['size_z']} Å")