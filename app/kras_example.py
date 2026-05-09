# app/gradio_app_kras.py
# KRAS G12D demo UI — pre-computed results
# Run with: python app/gradio_app_kras.py
# Opens at: http://localhost:7860

import gradio as gr
import pandas as pd
import os
import re
import glob

# Always resolve paths relative to project root
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

RESULTS_CSV    = "data/results/ranked_candidates.csv"
DOCKING_DIR    = "data/docking/results"
STRUCTURES_DIR = "data/structures/kras"
PDBQT_DIR      = "data/inputs/ligands/prepared/pdbqt"
PLOT_PATH      = "outputs/docking_results.png"
PLDDT_PATH     = "outputs/example_plddt_plot.png"


# ── Helpers ───────────────────────────────────────────────────────────────────
def load_results():
    if os.path.exists(RESULTS_CSV):
        return pd.read_csv(RESULTS_CSV)
    return None

def get_docked_pose(ligand_id):
    dlg_path = os.path.join(DOCKING_DIR, f"{ligand_id}.dlg")
    if not os.path.exists(dlg_path):
        return None
    pose_lines = []
    in_model = False
    with open(dlg_path) as f:
        for line in f:
            if "DOCKED: MODEL" in line:
                in_model = True
                pose_lines = []
                continue
            if "DOCKED: ENDMDL" in line and in_model:
                break
            if in_model and line.startswith("DOCKED:"):
                pose_lines.append(line[8:])
    return "".join(pose_lines) if pose_lines else None


# ── 3D Viewer ─────────────────────────────────────────────────────────────────
def make_3d_viewer(mol_name, ligand_id, energy):
    pose_pdb = get_docked_pose(ligand_id)
    receptor_pdb = ""
    receptor_pdb_path = os.path.join(STRUCTURES_DIR,
        "KRAS_G12D_Human_KRAS_proto-oncogene_GTPase_G12D_mutant_unrelaxed_rank_001_alphafold2_model_1_seed_000.pdb")
    if os.path.exists(receptor_pdb_path):
        with open(receptor_pdb_path) as f:
            receptor_pdb = f.read()

    if not pose_pdb:
        return f"<div style='padding:30px;color:#ff6b6b;font-family:monospace;text-align:center'>⚠ No docked pose found for {mol_name}</div>"

    receptor_escaped = receptor_pdb.replace('\\','\\\\').replace('`','\\`').replace('${','\\${')
    pose_escaped     = pose_pdb.replace('\\','\\\\').replace('`','\\`').replace('${','\\${')

    inner_html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<style>
  @import url('https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;600&family=Syne:wght@600;700&display=swap');
  * {{ margin:0; padding:0; box-sizing:border-box; }}
  body {{ background:#04080f; font-family:'JetBrains Mono',monospace; color:#8899aa; }}
  #header {{
    padding:10px 16px;
    font-size:11px;
    color:#00e5cc;
    border-bottom:1px solid #0d2030;
    display:flex;
    justify-content:space-between;
    align-items:center;
    letter-spacing:0.05em;
  }}
  #header .name {{ font-size:13px; color:#e0f0ff; font-weight:600; }}
  #header .energy {{
    background:linear-gradient(90deg,#00e5cc22,#00e5cc11);
    border:1px solid #00e5cc44;
    padding:3px 10px;
    border-radius:4px;
    color:#00e5cc;
    font-weight:600;
  }}
  #viewer {{ width:100%; height:370px; }}
  #footer {{
    padding:6px 16px;
    font-size:10px;
    color:#1a3040;
    border-top:1px solid #0d2030;
    letter-spacing:0.04em;
  }}
</style>
</head>
<body>
<div id="header">
  <span class="name">◈ {mol_name}</span>
  <span class="energy">{energy:.2f} kcal/mol</span>
</div>
<div id="viewer"></div>
<div id="footer">PROTEIN: cartoon ribbon &nbsp;·&nbsp; LIGAND: element sticks &nbsp;·&nbsp; WebGL via 3Dmol.js</div>
<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.6/3Dmol-min.js"></script>
<script>
window.addEventListener('load', function() {{
  var viewer = $3Dmol.createViewer(document.getElementById('viewer'), {{ backgroundColor: '#04080f' }});
  var receptorData = `{receptor_escaped}`;
  if (receptorData.trim().length > 0) {{
    viewer.addModel(receptorData, 'pdb');
    viewer.setStyle({{ model:0 }}, {{ cartoon:{{ color:'#00e5cc', opacity:0.85 }} }});
    viewer.addSurface($3Dmol.SurfaceType.VDW, {{ opacity:0.06, color:'#00e5cc' }}, {{ model:0 }});
  }}
  var ligandData = `{pose_escaped}`;
  viewer.addModel(ligandData, 'pdbqt');
  viewer.setStyle({{ model:1 }}, {{
    stick:  {{ colorscheme:'elementColors', radius:0.22 }},
    sphere: {{ colorscheme:'elementColors', scale:0.32 }}
  }});
  viewer.zoomTo({{ model:1 }});
  viewer.zoom(0.75);
  viewer.render();
}});
</script>
</body>
</html>"""

    srcdoc = inner_html.replace('"', '&quot;')
    return f'<iframe srcdoc="{srcdoc}" style="width:100%;height:430px;border:1px solid #0d2030;border-radius:8px;background:#04080f;" sandbox="allow-scripts"></iframe>'


# ── Event functions ───────────────────────────────────────────────────────────
def show_viewer(selection):
    if not selection:
        return "<div style='padding:40px;color:#1a3a4a;text-align:center;font-family:monospace;font-size:12px;'>SELECT A LIGAND TO VIEW ITS DOCKING POSE</div>"
    m = re.match(r'^(.+?)\s*\((-?[\d.]+)', selection)
    if not m:
        return "<p style='color:#ff6b6b'>Could not parse selection.</p>"
    mol_name = m.group(1).strip()
    energy   = float(m.group(2))
    ligand_id = None
    for pdbqt in glob.glob(os.path.join(PDBQT_DIR, "*.pdbqt")):
        with open(pdbqt) as f:
            for line in f:
                if "REMARK  Name" in line and mol_name in line:
                    ligand_id = os.path.splitext(os.path.basename(pdbqt))[0]
                    break
        if ligand_id:
            break
    if not ligand_id:
        return f"<p style='color:#ff6b6b;font-family:monospace'>PDBQT not found for {mol_name}</p>"
    return make_3d_viewer(mol_name, ligand_id, energy)


def show_info(selection):
    if not selection or load_results() is None:
        return ""
    mol_name = re.match(r'^(.+?)\s*\(', selection).group(1).strip()
    df  = load_results()
    row = df[df['mol_name'] == mol_name]
    if row.empty:
        return ""
    r = row.iloc[0]

    category_colors = {
        "excellent":  ("#ff4757", "#2d0a0a"),
        "strong hit": ("#ff9f43", "#2d1a00"),
        "hit":        ("#00e5cc", "#002a26"),
        "weak":       ("#546e7a", "#0d1a1f"),
    }
    col, bg = category_colors.get(r['category'], ("#546e7a", "#0d1a1f"))

    return f"""
    <div style="
        display:flex; gap:12px; flex-wrap:wrap; margin-top:12px;
        font-family:'JetBrains Mono',monospace; font-size:11px;
    ">
      <div style="
          background:{bg}; border:1px solid {col}44;
          border-left:3px solid {col};
          padding:10px 14px; border-radius:6px; min-width:160px;
      ">
        <div style="color:#8899aa;font-size:10px;letter-spacing:.08em;margin-bottom:4px;">COMPOUND</div>
        <div style="color:#e0f0ff;font-size:12px;font-weight:600;">{r['mol_name']}</div>
        <div style="color:{col};margin-top:4px;font-size:10px;letter-spacing:.06em;">● {r['category'].upper()}</div>
      </div>
      <div style="
          background:#060e18; border:1px solid #0d2030;
          padding:10px 14px; border-radius:6px; min-width:160px;
      ">
        <div style="color:#8899aa;font-size:10px;letter-spacing:.08em;margin-bottom:4px;">BINDING</div>
        <div style="color:#00e5cc;font-size:14px;font-weight:600;">{r['binding_energy_kcal_mol']:.2f} <span style="font-size:10px;color:#8899aa;">kcal/mol</span></div>
        <div style="color:#546e7a;margin-top:4px;">LE: {r['ligand_efficiency']:.3f} kcal/mol/atom</div>
      </div>
      <div style="
          background:#060e18; border:1px solid #0d2030;
          padding:10px 14px; border-radius:6px; min-width:200px;
      ">
        <div style="color:#8899aa;font-size:10px;letter-spacing:.08em;margin-bottom:6px;">PROPERTIES</div>
        <div style="display:grid;grid-template-columns:1fr 1fr;gap:4px 16px;">
          <div><span style="color:#546e7a;">MW</span> <span style="color:#c0d8e8;">{r['mol_weight']:.1f} Da</span></div>
          <div><span style="color:#546e7a;">LogP</span> <span style="color:#c0d8e8;">{r['logp']:.2f}</span></div>
          <div><span style="color:#546e7a;">TPSA</span> <span style="color:#c0d8e8;">{r['tpsa']:.1f} Å²</span></div>
          <div><span style="color:#546e7a;">HBD</span> <span style="color:#c0d8e8;">{int(r['hbd']) if 'hbd' in r else '—'}</span></div>
        </div>
      </div>
    </div>"""


def refresh_results():
    df = load_results()
    if df is None:
        return None
    cols = ['rank','mol_name','binding_energy_kcal_mol','mol_weight','logp','tpsa','ligand_efficiency','category']
    return df[[c for c in cols if c in df.columns]]


# ── CSS ───────────────────────────────────────────────────────────────────────
css = """
@import url('https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;500&family=Syne:wght@600;700;800&display=swap');

/* ── Global ── */
.gradio-container {
    max-width: 1280px !important;
    background: #04080f !important;
    font-family: 'JetBrains Mono', monospace !important;
}
body, .dark { background: #04080f !important; }

/* ── Header area ── */
#hero {
    text-align: center;
    padding: 40px 0 24px;
    position: relative;
}
#hero::before {
    content: '';
    position: absolute;
    top: 0; left: 50%; transform: translateX(-50%);
    width: 600px; height: 1px;
    background: linear-gradient(90deg, transparent, #00e5cc44, transparent);
}

/* ── Tabs ── */
.tab-nav button {
    font-family: 'JetBrains Mono', monospace !important;
    font-size: 11px !important;
    letter-spacing: 0.08em !important;
    color: #546e7a !important;
    border-bottom: 2px solid transparent !important;
    padding: 10px 20px !important;
    background: transparent !important;
    text-transform: uppercase !important;
}
.tab-nav button.selected {
    color: #00e5cc !important;
    border-bottom-color: #00e5cc !important;
}

/* ── Buttons ── */
button.primary {
    background: linear-gradient(135deg, #00e5cc, #0097a7) !important;
    color: #04080f !important;
    font-family: 'JetBrains Mono', monospace !important;
    font-weight: 600 !important;
    font-size: 11px !important;
    letter-spacing: 0.08em !important;
    border: none !important;
    border-radius: 6px !important;
    transition: opacity 0.2s !important;
}
button.primary:hover { opacity: 0.85 !important; }

button.secondary {
    background: transparent !important;
    color: #00e5cc !important;
    border: 1px solid #00e5cc44 !important;
    font-family: 'JetBrains Mono', monospace !important;
    font-size: 11px !important;
    letter-spacing: 0.06em !important;
    border-radius: 6px !important;
    transition: border-color 0.2s, background 0.2s !important;
}
button.secondary:hover {
    border-color: #00e5cc !important;
    background: #00e5cc11 !important;
}

/* ── Dataframe ── */
.svelte-1gfkn6j, table {
    font-family: 'JetBrains Mono', monospace !important;
    font-size: 11px !important;
}

/* ── Dropdown ── */
.wrap { border-color: #0d2030 !important; background: #060e18 !important; }
select, input {
    background: #060e18 !important;
    color: #c0d8e8 !important;
    font-family: 'JetBrains Mono', monospace !important;
    border-color: #0d2030 !important;
}

/* ── Labels ── */
label span, .label-wrap span {
    font-family: 'JetBrains Mono', monospace !important;
    font-size: 10px !important;
    letter-spacing: 0.08em !important;
    text-transform: uppercase !important;
    color: #546e7a !important;
}

/* ── Stat cards ── */
.stat-card {
    background: #060e18;
    border: 1px solid #0d2030;
    border-radius: 8px;
    padding: 16px 20px;
    text-align: center;
}
"""


# ── UI ────────────────────────────────────────────────────────────────────────
with gr.Blocks(title="KRAS Drug Finder") as demo:

    # ── Hero header ───────────────────────────────────────────────────────────
    gr.HTML("""
    <div id="hero">
      <div style="
          font-family:'Syne',sans-serif; font-size:28px; font-weight:800;
          color:#e0f0ff; letter-spacing:-0.01em; margin-bottom:6px;
      ">
        🧬 Protein Folding Drug Target Finder
      </div>
      <div style="
          font-family:'JetBrains Mono',monospace; font-size:11px;
          color:#546e7a; letter-spacing:0.12em; text-transform:uppercase;
      ">
        AlphaFold2 &nbsp;→&nbsp; fpocket &nbsp;→&nbsp; RDKit &nbsp;→&nbsp; AutoDock-GPU
        &nbsp;&nbsp;|&nbsp;&nbsp;
        Target: <span style="color:#00e5cc">KRAS G12D</span>
      </div>

      <!-- Stats row -->
      <div style="
          display:flex; justify-content:center; gap:24px; margin-top:28px;
          flex-wrap:wrap;
      ">
        <div style="background:#060e18;border:1px solid #0d2030;border-top:2px solid #00e5cc;
                    border-radius:8px;padding:14px 24px;text-align:center;min-width:110px;">
          <div style="font-family:'Syne',sans-serif;font-size:26px;font-weight:700;color:#00e5cc;">28</div>
          <div style="font-size:9px;color:#546e7a;letter-spacing:.1em;text-transform:uppercase;margin-top:2px;">Ligands Screened</div>
        </div>
        <div style="background:#060e18;border:1px solid #0d2030;border-top:2px solid #00e5cc;
                    border-radius:8px;padding:14px 24px;text-align:center;min-width:110px;">
          <div style="font-family:'Syne',sans-serif;font-size:26px;font-weight:700;color:#00e5cc;">16</div>
          <div style="font-size:9px;color:#546e7a;letter-spacing:.1em;text-transform:uppercase;margin-top:2px;">Hits Found</div>
        </div>
        <div style="background:#060e18;border:1px solid #0d2030;border-top:2px solid #00e5cc;
                    border-radius:8px;padding:14px 24px;text-align:center;min-width:110px;">
          <div style="font-family:'Syne',sans-serif;font-size:26px;font-weight:700;color:#00e5cc;">34s</div>
          <div style="font-size:9px;color:#546e7a;letter-spacing:.1em;text-transform:uppercase;margin-top:2px;">Docking Time</div>
        </div>
        <div style="background:#060e18;border:1px solid #0d2030;border-top:2px solid #ff9f43;
                    border-radius:8px;padding:14px 24px;text-align:center;min-width:130px;">
          <div style="font-family:'Syne',sans-serif;font-size:26px;font-weight:700;color:#ff9f43;">−8.72</div>
          <div style="font-size:9px;color:#546e7a;letter-spacing:.1em;text-transform:uppercase;margin-top:2px;">Best Energy (kcal/mol)</div>
        </div>
      </div>
    </div>
    """)

    # ── Tab 1: Results ────────────────────────────────────────────────────────
    with gr.Tab("📊  Results"):
        with gr.Row():
            refresh_btn = gr.Button("↻  Load Results", variant="secondary", scale=0)

        results_table = gr.Dataframe(
            label="Ranked Candidates",
            interactive=False,
            wrap=True,
        )

        gr.HTML("<div style='height:8px'></div>")

        with gr.Row():
            energy_plot = gr.Image(
                label="Binding Energy Analysis",
                value=PLOT_PATH if os.path.exists(PLOT_PATH) else None,
            )
            plddt_plot = gr.Image(
                label="AlphaFold2 pLDDT Confidence",
                value=PLDDT_PATH if os.path.exists(PLDDT_PATH) else None,
                
            )

        refresh_btn.click(fn=refresh_results, outputs=[results_table])

    # ── Tab 2: 3D Viewer ──────────────────────────────────────────────────────
    with gr.Tab("🔬  3D Viewer"):
        gr.HTML("""
        <div style="font-family:'JetBrains Mono',monospace;font-size:11px;
                    color:#546e7a;letter-spacing:.06em;padding:4px 0 12px;">
          SELECT A LIGAND TO VISUALIZE ITS DOCKED CONFORMATION INSIDE THE KRAS G12D POCKET
        </div>""")

        df_init = load_results()
        init_choices = (
            [f"{r['mol_name']} ({r['binding_energy_kcal_mol']:.2f} kcal/mol)"
             for _, r in df_init.iterrows()]
            if df_init is not None else []
        )

        with gr.Row():
            ligand_dropdown = gr.Dropdown(
                choices=init_choices,
                value=init_choices[0] if init_choices else None,
                label="Compound",
                scale=3,
            )
            view_btn = gr.Button("▶  Show 3D Pose", variant="primary", scale=1)

        info_card   = gr.HTML()
        viewer_html = gr.HTML(
            "<div style='padding:60px;color:#0d2030;text-align:center;"
            "font-family:JetBrains Mono,monospace;font-size:11px;letter-spacing:.1em;"
            "border:1px solid #0d2030;border-radius:8px;margin-top:8px;'>"
            "SELECT A COMPOUND AND CLICK SHOW 3D POSE</div>"
        )

        view_btn.click(fn=show_viewer, inputs=[ligand_dropdown], outputs=[viewer_html])
        ligand_dropdown.change(fn=show_info, inputs=[ligand_dropdown], outputs=[info_card])

    # ── Tab 3: About ──────────────────────────────────────────────────────────
    with gr.Tab("ℹ️  About"):
        gr.Markdown("""
## About This Project

A full end-to-end computational drug discovery pipeline running entirely on a local consumer GPU — no cloud, no subscription.

### Pipeline

| Stage | Tool | What it does |
|---|---|---|
| 1 | **AlphaFold2** (ColabFold) | Predicts 3D protein structure from amino acid sequence using Evoformer transformer |
| 2 | **fpocket** | Detects druggable binding cavities using Voronoi tessellation |
| 3 | **RDKit** | Filters ligands by Lipinski Rule of Five, generates 3D conformers |
| 4 | **AutoDock-GPU** | Docks all ligands in parallel across CUDA cores |
| 5 | **pandas + matplotlib** | Ranks, filters, and visualizes results |

### Target: KRAS G12D

KRAS G12D is the most common oncogenic mutation in human cancer (pancreatic, colorectal, lung).
Historically called "undruggable" until **Sotorasib (AMG-510)** received FDA approval in 2021.

The top hit — AMG_510_analog at -8.72 kcal/mol — belongs to the same scaffold class as Sotorasib,
validating the docking setup against known biology.

### Binding Energy Reference

| Energy | Interpretation |
|---|---|
| > −5 kcal/mol | Weak |
| −5 to −7 kcal/mol | Moderate |
| −7 to −9 kcal/mol | **Strong hit** |
| < −9 kcal/mol | Excellent |

### Hardware
Tested on NVIDIA GeForce RTX 3050 Laptop GPU (4GB VRAM) · WSL2 + Ubuntu 24
""")


# ── Launch ────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("Starting KRAS demo UI...")
    print("Open: http://localhost:7860")
    demo.launch(
    server_name="0.0.0.0",
    server_port=7860,
    share=False,
    show_error=True,
    css=css,                          # ← add
    theme=gr.themes.Base(             # ← add
        primary_hue="teal",
        neutral_hue="slate",
        font=gr.themes.GoogleFont("JetBrains Mono"),
    ),)
