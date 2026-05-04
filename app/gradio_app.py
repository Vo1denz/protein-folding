# app/gradio_app.py
# Run with: python app/gradio_app.py
# Opens at: http://localhost:7860

import gradio as gr
import pandas as pd
import subprocess
import os
import re
import glob
import base64
import shutil
import sys
import numpy as np
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────
RESULTS_CSV    = "data/results/ranked_candidates.csv"
DOCKING_DIR    = "data/docking/results"
STRUCTURES_DIR = "data/structures/kras"
PDBQT_DIR      = "data/inputs/ligands/prepared/pdbqt"
PLOT_PATH      = "outputs/docking_results.png"
PLDDT_PATH     = "outputs/example_plddt_plot.png"
AUTODOCK_BIN   = "/mnt/d/projects/Protein_Folding/protein-drug-finder/AutoDock-GPU/bin/autodock_gpu_128wi"
DAT_FILE       = "/mnt/d/projects/Protein_Folding/protein-drug-finder/AutoDock-GPU/input/1ac8/derived/AD4.1_bound.dat"

os.makedirs("app", exist_ok=True)


# ── WSL helpers ───────────────────────────────────────────────────────────────
def to_wsl_path(p):
    """Convert Windows path (C:\\foo\\bar) to WSL path (/mnt/c/foo/bar)."""
    p = str(p).replace("\\", "/")
    if sys.platform == "win32" and len(p) > 1 and p[1] == ":":
        p = "/mnt/" + p[0].lower() + p[2:]
    return p

def wsl_cmd(cmd):
    """Prefix command with wsl on Windows, pass through on Linux."""
    if sys.platform == "win32":
        return ["wsl"] + [
            to_wsl_path(c) if (os.sep in str(c) or (len(str(c)) > 1 and str(c)[1] == ':'))
            else str(c) for c in cmd
        ]
    return [str(c) for c in cmd]


# ── Streaming subprocess ───────────────────────────────────────────────────────
def stream_run(cmd, label, log_lines):
    """
    Run a command and stream output line by line via generator.

    Why this matters for AlphaFold2:
    - ColabFold takes 20-60 minutes. With capture_output=True the UI shows
      nothing until the process finishes — looks completely frozen.
    - subprocess.Popen with stdout=PIPE reads lines as they arrive.
    - stderr=STDOUT merges both streams into one so nothing is lost.
    - bufsize=1 = line-buffered: each newline flushes immediately.
    - This function is a generator — every yield sends a UI update instantly.
    """
    log_lines.append(f"\n{'='*55}")
    log_lines.append(f"▶  {label}")
    log_lines.append(f"{'='*55}")
    yield "\n".join(log_lines), None, None, None

    process = subprocess.Popen(
        wsl_cmd(cmd),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )

    for line in process.stdout:
        line = line.rstrip()
        if line:
            log_lines.append(line)
            yield "\n".join(log_lines[-60:]), None, None, None

    process.wait()
    log_lines.append(f"[exit code: {process.returncode}]")
    yield "\n".join(log_lines[-60:]), None, None, None


def run_silent(cmd):
    """Run a command silently, return (returncode, stderr)."""
    result = subprocess.run(wsl_cmd(cmd), capture_output=True, text=True)
    return result.returncode, result.stderr


# ── Helpers ───────────────────────────────────────────────────────────────────
def load_results():
    if os.path.exists(RESULTS_CSV):
        return pd.read_csv(RESULTS_CSV)
    return None

def get_docked_pose(ligand_id):
    """Extract best docked pose from DLG file as PDB-format string."""
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
    """
    Build iframe-based 3Dmol.js viewer.
    iframe is required because Gradio 4.x strips <script> tags from gr.HTML().
    An iframe srcdoc gets its own document — scripts run freely inside it.
    """
    pose_pdb = get_docked_pose(ligand_id)

    receptor_pdb = ""
    receptor_pdb_path = os.path.join(STRUCTURES_DIR,
        "KRAS_G12D_Human_KRAS_proto-oncogene_GTPase_G12D_mutant_unrelaxed_rank_001_alphafold2_model_1_seed_000.pdb")
    if os.path.exists(receptor_pdb_path):
        with open(receptor_pdb_path) as f:
            receptor_pdb = f.read()

    if not pose_pdb:
        return f"<div style='padding:20px;color:#e74c3c'>No docked pose found for {mol_name}</div>"

    receptor_escaped = receptor_pdb.replace('\\', '\\\\').replace('`', '\\`').replace('${', '\\${')
    pose_escaped     = pose_pdb.replace('\\', '\\\\').replace('`', '\\`').replace('${', '\\${')

    inner_html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<style>
  * {{ margin:0; padding:0; box-sizing:border-box; }}
  body {{ background:#0a0e1a; font-family:monospace; color:#ccc; }}
  #header {{ padding:8px 12px; font-size:12px; color:#00d4aa; border-bottom:1px solid #1a2a3a; }}
  #viewer {{ width:100%; height:380px; }}
  #footer {{ padding:6px 12px; font-size:11px; color:#556; }}
</style>
</head>
<body>
<div id="header">
  ◆ {mol_name} &nbsp;|&nbsp; Binding energy: <strong>{energy:.2f} kcal/mol</strong>
</div>
<div id="viewer"></div>
<div id="footer">Protein: teal cartoon + transparent surface &nbsp;|&nbsp; Ligand: element-colored sticks</div>
<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.6/3Dmol-min.js"></script>
<script>
window.addEventListener('load', function() {{
  var viewer = $3Dmol.createViewer(document.getElementById('viewer'), {{ backgroundColor: '#0a0e1a' }});
  var receptorData = `{receptor_escaped}`;
  if (receptorData.trim().length > 0) {{
    viewer.addModel(receptorData, 'pdb');
    viewer.setStyle({{ model: 0 }}, {{ cartoon: {{ color: '#00d4aa', opacity: 0.9 }} }});
    viewer.addSurface($3Dmol.SurfaceType.VDW, {{ opacity: 0.07, color: '#00d4aa' }}, {{ model: 0 }});
  }}
  var ligandData = `{pose_escaped}`;
  viewer.addModel(ligandData, 'pdbqt');
  viewer.setStyle({{ model: 1 }}, {{
    stick:  {{ colorscheme: 'elementColors', radius: 0.25 }},
    sphere: {{ colorscheme: 'elementColors', scale: 0.35 }}
  }});
  viewer.zoomTo({{ model: 1 }});
  viewer.zoom(0.75);
  viewer.render();
}});
</script>
</body>
</html>"""

    srcdoc = inner_html.replace('"', '&quot;')
    return f"""<iframe srcdoc="{srcdoc}"
      style="width:100%; height:440px; border:1px solid #1a2a3a;
             border-radius:8px; background:#0a0e1a;"
      sandbox="allow-scripts"></iframe>"""


# ── Pipeline (generator — streams output line by line) ────────────────────────
def run_pipeline(fasta_text, progress=gr.Progress()):
    """
    Generator function — every yield sends an update to Gradio immediately.
    Gradio detects generators automatically and streams each value to the UI.
    Yields: (status_text, results_df, energy_plot, plddt_plot)
    """
    log_lines = []

    def log(msg):
        log_lines.append(msg)
        return "\n".join(log_lines[-60:])

    if not fasta_text.strip():
        yield "Please enter a FASTA sequence.", None, None, None
        return

    lines        = fasta_text.strip().splitlines()
    header       = lines[0] if lines[0].startswith(">") else ">protein"
    protein_name = re.sub(r'[^\w]', '_', header[1:].split()[0])[:20]

    yield log(f"Starting pipeline for: {protein_name}"), None, None, None

    seq_dir     = Path("data/inputs/sequences")
    struct_dir  = Path(f"data/structures/{protein_name}")
    docking_dir = Path("data/docking")
    results_dir = Path("data/results")
    for d in [seq_dir, struct_dir, docking_dir, results_dir, Path(DOCKING_DIR)]:
        d.mkdir(parents=True, exist_ok=True)

    fasta_path = seq_dir / f"{protein_name}.fasta"
    fasta_path.write_text(fasta_text.strip() + "\n")
    yield log(f"FASTA saved → {fasta_path}"), None, None, None

    # ── Stage 1: ColabFold ────────────────────────────────────────────────────
    progress(0.05, desc="Stage 1/5: AlphaFold2 folding...")
    for status, *_ in stream_run([
        "colabfold_batch", "--num-recycle", "3",
        "--model-type", "alphafold2_ptm", "--use-gpu",
        to_wsl_path(fasta_path), to_wsl_path(struct_dir),
    ], "Stage 1/5 — AlphaFold2 / ColabFold", log_lines):
        yield status, None, None, None

    if "[exit code: " in log_lines[-1] and "[exit code: 0]" not in log_lines[-1]:
        yield log("❌ AlphaFold2 failed. See output above."), None, None, None
        return

    pdb_files = sorted(struct_dir.glob("*_relaxed_rank_001*.pdb"))
    if not pdb_files:
        pdb_files = sorted(struct_dir.glob("*rank_001*.pdb"))
    if not pdb_files:
        yield log("❌ No PDB file produced by AlphaFold2."), None, None, None
        return
    best_pdb = pdb_files[0]
    yield log(f"✓ Structure: {best_pdb.name}"), None, None, None

    # ── Stage 2: fpocket ──────────────────────────────────────────────────────
    progress(0.35, desc="Stage 2/5: Detecting binding pockets...")
    yield log("\n" + "="*55 + "\n▶  Stage 2/5 — fpocket\n" + "="*55), None, None, None

    work_pdb = struct_dir / "receptor_input.pdb"
    shutil.copy(best_pdb, work_pdb)

    rc, err = run_silent(["fpocket", "-f", to_wsl_path(work_pdb)])
    yield log(f"fpocket exit code: {rc}"), None, None, None
    if rc != 0:
        yield log(f"❌ fpocket failed:\n{err}"), None, None, None
        return

    out_dir = work_pdb.parent / (work_pdb.stem + "_out") / "pockets"
    pocket_files = sorted(out_dir.glob("pocket*_atm.pdb"),
                          key=lambda f: int(re.search(r'pocket(\d+)', f.name).group(1)))
    if not pocket_files:
        yield log("❌ fpocket found no pockets."), None, None, None
        return

    yield log(f"✓ Found {len(pocket_files)} pockets"), None, None, None

    best_center, best_count = None, 0
    for pf in pocket_files:
        coords = []
        for line in pf.read_text().splitlines():
            if line.startswith(("ATOM", "HETATM")):
                try:
                    coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                except ValueError:
                    continue
        if len(coords) > best_count:
            best_count = len(coords)
            best_center = np.mean(coords, axis=0)

    cx, cy, cz = best_center
    yield log(f"✓ Best pocket center: ({cx:.2f}, {cy:.2f}, {cz:.2f})"), None, None, None

    # ── Stage 3: Prepare receptor PDBQT ──────────────────────────────────────
    progress(0.45, desc="Stage 3/5: Preparing receptor PDBQT...")
    yield log("\n" + "="*55 + "\n▶  Stage 3/5 — Receptor PDBQT\n" + "="*55), None, None, None

    clean_pdb      = struct_dir / "receptor_clean.pdb"
    receptor_pdbqt = struct_dir / "receptor.pdbqt"

    clean_lines = [l for l in best_pdb.read_text().splitlines(keepends=True)
                   if l.startswith(("ATOM", "TER", "END"))]
    clean_pdb.write_text("".join(clean_lines))
    yield log(f"✓ Cleaned PDB: {len(clean_lines)} lines"), None, None, None

    run_silent(["obabel", to_wsl_path(clean_pdb), "-O", to_wsl_path(receptor_pdbqt),
                "--partialcharge", "gasteiger", "-xr", "-h"])
    atom_count = sum(1 for l in receptor_pdbqt.read_text().splitlines() if l.startswith("ATOM"))
    yield log(f"✓ PDBQT ready: {atom_count} ATOM records"), None, None, None

    # ── Stage 4: AutoGrid ─────────────────────────────────────────────────────
    progress(0.55, desc="Stage 4/5: AutoGrid + docking...")
    yield log("\n" + "="*55 + "\n▶  Stage 4/5 — AutoGrid4 + AutoDock-GPU\n" + "="*55), None, None, None

    receptor_base = to_wsl_path(receptor_pdbqt.resolve()).replace(".pdbqt", "")

    atom_types_receptor = set()
    for line in receptor_pdbqt.read_text().splitlines():
        if line.startswith("ATOM"):
            at = line[77:79].strip()
            if at and at != "Ty":
                atom_types_receptor.add(at)

    all_atom_types = list(atom_types_receptor) + [t for t in ["Cl", "F"] if t not in atom_types_receptor]
    map_lines = "\n".join(f"map {receptor_base}.{t}.map" for t in all_atom_types)

    gpf_content = f"""npts 40 40 40
parameter_file {DAT_FILE}
gridfld {receptor_base}.maps.fld
spacing 0.375
receptor_types {" ".join(atom_types_receptor)}
ligand_types {" ".join(all_atom_types)}
receptor {to_wsl_path(receptor_pdbqt.resolve())}
gridcenter {cx:.3f} {cy:.3f} {cz:.3f}
smooth 0.5
{map_lines}
elecmap {receptor_base}.e.map
dsolvmap {receptor_base}.d.map
dielectric -0.1465
"""
    gpf_path = docking_dir / "grid.gpf"
    glg_path = docking_dir / "grid.glg"
    gpf_path.write_text(gpf_content)
    yield log(f"✓ GPF written — box 40³, center ({cx:.2f},{cy:.2f},{cz:.2f})"), None, None, None

    for status, *_ in stream_run(
        ["autogrid4", "-p", to_wsl_path(gpf_path), "-l", to_wsl_path(glg_path)],
        "AutoGrid4 — generating energy maps", log_lines
    ):
        yield status, None, None, None

    if "[exit code: " in log_lines[-1] and "[exit code: 0]" not in log_lines[-1]:
        yield log("❌ AutoGrid4 failed."), None, None, None
        return

    # Fix doubled paths in fld
    fld_path = receptor_pdbqt.with_suffix("").resolve().parent / (receptor_pdbqt.stem + ".maps.fld")
    if fld_path.exists():
        fld_text = fld_path.read_text()
        fld_text = re.sub(
            r'(/mnt/\w+)?/mnt/\w+[^\s]*/data/structures/',
            lambda m: to_wsl_path(struct_dir) + "/",
            fld_text
        )
        fld_path.write_text(fld_text)

    map_count = len(list(receptor_pdbqt.parent.glob("*.map")))
    yield log(f"✓ {map_count} grid maps generated"), None, None, None

    # ── Docking ───────────────────────────────────────────────────────────────
    progress(0.65, desc="Stage 4/5: Docking ligands...")
    ligand_files    = sorted(Path(PDBQT_DIR).glob("*.pdbqt"))
    results_docking = Path(DOCKING_DIR)
    total   = len(ligand_files)
    success = 0

    yield log(f"\nDocking {total} ligands against {protein_name}..."), None, None, None

    for i, lf in enumerate(ligand_files, 1):
        out_base = results_docking / lf.stem
        run_silent([
            AUTODOCK_BIN,
            "--ffile", to_wsl_path(fld_path),
            "--lfile", to_wsl_path(lf),
            "--nrun", "20",
            "--resnam", to_wsl_path(out_base),
            "--gbest", "1", "--lsmet", "sw", "--seed", "1234",
        ])
        dlg = results_docking / (lf.stem + ".dlg")
        if dlg.exists() and dlg.stat().st_size > 0:
            energy_str = ""
            for line in dlg.read_text().splitlines():
                if "Estimated Free Energy of Binding" in line:
                    m = re.search(r'=\s*([-\d.]+)', line)
                    if m:
                        energy_str = f"{float(m.group(1)):.2f} kcal/mol"
                    break
            log_lines.append(f"  [{i:2d}/{total}] {lf.stem:20s} ✓  {energy_str}")
            success += 1
        else:
            log_lines.append(f"  [{i:2d}/{total}] {lf.stem:20s} ✗")
        yield "\n".join(log_lines[-60:]), None, None, None

    yield log(f"\n✓ Docking complete: {success}/{total} succeeded"), None, None, None

    if success == 0:
        yield log("❌ All docking runs failed."), None, None, None
        return

    # ── Stage 5: Parse results ────────────────────────────────────────────────
    progress(0.90, desc="Stage 5/5: Ranking results...")
    for status, *_ in stream_run(
        ["python", "scripts/05_parse_results.py"],
        "Stage 5/5 — Ranking candidates", log_lines
    ):
        yield status, None, None, None

    progress(1.0, desc="Done!")

    df = load_results()
    if df is None:
        yield log("⚠ Pipeline complete but no results CSV found."), None, None, None
        return

    display_df = df[['rank', 'mol_name', 'binding_energy_kcal_mol',
                      'mol_weight', 'logp', 'tpsa', 'ligand_efficiency', 'category']].head(20)
    n_hits = len(df[df['category'].isin(['hit', 'strong hit', 'excellent'])])
    best   = df.iloc[0]

    log_lines.append("\n" + "="*55)
    log_lines.append("✅ PIPELINE COMPLETE")
    log_lines.append(f"   Protein       : {protein_name}")
    log_lines.append(f"   Ligands docked: {success}/{total}")
    log_lines.append(f"   Hits (<-7)    : {n_hits}")
    log_lines.append(f"   Best hit      : {best['mol_name']}  {best['binding_energy_kcal_mol']:.2f} kcal/mol")
    log_lines.append(f"   Pocket center : ({cx:.2f}, {cy:.2f}, {cz:.2f})")
    log_lines.append("="*55)

    plot  = PLOT_PATH if os.path.exists(PLOT_PATH) else None
    plddt = PLDDT_PATH if os.path.exists(PLDDT_PATH) else None
    yield "\n".join(log_lines[-60:]), display_df, plot, plddt


# ── Viewer helpers ────────────────────────────────────────────────────────────
def show_viewer(selection):
    if not selection:
        return "<p style='color:gray'>Select a ligand above to view its docking pose.</p>"
    m = re.match(r'^(.+?)\s*\((-?[\d.]+)', selection)
    if not m:
        return "<p style='color:red'>Could not parse selection.</p>"
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
        return f"<p style='color:red'>Could not find PDBQT for {mol_name}</p>"
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
    hit_color = {"excellent": "#e74c3c", "strong hit": "#e67e22",
                 "hit": "#3498db", "weak": "#95a5a6"}.get(r['category'], 'gray')
    return f"""
    <div style='display:flex; gap:20px; flex-wrap:wrap; margin-top:10px;
                font-family:monospace; font-size:13px;'>
      <div style='background:#f8f9fa; padding:10px 16px; border-radius:6px;
                  border-left:4px solid {hit_color};'>
        <strong>{r['mol_name']}</strong><br>
        Category: <span style='color:{hit_color}'>{r['category']}</span>
      </div>
      <div style='background:#f8f9fa; padding:10px 16px; border-radius:6px;'>
        Energy: <strong>{r['binding_energy_kcal_mol']:.2f} kcal/mol</strong><br>
        Lig. Efficiency: {r['ligand_efficiency']:.3f}
      </div>
      <div style='background:#f8f9fa; padding:10px 16px; border-radius:6px;'>
        MW: {r['mol_weight']:.1f} Da &nbsp;|&nbsp; LogP: {r['logp']:.2f}<br>
        TPSA: {r['tpsa']:.1f} Å²
      </div>
    </div>"""


def refresh_results():
    df = load_results()
    if df is None:
        return None
    return df[['rank', 'mol_name', 'binding_energy_kcal_mol',
               'mol_weight', 'logp', 'tpsa', 'ligand_efficiency', 'category']]


# ── Gradio UI ─────────────────────────────────────────────────────────────────
css = """
.gradio-container { max-width: 1200px !important; }
#title { text-align: center; padding: 20px 0 10px; }
#title h1 { font-size: 2em; margin-bottom: 4px; }
#status_box textarea { font-family: monospace; font-size: 12px; }
"""

with gr.Blocks(title="Protein Drug Finder", css=css, theme=gr.themes.Soft()) as demo:

    gr.HTML("""
    <div id='title'>
      <h1>🧬 Protein Folding Drug Target Finder</h1>
      <p style='color:#666; font-size:15px;'>
        AlphaFold2 → fpocket → RDKit → AutoDock-GPU
      </p>
    </div>""")

    # ── Tab 1: Run Pipeline ───────────────────────────────────────────────────
    with gr.Tab("🚀 Run Pipeline"):
        gr.Markdown("""Enter a protein sequence in FASTA format and click **Run Pipeline**.
        AlphaFold2 takes 20–60 min — you'll see live output as it runs.
        For a quick demo use the **Results** and **3D Viewer** tabs.""")

        with gr.Row():
            with gr.Column(scale=1):
                fasta_input = gr.Textbox(
                    label="Protein sequence (FASTA format)",
                    placeholder=">BCL2_HUMAN\nMAHAGRTGYDNREIVMKYIHYKLSQRGYEWDAGDVGAAPPGAAP...",
                    lines=8,
                )
                run_btn = gr.Button("▶ Run Full Pipeline", variant="primary", size="lg")

            with gr.Column(scale=1):
                gr.Markdown("""**What each stage does:**
| Stage | Tool | Time |
|---|---|---|
| 1. Fold protein | AlphaFold2 (ColabFold) | 20–60 min |
| 2. Find pockets | fpocket | < 1 min |
| 3. Prep receptor | Open Babel | < 1 min |
| 4. Dock ligands | AutoDock-GPU | 1–5 min |
| 5. Rank results | pandas | < 1 min |""")

        # Live log box — streams ColabFold epochs, docking scores etc in real time
        status_box = gr.Textbox(
            label="Live pipeline output",
            interactive=False,
            lines=20,
            max_lines=40,
            elem_id="status_box",
        )

    # ── Tab 2: Results ────────────────────────────────────────────────────────
    with gr.Tab("📊 Results"):
        gr.Markdown("### Ranked drug candidates")
        refresh_btn   = gr.Button("↻ Load / Refresh Results", variant="secondary")
        results_table = gr.Dataframe(label="All ligands ranked by binding energy",
                                     interactive=False, wrap=True)
        with gr.Row():
            energy_plot = gr.Image(label="Binding energy distribution",
                                   value=PLOT_PATH if os.path.exists(PLOT_PATH) else None)
            plddt_plot  = gr.Image(label="AlphaFold2 pLDDT confidence",
                                   value=PLDDT_PATH if os.path.exists(PLDDT_PATH) else None)

        refresh_btn.click(fn=refresh_results, outputs=[results_table])

    # ── Tab 3: 3D Viewer ──────────────────────────────────────────────────────
    with gr.Tab("🔬 3D Viewer"):
        gr.Markdown("""### Docked binding pose viewer
Select a ligand to see how it fits inside the binding pocket.
Protein: teal cartoon + transparent surface. Ligand: element-colored sticks.""")

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
                label="Select ligand", scale=2,
            )
            view_btn = gr.Button("Show 3D pose", variant="primary", scale=1)

        viewer_html = gr.HTML(
            "<div style='padding:20px;color:#888;text-align:center'>"
            "Select a ligand and click Show 3D pose</div>")
        info_card = gr.HTML()

        view_btn.click(fn=show_viewer, inputs=[ligand_dropdown], outputs=[viewer_html])
        ligand_dropdown.change(fn=show_info, inputs=[ligand_dropdown], outputs=[info_card])

    # ── Tab 4: About ──────────────────────────────────────────────────────────
    with gr.Tab("ℹ️ About"):
        gr.Markdown("""
## Protein Folding Drug Target Finder

A full end-to-end computational drug discovery pipeline built on a consumer GPU.

### Pipeline
1. **AlphaFold2** (ColabFold) — predicts 3D protein structure from amino acid sequence
2. **fpocket** — detects druggable cavities on the protein surface using Voronoi tessellation
3. **RDKit** — filters ligands by Lipinski Rule of Five, generates 3D conformers
4. **AutoDock-GPU** — screens all ligands against the binding pocket in parallel on CUDA
5. **pandas + matplotlib** — ranks, filters, and visualizes results

### KRAS G12D Results
- 28 ligands screened in 34 seconds on RTX 3050
- 16 hits identified below -7 kcal/mol threshold
- Top hit: AMG_510_analog at -8.72 kcal/mol (Sotorasib scaffold class)

### Binding Energy Reference
| Energy | Interpretation |
|---|---|
| > -5 kcal/mol | Weak — unlikely drug candidate |
| -5 to -7 kcal/mol | Moderate — worth investigating |
| -7 to -9 kcal/mol | Strong — good drug candidate |
| < -9 kcal/mol | Excellent — high priority hit |""")

    # ── Wire run button → all output components ───────────────────────────────
    run_btn.click(
        fn=run_pipeline,
        inputs=[fasta_input],
        outputs=[status_box, results_table, energy_plot, plddt_plot],
    )

# ── Launch ────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("Starting Gradio app...")
    print("Open in browser: http://localhost:7860")
    demo.launch(server_name="0.0.0.0", server_port=7860, share=False, show_error=True)