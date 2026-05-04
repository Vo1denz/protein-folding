# Save as: verify_all.py (in your protein-drug-finder/ folder)
import subprocess
import sys
import os

results = {}

def check(name, fn):
    try:
        fn()
        results[name] = "✅ PASS"
    except Exception as e:
        results[name] = f"❌ FAIL — {e}"

# 1. JAX + GPU
def check_jax():
    import jax
    devices = jax.devices()
    gpu = [d for d in devices if 'cuda' in str(d).lower() or d.device_kind == 'gpu']
    assert gpu, f"No GPU in jax.devices(): {devices}"

# 2. RDKit
def check_rdkit():
    from rdkit import Chem
    mol = Chem.MolFromSmiles("CCO")  # ethanol
    assert mol is not None

# 3. fpocket
def check_fpocket():
    r = subprocess.run(["fpocket", "--version"], capture_output=True, text=True)
    assert r.returncode == 0 or "fpocket" in r.stderr.lower() or "fpocket" in r.stdout.lower()

# 4. AutoDock-GPU binary
def check_autodock():
    binary = "AutoDock-GPU/bin/autodock_gpu_128wi"
    assert os.path.exists(binary), f"Binary not found at {binary}"
    r = subprocess.run([binary, "--help"], capture_output=True, text=True)
    assert r.returncode in (0, 1)  # --help often returns 1

# 5. Open Babel
def check_obabel():
    r = subprocess.run(["obabel", "--version"], capture_output=True, text=True)
    assert "Open Babel" in r.stdout or "Open Babel" in r.stderr

# 6. Gradio
def check_gradio():
    import gradio as gr
    assert gr.__version__ >= "4.0"

# 7. Folder structure
def check_folders():
    required = [
        "data/inputs/sequences", "data/inputs/ligands",
        "data/structures", "data/pockets",
        "data/docking", "data/results",
        "scripts", "app", "outputs"
    ]
    missing = [f for f in required if not os.path.isdir(f)]
    assert not missing, f"Missing folders: {missing}"

# Run all checks
check("1. JAX + GPU",          check_jax)
check("2. RDKit",              check_rdkit)
check("3. fpocket",            check_fpocket)
check("4. AutoDock-GPU binary",check_autodock)
check("5. Open Babel",         check_obabel)
check("6. Gradio >= 4.0",      check_gradio)
check("7. Folder structure",   check_folders)

# Report
print("\n" + "="*45)
print("  Phase 0 — Toolchain Verification Report")
print("="*45)
for name, status in results.items():
    print(f"  {status}  {name}")

failed = [k for k, v in results.items() if "FAIL" in v]
print("="*45)
if not failed:
    print("  🎉 All checks passed! Ready for Phase 1.")
else:
    print(f"  ⚠️  {len(failed)} check(s) failed. Fix before proceeding.")
print("="*45 + "\n")