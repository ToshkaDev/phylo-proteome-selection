#!/usr/bin/env python3
"""Run the complete GTDB r232 selection workflow from the package root."""
from pathlib import Path
import subprocess, sys

ROOT = Path(__file__).resolve().parent
STEPS = [
    ROOT / "scripts" / "reconstruct_phylogeny.py",
    ROOT / "scripts" / "select_baseline.py",
    ROOT / "scripts" / "include_model_organisms.py",
    ROOT / "scripts" / "add_aliases_and_mk_d1.py",
]
for step in STEPS:
    print(f"\n=== {step.name} ===", flush=True)
    subprocess.run([sys.executable, str(step)], cwd=ROOT, check=True)
subprocess.run([sys.executable, str(ROOT / "validate.py")], cwd=ROOT, check=True)
