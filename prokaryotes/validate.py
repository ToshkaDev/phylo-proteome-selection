#!/usr/bin/env python3
"""Small set of release checks for the generated genome-selection tables."""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent
BASELINE = ROOT / "outputs" / "baseline" / "selected_genomes_3424_phylo_coverage.tsv"
MODEL = ROOT / "outputs" / "model_enriched" / "selected_genomes_model_organisms.tsv"
FINAL = ROOT / "outputs" / "model_enriched" / "final_selected_genomes.tsv"
AUDIT = ROOT / "outputs" / "model_enriched" / "model_organism_inclusion_audit.tsv"
EXPLICIT_AUDIT = ROOT / "outputs" / "model_enriched" / "explicit_model_additions_audit.tsv"

def read(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing expected output: {path.relative_to(ROOT)}")
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)

checks = []
def check(label, condition):
    if not condition:
        raise AssertionError(label)
    checks.append(label)

baseline, model, final = read(BASELINE), read(MODEL), read(FINAL)
audit, explicit = read(AUDIT), read(EXPLICIT_AUDIT)
check("baseline has exactly 3,424 rows and unique accessions", len(baseline) == 3424 and baseline.accession.nunique() == 3424)
check("baseline domain split is 3,182 Bacteria / 242 Archaea", (baseline.domain == "Bacteria").sum() == 3182 and (baseline.domain == "Archaea").sum() == 242)
check("model-aware table has unique accessions", model.accession.nunique() == len(model))
check("final table has unique accessions", final.accession.nunique() == len(final))
check("all four explicit model accessions are present", {"RS_GCF_000008685.2", "RS_GCF_000091665.1", "RS_GCF_900079115.1", "RS_GCF_008000775.2"}.issubset(set(final.accession)))
check("MK-D1 is marked as a quality exception", ((explicit.accession == "RS_GCF_008000775.2") & (explicit.quality_exception.str.lower() == "true")).any())
check("five intentionally ignored targets are absent from active audit", not any(audit.model_organism.str.contains(name, case=False, regex=False).any() for name in ["Caulobacter crescentus", "Mycoplasma genitalium", "Mycoplasma pneumoniae", "Streptomyces coelicolor", "Candidatus Carsonella ruddii"]))
print(f"PASS: {len(checks)}/{len(checks)} checks")
for i, label in enumerate(checks, 1): print(f"  {i}. {label}")
print(f"Counts: baseline={len(baseline):,}; model-aware={len(model):,}; final={len(final):,}")
