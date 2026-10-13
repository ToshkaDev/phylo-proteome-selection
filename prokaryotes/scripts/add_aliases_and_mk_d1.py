#!/usr/bin/env python3
"""Add GTDB r232 renamed model taxa and MK-D1 quality exception to model-aware selection.
Run after select_model_organisms.py. Five unavailable targets are intentionally ignored.
"""
from pathlib import Path
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
OUTDIR = BASE / 'outputs' / 'model_enriched'
META = BASE / 'inputs' / 'metadata_r232.tsv'
INPUT = OUTDIR / 'selected_genomes_model_organisms.tsv'
OUTPUT = OUTDIR / 'final_selected_genomes.tsv'
AUDIT = OUTDIR / 'explicit_model_additions_audit.tsv'
meta = pd.read_csv(META, sep='\t', dtype=str, low_memory=False).fillna('')
selected = pd.read_csv(INPUT, sep='\t', dtype=str, low_memory=False).fillna('')

# Current GTDB r232 names for user-requested historical model names.
targets = [
    ('Borrelia burgdorferi B31', 'RS_GCF_000008685.2', False, 'Borreliella burgdorferi; strain B31'),
    ('Methanococcus jannaschii', 'RS_GCF_000091665.1', False, 'Methanocaldococcus jannaschii'),
    ('Sulfolobus solfataricus', 'RS_GCF_900079115.1', False, 'Saccharolobus solfataricus'),
    ('Candidatus Prometheoarchaeum syntrophicum MK-D1', 'RS_GCF_008000775.2', True, 'Prometheoarchaeum syntrophicum; strain MK-D1'),
]
idx = meta.set_index('accession', drop=False)
records=[]
for label, acc, exception, match_note in targets:
    if acc not in idx.index:
        raise ValueError(f'Required accession absent from supplied metadata: {acc} ({label})')
    row=idx.loc[acc]
    if isinstance(row, pd.DataFrame): row=row.iloc[0]
    comp=float(row.get('checkm2_completeness','nan')); cont=float(row.get('checkm2_contamination','nan'))
    if comp < 80 or (cont > 5 and not exception):
        raise ValueError(f'Unexpected quality for {label}: completeness={comp}, contamination={cont}')
    already=(selected.accession == acc).any()
    action='already present; no duplicate added' if already else 'added as separate model-organism genome'
    if not already:
        rec={c: row[c] if c in row.index else '' for c in selected.columns}
        tax={}
        for part in str(row.get('gtdb_taxonomy','')).split(';'):
            if '__' in part:
                k,v=part.split('__',1); tax[k]=v
        rec.update({
            'domain': tax.get('d','').replace('d__',''),
            'order': tax.get('o','').replace('o__',''),
            'family': tax.get('f','').replace('f__',''),
            'selected_family': tax.get('f','').replace('f__',''),
            'genome_species': tax.get('s','').replace('s__',''),
            'genus': tax.get('g','').replace('g__',''),
            'selection_method_reconstructed': 'model-organism inclusion; explicit MK-D1 quality exception' if exception else 'model-organism inclusion; current GTDB name/alias match',
            'genome_choice_method': 'explicit model-organism inclusion',
            'genome_choice_reason': 'MK-D1 explicit quality exception' if exception else 'current GTDB taxonomy alias match',
            '_forced_model': 'true', '_model_label': label,
            '_match_reason': match_note,
            'quality_rank': 'quality exception' if exception else 'passes 80/5',
        })
        selected=pd.concat([selected, pd.DataFrame([rec], columns=selected.columns)], ignore_index=True)
    records.append({'model_organism':label,'accession':acc,'match':match_note,'completeness':comp,'contamination':cont,
                    'quality_exception':exception,'action':action})
selected=selected.drop_duplicates('accession', keep='first')
if selected.accession.nunique()!=len(selected): raise RuntimeError('Accession uniqueness validation failed')
selected.to_csv(OUTPUT, sep='\t', index=False)
pd.DataFrame(records).to_csv(AUDIT, sep='\t', index=False)
print(f'Wrote {len(selected)} unique accessions to {OUTPUT}')
print(f'Added {sum(r["action"].startswith("added") for r in records)} explicit alias/exception targets.')
print('MK-D1 is intentionally retained despite 6.22% contamination; consider it separately in quality-sensitive analyses.')
