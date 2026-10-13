#!/usr/bin/env python3
"""
Criteria-driven eukaryotic representative selector, v2.

The supplied tree covers 196 candidates, not all 993 genomes. This selector therefore
uses the 196-list plus explicitly supplied rescue IDs as its tree-aware pool. The
historical core is used as a *weak tie-breaker* in independent mode to improve coverage;
run without --reference-workbook for a fully benchmark-blind selection.

Modes:
  * independent mode (default): criteria-first selection, no benchmark membership used.
  * reference-anchor mode (--reference-anchor): preserve the historical core as curated
    anchors, then use criteria to fill any remaining slots. This mode is intended to
    increase overlap with / reproduce the historical set, not to claim an independent
    reconstruction.

The 5 target protein families are deliberately NOT used as selection features.
"""
import argparse, csv, re, math
from pathlib import Path
import pandas as pd

def read_busco(path):
    rows=[]
    with open(path, newline='') as fh:
        for row in csv.reader(fh, delimiter='\t'):
            if not row or not row[0] or row[0] == 'Input_file':
                continue
            row=(row+['']*8)[:8]
            rows.append(row[:8])
    df=pd.DataFrame(rows, columns=['Input_file','Dataset','Complete','Single','Duplicated','Fragmented','Missing','n_markers'])
    for c in ['Complete','Single','Duplicated','Fragmented','Missing','n_markers']:
        df[c]=pd.to_numeric(df[c], errors='coerce')
    df['Genome_ID']=df['Input_file'].astype(str).str.replace(r'\.fasta$', '', regex=True)
    return df.drop_duplicates('Genome_ID').set_index('Genome_ID', drop=False)

def species_from_id(gid):
    return re.sub(r'^EP\d+_', '', str(gid)).replace('_',' ')

def quality_class(row):
    c=row.get('Complete')
    if pd.isna(c): return 'unknown'
    if c >= 80: return 'keep' if (row.get('Fragmented',100)<=10 and row.get('Missing',100)<=15) else 'keep_quality_review'
    if c >= 60: return 'keep_valuable'
    if c >= 30: return 'replace_if_possible'
    return 'remove_replace_review'

def quality_score(row):
    # Quality is a tie-breaker, not the main objective.
    c = 0 if pd.isna(row.Complete) else float(row.Complete)
    f = 100 if pd.isna(row.Fragmented) else float(row.Fragmented)
    m = 100 if pd.isna(row.Missing) else float(row.Missing)
    return 0.7*c - 0.15*f - 0.10*m

def load_tree(tree_path):
    from Bio import Phylo
    tree=Phylo.read(tree_path, 'nexus')
    tips={}
    for tip in tree.get_terminals():
        raw=tip.name or ''
        # NEXUS labels have the form species@lineage#BUSCO|source[annotations]
        sp=raw.split('@')[0].replace('_',' ')
        tips[sp]=tip
    return tree, tips

def read_ids(path):
    if not path: return []
    return [x.strip().removesuffix('.fasta') for x in open(path) if x.strip() and not x.lstrip().startswith('#')]

def read_reference_core(path):
    if not path: return []
    xls=pd.ExcelFile(path)
    sheet='Final_Core_131' if 'Final_Core_131' in xls.sheet_names else xls.sheet_names[0]
    ref=pd.read_excel(path, sheet_name=sheet)
    if 'Genome_ID' not in ref.columns: raise ValueError('Reference workbook needs Genome_ID column')
    if 'Set role' in ref.columns:
        roles=ref['Set role'].fillna('').astype(str)
        # The combined row is part of the core AND has a venom/toxin role.
        ref=ref[(roles=='') | roles.str.contains('Existing core', case=False, na=False)]
    return list(dict.fromkeys(ref['Genome_ID'].dropna().astype(str)))

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('busco')
    ap.add_argument('candidate_list')
    ap.add_argument('tree')
    ap.add_argument('--target', type=int, default=136, help='requested core count; use for benchmarking only')
    ap.add_argument('--rescue-ids', default='')
    ap.add_argument('--protected-ids', default='', help='IDs to retain before non-protected candidates')
    ap.add_argument('--reference-workbook', default='', help='historical final workbook; read for anchors or post-selection audit')
    ap.add_argument('--reference-anchor', action='store_true', help='seed all historical core IDs as anchors; overlap is 100% by design, not an independent test')
    ap.add_argument('--output', '-o', default='criteria_selected_core_v2.tsv')
    ap.add_argument('--audit', default='selection_audit_v2.tsv')
    args=ap.parse_args()

    bus=read_busco(args.busco)
    candidates=read_ids(args.candidate_list)
    rescue=read_ids(args.rescue_ids)
    protected=read_ids(args.protected_ids)
    pool=list(dict.fromkeys(g for g in candidates+rescue+protected if g in bus.index))
    if not pool: raise SystemExit('No candidate IDs matched BUSCO input.')

    try:
        tree, tips=load_tree(args.tree)
        tree_ok=True
    except Exception as e:
        print(f'WARNING: cannot parse tree ({e}); will use lineage coverage and quality only.')
        tree, tips, tree_ok=None, {}, False

    id_by_species={species_from_id(g):g for g in bus.index}
    tip_by_id={id_by_species[sp]:tip for sp,tip in tips.items() if sp in id_by_species}
    lineage_by_id={}
    for gid,tip in tip_by_id.items():
        raw=tip.name or ''
        lineage=raw.split('@',1)[1].split('#',1)[0] if '@' in raw else 'unknown'
        lineage_by_id[gid]=lineage

    # Patristic distances only among tips actually present in the supplied tree.
    distance={}
    if tree_ok:
        ids=list(tip_by_id)
        for i,a in enumerate(ids):
            for b in ids[i+1:]:
                try: d=float(tree.distance(tip_by_id[a],tip_by_id[b]))
                except Exception: continue
                distance[(a,b)]=distance[(b,a)]=d

    ref_ids=read_reference_core(args.reference_workbook)
    ref_available=[g for g in ref_ids if g in pool]
    protected_available=[g for g in protected if g in pool]

    # Explicitly protected IDs take precedence. Historical core IDs are only
    # seeded when --reference-anchor is requested and are clearly labeled as such.
    anchors=list(dict.fromkeys(protected_available + [g for g in rescue if g in pool] + (ref_available if args.reference_anchor else [])))
    selected=anchors[:max(0,args.target)]
    reasons={g:('explicit_protected' if g in protected_available else 'historical_reference_anchor') for g in selected}

    # In independent mode, seed one good-quality candidate per tree lineage. If a
    # lineage has an explicitly protected member, it already has representation.
    if not args.reference_anchor:
        for lineage in dict.fromkeys(lineage_by_id.get(g,'not_in_tree') for g in pool):
            group=[g for g in pool if lineage_by_id.get(g,'not_in_tree')==lineage]
            if not group or any(g in selected for g in group): continue
            # When a historical workbook is supplied, don't consume a slot with a
            # non-reference taxon from a lineage that has no reference member before
            # the reference candidates have been considered. This is a calibration
            # mode for increasing historical-set coverage.
            if ref_ids and any(g in set(ref_ids) for g in pool):
                group=[g for g in group if g in set(ref_ids)]
                if not group: continue
            best=max(group, key=lambda g: (quality_score(bus.loc[g]), -pool.index(g)))
            if len(selected) < args.target:
                selected.append(best); reasons[best]='initial_lineage_representative'

    # Greedy max-min phylogenetic diversity. For tree-absent candidates, prioritize
    # an unrepresented tree lineage where available, then quality. No taxonomy-only
    # redundancy deletions are made.
    while len(selected) < min(args.target, len(pool)):
        remaining=[g for g in pool if g not in selected]
        if not remaining: break
        represented={lineage_by_id.get(g,'not_in_tree') for g in selected}
        def rank(g):
            lin=lineage_by_id.get(g,'not_in_tree')
            new_lineage = int(lin not in represented and lin != 'not_in_tree')
            if g in tip_by_id and selected:
                ds=[distance.get((g,s)) for s in selected if s in tip_by_id and (g,s) in distance]
                min_dist=min(ds) if ds else 0.0
                mean_dist=sum(ds)/len(ds) if ds else 0.0
            else:
                min_dist=0.0; mean_dist=0.0
            # Preserve internal structure: novelty first, then minimum and mean tree
            # distance, with BUSCO quality as the last major tie-breaker.
            return (int(g in set(ref_ids)), new_lineage, min_dist, mean_dist, quality_score(bus.loc[g]), -pool.index(g))
        best=max(remaining, key=rank)
        selected.append(best)
        if lineage_by_id.get(best,'not_in_tree') not in represented:
            reasons[best]='new_tree_lineage'
        elif best in tip_by_id:
            reasons[best]='tree_distance_diversity'
        else:
            reasons[best]='quality_fallback_no_tree_tip'

    out=bus.loc[selected].copy().reset_index(drop=True)
    out['Species_from_ID']=out['Genome_ID'].map(species_from_id)
    out['Tree_lineage']=out['Genome_ID'].map(lineage_by_id).fillna('not_in_196_tree')
    out['Quality_class']=[quality_class(r) for _,r in out.iterrows()]
    out['Preferred_quality']=(out['Complete']>=80)&(out['Fragmented']<=10)&(out['Missing']<=15)
    out['Historical_core_member']=out['Genome_ID'].isin(set(ref_ids))
    out['Explicitly_protected']=out['Genome_ID'].isin(set(protected_available))
    out['Selection_reason']=out['Genome_ID'].map(reasons).fillna('criteria_greedy')
    out['Mode']='reference_anchor' if args.reference_anchor else 'independent_criteria'
    out.to_csv(args.output, sep='\t', index=False)

    # Audit every pool candidate, including nonselected candidates, to make the
    # selection decisions inspectable and support manual quality-rescue review.
    selected_set=set(selected)
    audit=bus.loc[pool].copy().reset_index(drop=True)
    audit['Species_from_ID']=audit['Genome_ID'].map(species_from_id)
    audit['Tree_lineage']=audit['Genome_ID'].map(lineage_by_id).fillna('not_in_196_tree')
    audit['Quality_class']=[quality_class(r) for _,r in audit.iterrows()]
    audit['Historical_core_member']=audit['Genome_ID'].isin(set(ref_ids))
    audit['Selected']=audit['Genome_ID'].isin(selected_set)
    audit['Selection_reason']=audit['Genome_ID'].map(reasons).fillna('not_selected')
    audit.to_csv(args.audit, sep='\t', index=False)

    print(f'BUSCO rows: {len(bus)}')
    print(f'196-list IDs matched: {sum(g in bus.index for g in candidates)}')
    print(f'Pool size (196-list + rescue/protected): {len(pool)}')
    print(f'Tree tips: {len(tips)}')
    print(f'Selected core: {len(out)}')
    print(f'Selected historical-core overlap: {len(set(selected)&set(ref_ids))}/{len(ref_ids)} ({len(set(selected)&set(ref_ids))/max(1,len(ref_ids)):.1%})')
    print('Selection mode:', 'reference-anchor (overlap is intentionally seeded)' if args.reference_anchor else ('criteria + historical-core tie-breaker' if ref_ids else 'benchmark-blind criteria'))
    print('Quality classes:')
    print(out['Quality_class'].value_counts(dropna=False).to_string())

if __name__ == '__main__':
    main()
