#!/usr/bin/env python3
"""Select exactly 3,424 GTDB r232 genomes using reconstructed phylogeny.

The existing 3,424-genome selection is used ONLY to set per-order quotas and to
preserve five documented archaeal lineage-recovery exception genomes. Within each
order, families are selected for phylogenetic coverage from reconstructed family
pair distances: take one medoid per region when quota permits, then greedily add
families maximizing distance to the nearest already-selected family. If quota is
smaller than region count, choose region medoids by farthest-first.

Baseline quality filter is completeness >=80 and contamination
<=5, except for the five explicitly preserved legacy archaeal exception genomes.
"""
from pathlib import Path
import zipfile, pandas as pd, numpy as np

BASE=Path(__file__).resolve().parents[1]
WORK=BASE/'outputs'/'baseline'
META=BASE/'inputs'/'metadata_r232.tsv'
OLD=BASE/'inputs'/'legacy_selection_3424.tsv'
PHY=BASE/'inputs'/'phylogeny'
PAIRS=PHY/'family_pairwise_min_distances_reconstructed.tsv'
METRICS=PHY/'order_phylogeny_metrics_reconstructed.tsv'
REGIONS=PHY/'family_region_assignments_reconstructed.tsv'
OUT=WORK
OUT.mkdir(exist_ok=True)
FALSE={'','nan','none','na','missing','n/a'}

def truth(x): return str(x).strip().lower() in {'t','true','1','yes'}

def ref_score(r):
    cat=str(r.get('ncbi_refseq_category','')).lower()
    if 'reference genome' in cat or cat=='reference': return 5
    if 'representative genome' in cat or cat=='representative': return 4
    typ=str(r.get('ncbi_type_material_designation','')).strip()
    if typ.lower() not in FALSE: return 3
    if str(r.get('ncbi_isolate','')).strip().lower() not in FALSE: return 2
    if truth(r.get('gtdb_representative','')) or truth(r.get('gtdb_genome_representative','')): return 1
    return 0

def choose_genome(g):
    g=g.copy(); g['checkm2_completeness']=pd.to_numeric(g.checkm2_completeness,errors='coerce'); g['checkm2_contamination']=pd.to_numeric(g.checkm2_contamination,errors='coerce')
    g=g[(g.checkm2_completeness>=80)&(g.checkm2_contamination<=5)].copy()
    if g.empty: raise RuntimeError('Selected family has no eligible genome')
    g['_ref_score']=g.apply(ref_score,axis=1)
    mc=g.checkm2_completeness.max(); mn=g.checkm2_contamination.min()
    near=g[(mc-g.checkm2_completeness<=0.5)&(g.checkm2_contamination-mn<=0.5)].copy()
    if near.empty: near=g.copy()
    near=near.sort_values(['checkm2_completeness','checkm2_contamination','_ref_score','accession'],ascending=[False,True,False,True],kind='mergesort')
    return near.iloc[0]

print('Loading metadata and target quotas...',flush=True)

md=pd.read_csv(META,sep='\t',dtype=str,low_memory=False).fillna('')
md['checkm2_completeness']=pd.to_numeric(md.checkm2_completeness,errors='coerce'); md['checkm2_contamination']=pd.to_numeric(md.checkm2_contamination,errors='coerce')
tax=md.gtdb_taxonomy.str.split(';')
for i,c in [(0,'domain'),(3,'order'),(4,'family')]: 
    md[c]=tax.map(lambda x:x[i].split('__',1)[-1] if len(x)>i and '__' in x[i] else '')
md=md[(md.checkm2_completeness>=80)&(md.checkm2_contamination<=5)&md.domain.isin(['Bacteria','Archaea'])&md.order.ne('')&md.family.ne('')].copy()
old=pd.read_csv(OLD,sep='\t',dtype=str).fillna('')
quota=old.groupby(['domain','order']).size().to_dict()

# These five manually retained archaeal lineage-recovery exceptions fall outside the
# standard 80/5 filtered universe; retain them exactly to preserve the agreed baseline.
exc=old[old.selection_method.eq('archaeal lineage recovery exception')].copy()
exc_orders=set(zip(exc.domain,exc.order))
print(f'Eligible genomes: {len(md):,}; orders: {md.groupby(["domain","order"]).ngroups:,}; legacy exception rows: {len(exc)}',flush=True)
metrics=pd.read_csv(METRICS,sep='\t').fillna(''); regions=pd.read_csv(REGIONS,sep='\t').fillna(''); pairs=pd.read_csv(PAIRS,sep='\t')
# matrices and deterministic indexes
D={}

for (d,o),p in pairs.groupby(['domain','order'],sort=False):
    fs=sorted(set(p.family1.astype(str))|set(p.family2.astype(str))); ix={f:i for i,f in enumerate(fs)}; mat=np.zeros((len(fs),len(fs)),float)
    for r in p.itertuples(index=False):
        i,j=ix[str(r.family1)],ix[str(r.family2)]; mat[i,j]=mat[j,i]=float(r.patristic_distance)
    D[(d,o)]=(fs,mat,ix)
sel=[]; audits=[]

for (d,o),block in md.groupby(['domain','order'],sort=True):
    if (d,o) in exc_orders: continue
    fams=sorted(block.family.unique().tolist()); n=len(fams); target=int(quota.get((d,o),1)); target=min(target,n)
    if n==1:
        chosen=fams[:]; nreg=1; scale=np.nan; threshold=np.nan; regmap={f:1 for f in fams}; method='single eligible family'
    else:
        if (d,o) not in D: raise RuntimeError(f'missing pair distances {d}/{o}')
        fs,mat,ix=D[(d,o)]
        if set(fs)!=set(fams): raise RuntimeError(f'family-distance coverage mismatch {d}/{o}')
        mr=metrics[(metrics.domain==d)&(metrics.order==o)]
        if mr.empty: raise RuntimeError(f'missing metrics {d}/{o}')
        mr=mr.iloc[0]; scale=float(mr.median_nn_reconstructed); threshold=float(mr.threshold_reconstructed); nreg=int(mr.effective_families_reconstructed)
        rg=regions[(regions.domain==d)&(regions.order==o)]
        regmap=dict(zip(rg.family.astype(str),rg.effective_region_reconstructed.astype(int)))
        reg_members={r:sorted(f for f in fams if regmap.get(f)==r) for r in sorted(set(regmap.values()))}
        medoids=[]
        for r,members in reg_members.items():
            ids=[ix[f] for f in members]; sub=mat[np.ix_(ids,ids)]
            medoids.append((r,members[int(np.argmin(sub.mean(axis=1))) ]))
        reps=[f for _,f in medoids]
        def farthest_first(candidates, count, seeds=None):
            chosen=list(seeds or [])
            while len(chosen)<count:
                pool=[f for f in candidates if f not in chosen]
                if not pool: break
                if not chosen:
                    # first seed maximizes distance to its nearest other family
                    nxt=max(pool,key=lambda f:(min(mat[ix[f],ix[g]] for g in pool if g!=f) if len(pool)>1 else 0.0,f))
                else:
                    nxt=max(pool,key=lambda f:(min(mat[ix[f],ix[g]] for g in chosen),f))
                chosen.append(nxt)
            return chosen
        if target<=len(reps):
            chosen=farthest_first(reps,target)
            method='region medoids, farthest-first to per-order target'
        else:
            chosen=list(reps)
            chosen=farthest_first(fams,target,seeds=chosen)
            method='all region medoids + farthest-first family additions to per-order target'
    for rank,f in enumerate(chosen,1):
        sel.append({'domain':d,'order':o,'selected_family':f,'selection_rank':rank,'effective_region_reconstructed':regmap.get(f,1),'raw_families':n,'median_nn_reconstructed':scale,'threshold_reconstructed':threshold,'effective_families_reconstructed':nreg,'target_genomes_from_legacy_order_quota':int(quota.get((d,o),1)),'selection_method_reconstructed':method})
    audits.append({'domain':d,'order':o,'eligible_families':n,'legacy_order_quota':int(quota.get((d,o),1)),'selected_families':len(chosen),'regions':nreg,'median_nn_reconstructed':scale,'threshold_reconstructed':threshold,'selection_method':method,'families':';'.join(chosen)})

sf=pd.DataFrame(sel)

# Add the four out-of-filter lineage exceptions to genome rows, not to eligible-family optimization.
md_groups={(d,o,f):g.copy() for (d,o,f),g in md.groupby(['domain','order','family'],sort=False)}
rows=[]
for r in sf.itertuples(index=False):
    g=md_groups[(r.domain,r.order,r.selected_family)]
    best=choose_genome(g)
    rec=best.drop(labels=['_ref_score'],errors='ignore').to_dict()
    rec.update({k:getattr(r,k) for k in sf.columns if k not in {'domain','order','selected_family'}})
    rec['selected_family']=r.selected_family; rec['genome_choice_method']='near-optimal quality window + reference metadata tie-break'
    rows.append(rec)

# Keep exception records from the old table as-is; mark them explicitly.
for _,r in exc.iterrows():
    rec=r.to_dict(); rec['selection_method_reconstructed']='preserved legacy archaeal lineage-recovery exception'; rec['effective_region_reconstructed']='exception'; rec['genome_choice_method']='preserved legacy exception (outside 80/5 filter)'; rows.append(rec)
out=pd.DataFrame(rows).fillna('')
out.to_csv(OUT/'selected_genomes_3424_phylo_coverage.tsv',sep='\t',index=False)
sf.to_csv(OUT/'selected_families_3424_phylo_coverage.tsv',sep='\t',index=False)
pd.DataFrame(audits).to_csv(OUT/'order_selection_audit.tsv',sep='\t',index=False)

# Comparison summaries
oldfam=set(zip(old.domain,old.order,old.selected_family)); newfam=set(zip(out.domain.astype(str),out.order.astype(str),out.selected_family.astype(str)))
newacc=set(out.accession.astype(str)); oldacc=set(old.accession.astype(str))

# per-order phylogenetic coverage: nearest selected-family distance for every eligible family
coverage=[]
selmap={(d,o):set(g.selected_family.astype(str)) for (d,o),g in sf.groupby(['domain','order'])}
for (d,o),block in md.groupby(['domain','order'],sort=True):
    fs=sorted(block.family.unique().tolist()); chosen=selmap.get((d,o),set())
    if not chosen: continue
    if len(fs)==1:
        nearest=np.array([0.0])
        within=[]
    else:
        fs0,mat,ix=D[(d,o)]; chosen_ids=[ix[f] for f in sorted(chosen) if f in ix]
        allids=[ix[f] for f in fs]
        nearest=mat[np.ix_(allids,chosen_ids)].min(axis=1)
        ids=[ix[f] for f in sorted(chosen) if f in ix]
        within=mat[np.ix_(ids,ids)][np.triu_indices(len(ids),1)] if len(ids)>1 else []
    coverage.append({'domain':d,'order':o,'eligible_families':len(fs),'selected_families':len(chosen),'mean_distance_family_to_nearest_selected':float(np.mean(nearest)),'max_distance_family_to_nearest_selected':float(np.max(nearest)),'mean_pairwise_distance_selected_families':float(np.mean(within)) if len(within) else 0.0,'max_pairwise_distance_selected_families':float(np.max(within)) if len(within) else 0.0})
cov=pd.DataFrame(coverage); cov.to_csv(OUT/'per_order_phylogenetic_coverage.tsv',sep='\t',index=False)

# Compare coverage metrics on the exact same eligible family universe, using old selected family sets and reconstructed new set.
def coverage_for_family_set(fset):
    vals=[]
    for (d,o),block in md.groupby(['domain','order'],sort=True):
        fs=sorted(block.family.unique().tolist()); ch=sorted(f for f in fset if f[0]==d and f[1]==o and f[2] in set(fs))
        if not ch: continue
        if len(fs)==1: near=np.array([0.0]); within=[]
        else:
            fs0,mat,ix=D[(d,o)]; cids=[ix[x[2]] for x in ch]; aids=[ix[f] for f in fs]
            near=mat[np.ix_(aids,cids)].min(axis=1); sub=mat[np.ix_(cids,cids)]; within=sub[np.triu_indices(len(cids),1)] if len(cids)>1 else []
        vals.append((d,o,float(near.mean()),float(near.max()),float(np.mean(within)) if len(within) else 0.0))
    z=pd.DataFrame(vals,columns=['domain','order','mean_nearest_cover','max_nearest_cover','mean_pairwise_selected'])
    return {'orders_with_selected_families':len(z),'mean_order_mean_nearest_cover':z.mean_nearest_cover.mean(),'median_order_mean_nearest_cover':z.mean_nearest_cover.median(),'mean_order_max_nearest_cover':z.max_nearest_cover.mean(),'mean_order_mean_pairwise_selected':z.mean_pairwise_selected.mean()}

# Use the 2,579 reconstruction output from package.
base_new=pd.read_csv(BASE/'inputs'/'reference_selection_2579_families.tsv',sep='\t',dtype=str)
base_set=set(zip(base_new.domain,base_new.order,base_new.selected_family))
old_set={x for x in oldfam if x not in {(r.domain,r.order,r.selected_family) for r in exc.itertuples()}}
comparison=[]

for name,fsset in [('legacy_3424_family_set',old_set),('reconstructed_2579_family_set',base_set),('new_3419_eligible_family_set',set(zip(sf.domain,sf.order,sf.selected_family)))]:
    comparison.append({'selection':name,**coverage_for_family_set(fsset)})

pd.DataFrame(comparison).to_csv(OUT/'selection_phylogenetic_coverage_comparison.tsv',sep='\t',index=False)
summary=['GTDB r232 PHYLOGENETIC COVERAGE SELECTION TARGETING 3,424 GENOMES','',f'Eligible genomes under 80/5 filter: {len(md):,}',f'Eligible orders: {md.groupby(["domain","order"]).ngroups:,}',f'Genome rows output: {len(out):,}',f'Unique accessions: {out.accession.nunique():,}',f'Bacteria: {(out.domain=="Bacteria").sum():,}',f'Archaea: {(out.domain=="Archaea").sum():,}',f'Preserved out-of-filter archaeal exception rows: {len(exc)}',f'Families selected from eligible filtered universe: {len(sf):,}',f'Unique (domain,order,family) rows in output: {out.groupby(["domain","order","selected_family"]).ngroups:,}',f'Legacy family-set overlap (unique order/family): {len(oldfam & newfam):,} / {len(oldfam):,}',f'Legacy accession overlap: {len(oldacc & newacc):,} / {len(oldacc):,}', '', 'Selection logic:', '1. Use legacy per-order row counts as target quotas, not old family identities.', '2. Within each order, select reconstructed region medoids first when quota allows.', '3. Fill remaining quota by farthest-first maximin distance over all eligible families in that order.', '4. Select one genome per chosen eligible family using the quality/reference tie-break procedure.', '5. Preserve five manually added archaeal lineage-recovery exception genomes from the legacy set, which are outside the 80/5 filtered universe.', '', 'This is a new coverage-oriented reconstruction. It targets the agreed 3,424 genome count and 3,182/242 domain split.']
(OUT/'SUMMARY.txt').write_text('\n'.join(summary)+'\n')
print('\n'.join(summary))
print('\nCoverage comparison:')
print(pd.DataFrame(comparison).to_string(index=False))
