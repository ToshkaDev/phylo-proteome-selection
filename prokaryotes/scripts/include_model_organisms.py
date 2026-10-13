from pathlib import Path
import pandas as pd, numpy as np, re, zipfile, os
BASE=Path(__file__).resolve().parents[1]
META=BASE/'inputs'/'metadata_r232.tsv'
SEL=BASE/'outputs'/'baseline'/'selected_genomes_3424_phylo_coverage.tsv'
PAIRS=BASE/'inputs'/'phylogeny'/'family_pairwise_min_distances_reconstructed.tsv'
REGIONS=BASE/'inputs'/'phylogeny'/'family_region_assignments_reconstructed.tsv'
OUT=BASE/'outputs'/'model_enriched'; OUT.mkdir(parents=True, exist_ok=True)

# Names and strain preferences transcribed from the user-provided model-organism list.
# Aliases help account for genus transfers and older nomenclature; exact strain is preferred when available.
TARGETS = [
('Escherichia coli',['K-12','MG1655']),('Escherichia coli',['BL21']),('Escherichia coli',['O157:H7','EDL933']),
('Pseudomonas aeruginosa',['PAO1']),('Pseudomonas putida',['KT2440']),('Agrobacterium tumefaciens',['C58']),('Rhizobium leguminosarum',[]),('Sinorhizobium meliloti',['1021']),('Vibrio cholerae',['N16961']),('Myxococcus xanthus',['DK1622']),('Helicobacter pylori',['26695']),
('Bacillus subtilis',['168']),('Bacillus anthracis',['Ames']),('Staphylococcus aureus',['NCTC 8325','USA300']),('Streptococcus pneumoniae',['TIGR4']),('Streptococcus pyogenes',[]),('Lactococcus lactis',[]),('Listeria monocytogenes',['EGD-e']),('Clostridium acetobutylicum',[]),('Clostridioides difficile',[]),
('Mycobacterium tuberculosis',['H37Rv']),('Mycobacterium smegmatis',[]),('Corynebacterium glutamicum',[]),('Bifidobacterium longum',[]),('Synechocystis sp. PCC 6803',['PCC 6803']),('Synechococcus elongatus',['PCC 7942']),('Anabaena sp. PCC 7120',['PCC 7120']),('Bacteroides thetaiotaomicron',[]),('Chlamydia trachomatis',[]),('Deinococcus radiodurans',['R1']),('Thermus thermophilus',[]),
('Streptococcus agalactiae',[]),('Streptococcus mutans',[]),('Neisseria meningitidis',[]),('Neisseria gonorrhoeae',[]),('Legionella pneumophila',['Philadelphia 1']),('Salmonella enterica',['LT2']),('Yersinia pestis',['CO92']),('Bordetella pertussis',[]),('Borrelia burgdorferi',['B31']),('Treponema pallidum',[]),('Campylobacter jejuni',['NCTC 11168']),
('Candidatus Pelagibacter ubique',['HTCC1062']),('Prochlorococcus marinus',['MED4','MIT9313']),('Geobacter sulfurreducens',[]),('Geobacter metallireducens',[]),('Shewanella oneidensis',['MR-1']),('Nitrosomonas europaea',[]),('Nitrobacter winogradskyi',[]),('Desulfovibrio vulgaris',[]),('Azotobacter vinelandii',[]),('Aliivibrio fischeri',['ES114']),('Thermotoga maritima',[]),('Aquifex aeolicus',[]),('Buchnera aphidicola',[]),('Rickettsia prowazekii',[]),
('Methanococcus jannaschii',[]),('Halobacterium salinarum',['NRC-1']),('Haloferax volcanii',[]),('Methanosarcina acetivorans',[]),('Thermococcus kodakarensis',[]),('Sulfolobus solfataricus',[]),('Sulfolobus acidocaldarius',[]),('Pyrobaculum aerophilum',[]),('Candidatus Prometheoarchaeum syntrophicum',['MK-D1']),('Lokiarchaeum',[]),('Odinarchaeum',[]),('Methanothermobacter thermautotrophicus',['Delta H']),('Methanobrevibacter smithii',[]),('Methanococcoides burtonii',[]),('Nitrosopumilus maritimus',['SCM1']),('Pyrococcus furiosus',[]),('Metallosphaera sedula',[]),('Ignicoccus hospitalis',[]),('Nanoarchaeum equitans',['Kin4-M']),('Candidatus Micrarchaeum acidiphilum',[])
]

def norm(s): return re.sub(r'[^a-z0-9]+',' ',str(s).lower()).strip()
def truth(s): return str(s).lower() in {'t','true','1','yes'}
def refscore(r):
    c=str(r.get('ncbi_refseq_category','')).lower()
    if 'reference genome' in c or c=='reference': return 5
    if 'representative genome' in c or c=='representative': return 4
    if str(r.get('ncbi_type_material_designation','')).strip().lower() not in {'','na','nan','none'}: return 3
    if str(r.get('ncbi_isolate','')).strip().lower() not in {'','na','nan','none'}: return 2
    if truth(r.get('gtdb_representative','')): return 1
    return 0

def preferred_score(row, organism, strains):
    name=norm(row.ncbi_organism_name); taxonomy=norm(row.gtdb_taxonomy); ids=norm(row.ncbi_strain_identifiers+' '+row.ncbi_isolate)
    org=norm(organism)
    # candidate must have taxonomic or NCBI-name evidence for requested organism; avoid matches to wrong bins where name contradicts GTDB genus.
    base=org in name or org in taxonomy
    if not base:
        # For Candidatus/genus alias handling, allow distinctive organism token phrase.
        toks=org.replace('candidatus ','').split()
        base=bool(toks) and ' '.join(toks[:2]) in name
    if not base: 
        return None
    strain_score=0
    for i,s in enumerate(strains):
        if norm(s) and norm(s) in ids+' '+name:
            strain_score=max(strain_score,100-i)
    return (strain_score, refscore(row), float(row.checkm2_completeness), -float(row.checkm2_contamination), row.accession)

md=pd.read_csv(META,sep='\t',dtype=str).fillna('')
for c in ['checkm2_completeness','checkm2_contamination']:
    md[c]=pd.to_numeric(md[c],errors='coerce')
tax=md.gtdb_taxonomy.str.split(';')
for i,c in [(0,'domain'),(3,'order'),(4,'family')]:
    md[c]=tax.map(lambda x:x[i].split('__',1)[-1] if len(x)>i and '__' in x[i] else '')
md=md[(md.checkm2_completeness>=80)&(md.checkm2_contamination<=5)&md.domain.isin(['Bacteria','Archaea'])&md.order.ne('')&md.family.ne('')].copy()
md['_name_norm']=md.ncbi_organism_name.map(norm)
md['_tax_norm']=md.gtdb_taxonomy.map(norm)
md['_id_norm']=(md.ncbi_strain_identifiers+' '+md.ncbi_isolate).map(norm)
sel=pd.read_csv(SEL,sep='\t',dtype=str).fillna('')
regions=pd.read_csv(REGIONS,sep='\t',dtype=str).fillna('')
regions['effective_region_reconstructed']=pd.to_numeric(regions.effective_region_reconstructed,errors='coerce')
pairs=pd.read_csv(PAIRS,sep='\t')

# distance matrices for family swaps
D={}
for (d,o),p in pairs.groupby(['domain','order'],sort=False):
    fs=sorted(set(p.family1.astype(str))|set(p.family2.astype(str))); ix={f:i for i,f in enumerate(fs)}; mat=np.zeros((len(fs),len(fs)))
    for r in p.itertuples(index=False):
        i,j=ix[str(r.family1)],ix[str(r.family2)]; mat[i,j]=mat[j,i]=float(r.patristic_distance)
    D[(d,o)]=(ix,mat)

# normalize region lookup by domain/order/family
regmap={(r.domain,r.order,r.family):str(r.effective_region_reconstructed) for r in regions.itertuples(index=False)}

# baseline selected rows only; keep order quota and total count exactly constant
selected=sel.copy(); selected['_forced_model']='false'; selected['_model_label']=''; selected['_match_reason']=''
selected['checkm2_completeness']=pd.to_numeric(selected.checkm2_completeness,errors='coerce'); selected['checkm2_contamination']=pd.to_numeric(selected.checkm2_contamination,errors='coerce')

# Build matches; use GTDB genus/species taxonomy consistency to reduce false positive NCBI names.
matched=[]; unresolved=[]
for organism,strains in TARGETS:
    org_norm=norm(organism)
    mask=md['_name_norm'].str.contains(re.escape(org_norm),na=False) | md['_tax_norm'].str.contains(re.escape(org_norm),na=False)
    if not mask.any():
        toks=org_norm.replace('candidatus ','').split()
        if toks:
            mask=md['_name_norm'].str.contains(re.escape(' '.join(toks[:2])),na=False)
    cand=md[mask].copy()
    if cand.empty:
        unresolved.append({'model_organism':organism,'preferred_strains':';'.join(strains),'status':'not found in eligible 80/5 r232 metadata','candidate_count':0,'selected_accession':'','selected_family':'','action':'unresolved'})
        continue
    cand['_score']=cand.apply(lambda r: preferred_score(r,organism,strains),axis=1)
    # rank by strain exactness, then reference, completeness, contamination
    cand=cand.sort_values('_score',ascending=False,kind='mergesort')
    # Keep all candidates in audit but best one is the desired representative for this target.
    best=cand.iloc[0]
    matched.append({'model_organism':organism,'preferred_strains':';'.join(strains),'accession':best.accession,'ncbi_organism_name':best.ncbi_organism_name,'ncbi_strain_identifiers':best.ncbi_strain_identifiers,'domain':best.domain,'order':best.order,'family':best.family,'gtdb_taxonomy':best.gtdb_taxonomy,'completeness':best.checkm2_completeness,'contamination':best.checkm2_contamination,'candidate_count':len(cand),'strain_preference_score':best['_score'][0]})

# One representative per selected family; for multiple targets in same family, keep best strain-specific target if one exists.
matched_df=pd.DataFrame(matched)

# Resolve family targets in deterministic order: exact requested strain match first, then list order.
if not matched_df.empty:
    matched_df=matched_df.sort_values(['strain_preference_score','model_organism'],ascending=[False,True],kind='mergesort')
for t in matched_df.itertuples(index=False):
    idx=selected.index[(selected.domain==t.domain)&(selected.order==t.order)&(selected.selected_family==t.family)]
    # If family is already represented, replace only genome within that family.
    if len(idx):
        i=idx[0]; oldacc=selected.at[i,'accession']
        selected.at[i,'accession']=t.accession
        # copy core metadata from candidate
        rr=md[md.accession==t.accession].iloc[0]
        for c in md.columns:
            if c in selected.columns: selected.at[i,c]=rr[c]
        selected.at[i,'_forced_model']='true'; selected.at[i,'_model_label']=t.model_organism
        selected.at[i,'_match_reason']=f'model genome replaced genome within already-selected family; prior accession={oldacc}'
        continue
    # Family is absent. Can swap only if the target family is in an already represented region OR a selected region has redundant representation.
    key=(t.domain,t.order)
    if key not in D or (t.domain,t.order,t.family) not in regmap:
        unresolved.append({'model_organism':t.model_organism,'preferred_strains':t.preferred_strains,'status':'matched genome but no reconstructed region/distance data','candidate_count':t.candidate_count,'selected_accession':t.accession,'selected_family':t.family,'action':'unresolved'})
        continue
    ix,mat=D[key]; target_region=regmap[(t.domain,t.order,t.family)]
    order_idxs=selected.index[(selected.domain==t.domain)&(selected.order==t.order)]
    order_rows=selected.loc[order_idxs]
    # Count selected families per region; never remove last representative of any region.
    counts={}
    for _,r in order_rows.iterrows():
        rg=regmap.get((t.domain,t.order,r.selected_family),'NA'); counts[rg]=counts.get(rg,0)+1
    candidates=[]
    for i,r in order_rows.iterrows():
        oldfam=r.selected_family; oldreg=regmap.get((t.domain,t.order,oldfam),'NA')
        if counts.get(oldreg,0)<=1: continue
        if oldfam not in ix or t.family not in ix: continue
        # Replacement loss: increase in nearest-selected-family distance over all eligible families in this order.
        fs=sorted(set(ix)); current=[f for f in order_rows.selected_family.tolist() if f in ix]
        if oldfam not in current: continue
        proposed=[f for f in current if f!=oldfam]+[t.family]
        allids=[ix[f] for f in fs]; curids=[ix[f] for f in current]; newids=[ix[f] for f in proposed]
        before=mat[np.ix_(allids,curids)].min(axis=1).mean(); after=mat[np.ix_(allids,newids)].min(axis=1).mean()
        # Favor filling an unrepresented region, then minimal coverage loss; if same region, minimal loss.
        adds_region=target_region not in counts
        same_region=oldreg==target_region
        candidates.append((0 if adds_region else 1, after-before, 0 if same_region else 1, i, oldfam, oldreg, adds_region))
    if not candidates:
        # Preserve this model family by adding one explicitly documented exception genome,
        # because replacing any currently selected family would eliminate that region.
        rr=md[md.accession==t.accession].iloc[0].copy()
        rec={c: rr.get(c, '') for c in selected.columns if c in rr.index}
        rec.update({c: '' for c in selected.columns if c not in rec})
        rec['selected_family']=t.family
        rec['selection_method_reconstructed']='added model-organism exception to preserve phylogenetic region'
        rec['effective_region_reconstructed']=target_region
        rec['_forced_model']='true'; rec['_model_label']=t.model_organism
        rec['_match_reason']='added because no coverage-safe family swap was available; would otherwise remove an existing region'
        selected=pd.concat([selected,pd.DataFrame([rec])],ignore_index=True)
        unresolved.append({'model_organism':t.model_organism,'preferred_strains':t.preferred_strains,'status':'added as exception to preserve phylogenetic coverage','candidate_count':t.candidate_count,'selected_accession':t.accession,'selected_family':t.family,'action':'added; no safe swap available'})
        continue
    candidates.sort()
    _,delta,_,i,oldfam,oldreg,adds_region=candidates[0]
    # Replace genome row entirely with candidate metadata, while retaining selection/audit fields that are selection-level.
    rr=md[md.accession==t.accession].iloc[0]
    for c in md.columns:
        if c in selected.columns: selected.at[i,c]=rr[c]
    selected.at[i,'selected_family']=t.family
    selected.at[i,'_forced_model']='true'; selected.at[i,'_model_label']=t.model_organism
    selected.at[i,'_match_reason']=f'family swap {oldfam}->{t.family}; old region={oldreg}; new region={target_region}; mean-nearest-coverage delta={delta:.6g}; new_region={adds_region}'

# Ensure unique accessions. Family swaps preserve the baseline count; coverage-protecting additions intentionally increase it. Multiple model taxa in one family may collide.
selected['selected_family']=selected.apply(lambda r: r.get('selected_family') or '',axis=1)
if selected.accession.nunique()!=len(selected): 
    raise RuntimeError(f'accession uniqueness failed: {selected.accession.nunique()} / {len(selected)}')
if len(selected)<3424:
    raise RuntimeError(f'Genome count unexpectedly below baseline: {len(selected)}')

# The baseline contains a small number of legacy exception rows with repeated family labels;
# retain these rather than silently dropping legacy coverage.
# Audit selected target status, all requested organisms including not found / blocked.
sel_targets=[]
for _,t in matched_df.iterrows():
    found=selected[(selected.accession==t.accession)&(selected['_model_label']==t.model_organism)]
    if len(found): status='included'; action=found.iloc[0]['_match_reason']
    else:
        prior=selected[(selected.domain==t.domain)&(selected.order==t.order)&(selected.selected_family==t.family)]
        if len(prior): status='family represented, but requested genome not selected due to family-level one-genome rule/collision'; action=';'.join(prior['_model_label'].astype(str).tolist())
        else: status='matched candidate not retained'; action='candidate genome is absent from the final model-aware selection and its family is not represented there; review the order/family coverage constraints and target-collision decisions'
    sel_targets.append({'model_organism':t.model_organism,'preferred_strains':t.preferred_strains,'accession':t.accession,'ncbi_organism_name':t.ncbi_organism_name,'domain':t.domain,'order':t.order,'family':t.family,'candidate_count':t.candidate_count,'status':status,'action':action})
all_audit=pd.DataFrame(sel_targets)
unres=pd.DataFrame(unresolved)

if not unres.empty:
    # Enrich matching target rows with genuine resolution failures (e.g. missing distance/region data).
    # Do not duplicate matched targets that were already included as coverage exceptions.
    matched_keys={(str(r.model_organism), str(r.preferred_strains)) for r in matched_df.itertuples(index=False)}
    for u in unres.to_dict('records'):
        key=(str(u.get('model_organism','')), str(u.get('preferred_strains','')))
        if key in matched_keys and str(u.get('status','')).startswith('matched genome but'):
            mask=(all_audit.model_organism.astype(str)==key[0]) & (all_audit.preferred_strains.astype(str)==key[1])
            all_audit.loc[mask,'status']=u.get('status','unresolved')
            all_audit.loc[mask,'action']=u.get('action','unresolved')
    unres_only=unres[~unres.apply(lambda r: (str(r.get('model_organism','')), str(r.get('preferred_strains',''))) in matched_keys, axis=1)].copy()
    if not unres_only.empty: all_audit=pd.concat([all_audit,unres_only],ignore_index=True,sort=False)
    
selected.drop(columns=['_score'],errors='ignore').to_csv(OUT/'selected_genomes_model_organisms.tsv',sep='\t',index=False)
all_audit.to_csv(OUT/'model_organism_inclusion_audit.tsv',sep='\t',index=False)
matched_df.to_csv(OUT/'model_organism_candidates.tsv',sep='\t',index=False)
# family selection and order coverage audit
selected[['domain','order','selected_family','accession','ncbi_organism_name','_forced_model','_model_label','_match_reason']].to_csv(OUT/'selected_families_and_models.tsv',sep='\t',index=False)
summary=[f'Rows: {len(selected)}',f'Unique accessions: {selected.accession.nunique()}',f'Bacteria: {(selected.domain=="Bacteria").sum()}',f'Archaea: {(selected.domain=="Archaea").sum()}',f'Model target names in request: {len(TARGETS)}',f'Model taxa matched to at least one eligible metadata genome: {len(matched_df)}',f'Model genomes included as selected accession: {(all_audit.status=="included").sum() if not all_audit.empty else 0}',f'Target audit rows: {len(all_audit)}',f'Targets with no eligible match: {(all_audit.status=="not found in eligible 80/5 r232 metadata").sum()}',f'Targets blocked by family-level collision: {(all_audit.status=="family represented, but requested genome not selected due to family-level one-genome rule/collision").sum()}',f'Matched candidates not retained after coverage/collision decisions: {(all_audit.status=="matched candidate not retained").sum()}',f'Orders retained: {selected.groupby(["domain","order"]).ngroups}',f'Families represented: {selected.groupby(["domain","order","selected_family"]).ngroups}', '', 'Logic: start from the 3,424-genome baseline and preserve per-order quotas for safe swaps; prefer a model genome within an already-selected family; when the model family is not selected, replace a selected family only if another selected family still represents the outgoing family\'s region; among safe swaps, prioritize adding an unrepresented region and minimize change in mean nearest-selected-family distance.', 'Limit: one genome per selected family means different model strains/species in the same family can collide; unresolved targets are explicitly audited. Added model-organism exception rows only when a coverage-safe swap was unavailable.']
(OUT/'SUMMARY.txt').write_text('\n'.join(summary)+'\n')
print('\n'.join(summary)); print('\n',all_audit.status.value_counts(dropna=False).to_string())
