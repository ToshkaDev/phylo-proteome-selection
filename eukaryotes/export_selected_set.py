#!/usr/bin/env python3
"""Export selected core IDs plus the venom/toxin-only supplement.

Works with both workbook variants:
  * role-aware workbooks containing a ``Set role`` column; and
  * role-free workbooks containing the 146-row final set (137 core + 9 supplement-only).

When ``Set role`` is absent, the nine supplement-only entries are identified by
stable Genome_ID values / organism names from the curated supplement. Metadata and
column order are taken from the supplied workbook. Selected core IDs come from the
selector TSV; supplement-only entries are then appended without duplicating IDs.
"""
import argparse
import sys
import pandas as pd

# Stable identities of the nine curated supplement-only genomes. These are used
# only when the supplied workbook has no Set role column.
SUPPLEMENT_GENOME_IDS = {
    'NCBI_GCF_003254395.2_Apis_mellifera',
    'NCBI_GCF_000671375.1_Centruroides_sculpturatus',
    'NCBI_GCF_039654945.1_Dendrobates_tinctorius',
    'NCBI_GCF_016920785.2_Ixodes_scapularis',
    'NCBI_GCA_009733165.1_Naja_naja',
    'NCBI_GCF_001194135.2_Octopus_bimaculoides',
    'NCBI_GCF_004115215.2_Ornithorhynchus_anatinus',
    'EP00794_Sanderia_malayensis',
    'EP00083_Stegodyphus_mimosarum',
}
SUPPLEMENT_ORGANISMS = {
    'Apis mellifera', 'Centruroides sculpturatus', 'Dendrobates tinctorius',
    'Ixodes scapularis', 'Naja naja', 'Octopus bimaculoides',
    'Ornithorhynchus anatinus', 'Sanderia malayensis', 'Stegodyphus mimosarum',
}

def ensure_columns(frame, columns):
    """Return frame with all requested columns, in the requested order."""
    frame = frame.copy()
    for col in columns:
        if col not in frame.columns:
            frame[col] = pd.NA
    return frame[columns]

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--core-tsv', required=True, help='Selector output TSV containing Genome_ID')
    ap.add_argument('--reference-workbook', required=True, help='Curated workbook used as metadata source')
    ap.add_argument('--output-xlsx', default='selected_eukaryotic_set.xlsx')
    ap.add_argument('--output-tsv', default='selected_eukaryotic_set.tsv')
    ap.add_argument('--supplement-ids', help='Optional text file with one supplement Genome_ID per line; overrides built-in IDs when Set role is absent')
    args = ap.parse_args()

    core = pd.read_csv(args.core_tsv, sep='\t', dtype={'Genome_ID': str})
    if 'Genome_ID' not in core.columns:
        ap.error(f"Core TSV must contain a 'Genome_ID' column; found: {list(core.columns)}")
    selected = set(core['Genome_ID'].dropna().astype(str))

    ref = pd.read_excel(args.reference_workbook)
    required = {'Genome_ID', 'Organism'}
    missing_required = required - set(ref.columns)
    if missing_required:
        ap.error(f"Reference workbook is missing required columns: {sorted(missing_required)}")

    ref['Genome_ID'] = ref['Genome_ID'].astype(str)
    if 'Set role' in ref.columns:
        roles = ref['Set role'].fillna('').astype(str).str.strip()
        core_mask = roles.eq('') | roles.str.contains('Existing core', case=False, na=False)
        supplement_mask = roles.eq('Venom/toxin supplement')
    else:
        if args.supplement_ids:
            with open(args.supplement_ids, encoding='utf-8') as handle:
                supplement_ids = {line.strip() for line in handle if line.strip() and not line.lstrip().startswith('#')}
        else:
            supplement_ids = set(SUPPLEMENT_GENOME_IDS)
        supplement_mask = ref['Genome_ID'].isin(supplement_ids) | ref['Organism'].astype(str).isin(SUPPLEMENT_ORGANISMS)
        core_mask = ~supplement_mask
        print("Note: 'Set role' is absent; supplement rows identified using the curated nine-genome ID/name list.", file=sys.stderr)

    # Select core entries from workbook metadata where available.
    core_rows = ref.loc[core_mask & ref['Genome_ID'].isin(selected)].copy()
    found = set(core_rows['Genome_ID'])
    missing_ids = selected - found
    if missing_ids:
        # Use selector output to fill missing selected genomes, preserving workbook schema.
        extra = core.loc[core['Genome_ID'].astype(str).isin(missing_ids)].copy()
        extra = ensure_columns(extra, list(ref.columns))
        core_rows = pd.concat([core_rows, extra], ignore_index=True)

    # Supplement entries are always retained, except when already selected as core.
    supplement_rows = ref.loc[supplement_mask & ~ref['Genome_ID'].isin(selected)].copy()

    # Stable deduplication and readable sorting; avoid requiring optional columns.
    core_rows = core_rows.drop_duplicates('Genome_ID', keep='first')
    supplement_rows = supplement_rows.drop_duplicates('Genome_ID', keep='first')
    final = pd.concat([core_rows, supplement_rows], ignore_index=True)
    final = final.drop_duplicates('Genome_ID', keep='first')
    sort_cols = [c for c in ['Set role', 'Organism'] if c in final.columns]
    if sort_cols:
        final = final.sort_values(sort_cols, na_position='first')
    final = final.reset_index(drop=True)

    final.to_excel(args.output_xlsx, index=False)
    final.to_csv(args.output_tsv, sep='\t', index=False)

    print(f'Selected core IDs in TSV: {len(selected)}')
    print(f'Core rows exported: {len(core_rows)}')
    print(f'Venom/toxin-only rows added: {len(supplement_rows)}')
    print(f'Final rows: {len(final)}')
    print(f'Final unique genomes: {final["Genome_ID"].nunique()}')
    print(f'Workbook columns preserved: {len(ref.columns)}')
    print(f'Wrote {args.output_xlsx} and {args.output_tsv}')

if __name__ == '__main__':
    main()
