# Eukaryotic representative-set selection

This package implements the written criteria in `Eukaryotic_genome_selection_criteria.docx`
and the extended `Process.docx` procedure.

## Scientific principles encoded

1. Maximize phylogenetic/taxonomic coverage while maintaining useful proteome quality.
2. BUSCO is a quality-control criterion, not a hard exclusion filter:
   - Complete >=80%: keep; preferably Fragmented <=10% and Missing <=15%.
   - Complete 60-80%: valuable, especially when phylogenetically informative.
   - Complete 30-60%: replace if a better representative exists.
   - Complete <30%, especially Missing >60%: usually remove/replace, with exceptions for unusual/reduced or unique lineages.
3. Use tree distances where available; tree distance overrides simple taxonomy.
4. Preserve model/reference organisms; protect these explicitly using `protected_model_reference_ids.txt`.
5. Prefer replacement of a poor-quality representative by a better equivalent over unnecessary expansion.
6. Do not use trypsin, venom protein, GFP, ADA1, or MdNAT7 presence/absence as selection features.
7. Keep full NCBI taxonomy from the reference workbook for the exported set.

- `criteria_selector.py`: main algorithm. For a benchmark-blind comparison, rerun without --reference-workbook.
- `export_selected_set.py`: assembles the final set using the supplement manifest in the workbook. Adds venom/toxin supplement rows.

## Two modes

### A. Coverage-calibrated mode (recommended for increasing overlap with the previous set)

This mode uses historical core membership as a tie-breaker after explicit model/reference
protection. It is intentionally calibrated against the previous curated set, so the high
overlap is **not an independent validation** of the selection algorithm.

```bash
python -m pip install pandas openpyxl biopython

python criteria_selector.py \
  inputs/993_euks_full_stats.txt \
  inputs/196_euks_species.txt \
  inputs/tree.nex \
  --rescue-ids inputs/reference_rescue_ids_outside_196.txt \
  --protected-ids inputs/protected_model_reference_ids.txt \
  --reference-workbook inputs/final_eukaryotic_representative_set_146_with_venom_supplement_orig.xlsx \
  --target 137 \
  --output selected_core.tsv \
  --audit selection_audit.tsv

python export_selected_set.py \
  --core-tsv selected_core.tsv \
  --reference-workbook inputs/final_eukaryotic_representative_set_146_with_venom_supplement_orig.xlsx \
  --output-xlsx selected_eukaryotic_set.xlsx \
  --output-tsv selected_eukaryotic_set.tsv
```

The 137-core count is for comparison with the supplied final workbook: 136 rows have blank
`Set role`, and one row is labeled `Existing core + venom/toxin supplement`. There are nine
additional venom/toxin-only rows, for 146 unique genomes total.

## Outputs

- `selected_core.tsv`: selected core genomes with BUSCO metrics, tree lineage, quality class,
  historical-core flag, and reason.
- `selection_audit.tsv`: every candidate in the pool, selected or not, with the decision reason.
- `selected_eukaryotic_set.xlsx` / `.tsv`: selected core plus nine venom/toxin-only additions.

### B. Benchmark-blind mode

Omit `--reference-workbook`. This mode uses only the candidate list, tree, rescue IDs,
protected model/reference IDs, and BUSCO scores. It is a more honest test of what the criteria
alone can select, but is expected to overlap less with the historical manually curated set.

```bash
python criteria_selector.py \
  inputs/993_euks_full_stats.txt \
  inputs/196_euks_species.txt \
  inputs/tree.nex \
  --rescue-ids inputs/reference_rescue_ids_outside_196.txt \
  --protected-ids inputs/protected_model_reference_ids.txt \
  --target 137 \
  --output independent_core.tsv \
  --audit independent_audit.tsv
```

## Outputs

- `independent_core.tsv`: independently selected core genomes without calibration against the previous curated set
- `independent_audit.tsv`: every candidate in the pool, selected or not, with the decision reason.

## Interpretation

The selector is a transparent heuristic to prioritize review, not an automated substitute for
the documented biological decisions. Inspect `selection_audit.tsv` for low-BUSCO representatives,
tree lineages, and any candidate that would require manual review. The target count is a
benchmarking convenience; it is not a biological rule.
