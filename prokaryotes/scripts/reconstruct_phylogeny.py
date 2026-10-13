#!/usr/bin/env python3
"""Calculate GTDB r232 family-pair patristic distances and family regions.

Definitions:
- family-pair distance: minimum patristic distance between any eligible leaf
  in family A and any eligible leaf in family B, within the same GTDB order.
- median_nn: median across families of each family's nearest other-family distance.
- cutoff: 1.75 * median_nn.
- regions: complete-linkage clusters of families cut at that cutoff.
"""
from __future__ import annotations
import argparse
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_META = ROOT / "inputs/metadata_r232.tsv"
DEFAULT_BAC = ROOT / "inputs/trees/bac120_r232.sp_labels.tree"
DEFAULT_ARC = ROOT / "inputs/trees/ar53_r232.sp_labels.tree"
DEFAULT_OUT = ROOT / "inputs/phylogeny"


def parse_metadata(path: Path) -> pd.DataFrame:
    meta = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    required = {"checkm2_completeness", "checkm2_contamination", "gtdb_taxonomy"}
    missing = required - set(meta.columns)
    if missing:
        raise ValueError(f"Metadata is missing required columns: {sorted(missing)}")
    meta["completeness"] = pd.to_numeric(meta["checkm2_completeness"], errors="coerce")
    meta["contamination"] = pd.to_numeric(meta["checkm2_contamination"], errors="coerce")
    tax = meta["gtdb_taxonomy"].str.split(";")
    for idx, col in [(0, "domain"), (3, "order"), (4, "family")]:
        meta[col] = tax.map(lambda x: x[idx].split("__", 1)[-1] if len(x) > idx and "__" in x[idx] else "")
    meta["species"] = tax.map(lambda x: x[6] if len(x) > 6 and x[6].startswith("s__") else "")
    meta = meta[(meta.completeness >= 80) & (meta.contamination <= 5)
                & meta["order"].ne("") & meta["family"].ne("")
                & meta["species"].ne("")].copy()
    # A species label should map consistently to one order/family in GTDB metadata.
    inconsistent = meta.groupby("species")["order"].nunique().gt(1)
    if inconsistent.any():
        raise ValueError(f"{int(inconsistent.sum())} species labels map to multiple orders")
    inconsistent = meta.groupby("species")["family"].nunique().gt(1)
    if inconsistent.any():
        raise ValueError(f"{int(inconsistent.sum())} species labels map to multiple families")
    return meta


def compute_domain(domain: str, tree_path: Path, meta: pd.DataFrame):
    started = time.time()
    print(f"[{domain}] reading tree {tree_path}", flush=True)
    tree = Phylo.read(str(tree_path), "newick")
    d = meta[meta.domain == domain]
    species_map = {r.species: (r.order, r.family) for r in d.itertuples()}
    tipmap = {id(t): species_map[t.name] for t in tree.get_terminals() if t.name in species_map}
    mapped_species = {t.name for t in tree.get_terminals() if t.name in species_map}
    absent = sorted(set(d.species) - mapped_species)
    print(f"[{domain}] eligible rows={len(d):,}; mapped tree tips={len(tipmap):,}; absent species labels={len(absent):,}", flush=True)

    # Postorder dynamic programming: at each node keep the shortest distance
    # from that node to each family in each order. Combining different child
    # subtrees yields exact minimum cross-family tip-to-tip distances.
    pairdist: dict[tuple[str, str, str], float] = {}
    node_count = 0

    def walk(node):
        nonlocal node_count
        node_count += 1
        if node.is_terminal():
            x = tipmap.get(id(node))
            return {x[0]: {x[1]: 0.0}} if x else {}
        children = []
        for child in node.clades:
            summary = walk(child)
            branch = float(child.branch_length or 0.0)
            children.append({o: {f: dist + branch for f, dist in fs.items()}
                             for o, fs in summary.items()})
        for i, left in enumerate(children):
            for right in children[i + 1:]:
                for order in left.keys() & right.keys():
                    for fa, da in left[order].items():
                        for fb, db in right[order].items():
                            if fa == fb:
                                continue
                            f1, f2 = sorted((fa, fb))
                            key = (order, f1, f2)
                            dist = da + db
                            if dist < pairdist.get(key, float("inf")):
                                pairdist[key] = dist
        merged = {}
        for summary in children:
            for order, families in summary.items():
                dest = merged.setdefault(order, {})
                for family, dist in families.items():
                    if dist < dest.get(family, float("inf")):
                        dest[family] = dist
        return merged

    walk(tree.root)
    print(f"[{domain}] traversed {node_count:,} nodes; family pairs={len(pairdist):,}; {time.time()-started:.1f}s", flush=True)

    family_by_order = defaultdict(set)
    for row in d.itertuples():
        family_by_order[row.order].add(row.family)
    observed_by_order = defaultdict(set)
    for order, f1, f2 in pairdist:
        observed_by_order[order].update((f1, f2))

    coverage_rows, metric_rows, region_rows = [], [], []
    for order, famset in sorted(family_by_order.items()):
        families = sorted(famset)
        n = len(families)
        if n > 1:
            seen = observed_by_order[order]
            coverage_rows.append((domain, order, n, len(seen), len(seen) == n, n * (n - 1) // 2))
        if n == 1:
            metric_rows.append((domain, order, n, np.nan, np.nan, 1, 1.0))
            region_rows.append((domain, order, families[0], 1))
            continue
        matrix = np.zeros((n, n), dtype=float)
        ix = {f: i for i, f in enumerate(families)}
        for (o, f1, f2), dist in pairdist.items():
            if o == order and f1 in ix and f2 in ix:
                matrix[ix[f1], ix[f2]] = matrix[ix[f2], ix[f1]] = dist
        tri = matrix[np.triu_indices(n, 1)]
        if np.any(tri <= 0):
            raise ValueError(f"Missing or non-positive family-pair distance in {domain}/{order}")
        nearest = np.min(np.where(np.eye(n, dtype=bool), np.inf, matrix), axis=1)
        median_nn = float(np.median(nearest))
        cutoff = 1.75 * median_nn
        tree_linkage = linkage(squareform(matrix, checks=True), method="complete")
        labels = fcluster(tree_linkage, t=cutoff, criterion="distance")
        nregions = len(set(labels))
        target_k = min(nregions, int(np.ceil(1.75 * np.sqrt(n))))
        metric_rows.append((domain, order, n, median_nn, cutoff, nregions, target_k))
        for family, label in zip(families, labels):
            region_rows.append((domain, order, family, int(label)))
    return pairdist, coverage_rows, metric_rows, region_rows, len(absent)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--metadata", type=Path, default=DEFAULT_META)
    ap.add_argument("--bacteria-tree", type=Path, default=DEFAULT_BAC)
    ap.add_argument("--archaea-tree", type=Path, default=DEFAULT_ARC)
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()
    for path in (args.metadata, args.bacteria_tree, args.archaea_tree):
        if not path.is_file():
            ap.error(f"Required input not found: {path}")
    args.outdir.mkdir(parents=True, exist_ok=True)
    meta = parse_metadata(args.metadata)
    print(f"Quality-filtered eligible metadata rows: {len(meta):,}", flush=True)

    all_pairs, all_coverage, all_metrics, all_regions = [], [], [], []
    absent_counts = {}
    for domain, path in (("Bacteria", args.bacteria_tree), ("Archaea", args.archaea_tree)):
        pairs, coverage, metrics, regions, absent = compute_domain(domain, path, meta)
        all_pairs.extend((domain, o, f1, f2, dist) for (o, f1, f2), dist in pairs.items())
        all_coverage.extend(coverage)
        all_metrics.extend(metrics)
        all_regions.extend(regions)
        absent_counts[domain] = absent

    pd.DataFrame(all_pairs, columns=["domain", "order", "family1", "family2", "patristic_distance"]).to_csv(
        args.outdir / "family_pairwise_min_distances_reconstructed.tsv", sep="\t", index=False, float_format="%.12g")
    coverage = pd.DataFrame(all_coverage, columns=["domain", "order", "metadata_families", "families_in_pair_table", "all_families_covered", "expected_pair_count"])
    coverage.to_csv(args.outdir / "family_distance_coverage.tsv", sep="\t", index=False)
    metrics = pd.DataFrame(all_metrics, columns=["domain", "order", "raw_families", "median_nn_reconstructed", "threshold_reconstructed", "effective_families_reconstructed", "target_k_reconstructed"])
    metrics.to_csv(args.outdir / "order_phylogeny_metrics_reconstructed.tsv", sep="\t", index=False, float_format="%.12g")
    regions = pd.DataFrame(all_regions, columns=["domain", "order", "family", "effective_region_reconstructed"])
    regions.to_csv(args.outdir / "family_region_assignments_reconstructed.tsv", sep="\t", index=False)

    report = [
        "GTDB r232 PHYLOGENY RECONSTRUCTION",
        "==================================",
        f"Quality-filtered metadata rows: {len(meta):,}",
        f"Bacterial eligible species labels absent from tree: {absent_counts['Bacteria']:,}",
        f"Archaeal eligible species labels absent from tree: {absent_counts['Archaea']:,}",
        f"Family-pair distances: {len(all_pairs):,}",
        f"Orders with more than one family: {len(coverage):,}",
        f"Multi-family orders with all families covered: {int(coverage.all_families_covered.sum()):,}/{len(coverage):,}",
        "Family distance = minimum eligible-leaf patristic distance between two families in the same order.",
        "median_nn = median family nearest-neighbor distance; cutoff = 1.75 × median_nn.",
        "Family regions = complete-linkage clusters cut at that cutoff."
    ]
    (args.outdir / "RECONSTRUCTION_SUMMARY.txt").write_text("\n".join(report) + "\n")
    print("\n".join(report), flush=True)

if __name__ == "__main__":
    main()
