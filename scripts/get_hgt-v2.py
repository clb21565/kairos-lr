#!/usr/bin/env python3
"""
HGT (Horizontal Gene Transfer) Analysis Script

Criteria 1: Region has different taxonomy than its host contig.
Criteria 2: Internal regions shared between contigs with discordant taxonomy.

Union-Find grouping runs once across regions associated with HGT. Note that right now this step leads to very large clusters in some cases, need to explore how to stop the exploding connections. 

Can adjust mmseqs2 taxonomy settings minimum number of fragments, number of returned fragments, min agreed, and in general the minimum region length to consider for HGT prediction. 
"""

import argparse
import os
from collections import Counter

import pandas as pd
import networkx as nx

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

LINEAGE_COLS = ["d", "p", "c", "o", "f", "g", "s"]
REGION_TAXONOMY_COLS = [
    "region", "taxid", "level", "taxa",
    "fragments", "returned", "agreed", "percent_agreed", "lineage",
]
CONTIG_TAXONOMY_COLS = [
    "contig", "taxid", "level", "taxa",
    "fragments", "returned", "agreed", "percent_agreed", "lineage",
]
LEVEL_RANK = {
    "dHGT": "d", "pHGT": "p", "cHGT": "c",
    "oHGT": "o", "fHGT": "f", "gHGT": "g", "sHGT": "s",
}

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="HGT analysis: detect horizontal gene transfer events by comparing "
                    "region and contig taxonomy.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--path", "-p", default="../",
                        help="Base path where input TSV files are located.")
    parser.add_argument("--sample", "-s", required=True,
                        help="Sample prefix (e.g. '008-AR-23_assembly').")
    parser.add_argument("--min-fragments", type=int, default=1)
    parser.add_argument("--min-returned",  type=int, default=1)
    parser.add_argument("--min-agreed",    type=int, default=1)
    parser.add_argument("--min-region-len", type=int, default=1000,
                        help="Minimum length of internal regions.")
    parser.add_argument("--out-dir", "-o", default=None,
                        help="Output directory. Defaults to --path.")
    return parser.parse_args()

# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def read_taxonomy(filepath: str, id_col: str) -> pd.DataFrame:
    cols = REGION_TAXONOMY_COLS if id_col == "region" else CONTIG_TAXONOMY_COLS
    df = pd.read_csv(filepath, sep="\t", header=None, names=cols)
    df["lineage"] = df["lineage"].str.replace(r"[a-z]_", "", regex=True)
    lineage_split = df["lineage"].str.split(";", expand=True)
    lineage_split.columns = LINEAGE_COLS[: lineage_split.shape[1]]
    df = pd.concat([df.drop(columns=["lineage"]), lineage_split], axis=1)
    return df

# ---------------------------------------------------------------------------
# HGT classification
# ---------------------------------------------------------------------------
# connor notes: this function will silently fail if taxonomy strings are not consistent, as it expects for each increasing level of taxonomy that the previous ones are identical, but only up to 1 level below. so, if they have same genus names but belong to different families, it give an incorrect result at the genus species level. 
def hgt_levels(df: pd.DataFrame, suffix_a: str, suffix_b: str) -> dict:
    """
    Hierarchical HGT classification at each taxonomic rank.
    Each level only contains events where all higher ranks agree
    but this rank differs.
    """
    ranks  = ["d", "p", "c", "o", "f", "g", "s"]
    labels = ["dHGT", "pHGT", "cHGT", "oHGT", "fHGT", "gHGT", "sHGT"]
    results = {}
    for i, (label, rank) in enumerate(zip(labels, ranks)):
        a_col = f"{rank}_{suffix_a}"
        b_col = f"{rank}_{suffix_b}"
        if a_col not in df.columns or b_col not in df.columns:
            results[label] = pd.DataFrame() #stores an empty dataframe if either column is not in the df. this is a sanity check for a problem that doesn't exist. 
            continue
        valid = df[
            df[a_col].notna() & df[b_col].notna() &
            (df[a_col] != "None") & (df[b_col] != "None")
        ]
        if i == 0:
            subset = valid[valid[a_col] != valid[b_col]]
        else:
            prev = ranks[i - 1]
            prev_a, prev_b = f"{prev}_{suffix_a}", f"{prev}_{suffix_b}"
            subset = valid[
                (valid[prev_a] == valid[prev_b]) &
                (valid[a_col]  != valid[b_col])
            ]
        results[label] = subset.copy()
    return results

# ---------------------------------------------------------------------------
# Union-Find
# ---------------------------------------------------------------------------

def _union_find(pairs: list) -> dict:
    """Path-compressed union-find. Returns node → canonical root."""
    parent = {}

    def find(x):
        if x not in parent:
            parent[x] = x
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for a, b in pairs:
        union(a, b)
    return {x: find(x) for x in parent}

def make_hgt_groups(allHGT2: pd.DataFrame,
                    region_metadata: pd.DataFrame) -> pd.DataFrame:
    """
    Single global union-find over ALL HGT-implicated regions (any level).

    Inputs
    ------
    allHGT2        : union of all criteria-2 HGT overlap rows (any level)
    region_metadata: kairos_region_metadata.tsv

    Returns
    -------
    One row per region_idx that participates in any HGT, with:
      hgt_group_id, hgt_group_size, hgt_group_root
      + all region_metadata columns
    """
    id1, id2 = "region_idx_contig1", "region_idx_contig2"
    pairs = (
        allHGT2[[id1, id2]]
        .dropna()
        .drop_duplicates()
        .itertuples(index=False, name=None)
    )
    pairs = list(pairs)

    node_to_root = _union_find(pairs)
    root_counts  = Counter(node_to_root.values())
    root_rank    = {root: rank
                    for rank, (root, _) in enumerate(root_counts.most_common(), start=1)}

    groups_df = pd.DataFrame([
        {
            "region_idx":     node,
            "hgt_group_id":   root_rank[root],
            "hgt_group_size": root_counts[root],
            "hgt_group_root": root,
        }
        for node, root in node_to_root.items()
    ])

    # Join region metadata
    if not region_metadata.empty and "region_idx" in region_metadata.columns:
        groups_df = groups_df.merge(region_metadata, on="region_idx", how="left")

    # Recover contig from overlap table if region_metadata didn't supply it
    if "contig" not in groups_df.columns:
        idx_map = pd.concat([
            allHGT2[[id1, "contig1"]].rename(columns={id1: "region_idx", "contig1": "contig"}),
            allHGT2[[id2, "contig2"]].rename(columns={id2: "region_idx", "contig2": "contig"}),
        ]).drop_duplicates(subset="region_idx")
        groups_df = groups_df.merge(idx_map, on="region_idx", how="left")

    return groups_df.sort_values(["hgt_group_id", "region_idx"]).reset_index(drop=True)

# ---------------------------------------------------------------------------
# Donor-recipient annotation (Criteria 1)
# ---------------------------------------------------------------------------

def annotate_donor_recipient(hgt1_df: pd.DataFrame,
                             overlaps_annotated: pd.DataFrame,
                             rank: str) -> pd.DataFrame:
    """
    For a given criteria-1 HGT level, find cases where an overlapping region
    from a second contig has the SAME taxonomic assignment as the HGT region
    (i.e. matches the region taxonomy, not the host contig taxonomy).
    This second contig is the putative donor genome.

    Parameters
    ----------
    hgt1_df           : one level of criteria-1 HGT (e.g. oHGT1)
    overlaps_annotated: full overlap table with contig1/contig2 taxonomy columns
    rank              : taxonomic rank letter for this level (e.g. 'o')

    Returns
    -------
    DataFrame with columns:
        recipient, recipient_region_idx,
        donor, donor_region_idx,
        <rank>_region  (the shared taxon confirming the donor)
    """
    if hgt1_df.empty or overlaps_annotated.empty:
        return pd.DataFrame()

    rank_region = f"{rank}_region"
    rank_contig2 = f"{rank}_contig2"

    if "region_idx" not in hgt1_df.columns:
        return pd.DataFrame()
    if rank_region not in hgt1_df.columns or rank_contig2 not in overlaps_annotated.columns:
        return pd.DataFrame()

    merged = pd.merge(
        hgt1_df[["contig", "region_idx", rank_region]].drop_duplicates(),
        overlaps_annotated[["region_idx_contig1", "contig1",
                             "region_idx_contig2", "contig2", rank_contig2]].drop_duplicates(),
        left_on="region_idx", right_on="region_idx_contig1"
    )

    # The paired contig2 has the same taxon as the HGT region → putative donor
    has_donor = merged[merged[rank_contig2] == merged[rank_region]].copy()

    if has_donor.empty:
        return pd.DataFrame()

    result = has_donor[[
        "contig", "region_idx",
        "contig2", "region_idx_contig2",
        rank_region,
    ]].drop_duplicates()
    result.columns = [
        "recipient", "recipient_region_idx",
        "donor", "donor_region_idx",
        f"{rank}_shared_taxon",
    ]
    return result.reset_index(drop=True)

# ---------------------------------------------------------------------------
# Group report  (compact — one row per group)
# ---------------------------------------------------------------------------

def build_group_report(groups_df: pd.DataFrame) -> pd.DataFrame:
    """
    Compact summary: one row per hgt_group_id.

    Reported per group
    ------------------
    hgt_group_id, hgt_group_size,
    num_regions, num_contigs,
    total_region_len, avg_region_len, max_region_len,
    region_types  (value-counts string, e.g. _internal_:12,_end_edge_:3),
    num_distinct_domains, num_distinct_phyla, num_distinct_classes,
    num_distinct_orders,  num_distinct_families, num_distinct_genera,
    num_distinct_species,
    Full per-region detail stays in hgt_groups.tsv.
    """
    if groups_df.empty:
        return pd.DataFrame()

    def _vc_str(s):
        vc = s.value_counts()
        return ",".join(f"{k}:{v}" for k, v in vc.items())

    def _nunique(s):
        return s.dropna().nunique()

    agg_spec = {
        "hgt_group_size":   ("hgt_group_size", "first"),
        "num_regions":      ("region_idx",      "nunique"),
        "num_contigs":      ("contig",          "nunique"),
        "total_region_len": ("region_len",      "sum"),
        "avg_region_len":   ("region_len",      "mean"),
        "max_region_len":   ("region_len",      "max"),
        "region_types":     ("region_type",     _vc_str),
    }

    # Taxonomic diversity columns — only those present in groups_df
    for rank, alias in [
        ("d", "num_distinct_domains"),
        ("p", "num_distinct_phyla"),
        ("c", "num_distinct_classes"),
        ("o", "num_distinct_orders"),
        ("f", "num_distinct_families"),
        ("g", "num_distinct_genera"),
        ("s", "num_distinct_species"),
    ]:
        # taxonomy cols in groups_df come from region_metadata_and_taxa via join;
        # they may appear as d_region or just d depending on merge path
        for col in [f"{rank}_region", rank]:
            if col in groups_df.columns:
                agg_spec[alias] = (col, _nunique)
                break

    valid_spec = {
        k: pd.NamedAgg(column=v[0], aggfunc=v[1])
        for k, v in agg_spec.items()
        if v[0] in groups_df.columns
    }

    return (
        groups_df
        .groupby("hgt_group_id", as_index=True)
        .agg(**valid_spec)
        .reset_index()
        .sort_values("hgt_group_size", ascending=False)
        .reset_index(drop=True)
    )

# ---------------------------------------------------------------------------
# Graph stats - experimental right now 
# ---------------------------------------------------------------------------

def compute_graph_stats(hgt_df: pd.DataFrame,
                        donor_col: str, recipient_col: str,
                        agg_col: str = "region_len") -> pd.DataFrame:
    if hgt_df.empty or agg_col not in hgt_df.columns:
        return pd.DataFrame(columns=["node", "page_rank", "indegree", "outdegree"])
    fg = (
        hgt_df[[donor_col, recipient_col, agg_col]]
        .drop_duplicates()
        .groupby([donor_col, recipient_col], as_index=False)[agg_col]
        .mean()
    )
    g = nx.from_pandas_edgelist(fg, source=donor_col, target=recipient_col,
                                create_using=nx.DiGraph())
    pr = nx.pagerank(g)
    return pd.DataFrame({
        "node":      list(pr),
        "page_rank": list(pr.values()),
        "indegree":  [dict(g.in_degree())[n]  for n in pr],
        "outdegree": [dict(g.out_degree())[n] for n in pr],
    })

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    pth     = args.path
    sp      = args.sample
    out_dir = args.out_dir or pth
    os.makedirs(out_dir, exist_ok=True)

    mf, mr, ma = args.min_fragments, args.min_returned, args.min_agreed
    min_region_len = args.min_region_len

    # ------------------------------------------------------------------
    # 1. Load data
    # ------------------------------------------------------------------

    print(f"[1/4] Reading region taxonomy")
    region_taxonomy = read_taxonomy(
        os.path.join(pth, f"{sp}.potential_hgt_regions.taxonomy.tsv"), id_col="region")
    region_taxonomy["region"] = region_taxonomy["region"].str.replace(
        r"::(.*?)$", "", regex=True)

    print(f"[2/4] Reading contig taxonomy")
    contig_taxonomy = read_taxonomy(
        os.path.join(pth, f"{sp}.assignRes2.tsv"), id_col="contig")

    print(f"[3/4] Reading region metadata")
    region_metadata = pd.read_csv(
        os.path.join(pth, f"{sp}.kairos_region_metadata.tsv"), sep="\t")

    print(f"[4/4] Reading overlaps")
    overlaps = pd.read_csv(
        os.path.join(pth, f"{sp}.kairos_region_overlaps.tsv"), sep="\t")


    # ------------------------------------------------------------------
    # 2. Criteria 1 — region taxonomy != contig taxonomy
    # ------------------------------------------------------------------

    region_metadata_and_taxa = (
        region_metadata
        .merge(region_taxonomy, left_on="region_title", right_on="region")
        .merge(contig_taxonomy, on="contig", suffixes=("_region", "_contig"))
    )

    filtered = region_metadata_and_taxa[
        (region_metadata_and_taxa["taxid_region"] > 1) &
        (region_metadata_and_taxa["taxid_contig"] > 1) &
        (region_metadata_and_taxa["fragments_region"] >= mf) &
        (region_metadata_and_taxa["returned_region"]  >= mr) &
        (region_metadata_and_taxa["agreed_region"]    >= ma)
    ]

    print("\n[Criteria 1: region vs. contig HGT]")
    _hgt1 = hgt_levels(filtered, suffix_a="region", suffix_b="contig")

    dHGT1, pHGT1, cHGT1, oHGT1, fHGT1, gHGT1, sHGT1 = (
        _hgt1[k] for k in ["dHGT","pHGT","cHGT","oHGT","fHGT","gHGT","sHGT"])

    for label, df in _hgt1.items():
        print(f"  {label}: {len(df)} instances")

    allHGT1 = pd.concat(_hgt1.values(), ignore_index=True).drop_duplicates()
    print(f"  all criteria 1 HGT: {len(allHGT1)} predicted")

    # Donor-recipient annotation is wired in after overlaps_annotated is built (below).

    # ------------------------------------------------------------------
    # 3. Criteria 2 — overlapping internal regions between taxonomically
    #    discordant contigs
    # ------------------------------------------------------------------

    filtered_contig_taxonomy = contig_taxonomy[
        (contig_taxonomy["taxid"]     > 1) &
        (contig_taxonomy["fragments"] >= mf) &
        (contig_taxonomy["returned"]  >= mr) &
        (contig_taxonomy["agreed"]    >= ma)
    ]

    overlaps_annotated = (
        overlaps
        .merge(filtered_contig_taxonomy, left_on="contig1", right_on="contig")
        .merge(filtered_contig_taxonomy, left_on="contig2", right_on="contig",
               suffixes=("_contig1", "_contig2"))
    )

    print("\n[Criteria 2: contig vs. contig overlap HGT]")
    _hgt2 = hgt_levels(overlaps_annotated, suffix_a="contig1", suffix_b="contig2")

    dHGT2, pHGT2, cHGT2, oHGT2, fHGT2, gHGT2, sHGT2 = (
        _hgt2[k] for k in ["dHGT","pHGT","cHGT","oHGT","fHGT","gHGT","sHGT"])

    for label, df in _hgt2.items():
        print(f"  {label} (overlap): {len(df)} instances")

    allHGT2 = pd.concat(_hgt2.values(), ignore_index=True).drop_duplicates()
    print(f"  all criteria 2 HGT: {len(allHGT2)} instances")

    # ------------------------------------------------------------------
    # 3b. Donor-recipient annotation for each criteria-1 HGT level
    # ------------------------------------------------------------------

    print("\n[Criteria 1: donor-recipient annotation]")
    donor_recipient_frames = {}
    for label, df in _hgt1.items():
        rank = LEVEL_RANK[label]
        dr = annotate_donor_recipient(df, overlaps_annotated, rank)
        donor_recipient_frames[label] = dr
        print(f"  {label}: {len(dr)} putative donor-recipient pairs")

    allDR = pd.concat(donor_recipient_frames.values(), ignore_index=True).drop_duplicates()
    print(f"  all levels combined: {len(allDR)} pairs")

    # ------------------------------------------------------------------
    # 4. Global union-find over ALL HGT-implicated regions (any level)
    # ------------------------------------------------------------------

    print("\n[Union-Find grouping across all HGT levels]")
    hgt_groups = make_hgt_groups(allHGT2, region_metadata)
    n_groups   = hgt_groups["hgt_group_id"].nunique() if not hgt_groups.empty else 0
    print(f"  {n_groups} groups from {len(allHGT2)} overlap pairs")

    # Enrich hgt_groups with contig taxonomy so the report can count distinct taxa per group.
    # Select only contig + the 7 lineage rank columns (d,p,c,o,f,g,s) that are present.
    if "contig" in hgt_groups.columns:
        rank_cols_present = [r for r in LINEAGE_COLS if r in contig_taxonomy.columns]
        ctax_slim = (
            contig_taxonomy[["contig"] + rank_cols_present]
            .drop_duplicates(subset="contig")   # one row per contig, no duplicate col issue
        )
        hgt_groups = hgt_groups.merge(ctax_slim, on="contig", how="left", suffixes=("", "_ctax"))

    hgt_report = build_group_report(hgt_groups)

    # ------------------------------------------------------------------
    # 5. oHGT graph stats (criteria 1, order-level)
    # ------------------------------------------------------------------

    if not oHGT1.empty and "o_region" in oHGT1.columns:
        oHGT1_graph = compute_graph_stats(
            oHGT1[["o_region", "o_contig", "region_title", "region_len"]].drop_duplicates(),
            donor_col="o_region", recipient_col="o_contig", agg_col="region_len")
        oHGT1_graph.to_csv(
            os.path.join(out_dir, f"{sp}.oHGT1.graph_stats.tsv"), sep="\t", index=False)

    # ------------------------------------------------------------------
    # 6. Save
    # ------------------------------------------------------------------

    print()
    to_save = {
        # Criteria 1 — region vs contig
        f"{sp}.allHGT1": allHGT1,
        f"{sp}.dHGT1":   dHGT1,
        f"{sp}.pHGT1":   pHGT1,
        f"{sp}.cHGT1":   cHGT1,
        f"{sp}.oHGT1":   oHGT1,
        f"{sp}.fHGT1":   fHGT1,
        f"{sp}.gHGT1":   gHGT1,
        f"{sp}.sHGT1":   sHGT1,
        # Criteria 2 — contig vs contig overlap
        f"{sp}.allHGT2": allHGT2,
        f"{sp}.dHGT2":   dHGT2,
        f"{sp}.pHGT2":   pHGT2,
        f"{sp}.cHGT2":   cHGT2,
        f"{sp}.oHGT2":   oHGT2,
        f"{sp}.fHGT2":   fHGT2,
        f"{sp}.gHGT2":   gHGT2,
        f"{sp}.sHGT2":   sHGT2,
        # Global HGT groups (long-form, one row per region) and compact report
        f"{sp}.hgt_groups":  hgt_groups,
        f"{sp}.hgt_report":  hgt_report,
        # Donor-recipient pairs per criteria-1 HGT level
        f"{sp}.allDR":       allDR,
        f"{sp}.dHGT1_DR":    donor_recipient_frames.get("dHGT", pd.DataFrame()),
        f"{sp}.pHGT1_DR":    donor_recipient_frames.get("pHGT", pd.DataFrame()),
        f"{sp}.cHGT1_DR":    donor_recipient_frames.get("cHGT", pd.DataFrame()),
        f"{sp}.oHGT1_DR":    donor_recipient_frames.get("oHGT", pd.DataFrame()),
        f"{sp}.fHGT1_DR":    donor_recipient_frames.get("fHGT", pd.DataFrame()),
        f"{sp}.gHGT1_DR":    donor_recipient_frames.get("gHGT", pd.DataFrame()),
        f"{sp}.sHGT1_DR":    donor_recipient_frames.get("sHGT", pd.DataFrame()),
    }

    for stem, df in to_save.items():
        out = os.path.join(out_dir, f"{stem}.tsv")
        df.to_csv(out, sep="\t", index=False)
        print(f"Saved {stem} ({len(df)} rows) → {out}")

    print("\nDone.")

if __name__ == "__main__":
    main()
