# Features, counts and calls: what each file means

All coordinates are 0-based, half-open (BED style). A breakpoint at position *p* lies between
bases *p*-1 and *p*. `<name>` is the catalogue name.

## From regions to segments

kairos-dd finds ORFs that two catalogue contigs share, and regions-v2 merges nearby shared ORFs
into **regions** (`kairos_region_metadata.tsv`). Each region has a `region_type`:

- `_internal_`: the region lies inside the contig.
- `_start_edge_` / `_end_edge_`: the region lies within `edge_padding` bp of a contig end.

`kairos_region_overlaps.tsv` lists the pairs of regions that share ORFs.

`build_features.py` aligns the two regions of each pair (minimap2 asm5, ±`features.window` bp).
The **shared segment** then ends where the two contigs actually stop matching, not at ORF
boundaries. On each contig, a segment ending more than `features.min_flank` bp from the contig
end gets a **breakpoint** there, where shared sequence meets the contig's own sequence:

```
contig:  [ own flank ][======= shared segment =======][ own flank ]
                      ^ BP (side L)                   ^ BP (side R)
```

An internal region has two breakpoints. A region that reaches a contig end has one breakpoint
(or none).

## `results/detect/<name>.features/`

### `events.tsv`: one row per region pair

| column | meaning |
|---|---|
| `event_id` | `EV…` |
| `region_idx_1/2`, `contig_1/2` | the two regions and their contigs |
| `shared_start_i`, `shared_end_i` | the refined shared segment on contig *i* |
| `bp_left_i`, `bp_right_i` | breakpoint IDs on contig *i* (empty at a contig end) |
| `empty_site_i` | empty-site ID on contig *i* (only if it has both breakpoints) |
| `strand` | relative orientation of the two copies (+1 / -1; 0 if not refined) |
| `identity` | matches / alignment length of the shared segment |
| `shared_len` | segment length on contig 1 |
| `refined` | `False`: no alignment passed `min_identity`, so ORF coordinates were kept |
| `region_type_1/2` | region type from regions-v2 |

### `breakpoints.tsv`: one row per breakpoint

`bp_id`, `contig`, `pos` and `side` (L: the segment starts at `pos`; R: the segment ends at
`pos`).

- `window_start` / `window_end`: the junction window, ±`junction_flank` bp around `pos`.
- `flank_bp`: how much own flank is left after the other regions are removed.

If several events place a breakpoint within `features.snap` bp of the same spot, they share it.

### `flanks.bed`: contig, start, end, bp_id

This is the contig's own sequence next to each breakpoint (`flank_len` bp, with every other
region on the contig removed, because reads there may come from elsewhere).

Read depth here is the depth of the **host genome** at the junction. It sets the detection
limit: λ = flank depth × (L − 2a)/L spanning reads are expected if every host cell carries
the segment.

### `empty_sites.tsv`: one row per segment with both breakpoints

| column | meaning |
|---|---|
| `es_id` | `ES…` |
| `contig`, `shared_start`, `shared_end` | the segment that was cut out |
| `left_start`, `right_end` | the flanks joined in the empty-site sequence |
| `junction_pos` | position of the left/right join in that sequence |

An empty site is the contig as it would look **without** the segment. It is built from the
carrier contig. It does not need an "empty" sister contig in the catalogue.

### `junctions.fa`: sequences for the second-pass remapping

`junctions.fa` holds two kinds of sequence:

- **`BP…`**: the junction window of each breakpoint, copied straight from the catalogue. It is
  `junction_flank` bp of own flank plus `junction_flank` bp of shared segment (the *carrier*
  junction).
- **`ES…`**: the empty-site sequence, made by joining the left flank
  `[left_start, shared_start)` to the right flank `[shared_end, right_end)`. This is the
  *non-carrier* junction, a sequence that is in no assembly.

It is used only by `count_features.py --junctions`, for short reads when
`assess.second_pass: true`, because strobealign soft-clips reads instead of writing split
(SA) alignments. The second pass works like this:

1. Collect reads near each segment end that are soft-clipped by at least `overhang` bp, or
   are unmapped with their mate nearby.
2. Realign those reads (mappy, `sr` preset) to all BP and ES sequences at once.
3. Count a read under `empty_remap` only if its best hit is that segment's ES sequence and the
   alignment crosses `junction_pos` by at least `overhang` bp on each side.

The BP windows are competitive decoys. A clipped read that really belongs to the carrier
junction lands on its BP window and is not counted as empty. Long reads don't need this file,
because minimap2 reports the deletion or split alignment directly (`empty_split`).

### `contigs.bed`: contig, start, end

This is each event contig's **own sequence**: the whole contig minus all of its regions and
shared segments. Depth and breadth over it show whether the host contig itself is in a
sample, as opposed to only its mobile segment.

## `results/assess/<sample>/`

`<sample>.breakpoint_counts.tsv` has one row per breakpoint:

| column | meaning |
|---|---|
| `span` | reads aligned across `pos` with at least `overhang` bp on each side |
| `bridge` | read pairs with one mate in the own flank and the other in the segment (paired-end only) |
| `window_reads` | reads with ≥ `assess.window_min_aligned` (100) bp aligned inside ±`assess.window_half` (75) bp of `pos`: the original kairos assess window (75 bp of context + 75 bp of shared sequence), measured on the same catalogue alignments instead of a window-only salmon index. For short reads this forces ≥ 25 bp on each side, so it normally equals `span`. Readout only |
| `flank_depth` | mean depth over `flanks.bed`, at MAPQ ≥ `assess.min_mapq` |
| `mean_read_len`, `overhang` | used for λ |

`<sample>.empty_site_counts.tsv` has one row per empty site:

| column | meaning |
|---|---|
| `empty_split` | a deletion over the segment, or a split alignment that jumps across it |
| `empty_pairs` | read pairs whose insert is the segment length too long |
| `empty_remap` | reads counted in the second pass |
| `empty` | distinct templates with any of the above |

`<sample>.contig_counts.tsv` has one row per event contig, computed over `contigs.bed` at
MAPQ ≥ `assess.min_mapq`:

| column | meaning |
|---|---|
| `own_bp` | length of the contig's own sequence |
| `reads` | reads starting in it |
| `depth` | mean depth over it |
| `breadth` | fraction of its bases covered by at least one read |

Junction evidence uses `assess.junction_min_mapq` (default 1). Reads that lie partly in shared
sequence get modest MAPQ even when they are placed correctly.

## `results/call/`

An **observation** is one (group, timepoint, compartment). `call.pool` sets how its samples
are used:

- `observation` (default): the samples are pooled. Junction reads, empty-site reads, λ and
  contig depth are summed over them, and contig breadth combines as 1 − ∏(1 − bᵢ). This
  gives the lowest detection limit, and short + long reads of one sample give one call.
- `sample`: every sample is called on its own (one `states.tsv` row per sample × segment)
  and nothing is summed across libraries. `insitu.py` then gives each observation the
  strongest state among its samples (present > partial > absent > undetermined) and that
  sample's values. A cell is present if any one library shows the full region and contig,
  and absent if no library has *k* junction reads and at least one is deep enough to rule
  the segment out.

The `samples` column of `states.tsv` lists the samples behind each row.

### `segments.tsv`

A segment is one contig's copy of a shared region. Its row gives `segment_id`, `contig`,
breakpoints, empty site, `shared_start` / `shared_end`, `region_type`, and `events` (the
events it belongs to). If one segment takes part in several events, it appears once.

### `states.tsv`: one row per observation (or sample) × segment

| column | meaning |
|---|---|
| `n` | junction reads (spanning + bridging) summed over the segment's breakpoints |
| `n_left`, `n_right` | the same, per breakpoint (empty where there is no breakpoint) |
| `boundaries` | breakpoints with ≥ *k* (`call.min_reads`) junction reads |
| `n_breakpoints` | breakpoints the segment has (2 = internal) |
| `contig_state` | `present` if the pooled breadth of the contig's own sequence is ≥ `call.contig_min_breadth` and ≥ `call.contig_breadth_ratio` × (1 − e^−depth); else `not_detected` |
| `contig_depth`, `contig_breadth` | pooled values with `call.pool: observation`. Depths add; breadths combine as 1 − ∏(1 − bᵢ) |
| `lam`, `p_miss` | expected spanning reads if all hosts carry the segment, and P(N < *k* \| λ·f) |
| `empty` | empty-site templates |
| `carrier_frac_depth` | n_span / λ |
| `carrier_frac_empty` | n / (n + empty) |
| `n_window` | `window_reads` summed over the segment's breakpoints (and pooled samples) |
| `contig_window_prop` | fraction of the contig's breakpoint windows (over all its segments) with ≥ 1 read |
| `contig_window_pass` | `contig_window_prop` > `call.window_min_prop` (0.9): the original kairos assess detection rule (`get_props.py`) |

The last three reproduce kairos assess for comparison with the paper. They are not used to call
`state`. Two differences remain: windows sit on the mappy-refined breakpoints, not on the raw
asm5 alignment ends, and windows are not clustered at 99%, so reads tied between identical
junctions (MAPQ 0) are not counted.

`state` takes one of four values:

| state | rule |
|---|---|
| `present` | the segment has ≥ `call.min_boundaries` breakpoints, **every** breakpoint has ≥ *k* junction reads, **and** `contig_state` is present |
| `partial` | ≥ *k* junction reads in total, but not present: one boundary only, a contig-end segment, or the host contig not detected |
| `absent` | n < *k* and P(N < *k* \| λ·f) ≤ α |
| `undetermined` | anything else (too shallow), or no breakpoints |

### `insitu_events.tsv`

- A recipient or donor must be `present` at the timepoints that define the event.
- A segment that was `partial` earlier counts as **seen**, so it can't be called new later.
- A `partial` donor never gives a direction.
- With taxonomy, events are filtered on `lineage_relation`, the relation between the
  recipient and donor contigs' lineages:
  - `discordant`: both are named at some rank and the names differ there (`discordant_rank`).
  - `same`: identical at every rank down to species. **Never reported.**
  - `unresolved`: identical as far as both are named, but at least one stops before species
    (this includes unclassified contigs).

  `shared_rank` is the deepest rank at which both lineages agree.

`call.lineage_filter` sets which relations are kept:

- `discordant` (default) keeps only discordant events.
- `not_same` also keeps unresolved events.
- `none` keeps every event.

Removed events are written to `insitu_excluded.tsv`, with a `reason` column.

Other columns:

- `recipient_reads_left` / `_right` and `recipient_contig_breadth` show the evidence behind
  the recipient call.
- `recipient_samples` lists the samples behind the recipient call (with `call.pool: sample`,
  the single sample whose state was used). `samples_present` is the number of `states.tsv` rows
  that call the recipient present at t: always 1 when pooled, and the number of libraries
  that show it on their own with `call.pool: sample`.
- `persistence` can be `yes`, `partial`, `no`, `undetermined` or `NA`.
