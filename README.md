# kairos-lr

A Snakemake pipeline that detects horizontal gene transfer (HGT) in a metagenomic contig catalogue, tracks each transferred region in sequenced samples, and calls transfers that happened within the sampled system over time (in situ HGT).

> **Status:** pre-release / actively developed. Interfaces, config keys, and output schemas may still change before the first tagged version.

## Overview

The workflow has three stages; each can be run on its own (`snakemake detect | assess | call`).

1. **Detect** (once, on the catalogue): finds regions shared between contigs and classifies them as HGT by taxonomy.
   - **ORF calling:** Prodigal (meta mode), parallelized by splitting the input.
   - **Clustering:** mmseqs2 on nucleotide ORFs (99% identity) to find genes shared between contigs.
   - **derep-detect and regions:** shared ORFs are merged into regions, and region pairs are found.
   - **Taxonomy:** mmseqs2 against GTDB, per contig and per region.
   - **HGT classification:** taxonomic discordance at seven ranks, union-find grouping, and taxonomy-based donor–recipient inference.
   - **Annotation:** CARD, mobileOG-db and RIPDB (all optional).
   - **Junction features:**
     - The two copies of each shared region are aligned, so the shared segment ends where the contigs stop matching, not at ORF boundaries.
     - The ends of each segment are its breakpoints. For each one the workflow records the contig's own sequence next to it (the flank) and the contig with the segment removed (the empty site).
2. **Assess** (per sample): maps every sample's reads to the whole catalogue (competitive mapping) and counts reads at each feature.
   - **Junction reads:** reads spanning a breakpoint, and read pairs bridging it.
   - **Flank depth:** the read depth next to each breakpoint; it sets the detection limit.
   - **Empty-site evidence:** split or deletion reads, long-insert read pairs, and short reads realigned in a second pass.
3. **Call** (all samples):
   - **States:** samples sharing a group, timepoint and compartment are pooled into one observation. Each segment in each observation is called present, partial, absent or undetermined. Present means that both region boundaries have junction reads and the host contig's own sequence is detected. Absent uses a Poisson detection limit.
   - **In situ events:** a recipient whose segment appears at time *t* after being absent, with a donor that carried it earlier. Each event gets an evidence tier, a direction and a persistence flag.

The design (feature definitions, detection limit, in situ rules) is described in `docs/features.md` (every feature file and column) and in the docstrings of each script in `scripts/`.

### Detection criteria (detect)

- **Criteria 1:** a region's taxonomic assignment differs from its host contig's assignment.
- **Criteria 2:** two taxonomically discordant contigs share an overlapping internal region (via pairwise homology search).

For each pair/region, the pipeline records the most inclusive taxonomic rank at which discordance occurs (`domain` … `species`), then groups all implicated contigs globally using a union-find over every HGT event found at any rank. Where a paired contig shares the region's taxon rather than the host contig's, it is flagged as the putative donor.

### In situ calls (call)

Within a group, at timepoint *t* (`scripts/insitu.py`):

- **Present** means the full region is seen: every boundary (breakpoint) of the segment has junction reads, and the contig's own sequence (the contig minus its regions) is covered. Contig-end regions with one boundary, or junction reads without the host contig, are `partial`.
- **Recipient new:**
  - its segment is present in an observed compartment at *t*;
  - it was not seen (present or partial) in any prior compartment at an earlier timepoint;
  - it is not seen in another compartment at *t*, which rules out same-day import (e.g. influent → effluent).
- **Donor first:** the other contig's segment was present at an earlier timepoint; this is what gives the direction. If both appear at *t*, the event is reported with direction `unresolved`.
- **Lineage:** HGT is between lineages. Events whose recipient and donor contigs have the same lineage are never reported. By default (`call.lineage_filter: discordant`), the lineages must also be known to differ at some rank. Removed events go to `insitu_excluded.tsv`.
- **Tier:** `A` if the recipient's segment was confidently absent at the previous timepoint; `B` if that was undetermined (immigration or below detection).
- **Persistence:** whether the segment is still present at the next timepoint.

## Installation

```bash
git clone https://github.com/clb21565/kairos-lr.git
cd kairos-lr
conda env create -f environment.yml     # Snakemake launcher env
conda activate kairos-lr
```

Each rule's tools come from its own conda env in `workflow/envs/`, created on first use with `--use-conda` (on by default in the SLURM profile).

## Input

1. **A dereplicated contig catalogue** (FASTA). Build it beforehand:
   - Long-read assemblies of the study's samples are preferred, because they keep each transferred region and its flanks on one contig. Short-read assemblies work when a study has no long reads.
   - Prefix the contigs with their assembly ID, pool them, and dereplicate at 99%.
   - Contig names must be unique and contain no spaces.
2. **A sample sheet** (TSV; needed for assess and call; see `config/samples.example.tsv`). One row per sequenced sample:

   | column | needed for | example | notes |
   |---|---|---|---|
   | `sample` | always | `CB_2020-08-04_INF` | unique ID |
   | `reads` | assess | `r1.fq.gz,r2.fq.gz` | one path, or two comma-separated for paired-end |
   | `read_type` | assess | `short_pe`, `short_se`, `ont`, `hifi` | picks the aligner and the detection model |
   | `timepoint` | in situ | `2020-08-04` or `12` | only the order matters (numbers, dates or text) |
   | `group` | in situ | `CB` | in situ calls compare samples within one group (e.g. site) |
   | `compartment` | in situ | `INF`, `EFF` | free text; config says which compartments are observed / prior |
   | `assembly` | optional | `CB_0804_INF` | provenance only |

   A sample sequenced with both short and long reads gets two rows (different `sample` IDs, same group, timepoint and compartment). The rows are pooled in the call stage.

## Configuration

Edit `config/config.yaml`. The essentials:

| key | description |
|---|---|
| `catalogue` | catalogue FASTA |
| `samples` | sample sheet (empty: detect only) |
| `gtdb_db` | mmseqs-formatted GTDB database (required for detect taxonomy) |
| `ripdb`, `mobileog`, `card` | protein databases for annotation (each optional) |
| `assess.short_aligner` | `strobealign` (default) or `bwa-mem2` |
| `assess.short_read_length` | read length the strobealign index is built for |
| `assess.long_read_presets` | minimap2 preset per long `read_type` (`map-ont`, `map-hifi`, `lr:hq`) |
| `call.observed_compartments` / `call.prior_compartments` | where recipients are looked for / what counts as seen before (empty = all compartments) |
| `call.min_reads`, `call.alpha`, `call.absent_fraction` | detection limit: ≥ *k* junction reads per boundary; absent when P(< *k* reads \| carried) ≤ α |
| `call.min_boundaries`, `call.contig_min_breadth`, `call.contig_breadth_ratio` | present needs this many boundaries (2 = internal regions) and the host contig's own sequence covered |
| `call.lineage_filter` | `discordant` (default), `not_same` or `none`: which events to keep by recipient/donor lineage |

Every key is documented in the file and validated against `workflow/schemas/`.

## Usage

```bash
snakemake --use-conda -j 8                    # everything that the inputs allow
snakemake --use-conda -j 8 detect             # catalogue only
snakemake --use-conda -j 8 assess             # needs only the features, not taxonomy
snakemake --use-conda -j 8 call
snakemake --profile profiles/slurm            # on SLURM; edit partition/resources in the profile
```

On SLURM, submit from the repo root. Each script reads `$WORKDIR/config.yaml` (default: the repo's `config/config.yaml`):

```bash
WORKDIR=/path/to/run sbatch slurm/kairos_driver.sbatch        # 1-CPU driver job; each rule is its own job (profiles/slurm)
WORKDIR=/path/to/run sbatch slurm/kairos_single_node.sbatch   # alternative: everything in one 64-CPU / 650G job
```

Any config key can be overridden on the command line, e.g. `--config catalogue=cat.fa samples=samples.tsv`.

## Output

| path | contents |
|---|---|
| `results/detect/<name>.hgt/` | HGT tables per rank and criterion, groups, group report, donor–recipient pairs |
| `results/detect/<name>-{card,mobileog,ripdb}.tsv` | annotations |
| `results/detect/<name>.features/events.tsv` | one row per region pair: refined shared segment and breakpoints on each contig |
| `results/detect/<name>.features/breakpoints.tsv`, `flanks.bed`, `empty_sites.tsv`, `junctions.fa`, `contigs.bed` | the features assess counts reads at (see `docs/features.md`) |
| `results/assess/<sample>/<sample>.bam` | reads mapped to the catalogue |
| `results/assess/<sample>/<sample>.breakpoint_counts.tsv` | spanning reads, bridging pairs and flank depth per breakpoint |
| `results/assess/<sample>/<sample>.empty_site_counts.tsv` | empty-site evidence per empty site |
| `results/assess/<sample>/<sample>.contig_counts.tsv` | depth and breadth of each event contig's own sequence (contig minus its regions) |
| `results/assess/<sample>/<sample>.coverage.tsv` | contig coverage and breadth (`samtools coverage`) |
| `results/call/states.tsv` | present / partial / absent / undetermined per observation and segment, with reads per boundary, contig state, expected reads (λ) and carrier fraction |
| `results/call/segments.tsv`, `event_segments.tsv` | segment definitions and which segments make up each event |
| `results/call/insitu_events.tsv` | in situ HGT events: recipient, donor, direction, tier, persistence, lineages and their relation |
| `results/call/insitu_excluded.tsv` | events removed by the lineage filter, with the reason |

## Testing

See `docs/testing.md`:
- `python -m pytest tests/` runs the unit tests.
- `python .tests/run_sim.py` runs a simulated study end to end in about 10 s and checks the results against the simulation.
- Tiers 2–5 cover runs on SLURM with real data.

## Visualization

[see HGT_viz repo for current draft streamlit app](https://github.com/clb21565/hgt_viz)

## Repository structure

```
kairos-lr/
├── Snakefile                  # entry point; config, sample sheet, targets
├── config/
│   ├── config.yaml
│   └── samples.example.tsv
├── workflow/
│   ├── rules/                 # orf_calling, clustering, hgt_detection, taxonomy,
│   │                          # hgt_classification, annotation, features, assess, call
│   ├── envs/                  # per-rule conda envs
│   └── schemas/               # config and sample sheet validation
├── scripts/                   # kairos-dd-v2, regions-v2, get_hgt-v2, build_features,
│                              # count_features, detect_state, insitu
├── profiles/slurm/            # SLURM profile (per-rule resources)
├── slurm/                     # sbatch scripts: driver job, single-node job
├── tests/                     # unit tests (pytest)
├── .tests/                    # simulated study (make_sim_data.py, run_sim.py)
├── docs/
└── environment.yml            # Snakemake launcher env
```

## License

GNU General Public License v3.0 ([GPL-3.0](LICENSE)).
