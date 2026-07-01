# kairos-lr

A Snakemake pipeline for detecting and classifying horizontal gene transfer (HGT) events from long-read metagenomic assemblies.

> **Status:** pre-release / actively developed. Interfaces, config keys, and output schemas may still change before the first tagged version.

## Overview

kairos-lr takes a long-read metagenomic assembly and identifies contigs that carry regions of likely horizontally-transferred DNA, based on taxonomic discordance between a region and its host contig, and/or overlapping regions shared between taxonomically distinct contigs. It reports per-event, per-group, and donor–recipient summaries, and ships with a companion interactive visualization tool for inspecting HGT groups as clinker-style linear gene maps connected by Bezier arcs.

The pipeline covers:

1. **ORF calling** — Prodigal (meta mode), parallelized via input splitting
2. **Clustering / homology search** — identification of homologous gene pairs across contigs
3. **Taxonomic assignment** — per-contig lineage assignment
4. **HGT detection** — two complementary criteria (below)
5. **HGT classification** — hierarchical classification across seven taxonomic ranks (domain → species), union-find grouping of implicated contigs, and donor–recipient inference
6. **Annotation** — functional annotation of HGT-associated regions

### Detection criteria

- **Criteria 1:** a region's taxonomic assignment differs from its host contig's assignment.
- **Criteria 2:** two taxonomically discordant contigs share an overlapping internal region (via pairwise homology search).

For each pair/region, the pipeline records the most inclusive taxonomic rank at which discordance occurs (`domain`, `phylum`, `class`, `order`, `family`, `genus`, `species`), then groups all implicated contigs globally (not per-rank) using a union-find over every HGT event found at any rank. Where a paired contig in the overlap table shares the region's taxon rather than the host contig's, that paired contig is flagged as the putative donor.

## Installation

```bash
git clone https://github.com/<your-org>/kairos-lr.git
cd kairos-lr
conda env create -f environment.yml
conda activate kairos-lr
```

## Configuration

Copy the example config and fill in the required reference database paths:

```bash
cp config/config.yaml.example config/config.yaml
```

`config/config.yaml` requires (no defaults are provided — the pipeline validates these at startup):

| key | description |
|---|---|
| `ripdb` | path to RIP database |
| `mobileog` | path to mobileOG-db |
| `mobtyper` | path to MOB-typer database |
| `card` | path to CARD database |
| `gtdb_db` | path to GTDB reference |

A sample sheet mapping sample names to assembly FASTAs goes in `config/samples.tsv`.

## Usage

```bash
snakemake --cores <N> --configfile config/config.yaml \
  --config fasta="<sample>_assembly.fasta" threads=<N>
```

For cluster/HPC execution (SLURM), see `workflow/profiles/` (or your cluster's Snakemake profile) — each sample can be submitted as an independent job.

## Output

Results are written to `result/<sample>/`, including:

- Per-level HGT event tables (one per taxonomic rank)
- Global HGT group and group-report tables (`hgt_groups.tsv`, compact one-row-per-group summary with distinct-taxon counts per rank)
- Per-level donor–recipient pair tables
- Functional annotation of HGT-flanking regions

~25 TSVs are produced per sample in total; see `docs/output_schema.md` for full column definitions (WIP).

## Visualization

A companion interactive tool renders HGT groups as stacked linear gene maps (one row per contig) connected by Bezier arcs between homologous gene pairs, in the style of [clinker](https://github.com/gamcil/clinker). Arc color encodes divergence rank (e.g. red for domain-level discordance down to yellow for genus-level), with a rank filter to isolate, for example, only cross-phylum HGT events.

See `viz/README.md` for setup and usage (separate from the main Snakemake pipeline).

## Repository structure

```
kairos-lr/
├── Snakefile                  # entry point; includes rules/, defines rule all
├── config/
│   ├── config.yaml.example
│   └── samples.tsv
├── workflow/
│   ├── rules/
│   │   ├── orf_calling.smk
│   │   ├── clustering.smk
│   │   ├── taxonomy.smk
│   │   ├── hgt_detection.smk
│   │   ├── hgt_classification.smk
│   │   └── annotation.smk
│   └── envs/
├── scripts/                   # Python/R logic called by rules
├── viz/                       # interactive HGT group visualizer
├── test/
│   ├── data/
│   └── expected/
├── docs/
│   ├── hgt_methodology.md
│   └── output_schema.md
└── environment.yml
```

## Citing

If you use kairos-lr in your research, please cite: *(citation pending publication)*

## License

[MIT](LICENSE)
