# kairos-lr

A Snakemake pipeline for detecting and classifying horizontal gene transfer (HGT) events from long-read metagenomic assemblies.

> **Status:** pre-release / actively developed. Interfaces, config keys, and output schemas may still change before the first tagged version.

## Overview

kairos-lr takes a long-read metagenomic assembly and identifies contigs that carry regions of likely horizontally-transferred DNA, based on taxonomic discordance between a region and its host contig, and/or overlapping regions shared between taxonomically distinct contigs. It reports per-event, per-group, and donor–recipient summaries. Do quick assessment of results using [companion](https://github.com/clb21565/hgt_viz.git) streamlit tool. 

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
git clone https://github.com/clb21565/kairos-lr.git
cd kairos-lr
conda env create -f environment.yml
conda activate kairos-lr
```

## Configuration

Copy the example config and fill in the required reference database paths:

```bash
vi config.yaml
```

'config.yaml` requires (no defaults are provided — the pipeline validates these at startup): 

| key | description |
|---|---|
| `ripdb` | path to RIP database |
| `mobileog` | path to mobileOG-db |
| `mobtyper` | path to MOB-typer database |
| `card` | path to CARD database |
| `gtdb_db` | path to GTDB reference |



## Usage

```bash
snakemake --cores <N> --configfile config/config.yaml \
  --config fasta="<sample>_assembly.fasta" threads=<N>
```


## Output

Results are written to `result/<sample>/`, including:

- Per-level HGT event tables (one per taxonomic rank)
- Global HGT group and group-report tables (`hgt_groups.tsv`, compact one-row-per-group summary with distinct-taxon counts per rank)
- Per-level donor–recipient pair tables
- Functional annotation of HGT-flanking regions

~25 TSVs are produced per sample in total; see `docs/output_schema.md` for full column definitions (WIP).

## Visualization

[see HGT_viz repo for current draft streamlight app](https://github.com/clb21565/hgt_viz)

## Repository structure

```
kairos-lr/
├── Snakefile                  # entry point; includes rules/, defines rule all
├── config.yaml
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
├── test/
│   ├── data/
│   └── expected/
├── docs/
│   ├── hgt_methodology.md
│   └── output_schema.md
└── environment.yml
```


## License

[MIT](LICENSE)
