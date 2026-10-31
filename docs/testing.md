# Testing kairos-lr

Testing runs in tiers: fast checks on any laptop or login node first, and heavy runs on the
HPC only once those pass. Each tier answers a different question.

| Tier | Where | Data | Question | Time |
|---|---|---|---|---|
| 0 unit | local | hand-built contigs | does each script compute what it claims? | seconds |
| 1 simulation | local | `.tests/make_sim_data.py` (5 hosts, 5 samples) | does the workflow give the simulated answer end to end? | ~10 s |
| 2 conda / SLURM | HPC | tier 1 data, `--use-conda`, slurm profile | do the per-rule envs and the profile work? | ~30 min (env builds) |
| 3 aligner benchmark | HPC | simulated strains at 99.0–99.5% ANI | which short-read aligner places junction reads best? | hours |
| 4 real subset | HPC | one real catalogue + ~6 samples | runtime, memory, sanity of counts on real data | hours–1 day |
| 5 full study | HPC | the influent/effluent year | in situ calls vs the 240903 R results | days |

## Tier 0: unit tests (local)

```bash
python -m pytest tests/
```

`tests/test_build_features.py` checks how breakpoints are refined: a reverse-strand partner,
a segment that runs to the contig end, breakpoints merged across partners, empty-site
sequences, and flanks that overlap another region. Add a test file per script as scripts
are added (`detect_state.py`, `insitu.py`).

## Tier 1: simulated study (local)

```bash
python .tests/run_sim.py                          # strobealign
python .tests/run_sim.py --outdir .tests/sim_bwa \
    --configfile config/config.yaml bwa.yaml      # bwa-mem2 (bwa.yaml: assess: {short_aligner: bwa-mem2})
```

This simulates a ~4 kb gene cassette shared by host1–3. host1 also exists without the
cassette, so the empty site is present in the community. Five samples are made:

- t1: host1 does not carry the cassette.
- t2: half of host1 carries it.
- t3: host1 is at only 3×.

Both short (wgsim) and long (1% error) reads are simulated. The script runs detect up to the
features, then assess (taxonomy is skipped: it needs GTDB). It then checks the counts against
the simulation: junction reads, empty-site reads and pairs, flank depth, and absent hosts.

The tools must be on PATH, without conda. On this VM they are in a micromamba env:

```bash
export MAMBA_ROOT_PREFIX=~/.local/micromamba/root
~/.local/micromamba/bin/micromamba run -n kdev python .tests/run_sim.py
```

`cluster_method=easy-linclust` is used because nucleotide `easy-cluster` needs more than 3 GB.

Planned additions as stages land:
- **Stage 3:** check present/absent/undetermined per sample against the truth, and the
  expected in situ call: host1 gains the cassette at t2 (tier A, because t1 had empty-site
  evidence), and the t3 call is undetermined.
- **Second pass:** check empty-site split reads for strobealign.
- **Edge cases:** a donor first seen at the same time point (unresolved direction), and a
  segment at a contig end (one breakpoint only).

## Tier 2: environments and SLURM (HPC)

```bash
conda env create -f environment.yml && conda activate kairos-lr
python .tests/run_sim.py --use-conda                         # local executor, per-rule envs
python .tests/run_sim.py --use-conda --profile profiles/slurm --keep   # once envs exist: via SLURM
```

This checks that every `workflow/envs/*.yaml` solves and that the tool versions behave like
the dev env (for example, strobealign MAPQ and samtools flags). It also checks that the
profile's partition, account and resources are accepted. Run it once after any env change.

## Tier 3: short-read aligner benchmark (HPC)

Open question from the design: should strobealign or bwa-mem2 be the default?

1. Take real contigs from the catalogue. Make strain variants at 99.0, 99.5 and 99.8% ANI
   (substitutions and small indels), and insert a known cassette into some of them at
   known breakpoints.
2. Simulate 150 bp pairs at 2–30× depth.
3. Run assess with strobealign, bwa-mem2 and bowtie2.
4. Score each aligner on:
   - junction reads per breakpoint against the expected λ (sensitivity);
   - junction reads assigned to the wrong strain (false junctions);
   - empty-site evidence;
   - runtime and memory.

The choice of `junction_min_mapq` (currently 1) is part of this. Strobealign gives junction
reads with 50 bp of own flank a MAPQ of about 15, so a MAPQ ≥ 20 filter would discard most
correct junction reads.

## Tier 4: real subset (HPC)

One real catalogue (or one site's long-read assemblies) and about 6 samples covering
short-only, long-only and both read types:

- **Detect:** time and memory of each rule. Use these to set `profiles/slurm` resources,
  especially mmseqs clustering and GTDB taxonomy. Also compare `easy-cluster` with
  `easy-linclust` on cluster count and region count.
- **Features:** how many events and breakpoints there are, the share of events refined by
  alignment, and the distribution of shared-segment lengths and identities.
- **Assess:**
  - flank depth against contig coverage from `samtools coverage` (they should agree for
    single-copy contigs);
  - junction reads against flank depth: the ratio should sit near the expected
    (L − 2a)/L when carriers are 100%;
  - how often the empty site is seen.

## Tier 5: full study (HPC)

The influent/effluent year:

1. Run the whole workflow.
2. Compare the in situ events with the regenerated R results (`insitu_hgt_weekly.R`, with
   the fixes). They can't match exactly, because the old method used salmon on edge windows.
3. Report the overlap, then explain the differences by tier and detection state.
