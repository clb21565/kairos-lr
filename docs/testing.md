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

- `tests/test_build_features.py` checks how breakpoints are refined: a reverse-strand partner,
  a segment that runs to the contig end, breakpoints merged across partners, empty-site
  sequences, and flanks that overlap another region.
- `tests/test_call.py` checks the detection limit (Poisson threshold) and the in situ rules:
  - tier A vs B;
  - unresolved direction;
  - donor first seen later;
  - recipient already seen;
  - same-day import between compartments;
  - observed and prior compartments (the CIWARS setup);
  - date ordering;
  - discordant rank.

## Tier 1: simulated study (local)

```bash
python .tests/run_sim.py                          # strobealign
python .tests/run_sim.py --outdir .tests/sim_bwa --short-aligner bwa-mem2
```

This simulates a ~4 kb gene cassette shared by host1–3. host1 also exists without the
cassette, so the empty site is present in the community. Five samples are made:

- t1: host1 does not carry the cassette.
- t2: half of host1 carries it.
- t3: host1 is at only 3×.

Both short (wgsim) and long (1% error) reads are simulated. The script runs detect up to the
features, then assess and call (taxonomy is skipped: it needs GTDB). It then checks the
results against the simulation:

- counts: junction reads, empty-site reads and pairs, flank depth, and absent hosts;
- states: host1 absent at t1 and present at t2, with a carrier fraction of about 0.5;
- the in situ event: host1 is the recipient at t2, host2 the donor last seen at t1; it is
  tier A and persists.

The tools must be on PATH, without conda. On this VM they are in a micromamba env:

```bash
export MAMBA_ROOT_PREFIX=~/.local/micromamba/root
~/.local/micromamba/bin/micromamba run -n kdev python .tests/run_sim.py
```

`cluster_method=easy-linclust` is used because nucleotide `easy-cluster` needs more than 3 GB.

Planned additions:
- **Second pass:** check empty-site split reads for strobealign.
- **Edge cases:** a donor first seen at the same time point (unresolved direction), and a
  segment at a contig end (one breakpoint only).

## Tier 2: environments and SLURM (HPC)

```bash
conda env create -f environment.yml && conda activate kairos-lr
# 1. local executor (e.g. in an interactive job), per-rule conda envs
python .tests/run_sim.py --use-conda
# 2. the same through SLURM, in a separate output directory
python .tests/run_sim.py --use-conda --profile profiles/slurm --outdir .tests/sim_slurm
```

This checks that every `workflow/envs/*.yaml` solves and that the tool versions behave like
the dev env (for example, strobealign MAPQ and samtools flags). It also checks that the
profile's partition, account and resources are accepted. Run it once after any env change.

## Tier 3: short-read aligner benchmark (HPC)

Open question from the design: should strobealign or bwa-mem2 be the default, and which
`junction_min_mapq`? Script: `tools/aligner_benchmark.py` (env `tools/aligner_benchmark.yaml`,
job `slurm/aligner_benchmark.sbatch`).

1. Families from real catalogue contigs: a recipient backbone with a 5 kb cassette inserted
   (`F_R`), a donor carrying the same cassette (`F_D`) and, for sister families, the recipient
   backbone at 98.5% ANI without the cassette (`F_S`). Features are written from the known
   coordinates, so detection is not part of the test. `--background-mb` adds unused catalogue
   contigs as decoys.
2. Family states: carrier, empty, mixed (50% carriers), sister_carrier, sister_empty.
3. One sample per sample ANI (99.0, 99.5, 99.8, 100% to the catalogue, substitutions plus small
   indels; the cassette stays identical) and depth (2, 5, 10, 30×); 150 bp pairs with wgsim.
4. Map with strobealign, bwa-mem2, bowtie2 (end-to-end) and bowtie2 `--local`. Run
   `count_features.py` at junction MAPQ 0, 1, 5, 10 and 20.
5. Score each recipient segment as `detect_state.py` would (present / absent / undetermined,
   k = 2, alpha = 0.05) against the truth, plus:
   - spanning reads against the expected λ = c(L − 2a)/L (sensitivity);
   - junction reads on segments that should not be present (false junctions);
   - flank depth on sister_carrier against the true depth (reads lost to the sister contig);
   - empty-site templates against the non-carrier depth;
   - index and mapping time and peak memory.

```bash
conda env create -f tools/aligner_benchmark.yaml
CATALOGUE=/projects/MA/HGT/<catalogue>.fa OUTDIR=/projects/MA/HGT/bench \
    sbatch slurm/aligner_benchmark.sbatch
```

Results are written to `$OUTDIR/results/`: `summary.tsv`, `segments.tsv` and `resources.tsv`.
The job log ends with one table per MAPQ. Rerunning skips finished steps. Local smoke test
(random sequence, about 2 min):

```bash
python tools/aligner_benchmark.py --outdir .tests/bench --synthetic --families-per-state 2 \
    --backbone-len 20000 --anis 99.5,100 --depths 5,20 --threads 2
```

On the smoke test, all four aligners call every segment correctly, apart from empty families
at 5×, which are below the detection limit. On the sister_carrier families, though, bowtie2
puts only 0.2–0.4 of the true flank depth on the recipient at MAPQ ≥ 20 (bwa-mem2 1.0,
strobealign 0.9). The real-contig run will show whether that holds.

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

1. Start small: one site and about 10 consecutive weeks, influent and effluent. Use the
   CIWARS setup (`call.observed_compartments: [INF]`, `call.prior_compartments: [INF, EFF]`) and
   set `timepoint` in the sample sheet to the R `ordered_week`.
2. Regenerate the R results for the same weeks with the fixed `get_isHGTs_weekly`
   (`insitu_hgt_weekly.R`).
3. Compare:

   ```bash
   python tools/compare_insitu_r.py --r 240903-phgt-res.csv --rank p \
       --insitu results/call/insitu_events.tsv --call-dir results/call --out cmp_p.tsv
   ```

   This reports:
   - the edges found by both methods;
   - the new-only events, by tier;
   - the R-only edges, each with a reason (`no_event`, `recipient_not_new`,
     `recipient_not_seen`, `donor_not_before`, `other_timepoint`).

   They can't match exactly, because the old method used salmon on edge windows and
   contig-level presence. Use `--contig-map` if the catalogue renamed contigs.
4. If the reasons make sense, run the whole year.

## Long-read SV comparison (optional, Sniffles2)

kairos only sees an event when the catalogue holds both a carrier and a related contig. Long
reads can show insertions and deletions in catalogue contigs that kairos does not model.
`tools/sv_compare.py` runs on Sniffles2 calls against the assess BAMs and sorts each SV into:

- matches a kairos empty site or breakpoint;
- a novel insertion of a known kairos segment (a known mobile segment at a new site);
- novel catalogue sequence;
- no hit.

Many well-supported novel insertions would mean an SV evidence track is worth adding; mostly
matches means it is an independent check of the empty-site counts.

Local check (after `run_sim.py`). It adds a sample where host4 carries the cassette in 60% of
cells:

```bash
MAMBA_ROOT_PREFIX=~/.local/micromamba/root ~/.local/micromamba/bin/micromamba run -n ksv \
    python .tests/run_sv_sim.py
```

On the HPC, from the run directory after assess:

```bash
conda env create -f $REPO/tools/sv_compare.yaml
REPO=/path/to/kairos-lr SAMPLES=samples.tsv sbatch $REPO/slurm/sv_compare.sbatch
```

Notes:
- Sniffles2 skips contigs shorter than 1 Mb unless `--all-contigs` is given; the sbatch sets it.
- `--mosaic --mosaic-include-germline` keeps low-frequency SVs as well as fixed ones.
- Use `vaf_recount` (alt/ref reads recounted from the BAM), not the Sniffles2 FORMAT DR/DV. In
  the simulation, a 3.8 kb insertion in 60% of reads was written as DV=28, DR=0 (VAF 1.0); the
  recount gives 0.59. For deletions both agree with the kairos empty/spanning counts.
