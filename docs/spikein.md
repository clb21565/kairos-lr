# Spike-in benchmark: planted in situ events

`tools/spikein.py` plants HGT events into real metagenome reads with a known answer. Each event
is a donor genome carrying a mobile element and a recipient genome that is present first
without it (the empty site) and, from the transfer timepoint on, with it inserted at a known
position (with the element's target-site duplication). A **scenario** decides who is present
where and when, and `truth.tsv` records what kairos *should* report, including the cases where
the right answer is "no event".

| Scenario | Set-up (S = observed compartment, I = prior-only, k = transfer) | Expected |
|---|---|---|
| `insitu` | donor everywhere from t0; recipient in S, empty before k, carrier from k | called |
| `carrier_frac` | `insitu`, only a fraction of the recipients carry it | called |
| `immigration` | `insitu`, but the carrier already flows in through I before k | no event |
| `preexisting` | carrier present in S from t0 | no event |
| `rare_bloom` | `insitu`, but the carrier is in S before k below detection (`--rare-depth`) | not tier A |
| `donor_late` | the donor only appears at k | direction unresolved |
| `same_lineage` | donor and recipient of the same species | excluded by the lineage filter |
| `plasmid` | the element is a separate replicon: no junction in the host | no event (undetectable) |

`immigration` and `rare_bloom` are the ones that matter most: they test whether kairos tells a
transfer apart from immigration and from a rare carrier that blooms, which is the claim it
makes.

## Inputs

- `genomes.tsv`: `genome_id  fasta  lineage` (GTDB lineage). Complete, closed genomes not
  present in the background (check with `mash screen` / `sourmash`). Only the longest record
  (the chromosome) is used.
- `elements.tsv`: `element_id  fasta  type  tsd`. Types `is`, `transposon`, `cassette`,
  `prophage`, `ice` (inserted) and `plasmid` (separate replicon); `tsd` = target-site
  duplication length (bp). E.g. ISfinder IS elements and transposons, ICEberg, PLSDB.
- Background reads: one real metagenome per spike-in sample (timepoint x compartment), so the
  background changes over time as real data does.

Each genome and, by default, each element is used by one event only. A second carrier of the
same element is a genuine earlier donor as far as the sequence can tell, which changes the
right answer of `donor_late` and `same_lineage` (`--reuse-elements` allows it; the scorer
reports such donors as `same_element`).

## Steps

```bash
# 1. events: insitu over ranks x depths, the other scenarios once per rank
python tools/spikein.py design --genomes genomes.tsv --elements elements.tsv --out events.tsv \
    --ranks p,c,o,f,g --depths 0.5,1,2,5,10,30 --carrier-fracs 0.05,0.2,0.5,1

# 2. sequences, composition per sample, truth, oracle catalogue, kairos sample sheet
python tools/spikein.py build --events events.tsv --genomes genomes.tsv --elements elements.tsv \
    --out-dir build --timepoints 1,2,3,4,5 --transfer-index 2 \
    --compartments SYS,INF --observed SYS --window 200000

# 3. reads per sample (SLURM array), background prepended
n=$(( $(wc -l < build/samples.tsv) - 2 ))
BUILD=$PWD/build BACKGROUNDS=backgrounds.tsv SIM=iss sbatch --array=0-$n slurm/spikein_reads.sbatch

# 4. kairos on the oracle catalogue (or on a catalogue assembled from the spiked reads),
#    with call.observed_compartments: [SYS]

# 5. score (oracle catalogue: contig names are the planted sequence names)
python tools/spikein.py score --truth build/truth.tsv --events results/call/insitu_events.tsv \
    --excluded results/call/insitu_excluded.tsv --out score.tsv
```

`--window` uses a random window of each genome instead of the whole chromosome, which cuts
the simulated reads (and cost) by the same factor; the window ends become contig ends.

## Outputs of `build`

| File | Content |
|---|---|
| `truth.tsv` | one row per event: scenario, element, donor, recipient, discordant rank, depths, carrier fraction, transfer timepoint, `expect`, contig names, insertion coordinates in donor and recipient |
| `composition.tsv` | per sample: which sequence (donor, empty recipient, carrier, plasmid) at which depth |
| `variants.fa` | every planted sequence |
| `oracle_catalogue.fa` | carriers and donors (`--oracle-empty` adds the empty recipients, as a real catalogue would have them) |
| `oracle_lineages.tsv` | lineages of the oracle contigs in the format of `insitu.py --contig-taxonomy` |
| `junctions.fa` | 2 x `--junction-flank` bp around each recipient breakpoint, to find the planted events in an assembled catalogue |
| `samples.tsv` | kairos sample sheet for the spiked reads |

## Scoring

Each planted event is judged by its recipient: `called_A` / `called_B` (directed, right
recipient, right timepoint), `wrong_time`, `wrong_direction`, `unresolved`, `excluded` or
`none`, plus whether the donor was the planted one (`planted`), another carrier of the same
element (`same_element`) or something else. `correct` compares the outcome with `expect`
(`not_A` accepts anything but a tier A call). Reported events whose recipient was never
planted are false positives.

Report recall per stage when running on an assembled catalogue: assembled (junction sequences
found in the catalogue) -> region detected -> assessed -> called. Running the same reads on the
oracle catalogue separates assembly losses from calling losses.

## Local test

`python .tests/run_spikein.py` runs one event per scenario (insitu at 10x and 2x) on synthetic
genomes through detect features, assess and call on the oracle catalogue, and checks the
expected outcomes (needs the same tools as `.tests/run_sim.py`).
