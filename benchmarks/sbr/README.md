# SBR benchmark (Brown et al. 2024, Nat Commun)

Six deeply sequenced Illumina samples from BioProject PRJNA1020581: one 10% hospital-sewage
reactor (0.1 O1) sampled as mixed liquor before (2020-03-20) and after (2020-04-10), and
four feed samples in between. Timepoints and compartments come from the SRA BioSample
metadata (collection_date, distinguishing_features); the sample aliases were mixed up at
sequencing. See `config.yaml` for the run table.

An in situ event: a contig new in the 04-10 mixed liquor (absent 03-20) with its donor
seen in an earlier sample.

```bash
cd benchmarks/sbr
bash fetch_reads.sh                       # ~116 GB gzipped, md5-checked
sbatch --array=0-5 assemble.sbatch        # fastp + metaSPAdes per sample
ls $PWD/assemblies/*.fasta > assemblies.txt
cd ../..                                  # the rest runs from the repo root
LIST=benchmarks/sbr/assemblies.txt OUT=$PWD/benchmarks/sbr/catalogue \
    sbatch slurm/build_catalogue.sbatch
# then run kairos-lr with
#   --configfile config/config.yaml benchmarks/sbr/config.yaml
```
