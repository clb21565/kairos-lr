# DESCRIPTION: Snakemake workflow for identifying HGT regions in long-read assemblies
# AUTHOR: connor brown
# LICENSE: GNU General Public License

import os
import sys

from snakemake.utils import min_version
min_version("9.9.0")

wf_v = config.get("wf_version", "0.1.0")

# ── General parameters ──────────────────────────────────────────────────────
THREADS = config.get("threads", 32)

# ── Parameters for identifying HGT (mmseqs/orf clustering) ─────────────────
MINSEQID     = config.get("min_seq_id", 0.99)
MINQUERYCOVER = config.get("min_cov", 0.8)
MINLEN       = config.get("minlen", 1000)

# ── Annotation parameters ───────────────────────────────────────────────────
# RIP (replication initiation protein) database
RIPid = config.get("RIPid", 20)
RIPqc = config.get("RIPqc", 80)
RIPsc = config.get("RIPsc", 80)
# mobileOG-db
MGEid = config.get("MGEid", 20)
MGEqc = config.get("MGEqc", 80)
MGEsc = config.get("MGEsc", 80)
# CARD
ARGid = config.get("ARGid", 90)
ARGqc = config.get("ARGqc", 50)
ARGsc = config.get("ARGsc", 50)

# ── Region detection parameters ─────────────────────────────────────────────
# slop: extra basepairs added either side of a putative HGT region for context
SLOP = config.get("slop", 1000)
# padding: used to classify edge vs internal regions
PADDING = config.get("padding", 100)

# ── Reference database paths ────────────────────────────────────────────────
# These point to large reference data that lives outside the repo — every
# user must configure these for their own system, either via config.yaml or
# --config on the command line. There is no portable default for these.
RIPDB    = config.get("ripdb",    None)
MOBILEOG = config.get("mobileog", None)
MOBTYPER = config.get("mobtyper", None)
CARD     = config.get("card",     None)
GTDB_DB  = config.get("gtdb_db",  None)

# ── Script paths ────────────────────────────────────────────────────────────
# Default to the scripts/ directory shipped alongside this Snakefile, so a
# fresh `git clone` + `snakemake` works with zero config. Override scripts_dir
# (or individual script paths) in config.yaml if running a different version.
SCRIPTS_DIR = config.get("scripts_dir", os.path.join(workflow.basedir, "scripts"))

kairos_dd = config.get("kairos_dd", os.path.join(SCRIPTS_DIR, "kairos-dd-v2.py"))
regions   = config.get("regions",   os.path.join(SCRIPTS_DIR, "regions2.py"))
get_HGT   = config.get("get_HGT",   os.path.join(SCRIPTS_DIR, "get_hgt-v2.py"))

# ── Input ────────────────────────────────────────────────────────────────────
FASTA  = config.get("fasta", "input.fasta")
SAMPLE = os.path.basename(FASTA).rsplit(".fasta", 1)[0]

wildcard_constraints:
    sample   = r"[^/.]+",
    splitnum = r"\d+"


# ── Startup checks ───────────────────────────────────────────────────────────
onstart:
    required_paths = [
        (RIPDB,     "RIP database (config: ripdb)"),
        (MOBILEOG,  "mobileOG database (config: mobileog)"),
        (MOBTYPER,  "Mobtyper database (config: mobtyper)"),
        (CARD,      "CARD database (config: card)"),
        (GTDB_DB,   "GTDB taxonomy database (config: gtdb_db)"),
        (kairos_dd, "kairos-dd script (config: kairos_dd)"),
        (regions,   "regions script (config: regions)"),
        (get_HGT,   "get_HGT script (config: get_HGT)"),
    ]
    for path, name in required_paths:
        if path is None:
            sys.exit(f"Missing required config value for {name}. "
                     f"Set it in config.yaml or via --config.")
        if not os.path.exists(path):
            sys.exit(f"{name} not found: {path}")

    if not os.path.exists(f"{SAMPLE}.fasta"):
        sys.exit(f"Input fasta not found: {SAMPLE}.fasta")


onsuccess:
    from datetime import datetime
    now = datetime.now()
    print(f"kairos-lr {wf_v} completed at {now.strftime('%Y/%m/%d %H:%M:%S')}")
    print("Thank you for using kairos-lr")


# ── Final targets ────────────────────────────────────────────────────────────
rule all:
    input:
        expand("result/{sample}.group_metadata.tsv", sample=SAMPLE),
        expand("result/{sample}-ripdb.tsv", sample=SAMPLE),
        expand("result/{sample}-card.tsv", sample=SAMPLE),
        expand("result/{sample}.assignRes2.tsv", sample=SAMPLE),
        expand("result/{sample}.potential_hgt_regions.taxonomy.tsv", sample=SAMPLE),
        expand("result/{sample}.hgt/", sample=SAMPLE),
        # Uncomment if/when mobileOG annotation is re-enabled:
        # expand("result/{sample}-mobileog.tsv", sample=SAMPLE),


# ── Rule modules ─────────────────────────────────────────────────────────────
include: "workflow/rules/orf_calling.smk"
include: "workflow/rules/clustering.smk"
include: "workflow/rules/hgt_detection.smk"
include: "workflow/rules/taxonomy.smk"
include: "workflow/rules/hgt_classification.smk"
include: "workflow/rules/annotation.smk"
