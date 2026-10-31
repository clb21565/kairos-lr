# DESCRIPTION: Snakemake workflow for detecting HGT regions in a dereplicated contig catalogue
#              (detect) and scoring them in sequenced samples (assess)
# AUTHOR: connor brown
# LICENSE: GNU General Public License

import os
import sys

import pandas as pd
from snakemake.utils import min_version, validate

min_version("9.9.0")

configfile: "config/config.yaml"
validate(config, "workflow/schemas/config.schema.yaml")

wf_v = config.get("wf_version", "0.2.0-dev")

# ── General parameters ──────────────────────────────────────────────────────
THREADS = config["threads"]

# ── Catalogue (one dereplicated FASTA, built by the user) ───────────────────
# Everything from detect is written under results/detect/<catalogue_name>.*
CATALOGUE = config["catalogue"]
CAT = config["catalogue_name"]
DET = f"results/detect/{CAT}"

# ── Samples (only needed for assess / call) ─────────────────────────────────
SAMPLES = None
if config.get("samples"):
    SAMPLES = pd.read_csv(config["samples"], sep="\t", dtype=str, comment="#").fillna("")
    validate(SAMPLES, "workflow/schemas/samples.schema.yaml")
    if SAMPLES["sample"].duplicated().any():
        dups = SAMPLES.loc[SAMPLES["sample"].duplicated(), "sample"].unique()
        sys.exit(f"Duplicate sample IDs in {config['samples']}: {', '.join(dups)}")
    SAMPLES = SAMPLES.set_index("sample", drop=False)

# ── Parameters for identifying HGT (mmseqs ORF clustering, derep-detect) ────
MINSEQID      = config["min_seq_id"]
MINQUERYCOVER = config["min_cov"]
MINLEN        = config["minlen"]
MIN_CONTIG_ORFS = config["min_contig_orfs"]

# ── Region detection parameters ─────────────────────────────────────────────
# merge_distance: candidate HGT ORFs closer than this (bp) are merged into one region
# edge_padding:   regions within this distance (bp) of a contig end are classed as edge regions
# slop:           extra bp either side of a region when extracting its sequence (taxonomy)
MERGE_DISTANCE = config["merge_distance"]
EDGE_PADDING   = config["edge_padding"]
SLOP           = config["slop"]

# ── Annotation parameters ───────────────────────────────────────────────────
# Each database is optional: leave its path empty to skip that annotation.
ANNOTATION_DBS = {
    name: {"db": config[name], **config["annotation"][name]}
    for name in ("ripdb", "mobileog", "card")
    if config.get(name)
}
GTDB_DB = config["gtdb_db"]

# ── Script paths ────────────────────────────────────────────────────────────
SCRIPTS_DIR = config.get("scripts_dir") or os.path.join(workflow.basedir, "scripts")

kairos_dd = os.path.join(SCRIPTS_DIR, "kairos-dd-v2.py")
regions   = os.path.join(SCRIPTS_DIR, "regions-v2.py")
get_HGT   = os.path.join(SCRIPTS_DIR, "get_hgt-v2.py")


# ── Startup checks ───────────────────────────────────────────────────────────
onstart:
    required_paths = [
        (CATALOGUE, "catalogue FASTA (config: catalogue)"),
        (GTDB_DB,   "GTDB taxonomy database (config: gtdb_db)"),
        *[(v["db"], f"{k} database (config: {k})") for k, v in ANNOTATION_DBS.items()],
    ]
    for path, name in required_paths:
        if not os.path.exists(path):
            sys.exit(f"{name} not found: {path}")


onsuccess:
    from datetime import datetime
    now = datetime.now()
    print(f"kairos-lr {wf_v} completed at {now.strftime('%Y/%m/%d %H:%M:%S')}")
    print("Thank you for using kairos-lr")


# ── Final targets ────────────────────────────────────────────────────────────
def detect_targets():
    return [
        f"{DET}.assignRes2.tsv",
        f"{DET}.potential_hgt_regions.taxonomy.tsv",
        f"{DET}.hgt/",
        *[f"{DET}-{name}.tsv" for name in ANNOTATION_DBS],
    ]

rule all:
    input:
        detect_targets()

rule detect:
    input:
        detect_targets()


# ── Rule modules ─────────────────────────────────────────────────────────────
include: "workflow/rules/orf_calling.smk"
include: "workflow/rules/clustering.smk"
include: "workflow/rules/hgt_detection.smk"
include: "workflow/rules/taxonomy.smk"
include: "workflow/rules/hgt_classification.smk"
include: "workflow/rules/annotation.smk"
