#!/usr/bin/env python3
"""
from kairos, extract hgt regions
"""

import argparse
import pandas as pd
import sys


def load_data(merged_bed_file, overlaps_file, fai_file, edge_padding=1000):
    merged_bed = pd.read_csv(
        merged_bed_file, sep='\t', header=None,
        names=['contig', 'start', 'stop', 'orf_cts', 'orfs_encompassed']
    )
    overlaps = pd.read_csv(overlaps_file, sep='\t')
    overlaps['orf1'] = overlaps['contig1'] + '_' + overlaps['identical_orf_position_contig1'].astype(str)
    overlaps['orf2'] = overlaps['contig2'] + '_' + overlaps['identical_orf_position_contig2'].astype(str)
    fai = pd.read_csv(
        fai_file, sep='\t', header=None,
        names=['contig', 'contig_length', 'offset', 'linebases', 'qualoffset']
    )[['contig', 'contig_length']]
    return merged_bed, overlaps, fai, edge_padding


def get_regions_for_contig(contig, merged_bed, fai, edge_padding):
    idf = merged_bed[merged_bed['contig'] == contig].copy()
    idf['region_idx'] = [f"{contig}_region_{i+1}" for i in range(len(idf))]
    idf['region_len'] = idf['stop'] - idf['start']
    idf = idf.merge(fai, on='contig')
    idf['end_edge_cutoff'] = idf['contig_length'] - edge_padding
    idf['start_edge_cutoff'] = 1 + edge_padding

    def categorize_region(row):
        if row['stop'] > row['end_edge_cutoff']:
            return '_end_edge_'
        elif row['start'] < row['start_edge_cutoff']:
            return '_start_edge_'
        else:
            return '_internal_'

    idf['region_type'] = idf.apply(categorize_region, axis=1)
    return idf


def get_all_regions(merged_bed, fai, edge_padding):
    contigs = merged_bed['contig'].unique()
    region_metadata_list = []
    for contig in contigs:
        region_data = get_regions_for_contig(contig, merged_bed, fai, edge_padding)
        region_metadata_list.append(region_data)
    region_metadata = pd.concat(region_metadata_list, ignore_index=True)
    region_metadata['region_title'] = (
        region_metadata['region_idx'] +
        region_metadata['region_type'] +
        'length=' +
        region_metadata['region_len'].astype(str)
    )
    return region_metadata


def get_region_overlaps(region_metadata, overlaps):
    region_orfs = region_metadata[['contig', 'region_idx', 'orfs_encompassed']].copy()
    region_orfs = region_orfs.assign(
        orfs_encompassed=region_orfs['orfs_encompassed'].str.split(',')
    ).explode('orfs_encompassed')

    matched_orfs = overlaps[['orf1', 'orf2']].copy()

    region_overlaps = matched_orfs.merge(
        region_orfs, left_on='orf1', right_on='orfs_encompassed', how='inner'
    ).rename(columns={'orf1': 'orf_region1'})

    region_overlaps = region_overlaps.merge(
        region_orfs, left_on='orf2', right_on='orfs_encompassed',
        how='inner', suffixes=('_contig1', '_contig2')
    ).rename(columns={'orf2': 'orf_region2'})

    region_overlaps = region_overlaps[[
        'contig_contig1', 'contig_contig2',
        'region_idx_contig1', 'region_idx_contig2',
        'orf_region1', 'orf_region2'
    ]].rename(columns={
        'contig_contig1': 'contig1',
        'contig_contig2': 'contig2'
    })
    return region_overlaps


def filter_reciprocals(hgt_df):
    hgt_df = hgt_df.copy()
    hgt_df['region_key'] = hgt_df['region_idx_contig1'] + '___' + hgt_df['region_idx_contig2']
    pairs = hgt_df[['region_idx_contig1', 'region_idx_contig2']].drop_duplicates()
    pairs['sorted_pair'] = pairs.apply(
        lambda row: tuple(sorted([row['region_idx_contig1'], row['region_idx_contig2']])),
        axis=1
    )
    filtered_pairs = pairs[~pairs.duplicated(subset='sorted_pair', keep='first')].copy()
    filtered_pairs['region_key'] = (
        filtered_pairs['region_idx_contig1'] + '___' + filtered_pairs['region_idx_contig2']
    )
    filtered_df = hgt_df[hgt_df['region_key'].isin(filtered_pairs['region_key'])].copy()
    return filtered_df.drop(columns=['region_key'])


def create_bed_and_extract(region_metadata, fasta_file, bed_output, fasta_output,
                           slop=0, fai_file=None):
    import subprocess
    import os

    bed_data = region_metadata[['contig', 'start', 'stop', 'region_title']].copy()
    bed_data.to_csv(bed_output, sep='\t', header=False, index=False)
    print(f"\nCreated BED file: {bed_output}", file=sys.stderr)
    print(f"  Total regions: {len(bed_data)}", file=sys.stderr)

    try:
        subprocess.run(['bedtools', '--version'], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("WARNING: bedtools not found. Skipping sequence extraction.", file=sys.stderr)
        return

    if not os.path.exists(fasta_file):
        print(f"WARNING: FASTA file not found: {fasta_file}", file=sys.stderr)
        return

    if slop > 0:
        if not fai_file or not os.path.exists(fai_file):
            print(f"WARNING: FAI file required for slop but not found.", file=sys.stderr)
            return
        print(f"\nAdding {slop} bp padding (bedtools slop)...", file=sys.stderr)

    print(f"\nExtracting sequences with bedtools...", file=sys.stderr)
    try:
        if slop > 0:
            slop_proc = subprocess.Popen(
                ['bedtools', 'slop', '-i', bed_output, '-g', fai_file, '-b', str(slop)],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            getfasta_proc = subprocess.Popen(
                ['bedtools', 'getfasta', '-fi', fasta_file, '-bed', '-', '-name', '-fo', fasta_output],
                stdin=slop_proc.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            slop_proc.stdout.close()
            _, stderr = getfasta_proc.communicate()
            if getfasta_proc.returncode != 0:
                raise subprocess.CalledProcessError(getfasta_proc.returncode, 'getfasta', stderr=stderr)
        else:
            subprocess.run(
                ['bedtools', 'getfasta', '-fi', fasta_file, '-bed', bed_output,
                 '-name', '-fo', fasta_output],
                capture_output=True, text=True, check=True
            )

        with open(fasta_output, 'r') as f:
            seq_count = sum(1 for line in f if line.startswith('>'))
        print(f"  Extracted {seq_count} sequences → {fasta_output}", file=sys.stderr)

    except subprocess.CalledProcessError as e:
        print(f"ERROR running bedtools: {e.stderr}", file=sys.stderr)
        raise


def main():
    parser = argparse.ArgumentParser(
        description='Define putative HGT regions from kairos output.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-b', '--bed', required=True,
                        help='Merged windows BED file from bedtools')
    parser.add_argument('-o', '--overlaps', required=True,
                        help='Deduplicated overlaps TSV file from kairos')
    parser.add_argument('-f', '--fai', required=True,
                        help='FASTA index file (.fai)')
    parser.add_argument('-p', '--padding', type=int, default=1000,
                        help='Edge padding size for region classification')
    parser.add_argument('--region-metadata-output', default='region_metadata.tsv',
                        help='Output file for region metadata')
    parser.add_argument('--region-overlaps-output', default='region_overlaps.tsv',
                        help='Output file for region overlaps')
    parser.add_argument('--filtered-overlaps-output', default='region_overlaps_filtered.tsv',
                        help='Output file for filtered region overlaps')
    parser.add_argument('--filter-reciprocals', action='store_true',
                        help='Remove reciprocal overlaps, keeping only the first occurrence of each pair')
    parser.add_argument('--extract-regions', action='store_true',
                        help='Create BED file and extract sequences. Requires --input-fasta.')
    parser.add_argument('--input-fasta', default=None,
                        help='Input FASTA file for sequence extraction')
    parser.add_argument('--extracted-fasta', default='hgt_regions_extracted.fasta',
                        help='Output FASTA file with extracted regions')
    parser.add_argument('--slop', type=int, default=0,
                        help='Basepairs to add either side of each region (bedtools slop)')
    parser.add_argument('--bed-output', default='hgt_regions.bed',
                        help='Output BED file for bedtools extraction')

    args = parser.parse_args()

    # ── Load ────────────────────────────────────────────────────────────────
    print("Loading data files...", file=sys.stderr)
    merged_bed, overlaps, fai, edge_padding = load_data(
        args.bed, args.overlaps, args.fai, args.padding
    )
    print(f"  {len(merged_bed)} regions from BED", file=sys.stderr)
    print(f"  {len(overlaps)} overlaps", file=sys.stderr)
    print(f"  {len(fai)} contigs from FAI", file=sys.stderr)

    # ── Region metadata ─────────────────────────────────────────────────────
    print("\nProcessing regions...", file=sys.stderr)
    region_metadata = get_all_regions(merged_bed, fai, edge_padding)
    print(f"  {len(region_metadata)} regions identified", file=sys.stderr)

    # ── Region overlaps ─────────────────────────────────────────────────────
    print("\nIdentifying region overlaps...", file=sys.stderr)
    region_overlaps = get_region_overlaps(region_metadata, overlaps)
    print(f"  {len(region_overlaps)} region overlaps found", file=sys.stderr)

    # ── Save region metadata and overlaps ───────────────────────────────────
    region_metadata.to_csv(args.region_metadata_output, sep='\t', index=False)
    print(f"\nSaved region metadata → {args.region_metadata_output}", file=sys.stderr)

    region_overlaps.to_csv(args.region_overlaps_output, sep='\t', index=False)
    print(f"Saved region overlaps → {args.region_overlaps_output}", file=sys.stderr)

    # ── Filter reciprocals ──────────────────────────────────────────────────
    print("\nFiltering reciprocal overlaps...", file=sys.stderr)
    region_overlaps_filtered = filter_reciprocals(region_overlaps)
    removed = len(region_overlaps) - len(region_overlaps_filtered)
    print(f"  Removed {removed} reciprocals, retained {len(region_overlaps_filtered)} pairs",
          file=sys.stderr)
    region_overlaps_filtered.to_csv(args.filtered_overlaps_output, sep='\t', index=False)
    print(f"Saved filtered overlaps → {args.filtered_overlaps_output}", file=sys.stderr)

    # ── Sequence extraction ─────────────────────────────────────────────────
    if args.extract_regions and args.input_fasta:
        create_bed_and_extract(
            region_metadata,
            args.input_fasta,
            args.bed_output,
            args.extracted_fasta,
            slop=args.slop,
            fai_file=args.fai
        )

    print("\nDone!", file=sys.stderr)


if __name__ == '__main__':
    main()