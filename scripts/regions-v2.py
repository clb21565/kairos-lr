#!/usr/bin/env python3
"""
from kairos, extract hgt regions
"""

import argparse
import pandas as pd
import sys

# load the expanded bed files (so that orfs within the same region are grouped together), and the overlaps file (so that we can identify which regions share orfs)
def load_data(merged_bed_file, overlaps_file, fai_file, edge_padding=1000):
    """
    Load all input data files.
    
    Parameters:
    -----------  
    merged_bed_file : str
        Path to merged windows BED file from bedtools
    overlaps_file : str
        Path to deduplicated overlaps TSV from kairos
    fai_file : str
        Path to FASTA index file
    edge_padding : int
        Padding size for checking end region in contig (default: 1000)
    
    Returns:
    --------
    tuple : (merged_bed, overlaps, fai, edge_padding)
    """
    # Load merged BED file
    merged_bed = pd.read_csv(
        merged_bed_file,
        sep='\t',
        header=None,
        names=['contig', 'start', 'stop', 'orf_cts', 'orfs_encompassed']
    )
    
    # Load overlaps file
    overlaps = pd.read_csv(overlaps_file, sep='\t')
    overlaps['orf1'] = overlaps['contig1'] + '_' + overlaps['identical_orf_position_contig1'].astype(str)
    overlaps['orf2'] = overlaps['contig2'] + '_' + overlaps['identical_orf_position_contig2'].astype(str)
    
    # Load FASTA index
    fai = pd.read_csv(
        fai_file,
        sep='\t',
        header=None,
        names=['contig', 'contig_length', 'offset', 'linebases', 'qualoffset']
    )[['contig', 'contig_length']]
    
    return merged_bed, overlaps, fai, edge_padding

# define regions and categorize them based on the distance to the contig ends (i.e., are they internal or edge regions)
def get_regions_for_contig(contig, merged_bed, fai, edge_padding):
    """
    Get region metadata for a specific contig.
    
    Parameters:
    -----------
    contig : str
        Contig name
    merged_bed : pd.DataFrame
        Merged BED data
    fai : pd.DataFrame
        FASTA index data
    edge_padding : int
        Padding size for edge determination
    
    Returns:
    --------
    pd.DataFrame : Region metadata for the contig
    """
    # Filter for specific contig
    idf = merged_bed[merged_bed['contig'] == contig].copy()
    
    # Assign unique region IDs
    idf['region_idx'] = [f"{contig}_region_{i+1}" for i in range(len(idf))]
    
    # Calculate region length
    idf['region_len'] = idf['stop'] - idf['start']
    
    # Merge with contig metadata
    idf = idf.merge(fai, on='contig')
    
    # Calculate edge cutoffs
    idf['end_edge_cutoff'] = idf['contig_length'] - edge_padding
    idf['start_edge_cutoff'] = 1 + edge_padding
    
    # Determine region type
    def categorize_region(row):
        if row['stop'] > row['end_edge_cutoff']:
            return '_end_edge_'
        elif row['start'] < row['start_edge_cutoff']:
            return '_start_edge_'
        else:
            return '_internal_'
    
    idf['region_type'] = idf.apply(categorize_region, axis=1)
    
    return idf

# process all contigs to get complete region metadata, including region type and length, and create a region title for each region (e.g., "contig1_region_1_internal_length=5000")
def get_all_regions(merged_bed, fai, edge_padding):
    """
    Get region metadata for all contigs.
    
    Parameters:
    -----------
    merged_bed : pd.DataFrame
        Merged BED data
    fai : pd.DataFrame
        FASTA index data
    edge_padding : int
        Padding size for edge determination
    
    Returns:
    --------
    pd.DataFrame : Complete region metadata
    """
    contigs = merged_bed['contig'].unique()
    
    region_metadata_list = []
    for contig in contigs:
        region_data = get_regions_for_contig(contig, merged_bed, fai, edge_padding)
        region_metadata_list.append(region_data)
    
    region_metadata = pd.concat(region_metadata_list, ignore_index=True)
    
    # Create region title
    region_metadata['region_title'] = (
        region_metadata['region_idx'] + 
        region_metadata['region_type'] + 
        'length=' + 
        region_metadata['region_len'].astype(str)
    )
    
    return region_metadata

# identify overlaps between regions based on ORF matches, and create a dataframe with region pairs that share ORFs (i.e., if region1 and region2 both contain orfs that match in the overlaps file, then region1 and region2 are considered to overlap)
def get_region_overlaps(region_metadata, overlaps):
    """
    Identify overlaps between regions based on ORF matches.
    
    Parameters:
    -----------
    region_metadata : pd.DataFrame
        Region metadata with ORF information
    overlaps : pd.DataFrame
        ORF overlap data
    
    Returns:
    --------
    pd.DataFrame : Region overlap information
    """
    
    # Prepare region-ORF mapping
    region_orfs = region_metadata[['contig', 'region_idx', 'orfs_encompassed']].copy()
    
    # Split comma-separated ORFs into separate rows
    region_orfs = region_orfs.assign(
        orfs_encompassed=region_orfs['orfs_encompassed'].str.split(',')
    ).explode('orfs_encompassed')
    
    # Select relevant columns from overlaps
    matched_orfs = overlaps[['orf1', 'orf2']].copy()
    
    # Merge to get region associations
    region_overlaps = matched_orfs.merge(
        region_orfs,
        left_on='orf1',
        right_on='orfs_encompassed',
        how='inner'
    ).rename(columns={'orf1': 'orf_region1'})
    
    region_overlaps = region_overlaps.merge(
        region_orfs,
        left_on='orf2',
        right_on='orfs_encompassed',
        how='inner',
        suffixes=('_contig1', '_contig2')
    ).rename(columns={'orf2': 'orf_region2'})
    
    # Select and rename final columns
    region_overlaps = region_overlaps[[
        'contig_contig1', 'contig_contig2',
        'region_idx_contig1', 'region_idx_contig2',
        'orf_region1', 'orf_region2'
    ]].rename(columns={
        'contig_contig1': 'contig1',
        'contig_contig2': 'contig2'
    })
    
    return region_overlaps

# remove reciprocal overlaps where region1+region2 also appears as region2+region1. Keeps only the first occurrence of each unique region pair.
def filter_reciprocals(hgt_df):
    """
    Remove reciprocal overlaps where region1+region2 also appears as region2+region1.
    Keeps only the first occurrence of each unique region pair.
    
    Parameters:
    -----------
    hgt_df : pd.DataFrame
        DataFrame with region_idx_contig1 and region_idx_contig2 columns
    
    Returns:
    --------
    pd.DataFrame : Filtered dataframe without reciprocal duplicates
    """
    # Create a key for each pair
    hgt_df = hgt_df.copy()
    hgt_df['region_key'] = hgt_df['region_idx_contig1'] + '___' + hgt_df['region_idx_contig2']
    
    # Get unique pairs
    pairs = hgt_df[['region_idx_contig1', 'region_idx_contig2']].drop_duplicates()
    
    # Sort each pair to identify reciprocals
    # Create a sorted pair key where the alphabetically first region comes first
    pairs['sorted_pair'] = pairs.apply(
        lambda row: tuple(sorted([row['region_idx_contig1'], row['region_idx_contig2']])),
        axis=1
    )
    
    # Keep only first occurrence of each sorted pair
    filtered_pairs = pairs[~pairs.duplicated(subset='sorted_pair', keep='first')].copy()
    filtered_pairs['region_key'] = (
        filtered_pairs['region_idx_contig1'] + '___' + filtered_pairs['region_idx_contig2']
    )
    
    # Filter original dataframe
    filtered_df = hgt_df[hgt_df['region_key'].isin(filtered_pairs['region_key'])].copy()
    
    # Remove the temporary key column
    filtered_df = filtered_df.drop(columns=['region_key'])
    
    return filtered_df

# find groups of connected regions (regions that share overlaps) using union-find algorithm to identify connected components. This will allow us to group regions that are connected through shared ORFs, even if they don't directly overlap with each other.
def find_connected_components(region_overlaps):
    """
    Find groups of connected regions (regions that share overlaps).
    Uses union-find algorithm to identify connected components.
    
    Parameters:
    -----------
    region_overlaps : pd.DataFrame
        DataFrame with region overlap information
    
    Returns:
    --------
    dict : Mapping of region_idx to group_id
    """
    # Build a graph of region connections
    from collections import defaultdict
    
    # Union-Find data structure
    parent = {}
    
    def find(x):
        if x not in parent:
            parent[x] = x
        if parent[x] != x:
            parent[x] = find(parent[x])  # Path compression
        return parent[x]
    
    def union(x, y):
        root_x = find(x)
        root_y = find(y)
        if root_x != root_y:
            parent[root_x] = root_y
    
    # Union all connected regions
    for _, row in region_overlaps.iterrows():
        union(row['region_idx_contig1'], row['region_idx_contig2'])
    
    # Group regions by their root
    groups = defaultdict(list)
    for region in parent.keys():
        root = find(region)
        groups[root].append(region)
    
    # Create region to group_id mapping
    region_to_group = {}
    for group_id, (root, regions) in enumerate(groups.items(), 1):
        for region in regions:
            region_to_group[region] = group_id
    
    return region_to_group, groups

def extract_region_groups(region_metadata, region_overlaps, output_file):
    """
    Identify groups of overlapping regions and write group metadata.
    Each group is named after its longest region.
    
    Parameters:
    -----------
    region_metadata : pd.DataFrame
        Region metadata with coordinates
    region_overlaps : pd.DataFrame
        Region overlap information
    output_file : str
        Path for output metadata file
    """
    import os
    
       
    print("\nIdentifying groups of overlapping regions...", file=sys.stderr)
    region_to_group, groups = find_connected_components(region_overlaps)
    print(f"Found {len(groups)} groups of overlapping regions", file=sys.stderr)
    
    region_meta_dict = region_metadata.set_index('region_idx').to_dict('index')
    
    group_info_list = []
    
    for group_id, (root, regions) in enumerate(groups.items(), 1):
        group_regions_meta = []
        for region in regions:
            if region in region_meta_dict:
                meta = region_meta_dict[region]
                meta['region_idx'] = region
                group_regions_meta.append(meta)
        
        if not group_regions_meta:
            continue
        
        longest_region = max(group_regions_meta, key=lambda x: x['region_len'])
        group_name = f"group_{group_id:04d}_{longest_region['region_idx']}_len{longest_region['region_len']}"
        
        group_info_list.append({
            'group_id': group_id,
            'group_name': group_name,
            'num_regions': len(group_regions_meta),
            'longest_region': longest_region['region_idx'],
            'longest_region_length': longest_region['region_len'],
            'regions': ','.join(sorted(regions))
        })
    
    group_info_df = pd.DataFrame(group_info_list)
    group_info_df.to_csv(output_file, sep='\t', index=False)
    
    print(f"\nIdentified {len(group_info_list)} groups", file=sys.stderr)
    print(f"Group metadata saved to {output_file}", file=sys.stderr)
    
    return group_info_df


def create_bed_and_extract(region_metadata, fasta_file, bed_output, fasta_output, slop=0, fai_file=None):
    """
    Create a BED file from region metadata and extract sequences using bedtools.
    
    Parameters:
    -----------
    region_metadata : pd.DataFrame
        Region metadata with coordinates and region_title
    fasta_file : str
        Path to input FASTA file
    bed_output : str
        Path for output BED file
    fasta_output : str
        Path for output extracted FASTA file
    slop : int
        Number of basepairs to add to left and right of each region (default: 0)
    fai_file : str
        Path to FASTA index file (required if slop > 0)
    """
    import subprocess
    import os
    
    # Create BED file with region_title as the name field
    bed_data = region_metadata[['contig', 'start', 'stop', 'region_title']].copy()
    
    # BED format: chrom, start, end, name
    bed_data.to_csv(bed_output, sep='\t', header=False, index=False)
    
    print(f"\nCreated BED file: {bed_output}", file=sys.stderr)
    print(f"  Total regions: {len(bed_data)}", file=sys.stderr)
    
    # Check if bedtools is available
    try:
        subprocess.run(['bedtools', '--version'], 
                      capture_output=True, 
                      check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("WARNING: bedtools not found. Please install bedtools to extract sequences.", 
              file=sys.stderr)
        print(f"You can manually extract sequences with:", file=sys.stderr)
        if slop > 0:
            print(f"  bedtools slop -i {bed_output} -g {fai_file} -b {slop} | bedtools getfasta -fi {fasta_file} -bed - -name -fo {fasta_output}", 
                  file=sys.stderr)
        else:
            print(f"  bedtools getfasta -fi {fasta_file} -bed {bed_output} -name -fo {fasta_output}", 
                  file=sys.stderr)
        return
    
    # Check if FASTA file exists
    if not os.path.exists(fasta_file):
        print(f"WARNING: FASTA file not found: {fasta_file}", file=sys.stderr)
        print(f"BED file created, but skipping sequence extraction.", file=sys.stderr)
        return
    
    # If slop > 0, check for FAI file
    if slop > 0:
        if not fai_file or not os.path.exists(fai_file):
            print(f"WARNING: FAI file required for slop but not found: {fai_file}", file=sys.stderr)
            print(f"BED file created, but skipping sequence extraction.", file=sys.stderr)
            print(f"You can create an FAI file with: samtools faidx {fasta_file}", file=sys.stderr)
            return
        print(f"\nAdding {slop} bp padding to each side of regions (bedtools slop)...", file=sys.stderr)
    
    # Run bedtools getfasta (with or without slop)
    print(f"\nExtracting sequences with bedtools...", file=sys.stderr)
    try:
        if slop > 0:
            # First run bedtools slop, then pipe to getfasta
            slop_cmd = [
                'bedtools', 'slop',
                '-i', bed_output,
                '-g', fai_file,
                '-b', str(slop)
            ]
            
            getfasta_cmd = [
                'bedtools', 'getfasta',
                '-fi', fasta_file,
                '-bed', '-',  # Read from stdin
                '-name',
                '-fo', fasta_output
            ]
            
            # Run slop and pipe to getfasta
            slop_proc = subprocess.Popen(slop_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            getfasta_proc = subprocess.Popen(getfasta_cmd, stdin=slop_proc.stdout, 
                                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            slop_proc.stdout.close()  # Allow slop_proc to receive SIGPIPE
            
            stdout, stderr = getfasta_proc.communicate()
            
            if getfasta_proc.returncode != 0:
                raise subprocess.CalledProcessError(getfasta_proc.returncode, getfasta_cmd, stderr=stderr)
            
            print(f"Successfully extracted sequences with {slop} bp padding to: {fasta_output}", file=sys.stderr)
        else:
            # No slop, just run getfasta
            cmd = [
                'bedtools', 'getfasta',
                '-fi', fasta_file,
                '-bed', bed_output,
                '-name',  # Use the name field (region_title) in output
                '-fo', fasta_output
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            print(f"Successfully extracted sequences to: {fasta_output}", file=sys.stderr)
        
        # Count sequences in output
        with open(fasta_output, 'r') as f:
            seq_count = sum(1 for line in f if line.startswith('>'))
        print(f"  Extracted {seq_count} sequences", file=sys.stderr)
        
    except subprocess.CalledProcessError as e:
        print(f"ERROR running bedtools: {e.stderr}", file=sys.stderr)
        raise


def main():
    parser = argparse.ArgumentParser(
        description='define putative HGT regions and groups',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '-b', '--bed',
        required=True,
        help='Merged windows BED file from bedtools'
    )
    parser.add_argument(
        '-o', '--overlaps',
        required=True,
        help='Deduplicated overlaps TSV file from kairos'
    )
    parser.add_argument(
        '-f', '--fai',
        required=True,
        help='FASTA index file (.fai)'
    )
    parser.add_argument(
        '-p', '--padding',
        type=int,
        default=1000,
        help='Edge padding size for region classification'
    )
    parser.add_argument(
        '--region-metadata-output',
        default='region_metadata.tsv',
        help='Output file for region metadata'
    )
    parser.add_argument(
        '--region-overlaps-output',
        default='region_overlaps.tsv',
        help='Output file for region overlaps'
    )
    parser.add_argument(
        '--group-metadata-output',
        default='group_metadata.tsv',
        help='Output file for group metadata'
    )
    parser.add_argument(
        '--filtered-overlaps-output',
        default='region_overlaps_filtered.tsv',
        help='Output file for filtered region overlaps (only used if --filter-reciprocals is set)'
    )

    parser.add_argument(
        '--filter-reciprocals',
        action='store_true',
        help='Whether to filter out reciprocal overlaps (region1-region2 and region2-region1). If set, only the first occurrence of each unique region pair will be retained for grouping.'
    )

    parser.add_argument(
        '--extract-regions',
        action='store_true',
        help='Whether to create BED file and extract sequences for identified regions. Requires --input-fasta to be set.'
    )

    parser.add_argument(
        '--input-fasta',
        default=None, 
        help='Input FASTA file for sequence extraction (required if --extract-regions is set)'
    )

    parser.add_argument(
        '--extracted-fasta',
        default='hgt_regions_extracted.fasta',
        help='Output FASTA file with extracted regions'
    )
    
    parser.add_argument(
        '--slop',
        type=int,
        default=0,
        help='Number of basepairs to add to left and right of each region (bedtools slop)'
    )
    parser.add_argument(
        '--bed-output',
        default='hgt_regions.bed',
        help='Output BED file for bedtools extraction'
    )
    
    args = parser.parse_args()
    
    # Load data
    print("Loading data files...", file=sys.stderr)
    merged_bed, overlaps, fai, edge_padding = load_data(
        args.bed,
        args.overlaps,
        args.fai,
        args.padding
    )
    
    print(f"Loaded {len(merged_bed)} regions from BED file", file=sys.stderr)
    print(f"Loaded {len(overlaps)} overlaps", file=sys.stderr)
    print(f"Loaded {len(fai)} contigs from FAI", file=sys.stderr)
    
    # Get region metadata
    print("\nProcessing regions...", file=sys.stderr)
    region_metadata = get_all_regions(merged_bed, fai, edge_padding)
    print(f"Identified {len(region_metadata)} total regions", file=sys.stderr)
    
    # Get region overlaps
    print("\nIdentifying region overlaps...", file=sys.stderr)
    region_overlaps = get_region_overlaps(region_metadata, overlaps)
    print(f"Found {len(region_overlaps)} region overlaps", file=sys.stderr)
    
    # Save output
    print(f"\nSaving region metadata to {args.region_metadata_output}", file=sys.stderr)
    region_metadata.to_csv(args.region_metadata_output, sep='\t', index=False)
    
    print(f"Saving region overlaps to {args.region_overlaps_output}", file=sys.stderr)
    region_overlaps.to_csv(args.region_overlaps_output, sep='\t', index=False)
    

    print("\nFiltering reciprocal overlaps...", file=sys.stderr)
    region_overlaps_filtered = filter_reciprocals(region_overlaps)
    removed_count = len(region_overlaps) - len(region_overlaps_filtered)
    print(f"Removed {removed_count} reciprocal overlaps", file=sys.stderr)
    print(f"Retained {len(region_overlaps_filtered)} unique region pairs", file=sys.stderr)
    print(f"Saving filtered region overlaps to {args.filtered_overlaps_output}", file=sys.stderr)
    region_overlaps_filtered.to_csv(args.filtered_overlaps_output, sep='\t', index=False)


    # Use filtered overlaps if available, otherwise use all overlaps
    overlaps_for_grouping = (region_overlaps_filtered
                            if args.filter_reciprocals 
                            else region_overlaps)
        
    group_summary = extract_region_groups(
        region_metadata,
        overlaps_for_grouping,
        args.group_metadata_output)


    # Extract sequences if FASTA file provided
    if args.input_fasta:
        create_bed_and_extract(
            region_metadata,
            args.input_fasta,
            args.bed_output,
            args.extracted_fasta,
            slop=args.slop,
            fai_file=args.fai
        )
    

    print(f"Saving group summary to {args.group_metadata_output}", file=sys.stderr)
    group_summary.to_csv(args.group_metadata_output, sep='\t', index=False)
    print("\nDone!", file=sys.stderr)


if __name__ == '__main__':
    main()
