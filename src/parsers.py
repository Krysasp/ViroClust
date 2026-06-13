#!/usr/bin/env python3
"""
Parsers module for ViroClust.

Handles FASTA and cluster file parsing, including extraction of sequence 
length information and reference sequence markers ('*').
"""

import re


def parse_fasta(fasta_path):
    """
    Parse FASTA file and return dictionary of sequences with metadata extraction.
    
    Args:
        fasta_path: path to FASTA file
    
    Returns:
        dict mapping headers to sequences
    """
    sequences = {}
    current_header = None
    current_seq = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_header:
                    sequences[current_header] = ''.join(current_seq)
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line.upper())
        
        if current_header:
            sequences[current_header] = ''.join(current_seq)
    
    return sequences


def extract_fasta_metadata(fasta_sequences):
    """
    Extract metadata from FASTA headers including taxID, species, and other info.
    
    Expected FASTA header format: ACCESSION|Species|taxID=1234|other_info...
    
    Args:
        fasta_sequences: dict mapping FASTA headers to sequences
    
    Returns:
        dict mapping headers to metadata dicts with keys:
        - accession: accession ID (before first '|')
        - taxID: taxonomic ID if present
        - species: species name if present
        - raw_header: original header string
        - additional_info: dict of other key=value pairs
    """
    metadata = {}
    
    for header, seq in fasta_sequences.items():
        meta = {
            'accession': header.split('|')[0],
            'taxID': None,
            'species': None,
            'raw_header': header,
            'additional_info': {}
        }
        
        # Parse header components
        parts = header.split('|')
        
        # Extract species (second part after accession)
        if len(parts) > 1:
            meta['species'] = parts[1].strip()
        
        # Extract taxID and other key=value pairs
        for part in parts[2:]:
            part = part.strip()
            if '=' in part:
                key, value = part.split('=', 1)
                if key.strip() == 'taxID':
                    meta['taxID'] = value.strip()
                else:
                    meta['additional_info'][key.strip()] = value.strip()
            else:
                # Store as additional info if not key=value format
                meta['additional_info'][f'field_{len(meta["additional_info"])}'] = part
        
        # If taxID not found in key=value format, try to find it as standalone number
        if meta['taxID'] is None:
            for part in parts:
                if part.isdigit():
                    meta['taxID'] = part
                    break
        
        metadata[header] = meta
    
    return metadata


def parse_fasta_with_metadata(fasta_path):
    """
    Parse FASTA file and return sequences with extracted metadata.
    
    Args:
        fasta_path: path to FASTA file
    
    Returns:
        tuple of (sequences, metadata)
        - sequences: dict mapping headers to sequences
        - metadata: dict mapping headers to metadata dicts
    """
    sequences = parse_fasta(fasta_path)
    metadata = extract_fasta_metadata(sequences)
    return sequences, metadata


def parse_cluster_file(cluster_path, max_clusters=None):
    """
    Parse cluster file and return dictionary mapping cluster numbers to list of accession IDs.
    Also extracts sequence length information and reference sequence markers ('*').
    
    Cluster file format (from CD-HIT):
    >Cluster 0
    0   36511aa, >MN486048|Human... at 97.11%
    1   36510aa, >MN486049|Human... at 97.10% *
    
    The '*' at the end of a line indicates the reference sequence.
    
    Args:
        cluster_path: path to cluster file
        max_clusters: if set, only read first N clusters (for testing)
    
    Returns:
        tuple of (clusters, length_info, reference_info)
        - clusters: dict mapping cluster_num to list of accession IDs
        - length_info: dict mapping cluster_num to list of (accession, length, is_reference) tuples
        - reference_info: dict mapping cluster_num to reference sequence info
    """
    clusters = {}
    length_info = {}
    reference_info = {}
    current_cluster = None
    cluster_count = 0
    
    with open(cluster_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('>Cluster '):
                cluster_num = line.split()[1]
                current_cluster = cluster_num
                clusters[cluster_num] = []
                length_info[cluster_num] = []
                cluster_count += 1
                if max_clusters and cluster_count >= max_clusters:
                    break
            else:
                # Parse cluster member line
                # Format: "0   36511aa, >MN486048|Human... at 97.11%" or "1   36510aa, >MN486049|Human... * at 97.10%"
                # The '*' marker appears before "at" or at the end
                
                # Extract accession ID
                match = re.search(r'>([A-Za-z0-9_.]+)\|', line)
                if match:
                    accession = match.group(1)
                    clusters[current_cluster].append(accession)
                
                # Extract sequence length (e.g., "36511aa")
                length_match = re.search(r'(\d+)aa', line)
                if length_match:
                    length = int(length_match.group(1))
                    
                    # Check for reference marker '*'
                    # The '*' can appear in two formats:
                    # 1. " >ACCESSION|... at 97.11% *" (at the end)
                    # 2. " >ACCESSION|... * at 97.11%" (before "at")
                    is_reference = False
                    if ' *' in line or line.endswith('*'):
                        is_reference = True
                    
                    length_info[current_cluster].append({
                        'accession': accession,
                        'length': length,
                        'is_reference': is_reference,
                        'index': len(length_info[current_cluster])
                    })
                    
                    # Track reference sequence info
                    if is_reference:
                        reference_info[current_cluster] = {
                            'reference_id': accession,
                            'reference_idx': len(length_info[current_cluster]) - 1,
                            'reference_length': length
                        }
    
    return clusters, length_info, reference_info


def extract_cluster_sequences(fasta_sequences, cluster_sequences):
    """
    Extract sequences from FASTA based on cluster file listings.
    Returns dict mapping cluster_num to list of (full_header, sequence) tuples.
    Only includes sequences that are found in the FASTA file (excludes missing sequences).
    
    Args:
        fasta_sequences: dict mapping FASTA headers to sequences
        cluster_sequences: dict mapping cluster_num to list of accession IDs
    
    Returns:
        dict mapping cluster_num to list of (header, sequence) tuples (only found sequences)
    """
    cluster_data = {}
    
    for cluster_num, accessions in cluster_sequences.items():
        cluster_seqs = []
        for accession in accessions:
            # Try exact match first
            if accession in fasta_sequences:
                cluster_seqs.append((accession, fasta_sequences[accession]))
            else:
                # Try matching by accession prefix (before first '|')
                matched = False
                for header, seq in fasta_sequences.items():
                    header_accession = header.split('|')[0]
                    if header_accession == accession:
                        cluster_seqs.append((header, seq))
                        matched = True
                        break
                
                if not matched:
                    # Sequence not found in FASTA - skip it (will be reported in missing sequence report)
                    pass
        
        cluster_data[cluster_num] = cluster_seqs
    
    return cluster_data


def extract_cluster_sequences_with_info(fasta_sequences, cluster_sequences, length_info):
    """
    Extract sequences from FASTA with additional length and reference information.
    Returns enhanced cluster data with sequence metadata.
    Only includes sequences that are found in the FASTA file (excludes missing sequences).
    
    Args:
        fasta_sequences: dict mapping FASTA headers to sequences
        cluster_sequences: dict mapping cluster_num to list of accession IDs
        length_info: dict mapping cluster_num to list of (accession, length, is_reference) tuples
    
    Returns:
        tuple of (cluster_data, reference_sequence_list, missing_sequences)
        - cluster_data: dict mapping cluster_num to list of (header, sequence, length, is_reference) tuples
        - reference_sequence_list: list of (cluster_num, reference_header, reference_length) tuples
        - missing_sequences: dict mapping cluster_num to list of missing accession IDs
    """
    cluster_data = {}
    reference_sequence_list = []
    missing_sequences = {}
    
    for cluster_num, accessions in cluster_sequences.items():
        cluster_seqs = []
        cluster_ref = None
        cluster_missing = []
        
        for i, accession in enumerate(accessions):
            # Get length and reference info from length_info
            len_info = length_info.get(cluster_num, [{}])[i] if i < len(length_info.get(cluster_num, [])) else {}
            seq_length = len_info.get('length', 0)
            is_reference = len_info.get('is_reference', False)
            
            # Try exact match first
            if accession in fasta_sequences:
                seq = fasta_sequences[accession]
            else:
                # Try matching by accession prefix
                matched = False
                for header, s in fasta_sequences.items():
                    header_accession = header.split('|')[0]
                    if header_accession == accession:
                        seq = s
                        matched = True
                        break
                
                if not matched:
                    # Sequence not found in FASTA - track it but skip it
                    cluster_missing.append(accession)
                    continue
            
            cluster_seqs.append((accession, seq, seq_length, is_reference))
            
            if is_reference:
                cluster_ref = (cluster_num, accession, seq_length)
        
        cluster_data[cluster_num] = cluster_seqs
        if cluster_ref:
            reference_sequence_list.append(cluster_ref)
        if cluster_missing:
            missing_sequences[cluster_num] = cluster_missing
    
    return cluster_data, reference_sequence_list, missing_sequences


def build_sequence_metadata_index(fasta_metadata, cluster_data):
    """
    Build a comprehensive metadata index linking cluster sequences to their metadata.
    
    Args:
        fasta_metadata: dict mapping FASTA headers to metadata dicts
        cluster_data: dict mapping cluster_num to list of (header, seq, length, is_ref) tuples
    
    Returns:
        dict mapping cluster_num to list of enriched metadata tuples:
        (header, seq, length, is_ref, metadata_dict, taxID, species, renamed_header)
    """
    metadata_index = {}
    
    for cluster_num, seq_list in cluster_data.items():
        enriched_seqs = []
        for item in seq_list:
            header, seq, length, is_ref = item
            
            # Get metadata from FASTA metadata
            meta = fasta_metadata.get(header, {
                'accession': header.split('|')[0],
                'taxID': None,
                'species': None,
                'raw_header': header,
                'additional_info': {}
            })
            
            taxID = meta.get('taxID')
            species = meta.get('species')
            
            # Generate renamed header format: CLUST{cluster_num}_TAX{taxID}_{accession}
            renamed_header = f"CLUST{cluster_num.zfill(4)}_TAX{taxID or 'unknown'}_{meta['accession']}"
            
            enriched_seqs.append((header, seq, length, is_ref, meta, taxID, species, renamed_header))
        
        metadata_index[cluster_num] = enriched_seqs
    
    return metadata_index
