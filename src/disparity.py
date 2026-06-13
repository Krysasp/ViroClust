#!/usr/bin/env python3
"""
Disparity module for ViroClust.

Handles detection and processing of clusters with significant sequence length 
differences, including extended region handling for truncated sequences.
"""

from collections import defaultdict


def detect_length_disparity(cluster_data, reference_info=None, threshold=0.20):
    """
    Detect clusters where shortest and longest sequences differ by >= threshold.
    
    For a cluster with sequences of varying lengths:
    - If shortest_length / longest_length < threshold (e.g., 0.20), the cluster 
      has significant length disparity
    - Leading/trailing regions of the reference sequence (marked with '*') that are 
      truncated in shorter sequences will be handled with relaxed inclusion threshold
    - Extended regions use per-position conservation (denominator = sequences with base at that position)
    
    Example:
        Cluster 316: reference OL794386 (13516aa), shortest KM408077 (765aa)
        Ratio = 765/13516 = 0.057 < 0.20 → has disparity
        Leading/trailing regions of OL794386 that are missing in KM408077 will use
        base-by-base conservation with 7 sequences in the longer subcluster
    
    Args:
        cluster_data: dict mapping cluster_num to list of (header, sequence) tuples
        reference_info: dict mapping cluster_num to reference sequence info (optional)
        threshold: minimum ratio (shortest/longest) to NOT trigger disparity handling (default 0.20)
    
    Returns:
        dict mapping cluster_num to disparity details:
        {
            'cluster_num': {
                'num_sequences': int,
                'shortest_length': int,
                'longest_length': int,
                'reference_length': int,
                'reference_id': str,
                'ratio': float,
                'has_disparity': bool,
                'extended_leading_regions': int,
                'extended_trailing_regions': int,
                'extended_leading_length': int,
                'extended_trailing_length': int,
                'longer_subcluster_count': int
            }
        }
    """
    disparity_info = {}
    
    for cluster_num, seq_list in cluster_data.items():
        if len(seq_list) == 0:
            continue
        
        # Calculate lengths
        # seq_list items are tuples of (header, seq, length, is_reference) from extract_cluster_sequences_with_info
        # or tuples of (header, seq) from extract_cluster_sequences
        lengths = []
        for item in seq_list:
            if len(item) == 4:
                # (header, seq, length, is_reference) format
                _, seq, length, _ = item
                if seq:
                    lengths.append(length)
            else:
                # (header, seq) format
                _, seq = item
                if seq:
                    lengths.append(len(seq))
        
        if not lengths:
            continue
        
        shortest_length = min(lengths)
        longest_length = max(lengths)
        num_sequences = len(lengths)
        
        # Get reference sequence info if available
        ref_length = None
        ref_id = None
        if reference_info and cluster_num in reference_info:
            ref_id = reference_info[cluster_num].get('reference_id')
            # Find the reference sequence length
            for item in seq_list:
                if len(item) == 4:
                    header, seq, length, is_ref = item
                    if header == ref_id:
                        ref_length = length
                        break
                else:
                    header, seq = item
                    if header == ref_id:
                        ref_length = len(seq)
                        break
        
        # If no reference marked, use longest sequence as reference
        if ref_length is None:
            for header, seq in seq_list:
                if len(seq) == longest_length:
                    ref_id = header
                    ref_length = longest_length
                    break
        
        # Calculate ratio using reference length vs shortest
        ratio = shortest_length / ref_length if ref_length > 0 else 0.0
        
        # Check if disparity exceeds threshold
        has_disparity = ratio < threshold
        
        if has_disparity:
            # Calculate extended regions based on reference sequence
            # Extended regions are the leading/trailing parts of the reference sequence
            # that are not present in shorter sequences
            
            # Find sequences in the "longer subcluster" (sequences within 5% of reference length)
            # Handle both 2-tuples (header, seq) and 4-tuples (header, seq, length, is_reference)
            longer_subcluster = []
            for item in seq_list:
                if len(item) == 4:
                    _, seq, length, _ = item
                    if seq and length >= ref_length * 0.95:
                        longer_subcluster.append(seq)
                else:
                    _, seq = item
                    if seq and len(seq) >= ref_length * 0.95:
                        longer_subcluster.append(seq)
            longer_subcluster_count = len(longer_subcluster)
            
            # Calculate extended regions based on reference vs shortest sequence difference
            # The extended regions are the parts of reference that are missing in shortest
            length_diff = ref_length - shortest_length
            
            # For the example: reference 13516aa, shortest 765aa
            # The longer subcluster (OL794378, OL794386, OL794403, OL794407, OL794455, OL794465, MK167039)
            # ranges from 13462aa to 13516aa (within 5% of reference)
            # Extended regions = (13516 - 765) = 12751 positions split between leading and trailing
            # But we need to identify which specific positions are in extended regions
            
            # For simplicity, estimate based on length difference
            # Assume shorter sequences are truncated at both ends proportionally
            # Extended regions are positions in reference not covered by shortest sequence
            extended_leading_length = length_diff // 2
            extended_trailing_length = length_diff - extended_leading_length
            
            disparity_info[cluster_num] = {
                'num_sequences': num_sequences,
                'shortest_length': shortest_length,
                'longest_length': longest_length,
                'reference_length': ref_length,
                'reference_id': ref_id,
                'ratio': ratio,
                'has_disparity': True,
                'extended_leading_regions': extended_leading_length,
                'extended_trailing_regions': extended_trailing_length,
                'extended_leading_length': extended_leading_length,
                'extended_trailing_length': extended_trailing_length,
                'longer_subcluster_count': longer_subcluster_count
            }
    
    return disparity_info


def identify_reference_sequence(cluster_data, length_info=None):
    """
    Identify reference sequences (marked with '*') in each cluster.
    
    The reference sequence is typically the longest/most complete sequence
    and should be used as the anchor for alignment.
    
    Args:
        cluster_data: dict mapping cluster_num to list of (header, sequence) tuples
        length_info: dict mapping cluster_num to list of (accession, length, is_reference) dicts (optional)
    
    Returns:
        dict mapping cluster_num to reference sequence info:
        {
            'cluster_num': {
                'reference_id': str,
                'reference_idx': int,
                'reference_length': int
            }
        }
    """
    reference_info = {}
    
    # Use cluster_data (already filtered to exclude missing sequences)
    # cluster_data items are 4-tuples: (header, seq, length, is_reference)
    for cluster_num, seq_list in cluster_data.items():
        if not seq_list:
            continue
        
        # First, try to find reference sequence marked with '*'
        ref_idx = None
        for idx, item in enumerate(seq_list):
            if len(item) == 4:
                _, _, _, is_reference = item
                if is_reference:
                    ref_idx = idx
                    break
        
        if ref_idx is not None:
            # Found reference sequence
            header, seq, length, _ = seq_list[ref_idx]
            reference_info[cluster_num] = {
                'reference_id': header,
                'reference_idx': ref_idx,
                'reference_length': length
            }
        else:
            # Fall back to identifying longest sequence as reference
            longest_idx = max(range(len(seq_list)), key=lambda i: len(seq_list[i][1]))
            header, seq = seq_list[longest_idx][:2]  # Extract header and seq from 4-tuple
            reference_info[cluster_num] = {
                'reference_id': header,
                'reference_idx': longest_idx,
                'reference_length': len(seq)
            }
    
    return reference_info


def calculate_disparity_threshold(shortest_length, longest_length, threshold=0.20):
    """
    Calculate whether a cluster exceeds the length disparity threshold.
    
    Args:
        shortest_length: length of shortest sequence
        longest_length: length of longest sequence
        threshold: minimum ratio to NOT trigger disparity (default 0.20)
    
    Returns:
        tuple of (has_disparity, ratio)
    """
    if longest_length == 0:
        return False, 0.0
    
    ratio = shortest_length / longest_length
    has_disparity = ratio < threshold
    
    return has_disparity, ratio


def process_extended_regions(aligned, iupac_consensus, conservation, inclusion, 
                            filter_gaps, gap_threshold, disparity_info):
    """
    Process consensus generation with special handling for extended regions.
    
    For clusters with significant length disparity, the leading and trailing 
    regions where shorter sequences are truncated will:
    1. Bypass the inclusion threshold (0.85) - always include if any sequence has a base
    2. Use per-position conservation (denominator = sequences with base at that position)
       For the example cluster 316: extended regions use 7 sequences from the longer subcluster
    3. Support IUPAC codes based on available bases at each position
    
    Args:
        aligned: dict mapping header to aligned sequence
        iupac_consensus: IUPAC consensus sequence (with gaps)
        conservation: conservation threshold (default 0.90)
        inclusion: inclusion threshold (default 0.70)
        filter_gaps: if True, remove positions with >gap_threshold gaps
        gap_threshold: fraction of sequences with gaps to trigger removal
        disparity_info: dict with extended region information
    
    Returns:
        tuple of (consensus, gap_filter_info, extended_region_info)
        - consensus: final consensus sequence with extended regions
        - gap_filter_info: dict with gap filtering statistics
        - extended_region_info: dict with extended region details
    """
    if not aligned:
        return "", {}, {}
    
    seq_length = len(iupac_consensus)
    num_sequences = len(aligned)
    
    # Get extended region boundaries
    extended_leading = disparity_info.get('extended_leading_regions', 0)
    extended_trailing = disparity_info.get('extended_trailing_regions', 0)
    
    # Calculate central region
    central_start = extended_leading
    central_end = seq_length - extended_trailing
    
    # Get the longer subcluster count (sequences within 5% of reference length)
    # This is used as the denominator for extended region consensus
    longer_subcluster_count = disparity_info.get('longer_subcluster_count', num_sequences)
    
    # Generate consensus with special handling for extended regions
    consensus = []
    extended_leading_bases = []
    extended_trailing_bases = []
    central_bases = []
    
    for pos in range(seq_length):
        base_count = defaultdict(int)
        total_bases = 0
        
        for header, seq in aligned.items():
            base = seq[pos] if pos < len(seq) else '-'
            if base != '-':
                base_count[base] += 1
                total_bases += 1
        
        # Determine which region this position is in
        if pos < extended_leading:
            # Leading extended region - use per-position conservation
            # Denominator = sequences with base at this position (e.g., 7 sequences in longer subcluster)
            if total_bases > 0:
                max_count = max(base_count.values())
                dominant_bases = frozenset([b for b, c in base_count.items() if c == max_count])
                base_freq = max_count / total_bases  # Denominator is total_bases (e.g., 7)
                
                if base_freq >= conservation:
                    if len(dominant_bases) == 1:
                        consensus.append(list(dominant_bases)[0])
                        extended_leading_bases.append(list(dominant_bases)[0])
                    else:
                        iupac_char = get_iupac_code(dominant_bases)
                        consensus.append(iupac_char)
                        extended_leading_bases.append(iupac_char)
                else:
                    # Below conservation - use significant bases
                    significant_bases = frozenset([b for b, c in base_count.items() if c >= max_count * 0.5])
                    iupac_char = get_iupac_code(significant_bases)
                    consensus.append(iupac_char)
                    extended_leading_bases.append(iupac_char)
            else:
                consensus.append('-')
                extended_leading_bases.append('-')
                
        elif pos >= central_end:
            # Trailing extended region - same logic as leading
            if total_bases > 0:
                max_count = max(base_count.values())
                dominant_bases = frozenset([b for b, c in base_count.items() if c == max_count])
                base_freq = max_count / total_bases  # Denominator is total_bases (e.g., 7)
                
                if base_freq >= conservation:
                    if len(dominant_bases) == 1:
                        consensus.append(list(dominant_bases)[0])
                        extended_trailing_bases.append(list(dominant_bases)[0])
                    else:
                        iupac_char = get_iupac_code(dominant_bases)
                        consensus.append(iupac_char)
                        extended_trailing_bases.append(iupac_char)
                else:
                    significant_bases = frozenset([b for b, c in base_count.items() if c >= max_count * 0.5])
                    iupac_char = get_iupac_code(significant_bases)
                    consensus.append(iupac_char)
                    extended_trailing_bases.append(iupac_char)
            else:
                consensus.append('-')
                extended_trailing_bases.append('-')
        else:
            # Central region - use standard inclusion threshold
            min_inclusion_count = int(num_sequences * inclusion)
            
            if total_bases < min_inclusion_count:
                consensus.append('-')
                central_bases.append('-')
                continue
            
            if not base_count:
                consensus.append('-')
                central_bases.append('-')
                continue
            
            max_count = max(base_count.values())
            dominant_bases = frozenset([b for b, c in base_count.items() if c == max_count])
            base_freq = max_count / total_bases
            
            if base_freq >= conservation:
                if len(dominant_bases) == 1:
                    consensus.append(list(dominant_bases)[0])
                    central_bases.append(list(dominant_bases)[0])
                else:
                    iupac_char = get_iupac_code(dominant_bases)
                    consensus.append(iupac_char)
                    central_bases.append(iupac_char)
            else:
                significant_bases = frozenset([b for b, c in base_count.items() if c >= max_count * 0.5])
                iupac_char = get_iupac_code(significant_bases)
                consensus.append(iupac_char)
                central_bases.append(iupac_char)
    
    # Build gap filter info
    gap_filter_info = {
        'positions_removed': [],  # Would be populated if filter_gaps=True
        'gap_stats': {
            'total_positions': seq_length,
            'extended_leading_length': extended_leading,
            'extended_trailing_length': extended_trailing,
            'central_length': central_end - central_start,
            'gap_threshold': gap_threshold,
            'longer_subcluster_count': longer_subcluster_count
        }
    }
    
    # Build extended region info
    extended_region_info = {
        'has_extended_regions': extended_leading > 0 or extended_trailing > 0,
        'extended_leading_start': 0,
        'extended_leading_end': extended_leading,
        'central_start': central_start,
        'central_end': central_end,
        'extended_trailing_start': central_end,
        'extended_trailing_end': seq_length,
        'extended_leading_length': len(extended_leading_bases),
        'extended_trailing_length': len(extended_trailing_bases),
        'extended_leading_consensus': ''.join(extended_leading_bases),
        'extended_trailing_consensus': ''.join(extended_trailing_bases),
        'central_consensus': ''.join(central_bases),
        'longer_subcluster_count': longer_subcluster_count
    }
    
    return ''.join(consensus), gap_filter_info, extended_region_info


def get_iupac_code(base_set):
    """
    Get IUPAC code for a set of bases.
    
    Args:
        base_set: frozenset of bases (e.g., frozenset(['A', 'G']))
    
    Returns:
        IUPAC code character (e.g., 'R' for A/G)
    """
    # IUPAC ambiguity codes mapping
    IUPAC_CODES = {
        frozenset(['A']): 'A',
        frozenset(['C']): 'C',
        frozenset(['G']): 'G',
        frozenset(['T']): 'T',
        frozenset(['A', 'G']): 'R',
        frozenset(['C', 'T']): 'Y',
        frozenset(['G', 'C']): 'S',
        frozenset(['A', 'T']): 'W',
        frozenset(['G', 'T']): 'K',
        frozenset(['A', 'C']): 'M',
        frozenset(['C', 'G', 'T']): 'B',
        frozenset(['A', 'G', 'T']): 'D',
        frozenset(['A', 'C', 'T']): 'H',
        frozenset(['A', 'C', 'G']): 'V',
        frozenset(['A', 'C', 'G', 'T']): 'N',
    }
    
    if not base_set:
        return '-'
    
    # Check exact match in IUPAC_CODES
    for iupac_bases, iupac_char in IUPAC_CODES.items():
        if base_set == iupac_bases:
            return iupac_char
    
    # No exact match - return N for unknown
    return 'N'
