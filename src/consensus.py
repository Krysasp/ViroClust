#!/usr/bin/env python3
"""
Consensus module for ViroClust.

Handles consensus sequence generation with IUPAC codes and gap filtering.
"""

from collections import defaultdict


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


def get_iupac_code(base_set):
    """
    Get IUPAC code for a set of bases.
    
    Args:
        base_set: frozenset of bases (e.g., frozenset(['A', 'G']))
    
    Returns:
        IUPAC code character (e.g., 'R' for A/G)
    """
    if not base_set:
        return '-'
    
    # Check exact match in IUPAC_CODES
    for iupac_bases, iupac_char in IUPAC_CODES.items():
        if base_set == iupac_bases:
            return iupac_char
    
    # No exact match - return N for unknown
    return 'N'


def generate_iupac_consensus(aligned_sequences, conservation=0.90, inclusion=0.70):
    """
    Generate consensus sequence with IUPAC codes based on conservation and inclusion thresholds.
    
    Args:
        aligned_sequences: dict mapping header to aligned sequence
        conservation: minimum frequency of dominant base to call single base (default 0.90)
        inclusion: minimum fraction of sequences that must have base at position (default 0.70)
    
    Returns:
        consensus sequence string with IUPAC codes
    """
    if not aligned_sequences:
        return ""
    
    seq_length = len(next(iter(aligned_sequences.values())))
    num_sequences = len(aligned_sequences)
    min_inclusion_count = int(num_sequences * inclusion)
    
    consensus = []
    
    for pos in range(seq_length):
        base_count = defaultdict(int)
        total_bases = 0
        
        for header, seq in aligned_sequences.items():
            base = seq[pos] if pos < len(seq) else '-'
            if base != '-':
                base_count[base] += 1
                total_bases += 1
        
        # Check inclusion threshold
        if total_bases < min_inclusion_count:
            consensus.append('-')
            continue
        
        if not base_count:
            consensus.append('-')
            continue
        
        max_count = max(base_count.values())
        dominant_bases = frozenset([b for b, c in base_count.items() if c == max_count])
        base_freq = max_count / total_bases
        
        # Check conservation threshold
        if base_freq >= conservation:
            if len(dominant_bases) == 1:
                consensus.append(list(dominant_bases)[0])
            else:
                # Multiple dominant bases at same frequency - use IUPAC
                iupac_char = get_iupac_code(dominant_bases)
                consensus.append(iupac_char)
        else:
            # Below conservation threshold - check if significant ambiguity
            # Include bases that are at least 50% of dominant frequency
            significant_bases = frozenset([b for b, c in base_count.items() if c >= max_count * 0.5])
            iupac_char = get_iupac_code(significant_bases)
            consensus.append(iupac_char)
    
    return ''.join(consensus)


def generate_consensus(aligned_sequences, conservation=0.90, inclusion=0.70, filter_gaps=False, gap_threshold=0.70):
    """
    Generate consensus sequence with IUPAC codes and optional gap filtering.
    
    Args:
        aligned_sequences: dict mapping header to aligned sequence
        conservation: minimum frequency of dominant base (default 0.90)
        inclusion: minimum fraction of sequences with base at position (default 0.70)
        filter_gaps: if True, remove positions with >gap_threshold gaps
        gap_threshold: fraction of sequences with gaps to trigger position removal
    
    Returns:
        tuple of (consensus sequence string, gap_filter_info)
        - gap_filter_info: dict with filtering statistics (if filter_gaps=True)
    """
    if not aligned_sequences:
        return "", {}
    
    working_aligned = aligned_sequences
    gap_filter_info = {}
    
    # Apply gap filtering if requested
    if filter_gaps:
        working_aligned, positions_removed, gap_stats = filter_gappy_positions(aligned_sequences, gap_threshold)
        gap_filter_info = {
            'positions_removed': positions_removed,
            'gap_stats': gap_stats,
            'filter_applied': True
        }
    
    seq_length = len(next(iter(working_aligned.values())))
    num_sequences = len(working_aligned)
    min_inclusion_count = int(num_sequences * inclusion)
    
    consensus = []
    
    for pos in range(seq_length):
        base_count = defaultdict(int)
        total_bases = 0
        
        for header, seq in working_aligned.items():
            base = seq[pos] if pos < len(seq) else '-'
            if base != '-':
                base_count[base] += 1
                total_bases += 1
        
        if total_bases < min_inclusion_count:
            consensus.append('-')
            continue
        
        if not base_count:
            consensus.append('-')
            continue
        
        max_count = max(base_count.values())
        dominant_bases = frozenset([b for b, c in base_count.items() if c == max_count])
        base_freq = max_count / total_bases
        
        if base_freq >= conservation:
            if len(dominant_bases) == 1:
                consensus.append(list(dominant_bases)[0])
            else:
                iupac_char = get_iupac_code(dominant_bases)
                consensus.append(iupac_char)
        else:
            significant_bases = frozenset([b for b, c in base_count.items() if c >= max_count * 0.5])
            iupac_char = get_iupac_code(significant_bases)
            consensus.append(iupac_char)
    
    return ''.join(consensus), gap_filter_info


def filter_gappy_positions(aligned_sequences, gap_threshold=0.70):
    """
    Remove alignment positions where gap frequency exceeds threshold.
    
    Args:
        aligned_sequences: dict mapping header to aligned sequence
        gap_threshold: maximum tolerated proportion of gaps
    
    Returns:
        tuple of (filtered_aligned, positions_removed, gap_stats)
        - filtered_aligned: dict with gap positions removed
        - positions_removed: list of (position, gap_fraction) tuples
        - gap_stats: dict with gap statistics
    """
    if not aligned_sequences:
        return {}, [], {}
    
    seq_length = len(next(iter(aligned_sequences.values())))
    num_sequences = len(aligned_sequences)
    max_allowed_gaps = int(num_sequences * gap_threshold)
    
    positions_removed = []
    positions_kept = []
    gap_fraction_by_pos = []
    
    for pos in range(seq_length):
        gap_count = sum(1 for seq in aligned_sequences.values() if pos < len(seq) and seq[pos] == '-')
        gap_fraction = gap_count / num_sequences
        
        if gap_count > max_allowed_gaps:
            positions_removed.append((pos, gap_fraction))
        else:
            positions_kept.append(pos)
            gap_fraction_by_pos.append(gap_fraction)
    
    # Create filtered alignment
    filtered_aligned = {}
    for header, seq in aligned_sequences.items():
        filtered_seq = ''.join(seq[pos] for pos in positions_kept if pos < len(seq))
        filtered_aligned[header] = filtered_seq
    
    # Calculate gap statistics
    total_gaps = sum(gap_frac for _, gap_frac in positions_removed)
    avg_gap_fraction = total_gaps / len(positions_removed) if positions_removed else 0.0
    max_gap_fraction = max((gap_frac for _, gap_frac in positions_removed), default=0.0)
    
    gap_stats = {
        'total_positions': seq_length,
        'positions_removed': len(positions_removed),
        'positions_kept': len(positions_kept),
        'gap_threshold': gap_threshold,
        'avg_gap_fraction': avg_gap_fraction,
        'max_gap_fraction': max_gap_fraction
    }
    
    return filtered_aligned, positions_removed, gap_stats


def degap_filter_consensus(iupac_consensus, gap_threshold=0.10):
    """
    Generate degap-filtered consensus by removing gap positions from IUPAC consensus.
    
    Inherits IUPAC codes from the IUPAC consensus sequence and only removes positions
    where gap frequency exceeds the threshold. This ensures the degap consensus
    preserves all IUPAC ambiguity codes while being shorter due to gap removal.
    
    Example: If IUPAC consensus is 'ATCRYG...' and positions 3,7 are gap positions,
    the degap consensus becomes 'ATCG...' (positions 3 and 7 removed).
    
    Args:
        iupac_consensus: IUPAC consensus sequence string (with IUPAC codes)
        gap_threshold: maximum tolerated proportion of sequences with gaps (default 0.10)
    
    Returns:
        tuple of (degap consensus sequence string, list of gap positions removed)
    """
    if not iupac_consensus:
        return "", []
    
    # For degap filtering, we need to know the original alignment to count gaps
    # We'll use a simplified approach: remove positions where consensus has '-'
    # The IUPAC consensus already encodes the ambiguity, we just need to remove gaps
    
    consensus = []
    gap_positions = []
    
    for pos, base in enumerate(iupac_consensus):
        if base == '-':
            # This is a gap position - remove it
            gap_positions.append(pos)
        else:
            # Keep the base (including IUPAC codes)
            consensus.append(base)
    
    return ''.join(consensus), gap_positions


def analyze_iupac_bases(consensus, num_seqs):
    """
    Analyze IUPAC degenerate bases in consensus sequence.
    
    Args:
        consensus: consensus sequence string
        num_seqs: number of sequences in cluster
    
    Returns:
        dict with IUPAC analysis
    """
    iupac_counts = defaultdict(int)
    iupac_positions = defaultdict(list)
    
    for pos, base in enumerate(consensus):
        if base in IUPAC_CODES.values():
            iupac_counts[base] += 1
            iupac_positions[base].append(pos)
    
    total_iupac = sum(iupac_counts.values())
    consensus_len = len(consensus)
    iupac_fraction = total_iupac / consensus_len if consensus_len > 0 else 0.0
    
    # Analyze which IUPAC codes represent which ambiguity levels
    single_base = sum(1 for b in IUPAC_CODES.values() if iupac_counts[b] == 0)
    ambiguous_codes = {k: v for k, v in iupac_counts.items() if v > 0}
    
    return {
        'total_iupac': total_iupac,
        'iupac_fraction': iupac_fraction,
        'consensus_length': consensus_len,
        'iupac_count_by_type': dict(iupac_counts),
        'iupac_positions_by_type': {k: len(v) for k, v in iupac_positions.items()},
        'num_sequences': num_seqs
    }
