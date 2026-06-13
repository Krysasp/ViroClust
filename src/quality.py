#!/usr/bin/env python3
"""
Quality module for ViroClust.

Handles quality scoring for consensus positions and consensus validation.
"""

from collections import defaultdict
from typing import Dict, List, Tuple, Optional, Any


# IUPAC ambiguity codes and their quality scores
# Quality score reflects confidence based on agreement percentage
IUPAC_QUALITY = {
    'A': 1.0, 'C': 1.0, 'G': 1.0, 'T': 1.0,  # Single bases - highest confidence
    'R': 0.9, 'Y': 0.9, 'S': 0.9, 'W': 0.9, 'K': 0.9, 'M': 0.9,  # 2-way ambiguity
    'B': 0.7, 'D': 0.7, 'H': 0.7, 'V': 0.7,  # 3-way ambiguity
    'N': 0.5,  # Any base - lowest confidence
    '-': 0.0,  # Gap
}


def calculate_position_quality(aligned_seqs: Dict[str, str], pos: int, 
                               consensus_base: str) -> float:
    """
    Calculate quality score for a consensus position.
    
    Quality score is based on:
    1. Agreement percentage (how many sequences match the consensus base)
    2. IUPAC ambiguity (lower for ambiguous codes)
    3. Gap frequency (lower if many gaps)
    
    Args:
        aligned_seqs: dict mapping header to aligned sequence
        pos: position in alignment
        consensus_base: consensus base at this position
    
    Returns:
        quality score between 0.0 and 1.0
    """
    if not aligned_seqs:
        return 0.0
    
    num_seqs = len(aligned_seqs)
    
    # Count bases at this position
    base_counts = defaultdict(int)
    gap_count = 0
    
    for header, seq in aligned_seqs.items():
        if pos < len(seq):
            base = seq[pos]
            if base == '-':
                gap_count += 1
            else:
                base_counts[base] += 1
        else:
            gap_count += 1  # Missing sequence treated as gap
    
    # Calculate agreement percentage
    if consensus_base == '-':
        # Gap position - quality based on gap frequency
        agreement = gap_count / num_seqs
        quality = agreement * IUPAC_QUALITY.get('-', 0.0)
    else:
        # Non-gap position - quality based on agreement and IUPAC
        if consensus_base in base_counts:
            agreement = base_counts[consensus_base] / num_seqs
        else:
            # Consensus is IUPAC code - calculate based on ambiguity
            agreement = sum(base_counts.get(b, 0) for b in get_iupac_bases(consensus_base)) / num_seqs
        
        iupac_quality = IUPAC_QUALITY.get(consensus_base, 0.5)
        quality = agreement * iupac_quality
    
    return quality


def get_iupac_bases(iupac_code: str) -> List[str]:
    """
    Get the list of bases represented by an IUPAC code.
    
    Args:
        iupac_code: IUPAC ambiguity code
    
    Returns:
        list of base characters
    """
    iupac_mapping = {
        'R': ['A', 'G'],
        'Y': ['C', 'T'],
        'S': ['G', 'C'],
        'W': ['A', 'T'],
        'K': ['G', 'T'],
        'M': ['A', 'C'],
        'B': ['C', 'G', 'T'],
        'D': ['A', 'G', 'T'],
        'H': ['A', 'C', 'T'],
        'V': ['A', 'C', 'G'],
        'N': ['A', 'C', 'G', 'T'],
        'A': ['A'],
        'C': ['C'],
        'G': ['G'],
        'T': ['T'],
        '-': ['-'],
    }
    return iupac_mapping.get(iupac_code, [])


def calculate_consensus_quality(consensus: str, aligned_seqs: Dict[str, str],
                               threshold: float = 0.5) -> Dict[str, float]:
    """
    Calculate quality scores for all positions in a consensus sequence.
    
    Args:
        consensus: consensus sequence
        aligned_seqs: dict mapping header to aligned sequence
        threshold: minimum quality score for acceptable positions
    
    Returns:
        dict with quality metrics:
        {
            'position_scores': list of quality scores per position,
            'mean_quality': mean quality score,
            'min_quality': minimum quality score,
            'max_quality': maximum quality score,
            'high_quality_positions': number of positions above threshold,
            'low_quality_positions': number of positions below threshold,
            'quality_distribution': dict of quality ranges with counts
        }
    """
    if not consensus:
        return {
            'position_scores': [],
            'mean_quality': 0.0,
            'min_quality': 0.0,
            'max_quality': 0.0,
            'high_quality_positions': 0,
            'low_quality_positions': 0,
            'quality_distribution': {}
        }
    
    position_scores = []
    quality_distribution = {
        'high': 0,      # >= 0.9
        'medium': 0,    # >= 0.7 and < 0.9
        'low': 0,       # >= 0.5 and < 0.7
        'very_low': 0,  # < 0.5
    }
    
    for pos, base in enumerate(consensus):
        score = calculate_position_quality(aligned_seqs, pos, base)
        position_scores.append(score)
        
        # Categorize by quality
        if score >= 0.9:
            quality_distribution['high'] += 1
        elif score >= 0.7:
            quality_distribution['medium'] += 1
        elif score >= 0.5:
            quality_distribution['low'] += 1
        else:
            quality_distribution['very_low'] += 1
    
    mean_quality = sum(position_scores) / len(position_scores) if position_scores else 0.0
    min_quality = min(position_scores) if position_scores else 0.0
    max_quality = max(position_scores) if position_scores else 0.0
    
    return {
        'position_scores': position_scores,
        'mean_quality': mean_quality,
        'min_quality': min_quality,
        'max_quality': max_quality,
        'high_quality_positions': quality_distribution['high'] + quality_distribution['medium'],
        'low_quality_positions': quality_distribution['low'] + quality_distribution['very_low'],
        'quality_distribution': quality_distribution
    }


def validate_consensus_against_reference(consensus: str, reference: str,
                                         gap_threshold: float = 0.1) -> Dict[str, Any]:
    """
    Validate consensus sequence against a known reference sequence.
    
    Args:
        consensus: consensus sequence to validate
        reference: reference sequence
        gap_threshold: maximum tolerated gap fraction
    
    Returns:
        dict with validation metrics:
        {
            'alignment_length': int,
            'matches': int,
            'mismatches': int,
            'gaps': int,
            'identity': float,
            'gap_fraction': float,
            'mismatch_positions': list of mismatch positions,
            'gap_positions': list of gap positions
        }
    """
    min_len = min(len(consensus), len(reference))
    
    matches = 0
    mismatches = 0
    gaps = 0
    mismatch_positions = []
    gap_positions = []
    
    for pos in range(min_len):
        cons_base = consensus[pos]
        ref_base = reference[pos]
        
        if cons_base == '-':
            gaps += 1
            gap_positions.append(pos)
        elif ref_base == '-':
            gaps += 1
            gap_positions.append(pos)
        elif cons_base == ref_base:
            matches += 1
        else:
            # Check if it's an IUPAC ambiguity match
            ref_bases = get_iupac_bases(ref_base) if len(ref_base) > 1 else [ref_base]
            cons_bases = get_iupac_bases(cons_base) if len(cons_base) > 1 else [cons_base]
            
            # Check for overlap
            if set(ref_bases) & set(cons_bases):
                matches += 1  # Ambiguous match
            else:
                mismatches += 1
                mismatch_positions.append(pos)
    
    identity = matches / min_len if min_len > 0 else 0.0
    gap_fraction = gaps / min_len if min_len > 0 else 0.0
    
    return {
        'alignment_length': min_len,
        'matches': matches,
        'mismatches': mismatches,
        'gaps': gaps,
        'identity': identity,
        'gap_fraction': gap_fraction,
        'mismatch_positions': mismatch_positions,
        'gap_positions': gap_positions
    }


def get_quality_report(consensus: str, aligned_seqs: Dict[str, str],
                       reference: Optional[str] = None) -> str:
    """
    Generate a quality report for a consensus sequence.
    
    Args:
        consensus: consensus sequence
        aligned_seqs: dict mapping header to aligned sequence
        reference: optional reference sequence for validation
    
    Returns:
        formatted quality report string
    """
    quality_data = calculate_consensus_quality(consensus, aligned_seqs)
    
    lines = [
        "=" * 60,
        "CONSENSUS QUALITY REPORT",
        "=" * 60,
        "",
        "Overall Metrics:",
        f"  Mean quality: {quality_data['mean_quality']:.3f}",
        f"  Min quality: {quality_data['min_quality']:.3f}",
        f"  Max quality: {quality_data['max_quality']:.3f}",
        "",
        "Quality Distribution:",
        f"  High (>=0.9): {quality_data['quality_distribution']['high']} positions",
        f"  Medium (>=0.7): {quality_data['quality_distribution']['medium']} positions",
        f"  Low (>=0.5): {quality_data['quality_distribution']['low']} positions",
        f"  Very low (<0.5): {quality_data['quality_distribution']['very_low']} positions",
        "",
        f"High quality positions: {quality_data['high_quality_positions']}",
        f"Low quality positions: {quality_data['low_quality_positions']}",
    ]
    
    if reference:
        validation = validate_consensus_against_reference(consensus, reference)
        lines.extend([
            "",
            "Reference Validation:",
            f"  Identity: {validation['identity']*100:.2f}%",
            f"  Matches: {validation['matches']}/{validation['alignment_length']}",
            f"  Mismatches: {validation['mismatches']}",
            f"  Gaps: {validation['gaps']}",
        ])
    
    lines.append("=" * 60)
    
    return "\n".join(lines)
