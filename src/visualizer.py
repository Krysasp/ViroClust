#!/usr/bin/env python3
"""
Visualizer module for ViroClust.

Handles visualization of extended regions and sequence alignments.
"""

from typing import Dict, List, Optional, Tuple


def visualize_extended_regions(consensus: str, extended_region_info: Dict,
                               max_width: int = 80) -> str:
    """
    Create a visual representation of extended regions in a consensus.
    
    Args:
        consensus: consensus sequence
        extended_region_info: dict with extended region boundaries
        max_width: maximum width for display
    
    Returns:
        formatted visualization string
    """
    if not extended_region_info.get('has_extended_regions'):
        return None
    
    lines = [
        "=" * max_width,
        "EXTENDED REGION VISUALIZATION",
        "=" * max_width,
        "",
        f"Consensus length: {len(consensus)}",
        f"Extended leading: {extended_region_info['extended_leading_length']} positions "
        f"(0 - {extended_region_info['extended_leading_end']-1})",
        f"Central region: {extended_region_info['central_end'] - extended_region_info['central_start']} positions "
        f"({extended_region_info['central_start']} - {extended_region_info['central_end']-1})",
        f"Extended trailing: {extended_region_info['extended_trailing_length']} positions "
        f"({extended_region_info['extended_trailing_start']} - {extended_region_info['extended_trailing_end']-1})",
        "",
    ]
    
    # Create visual bar
    full_line = []
    label_line = []
    
    for pos in range(min(len(consensus), max_width)):
        base = consensus[pos]
        
        if pos < extended_region_info['extended_leading_end']:
            # Extended leading region - use uppercase
            full_line.append(base.upper())
            label_line.append('^')
        elif pos >= extended_region_info['central_start'] and pos < extended_region_info['central_end']:
            # Central region - use lowercase
            full_line.append(base.lower())
            label_line.append('.')
        else:
            # Extended trailing region - use uppercase
            full_line.append(base.upper())
            label_line.append('v')
    
    # Join and format
    lines.append("Consensus (showing first {} positions):".format(min(len(consensus), max_width)))
    lines.append("".join(full_line))
    lines.append("".join(label_line))
    lines.append("")
    lines.append("Legend: ^ = extended leading, . = central, v = extended trailing")
    lines.append("=" * max_width)
    
    return "\n".join(lines)


def visualize_length_distribution(cluster_data: Dict[str, List[Tuple]]) -> str:
    """
    Create a visual representation of sequence length distribution in a cluster.
    
    Args:
        cluster_data: dict mapping cluster_num to list of (header, sequence) tuples
    
    Returns:
        formatted visualization string
    """
    lines = [
        "=" * 80,
        "SEQUENCE LENGTH DISTRIBUTION",
        "=" * 80,
        "",
    ]
    
    for cluster_num, seq_list in sorted(cluster_data.items(), key=lambda x: int(x[0])):
        if not seq_list:
            continue
        
        lengths = [len(seq) for _, seq in seq_list if seq]
        if not lengths:
            continue
        
        min_len = min(lengths)
        max_len = max(lengths)
        avg_len = sum(lengths) / len(lengths)
        
        lines.append(f"Cluster {cluster_num}: {len(lengths)} sequences")
        lines.append(f"  Min: {min_len}, Max: {max_len}, Avg: {avg_len:.1f}")
        lines.append(f"  Ratio: {min_len/max_len:.3f}")
        
        # Create histogram-style visualization
        normalized = [l / max_len * 40 for l in lengths]  # Scale to 40 chars max
        for i, (header, _) in enumerate(seq_list):
            if i < len(lengths):
                bar_len = int(normalized[i])
                bar = "#" * bar_len
                lines.append(f"  {header[:20]:20s} |{bar:<40s}| {lengths[i]}")
        
        lines.append("")
    
    return "\n".join(lines)


def visualize_alignment_summary(aligned_seqs: Dict[str, str], 
                                consensus: str,
                                gap_threshold: float = 0.1) -> str:
    """
    Create a visual summary of alignment quality.
    
    Args:
        aligned_seqs: dict mapping header to aligned sequence
        consensus: consensus sequence
        gap_threshold: threshold for identifying high-gap positions
    
    Returns:
        formatted visualization string
    """
    if not aligned_seqs or not consensus:
        return "No alignment data available"
    
    lines = [
        "=" * 80,
        "ALIGNMENT QUALITY SUMMARY",
        "=" * 80,
        "",
        f"Total sequences: {len(aligned_seqs)}",
        f"Consensus length: {len(consensus)}",
        "",
    ]
    
    # Calculate gap frequency per position
    gap_freq = []
    for pos in range(len(consensus)):
        gaps = sum(1 for seq in aligned_seqs.values() if pos < len(seq) and seq[pos] == '-')
        gap_freq.append(gaps / len(aligned_seqs))
    
    # High gap positions
    high_gap_positions = [i for i, gf in enumerate(gap_freq) if gf > gap_threshold]
    
    if high_gap_positions:
        lines.append(f"High gap positions (>{gap_threshold*100:.0f}% gaps): {len(high_gap_positions)}")
        
        # Show first 10 high gap regions
        regions = []
        start = high_gap_positions[0] if high_gap_positions else None
        for pos in high_gap_positions[1:]:
            if pos != start + len(regions[-1] if regions else []):
                if start == high_gap_positions[0]:
                    regions.append((start, pos))
                else:
                    regions.append((regions[-1][1], pos))
            if start is None:
                start = pos
        
        # Simplified region detection
        if len(high_gap_positions) > 1:
            gap_regions = []
            start = high_gap_positions[0]
            end = high_gap_positions[0]
            
            for pos in high_gap_positions[1:]:
                if pos <= end + 1:
                    end = pos
                else:
                    gap_regions.append((start, end))
                    start = pos
                    end = pos
            gap_regions.append((start, end))
            
            lines.append("High gap regions:")
            for i, (start, end) in enumerate(gap_regions[:5]):
                lines.append(f"  Positions {start}-{end} ({end-start+1} bp)")
    else:
        lines.append(f"No high gap positions detected (threshold: {gap_threshold*100:.0f}%)")
    
    lines.append("")
    lines.append("=" * 80)
    
    return "\n".join(lines)


def get_visualization_report(cluster_num: str, consensus: str, 
                            extended_region_info: Dict,
                            aligned_seqs: Dict[str, str],
                            gap_threshold: float = 0.1) -> str:
    """
    Generate a comprehensive visualization report for a cluster.
    
    Args:
        cluster_num: cluster number
        consensus: consensus sequence
        extended_region_info: dict with extended region information
        aligned_seqs: dict mapping header to aligned sequence
        gap_threshold: threshold for high-gap positions
    
    Returns:
        formatted visualization report
    """
    reports = [
        f"Cluster {cluster_num} Visualization Report",
        "=" * 80,
        "",
    ]
    
    # Extended region visualization
    if extended_region_info.get('has_extended_regions'):
        reports.append(visualize_extended_regions(consensus, extended_region_info))
        reports.append("")
    
    # Alignment summary
    reports.append(visualize_alignment_summary(aligned_seqs, consensus, gap_threshold))
    
    return "\n".join(reports)
