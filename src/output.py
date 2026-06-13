#!/usr/bin/env python3
"""
Output module for ViroClust.

Handles writing of various output files including FASTA, CSV, and summary reports.
"""

import os
from collections import defaultdict


def write_gap_analysis_summary(report_path, all_stats, cluster_data):
    """
    Write gap analysis summary to text file.
    
    Args:
        report_path: output report file path
        all_stats: dict mapping cluster_num to stats_dict
        cluster_data: dict mapping cluster_num to list of (header, seq) tuples
    """
    with open(report_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("GAP ANALYSIS SUMMARY\n")
        f.write("=" * 80 + "\n\n")
        
        f.write(f"{'Cluster':>8} {'Sequences':>10} {'Status':>12} {'OrigLen':>10} {'AlnLen':>10} {'Gaps':>8} {'Gap%':>8} {'Ins':>6}\n")
        f.write("-" * 80 + "\n")
        
        for cluster_num in sorted(all_stats.keys(), key=int):
            stats = all_stats.get(cluster_num, {})
            
            orig_len = stats.get('original_length', 0)
            aln_len = stats.get('alignment_length', 0)
            num_seqs = len(cluster_data.get(cluster_num, []))
            gaps = stats.get('gaps_in_consensus', 0)
            gap_pct = stats.get('gap_fraction', 0.0) * 100
            insertions = stats.get('insertions', 0)
            
            if orig_len == aln_len:
                status = 'SKIPPED'
            else:
                status = 'ALIGNED'
            
            f.write(f"{cluster_num:>8} {num_seqs:>10} {status:>12} {orig_len:>10} {aln_len:>10} {gaps:>8} {gap_pct:>8.2f} {insertions:>6}\n")
        
        f.write("-" * 80 + "\n")
        f.write("=" * 80 + "\n")
    
    print(f"\nWrote gap analysis summary to: {report_path}")


def write_gap_filter_csv(csv_path, gap_filter_data):
    """
    Write gap filtering analysis to CSV file.
    Reports on positions removed due to gap threshold exceeded.
    
    Args:
        csv_path: output CSV file path
        gap_filter_data: dict mapping cluster_num to (positions_removed, gap_stats) tuple
    """
    import csv
    
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Header row
        header = ['cluster', 'total_positions', 'positions_removed', 'positions_kept',
                  'gap_threshold', 'avg_gap_fraction', 'max_gap_fraction',
                  'removed_positions', 'removed_gap_fractions', 'gap_filter_impact']
        writer.writerow(header)
        
        # Data rows
        for cluster_num in sorted(gap_filter_data.keys(), key=int):
            positions_removed, gap_stats = gap_filter_data[cluster_num]
            
            # Format removed positions and fractions as semicolon-separated lists
            removed_positions = ";".join(str(pos) for pos, _ in positions_removed)
            removed_fractions = ";".join(f"{frac:.4f}" for _, frac in positions_removed)
            
            # Calculate impact
            total_pos = gap_stats.get('total_positions', 0)
            removed_count = len(positions_removed)
            if total_pos > 0:
                impact_pct = (removed_count / total_pos) * 100
                if impact_pct < 1:
                    impact = 'MINIMAL'
                elif impact_pct < 5:
                    impact = 'LOW'
                elif impact_pct < 20:
                    impact = 'MODERATE'
                else:
                    impact = 'HIGH'
            else:
                impact_pct = 0
                impact = 'NONE'
            
            row = [
                cluster_num,
                gap_stats.get('total_positions', 0),
                removed_count,
                gap_stats.get('positions_kept', 0),
                f"{gap_stats.get('gap_threshold', 0.70):.2f}",
                f"{gap_stats.get('avg_gap_fraction', 0.0):.4f}",
                f"{gap_stats.get('max_gap_fraction', 0.0):.4f}",
                removed_positions,
                removed_fractions,
                impact
            ]
            
            writer.writerow(row)
    
    print(f"Wrote gap filter analysis to: {csv_path}")


def write_iupac_csv(csv_path, all_iupac, all_stats, cluster_data):
    """
    Write IUPAC analysis to CSV file.
    
    Args:
        csv_path: output CSV file path
        all_iupac: dict mapping cluster_num to iupac_analysis dict
        cluster_data: dict mapping cluster_num to list of (header, seq) tuples
    """
    import csv
    
    # IUPAC codes in order
    iupac_order = ['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N']
    
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Header row
        header = ['cluster', 'status', 'num_sequences', 'consensus_length',
                  'total_iupac', 'iupac_fraction']
        # Add IUPAC code counts and positions
        for base in iupac_order:
            header.append(f'{base}_count')
            header.append(f'{base}_positions')
        writer.writerow(header)
        
        # Data rows
        for cluster_num in sorted(all_iupac.keys(), key=int):
            iupac_data = all_iupac.get(cluster_num, {})
            stats = all_stats.get(cluster_num, {})
            
            orig_len = stats.get('original_length', 0)
            aln_len = stats.get('alignment_length', 0)
            num_seqs = len(cluster_data.get(cluster_num, []))
            
            if orig_len == aln_len:
                status = 'SKIPPED'
            else:
                status = 'ALIGNED'
            
            row = [
                cluster_num,
                status,
                num_seqs,
                iupac_data['consensus_length'],
                iupac_data['total_iupac'],
                f"{iupac_data['iupac_fraction']*100:.2f}%"
            ]
            
            # Add IUPAC counts and positions
            iupac_count_str = ''
            for base in iupac_order:
                count = iupac_data['iupac_count_by_type'].get(base, 0)
                positions = iupac_data['iupac_positions_by_type'].get(base, 0)
                row.append(count)
                row.append(positions)
            
            writer.writerow(row)
    
    print(f"Wrote IUPAC analysis to: {csv_path}")


def write_oligo_analysis_csv(csv_path, all_consensus, all_iupac_consensus, all_gap_filter):
    """
    Write oligo analysis CSV comparing gap-filtered vs IUPAC consensus.
    Identifies positions where gap removal may affect oligo design.
    
    Args:
        csv_path: output CSV file path
        all_consensus: dict mapping consensus header to gap-filtered sequence
        all_iupac_consensus: dict mapping consensus header to IUPAC sequence
        all_gap_filter: dict mapping cluster_num to (positions_removed, gap_stats)
    """
    import csv
    
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Header row
        header = ['cluster', 'gap_filtered_length', 'iupac_length', 'length_difference',
                  'positions_removed_by_gap_filter', 'bases_lost_in_gap_regions', 
                  'iupac_bases_in_removed_positions', 'percent_length_change', 
                  'oligo_impact_assessment']
        writer.writerow(header)
        
        # Data rows
        for cluster_num in sorted(all_gap_filter.keys(), key=int):
            positions_removed, gap_stats = all_gap_filter[cluster_num]
            
            # Get consensus sequences
            gap_filtered_seq = None
            iupac_seq = None
            
            for header, seq in all_consensus.items():
                clean_header = header.lstrip('>')
                if f"|{cluster_num}|" in clean_header and "_IUPAC" not in clean_header:
                    gap_filtered_seq = seq
                    break
            
            for header, seq in all_iupac_consensus.items():
                clean_header = header.lstrip('>')
                if f"|{cluster_num}|" in clean_header and "_IUPAC" in clean_header:
                    iupac_seq = seq
                    break
            
            if gap_filtered_seq is None or iupac_seq is None:
                writer.writerow([cluster_num, 0, 0, 0, 0, 0, 0, 0.00, 'MISSING'])
                continue
            
            # Analyze bases lost in gap regions
            bases_lost = 0
            iupac_bases_in_removed = defaultdict(int)
            
            for pos, gap_frac in positions_removed:
                if pos < len(iupac_seq):
                    base = iupac_seq[pos]
                    if base != '-':
                        bases_lost += 1
                        iupac_bases_in_removed[base] += 1
            
            length_diff = len(iupac_seq) - len(gap_filtered_seq)
            percent_change = (length_diff / len(iupac_seq) * 100) if len(iupac_seq) > 0 else 0
            
            # Assess oligo impact
            if bases_lost == 0:
                impact = 'NONE'
            elif bases_lost < 10:
                impact = 'MINIMAL'
            elif bases_lost < 50:
                impact = 'LOW'
            elif bases_lost < 100:
                impact = 'MODERATE'
            elif bases_lost < 200:
                impact = 'HIGH'
            else:
                impact = 'SEVERE'
            
            # Format IUPAC bases as string
            iupac_base_str = ''.join(f"{b}:{c}" for b, c in sorted(iupac_bases_in_removed.items()))
            
            writer.writerow([
                cluster_num,
                len(gap_filtered_seq),
                len(iupac_seq),
                length_diff,
                len(positions_removed),
                bases_lost,
                iupac_base_str,
                f"{percent_change:.2f}",
                impact
            ])
    
    print(f"Wrote oligo analysis to: {csv_path}")


def write_degap_analysis_csv(csv_path, all_degap_consensus, all_iupac_consensus, all_gap_filter):
    """
    Write degap analysis CSV comparing degap-filtered vs IUPAC consensus.
    Shows which positions were removed due to gap filtering.
    
    Args:
        csv_path: output CSV file path
        all_degap_consensus: dict mapping consensus header to degap-filtered sequence
        all_iupac_consensus: dict mapping consensus header to IUPAC sequence
        all_gap_filter: dict mapping cluster_num to (positions_removed, gap_stats)
    """
    import csv
    
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Header row
        header = ['cluster', 'iupac_length', 'degap_length', 'length_difference',
                  'gap_positions_removed', 'gap_threshold', 'avg_gap_fraction',
                  'oligo_length_change', 'degap_impact']
        writer.writerow(header)
        
        # Data rows
        for cluster_num in sorted(all_gap_filter.keys(), key=int):
            positions_removed, gap_stats = all_gap_filter[cluster_num]
            
            # Get consensus sequences
            degap_seq = None
            iupac_seq = None
            
            for header, seq in all_degap_consensus.items():
                clean_header = header.lstrip('>')
                if f"|{cluster_num}|" in clean_header and "_degap" in clean_header:
                    degap_seq = seq
                    break
            
            for header, seq in all_iupac_consensus.items():
                clean_header = header.lstrip('>')
                if f"|{cluster_num}|" in clean_header and "_IUPAC" in clean_header:
                    iupac_seq = seq
                    break
            
            if degap_seq is None or iupac_seq is None:
                writer.writerow([cluster_num, 0, 0, 0, 0, 0.00, 0.00, 0, 'MISSING'])
                continue
            
            # Calculate metrics
            length_diff = len(iupac_seq) - len(degap_seq)
            gap_positions_str = ";".join(str(pos) for pos, _ in positions_removed)
            avg_gap_frac = gap_stats.get('avg_gap_fraction', 0.0)
            
            # Assess impact based on positions removed
            if len(positions_removed) == 0:
                impact = 'NONE'
            elif len(positions_removed) < 5:
                impact = 'MINIMAL'
            elif len(positions_removed) < 20:
                impact = 'LOW'
            elif len(positions_removed) < 50:
                impact = 'MODERATE'
            else:
                impact = 'HIGH'
            
            writer.writerow([
                cluster_num,
                len(iupac_seq),
                len(degap_seq),
                length_diff,
                len(positions_removed),
                f"{gap_stats.get('gap_threshold', 0.70):.2f}",
                f"{avg_gap_frac:.4f}",
                length_diff,
                impact
            ])
    
    print(f"Wrote degap analysis to: {csv_path}")


def write_missing_sequence_report(report_path, fasta_sequences, missing_sequences):
    """
    Report sequences that are listed in cluster file but not found in input FASTA.
    
    Args:
        report_path: output report file path
        fasta_sequences: dict mapping FASTA headers to sequences
        missing_sequences: dict mapping cluster_num to list of missing accession IDs
    
    Returns:
        dict with missing sequence information
    """
    # Build set of accessions from FASTA file
    fasta_accessions = set()
    for header in fasta_sequences.keys():
        accession = header.split('|')[0]
        fasta_accessions.add(accession)
    
    # Calculate total missing sequences
    total_missing = sum(len(missing_list) for missing_list in missing_sequences.values())
    
    # Write report
    with open(report_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("MISSING SEQUENCE REPORT\n")
        f.write("Sequences in cluster file but not found in input FASTA\n")
        f.write("=" * 80 + "\n\n")
        
        f.write("SUMMARY\n")
        f.write("-" * 80 + "\n")
        f.write(f"Total sequences in input FASTA:     {len(fasta_sequences)}\n")
        f.write(f"Total missing sequences:            {total_missing}\n\n")
        
        if total_missing > 0:
            f.write("MISSING SEQUENCES BY CLUSTER\n")
            f.write("-" * 80 + "\n")
            
            for cluster_num in sorted(missing_sequences.keys(), key=int):
                missing = missing_sequences[cluster_num]
                f.write(f"Cluster {cluster_num}: {len(missing)} missing sequence(s)\n")
                for acc in missing:
                    f.write(f"  - {acc}\n")
                f.write("\n")
                f.write("\n")
    
    print(f"Wrote missing sequence report to: {report_path}")


def write_length_disparity_csv(csv_path, disparity_info, reference_info, cluster_data):
    """
    Write length disparity analysis to CSV file.
    Reports on clusters with significant sequence length differences.
    
    For clusters with length disparity >= 20% (configurable):
    - Reference sequence (marked with '*') is identified
    - Extended regions (leading/trailing) are calculated
    - Longer subcluster count is tracked (sequences within 5% of reference length)
    
    Args:
        csv_path: output CSV file path
        disparity_info: dict mapping cluster_num to disparity details
        reference_info: dict mapping cluster_num to reference sequence info
        cluster_data: dict mapping cluster_num to list of (header, seq) tuples
    """
    import csv
    
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Header row
        header = ['cluster', 'num_sequences', 'shortest_length', 'longest_length', 
                  'reference_length', 'reference_id', 'length_ratio', 'has_disparity',
                  'extended_leading_regions', 'extended_trailing_regions', 
                  'extended_leading_length', 'extended_trailing_length',
                  'longer_subcluster_count']
        writer.writerow(header)
        
        # Data rows - only clusters with length disparity
        for cluster_num in sorted(disparity_info.keys(), key=int):
            info = disparity_info[cluster_num]
            
            writer.writerow([
                cluster_num,
                info['num_sequences'],
                info['shortest_length'],
                info['longest_length'],
                info.get('reference_length', info['longest_length']),
                info.get('reference_id', ''),
                f"{info['ratio']:.4f}",
                info['has_disparity'],
                info.get('extended_leading_regions', 0),
                info.get('extended_trailing_regions', 0),
                info.get('extended_leading_length', 0),
                info.get('extended_trailing_length', 0),
                info.get('longer_subcluster_count', info['num_sequences'])
            ])
    
    print(f"Wrote length disparity analysis to: {csv_path}")


def write_quality_report(csv_path, all_quality, all_stats):
    """
    Write quality summary report to CSV file.
    
    Args:
        csv_path: output CSV file path
        all_quality: dict mapping cluster_num to quality_data dict
        all_stats: dict mapping cluster_num to stats_dict
    """
    import csv
    
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Header row
        header = ['cluster', 'consensus_length', 'mean_quality', 'min_quality', 
                  'max_quality', 'high_quality_positions', 'low_quality_positions',
                  'high_quality_pct', 'quality_distribution']
        writer.writerow(header)
        
        # Data rows
        for cluster_num in sorted(all_quality.keys(), key=int):
            quality_data = all_quality.get(cluster_num, {})
            stats = all_stats.get(cluster_num, {})
            
            if not quality_data:
                writer.writerow([cluster_num, 0, 0.0, 0.0, 0.0, 0, 0, 0.0, 'EMPTY'])
                continue
            
            consensus_len = stats.get('alignment_length', 0)
            mean_quality = quality_data.get('mean_quality', 0.0)
            min_quality = quality_data.get('min_quality', 0.0)
            max_quality = quality_data.get('max_quality', 0.0)
            high_quality = quality_data.get('high_quality_positions', 0)
            low_quality = quality_data.get('low_quality_positions', 0)
            high_quality_pct = (high_quality / consensus_len * 100) if consensus_len > 0 else 0.0
            
            # Quality distribution string
            quality_dist = quality_data.get('quality_distribution', {})
            dist_str = f"H:{quality_dist.get('high', 0)} M:{quality_dist.get('medium', 0)} L:{quality_dist.get('low', 0)} V:{quality_dist.get('very_low', 0)}"
            
            writer.writerow([
                cluster_num,
                consensus_len,
                f"{mean_quality:.4f}",
                f"{min_quality:.4f}",
                f"{max_quality:.4f}",
                high_quality,
                low_quality,
                f"{high_quality_pct:.2f}",
                dist_str
            ])
    
    print(f"Wrote quality summary to: {csv_path}")


def write_metadata_csv(csv_path, metadata_index, cluster_data):
    """
    Write comprehensive sequence metadata to CSV file.
    
    Args:
        csv_path: output CSV file path
        metadata_index: dict mapping cluster_num to list of enriched metadata tuples
        cluster_data: dict mapping cluster_num to list of (header, seq, length, is_ref) tuples
    """
    import csv
    
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Header row
        header = ['cluster', 'accession', 'original_header', 'taxID', 'species', 
                  'sequence_length', 'is_reference', 'renamed_header', 
                  'additional_info']
        writer.writerow(header)
        
        # Data rows
        for cluster_num in sorted(metadata_index.keys(), key=int):
            enriched_seqs = metadata_index.get(cluster_num, [])
            
            for item in enriched_seqs:
                header, seq, length, is_ref, meta, taxID, species, renamed_header = item
                
                # Format additional info as semicolon-separated key=value pairs
                additional_str = ";".join(f"{k}={v}" for k, v in meta.get('additional_info', {}).items())
                
                writer.writerow([
                    cluster_num,
                    meta['accession'],
                    header,
                    taxID or '',
                    species or '',
                    length,
                    is_ref,
                    renamed_header,
                    additional_str
                ])
    
    print(f"Wrote sequence metadata to: {csv_path}")


def write_renamed_sequences_fasta(fasta_path, metadata_index, cluster_data):
    """
    Write FASTA file with renamed sequences using standardized naming convention.
    
    Naming format: CLUST{cluster_num}_TAX{taxID}_{accession}
    Example: CLUST0001_TAX9606_MN486048
    
    Args:
        fasta_path: output FASTA file path
        metadata_index: dict mapping cluster_num to list of enriched metadata tuples
        cluster_data: dict mapping cluster_num to list of (header, seq, length, is_ref) tuples
    """
    with open(fasta_path, 'w') as f:
        for cluster_num in sorted(metadata_index.keys(), key=int):
            enriched_seqs = metadata_index.get(cluster_num, [])
            
            for item in enriched_seqs:
                orig_header, seq, length, is_ref, meta, taxID, species, renamed_header = item
                
                # Write sequence with renamed header
                f.write(f">{renamed_header}\n{seq}\n")
    
    print(f"Wrote renamed sequences to: {fasta_path}")


def write_sequence_rename_map(map_path, metadata_index):
    """
    Write mapping file showing original header to renamed header correspondence.
    
    Args:
        map_path: output mapping file path
        metadata_index: dict mapping cluster_num to list of enriched metadata tuples
    """
    with open(map_path, 'w') as f:
        f.write("# Sequence Renaming Map\n")
        f.write("# Format: original_header -> renamed_header\n")
        f.write("#" + "=" * 80 + "\n\n")
        
        for cluster_num in sorted(metadata_index.keys(), key=int):
            enriched_seqs = metadata_index.get(cluster_num, [])
            
            f.write(f"\n# Cluster {cluster_num}\n")
            f.write("-" * 80 + "\n")
            
            for item in enriched_seqs:
                orig_header, seq, length, is_ref, meta, taxID, species, renamed_header = item
                f.write(f"{orig_header} -> {renamed_header}\n")
    
    print(f"Wrote sequence rename map to: {map_path}")
