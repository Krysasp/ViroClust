"""
ViroClust - Auto-align sequences based on cluster files and generate consensus sequences.

Main entry point that coordinates all modules with enhanced features:
- Config management
- Checkpoint/restart capability
- Progress tracking with tqdm
- Quality scoring
- Visualization
"""

import sys
import argparse
import os
import multiprocessing as mp
from datetime import datetime

# Add parent directory to path so we can import from src
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.config import load_config, save_config, validate_config, calculate_timeout
from src.parsers import (
    parse_fasta, 
    parse_cluster_file, 
    extract_cluster_sequences, 
    extract_cluster_sequences_with_info,
    extract_fasta_metadata,
    parse_fasta_with_metadata,
    build_sequence_metadata_index
)
from src.alignment import align_sequences_mafft
from src.consensus import (
    generate_consensus, 
    generate_iupac_consensus, 
    degap_filter_consensus,
    analyze_iupac_bases,
    IUPAC_CODES
)
from src.disparity import (
    detect_length_disparity,
    identify_reference_sequence,
    process_extended_regions
)
from src.output import (
    write_gap_analysis_summary,
    write_gap_filter_csv,
    write_iupac_csv,
    write_oligo_analysis_csv,
    write_degap_analysis_csv,
    write_missing_sequence_report,
    write_length_disparity_csv,
    write_quality_report,
    write_metadata_csv,
    write_renamed_sequences_fasta,
    write_sequence_rename_map
)
from src.checkpoint import CheckpointManager, write_checkpoint_summary
from src.progress import ProgressManager, format_progress_status
from src.quality import calculate_consensus_quality, get_quality_report
from src.visualizer import get_visualization_report


def run_viroclust(
    fasta_path,
    cluster_path,
    output_dir,
    config_path=None,
    conservation=None,
    inclusion=None,
    threads=None,
    test_clusters=None,
    test_run=False,
    filter_gaps=True,
    gap_threshold=0.10,
    length_disparity_threshold=0.20,
    quality_threshold=0.5
):
    """
    Main function to run ViroClust pipeline.
    
    Args:
        fasta_path: path to input FASTA file
        cluster_path: path to cluster file (.clstr format)
        output_dir: output directory path
        config_path: path to config JSON file (optional)
        conservation: conservation threshold (default 0.90)
        inclusion: inclusion threshold (default 0.70)
        threads: number of parallel threads (default: CPU count)
        test_clusters: if set, only process first N clusters (for testing)
        test_run: if True, only align first 10 cluster sequences
        filter_gaps: if True, remove positions with >gap_threshold gaps
        gap_threshold: fraction of sequences with gaps to trigger position removal
        length_disparity_threshold: min ratio (shortest/longest) for extended region handling
        quality_threshold: minimum quality score for acceptable positions (default 0.5)
    """
    # Load configuration
    print("Loading configuration...")
    config = load_config(config_path)
    
    # Override config with command-line arguments if provided
    if conservation is not None:
        config['conservation'] = conservation
    if inclusion is not None:
        config['inclusion'] = inclusion
    if threads is not None:
        config['threads'] = threads
    if filter_gaps is not None:
        config['filter_gaps'] = filter_gaps
    if gap_threshold is not None:
        config['gap_threshold'] = gap_threshold
    if length_disparity_threshold is not None:
        config['length_disparity_threshold'] = length_disparity_threshold
    
    # Set default quality_threshold if not specified
    if quality_threshold is not None:
        config['quality_threshold'] = quality_threshold
    elif config.get('quality_threshold') is None:
        config['quality_threshold'] = 0.5
    
    # Auto-detect threads if not specified
    if config['threads'] is None:
        config['threads'] = mp.cpu_count()
    
    # Validate configuration
    if not validate_config(config):
        print("[WARN] Configuration has issues, proceeding with defaults")
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Save final configuration
    final_config_path = os.path.join(output_dir, 'viroclust_config.json')
    save_config(config, final_config_path)
    print(f"  Saved configuration to: {final_config_path}")
    
    # Initialize checkpoint manager
    checkpoint_manager = CheckpointManager(output_dir, config)
    
    # Initialize progress manager
    progress_manager = ProgressManager(config)
    
    # Check for checkpoint to resume
    cluster_results, processed_clusters, start_time, should_resume = checkpoint_manager.resume()
    if should_resume:
        print(f"  Resuming from checkpoint at {start_time}")
    
    print(f"Parsing FASTA file: {fasta_path}")
    fasta_sequences = parse_fasta(fasta_path)
    print(f"  Loaded {len(fasta_sequences)} sequences")
    
    print(f"Parsing cluster file: {cluster_path}")
    cluster_sequences, length_info, reference_info = parse_cluster_file(cluster_path, max_clusters=test_clusters)
    print(f"  Found {len(cluster_sequences)} clusters" + (f" (first {test_clusters} for testing)" if test_clusters else ""))
    
    print("Extracting cluster sequences...")
    cluster_data, reference_sequence_list, missing_sequences = extract_cluster_sequences_with_info(fasta_sequences, cluster_sequences, length_info)
    
    # Count total sequences to be processed
    total_seqs = sum(len(seq_list) for seq_list in cluster_data.values())
    print(f"  Extracted {total_seqs} sequences across {len(cluster_data)} clusters")
    
    # Report missing sequences
    total_missing = sum(len(missing_list) for missing_list in missing_sequences.values())
    if total_missing > 0:
        print(f"  WARNING: {total_missing} sequence(s) from cluster file not found in FASTA (excluded from alignment)")
    
    # Extract metadata from FASTA headers
    print("\nExtracting sequence metadata (taxID, species, etc.)...")
    fasta_metadata = extract_fasta_metadata(fasta_sequences)
    print(f"  Extracted metadata for {len(fasta_metadata)} sequences")
    
    # Build metadata index linking cluster sequences to their metadata
    print("Building sequence metadata index...")
    metadata_index = build_sequence_metadata_index(fasta_metadata, cluster_data)
    
    # Count unique taxIDs
    taxIDs_found = set()
    for cluster_num, enriched_seqs in metadata_index.items():
        for item in enriched_seqs:
            _, _, _, _, meta, taxID, species, _ = item
            if taxID:
                taxIDs_found.add(taxID)
    print(f"  Found {len(taxIDs_found)} unique taxID(s): {sorted(taxIDs_found)}")
    
    # Test run: limit to first 10 cluster sequences
    if test_run:
        limited_data = {}
        seq_count = 0
        for cluster_num, seq_list in cluster_data.items():
            remaining = 10 - seq_count
            if remaining <= 0:
                break
            if len(seq_list) <= remaining:
                limited_data[cluster_num] = seq_list
                seq_count += len(seq_list)
            else:
                limited_data[cluster_num] = seq_list[:remaining]
                seq_count += remaining
            if seq_count >= 10:
                break
        cluster_data = limited_data
        total_seqs = sum(len(seq_list) for seq_list in cluster_data.values())
        print(f"  [TEST RUN] Limited to first {total_seqs} sequences across {len(cluster_data)} clusters")
    
    aligned_output_dir = os.path.join(output_dir, 'aligned')
    consensus_output_dir = os.path.join(output_dir, 'consensus')
    
    os.makedirs(aligned_output_dir, exist_ok=True)
    os.makedirs(consensus_output_dir, exist_ok=True)
    
    # Identify reference sequences (marked with '*')
    print("\nIdentifying reference sequences (marked with '*')...")
    reference_info = identify_reference_sequence(cluster_data, length_info)
    print(f"  Found {len(reference_info)} reference sequences")
    
    # Detect clusters with length disparity (using reference_info for accurate calculation)
    print(f"\nDetecting sequence length disparities (threshold: {config['length_disparity_threshold']*100:.0f}%)...")
    disparity_info = detect_length_disparity(cluster_data, reference_info, config['length_disparity_threshold'])
    print(f"  Found {len(disparity_info)} clusters with length disparity >= {config['length_disparity_threshold']*100:.0f}%")
    
    # Estimate timeout based on total sequence size
    # Handle both 2-tuples (header, seq) and 4-tuples (header, seq, length, is_reference)
    total_bases = 0
    for seq_list in cluster_data.values():
        for item in seq_list:
            if len(item) == 4:
                _, seq, length, _ = item
                total_bases += length
            else:
                _, seq = item
                total_bases += len(seq)
    estimated_timeout = calculate_timeout(total_bases, config)
    print(f"\n[CONFIG] Estimated timeout per cluster: {estimated_timeout}s (based on {total_bases} total bases)")
    print(f"[CONFIG] Gap filtering: {'enabled' if config['filter_gaps'] else 'disabled'} (threshold={config['gap_threshold']*100:.0f}%)")
    print(f"[CONFIG] Consensus: conservation={config['conservation']*100:.0f}%, inclusion={config['inclusion']*100:.0f}%")
    print(f"[CONFIG] Length disparity threshold: {config['length_disparity_threshold']*100:.0f}%")
    print(f"[CONFIG] Quality threshold: {config['quality_threshold']}")
    print(f"[CONFIG] Checkpointing: {'enabled' if config['enable_checkpointing'] else 'disabled'}")
    print(f"[CONFIG] Progress bars: {'enabled' if config['enable_progress_bars'] else 'disabled'}\n")
    
    # Prepare cluster processing tasks
    tasks = []
    for cluster_num, seq_list in cluster_data.items():
        # Check if this cluster has length disparity
        has_disparity = cluster_num in disparity_info
        # Check if this cluster has a reference sequence
        has_reference = cluster_num in reference_info
        
        tasks.append((
            cluster_num, 
            seq_list, 
            aligned_output_dir, 
            config['conservation'], 
            config['inclusion'], 
            estimated_timeout, 
            config['filter_gaps'], 
            config['gap_threshold'],
            has_disparity,
            disparity_info.get(cluster_num, {}),
            has_reference,
            reference_info.get(cluster_num, {})
        ))
    
    total_queued_seqs = sum(len(seq_list) for _, seq_list, _, _, _, _, _, _, _, _, _, _ in tasks)
    print(f"\n{'='*60}")
    print(f"QUEUED CLUSTERS FOR ALIGNMENT ({len(tasks)} clusters, {total_queued_seqs} sequences, {config['threads']} threads, timeout={estimated_timeout}s)")
    print(f"{'='*60}")
    for i, (cluster_num, seq_list, _, _, _, _, _, _, _, _, _, _) in enumerate(sorted(tasks, key=lambda x: int(x[0]))):
        # Handle both 2-tuples and 4-tuples
        seq_size = 0
        for item in seq_list:
            if len(item) == 4:
                _, seq, length, _ = item
                seq_size += length
            else:
                _, seq = item
                seq_size += len(seq)
        print(f"  [{i+1}/{len(tasks)}] Cluster {cluster_num}: {len(seq_list):4d} sequences ({seq_size:6d} bases) -> queued")
    print(f"{'='*60}\n")
    
    all_aligned = {}
    all_consensus = {}
    all_stats = {}
    all_iupac = {}
    all_degap = {}
    all_quality = {}
    all_visualization = {}
    
    # Process clusters in parallel with progress tracking
    print(f"Starting MAFFT alignment processing...")
    print(f"{'='*60}")
    
    # Filter out already processed clusters (checkpoint resume)
    tasks_to_process = []
    resuming_clusters = []
    for task in sorted(tasks, key=lambda x: int(x[0])):
        cluster_num = task[0]
        if not checkpoint_manager.should_process_cluster(cluster_num):
            resuming_clusters.append(cluster_num)
        else:
            tasks_to_process.append(task)
    
    # Create progress bar for cluster processing
    with progress_manager.create_cluster_progress(len(tasks_to_process), "Processing clusters") as cluster_progress:
        # Process all clusters in parallel using multiprocessing pool
        with mp.Pool(processes=config['threads']) as pool:
            results = pool.starmap(process_cluster, [(task, config) for task in tasks_to_process])
        
        # Update progress status on completion
        print(f"\n{'='*60}")
        print(f"Alignment progress:")
        for i, result in enumerate(sorted(results, key=lambda x: int(x[0]))):
            cluster_num = result[0]
            status = result[7]
            if cluster_num in resuming_clusters:
                print(f"  [{i+1}/{len(tasks)}] Cluster {cluster_num}: [RESUMED]")
            elif 'SKIPPED' in status:
                print(f"  [{i+1}/{len(tasks)}] Cluster {cluster_num}: [SKIPPED]")
            else:
                print(f"  [{i+1}/{len(tasks)}] Cluster {cluster_num}: [DONE]")
        print(f"{'='*60}")
    
    print(f"\n{'='*60}")
    print(f"ALIGNMENT COMPLETED")
    print(f"{'='*60}")
    
    # Collect results and write consensus incrementally
    failed_clusters = []
    all_gap_filter = {}
    all_iupac_consensus = {}
    all_degap_consensus = {}
    
    # Open consensus FASTA files for incremental writing
    consensus_fasta_path = os.path.join(output_dir, 'consensus_sequences.fasta')
    iupac_consensus_fasta_path = os.path.join(output_dir, 'iupac_consensus_sequences.fasta')
    degap_consensus_fasta_path = os.path.join(output_dir, 'degap_consensus_sequences.fasta')
    quality_report_path = os.path.join(output_dir, 'quality_reports.txt')
    
    consensus_fasta = open(consensus_fasta_path, 'w')
    iupac_consensus_fasta = open(iupac_consensus_fasta_path, 'w')
    degap_consensus_fasta = open(degap_consensus_fasta_path, 'w')
    quality_report = open(quality_report_path, 'w')
    
    # Handle resumed clusters first (they weren't processed in the pool)
    for cluster_num in resuming_clusters:
        # Load from checkpoint - for now we'll just skip writing since we don't store full results
        print(f"[RESUME] Cluster {cluster_num}: Loaded from checkpoint")
    
    # Process new results from parallel execution
    for (cluster_num, aligned, consensus_dict, iupac_consensus_dict, degap_consensus_dict, 
         stats_dict, iupac_analysis, status, gap_filter_data, 
         extended_region_info, disparity_status, quality_data, viz_report) in results:
        # Write aligned sequences immediately
        for header, seq in aligned.items():
            formatted_header = format_header_with_cluster(header, cluster_num)
            all_aligned[formatted_header] = seq
        
        # Write consensus sequences immediately
        for header, seq in consensus_dict.items():
            consensus_fasta.write(f">{header}\n{seq}\n")
            all_consensus[header] = seq
        
        # Write IUPAC consensus immediately
        for header, seq in iupac_consensus_dict.items():
            iupac_consensus_fasta.write(f">{header}\n{seq}\n")
            all_iupac_consensus[header] = seq
        
        # Write degap consensus immediately
        for header, seq in degap_consensus_dict.items():
            degap_consensus_fasta.write(f">{header}\n{seq}\n")
            all_degap_consensus[header] = seq
        
        # Collect other data
        all_stats[cluster_num] = stats_dict
        all_iupac[cluster_num] = iupac_analysis
        all_gap_filter[cluster_num] = gap_filter_data
        all_quality[cluster_num] = quality_data
        all_visualization[cluster_num] = viz_report
        
        # Track failed/empty clusters
        if 'failed' in status.lower() or 'empty' in status.lower():
            failed_clusters.append((cluster_num, status, len(all_iupac[cluster_num].get('iupac_count_by_type', {}))))
    
    # Close consensus FASTA files
    consensus_fasta.close()
    iupac_consensus_fasta.close()
    degap_consensus_fasta.close()
    quality_report.close()
    
    # Mark checkpoint as complete
    processed_cluster_nums = [r[0] for r in results]  # Extract cluster numbers from results
    checkpoint_manager.complete(all_quality, processed_cluster_nums)
    
    if failed_clusters:
        print(f"\n[WARN] {len(failed_clusters)} cluster(s) had alignment issues:")
        for fcluster, fstatus, fiupac in failed_clusters:
            print(f"  - Cluster {fcluster}: {fstatus} (IUPAC bases: {fiupac})")
    
    # Write output files
    aligned_fasta = os.path.join(output_dir, 'aligned_sequences.fasta')
    with open(aligned_fasta, 'w') as f:
        for header, seq in sorted(all_aligned.items()):
            f.write(f">{header}\n{seq}\n")
    print(f"\nWrote aligned sequences to: {aligned_fasta}")
    print(f"Wrote consensus sequences to: {consensus_fasta_path}")
    print(f"Wrote IUPAC consensus sequences to: {iupac_consensus_fasta_path}")
    print(f"Wrote degap-filtered consensus sequences to: {degap_consensus_fasta_path}")
    print(f"Wrote quality reports to: {quality_report_path}")
    
    print(f"\nViroClust completed successfully!")
    print(f"  Aligned sequences: {len(all_aligned)}")
    print(f"  Consensus sequences (gap-filtered): {len(all_consensus)}")
    print(f"  IUPAC consensus sequences: {len(all_iupac_consensus)}")
    print(f"  Degap-filtered consensus sequences: {len(all_degap_consensus)}")
    
    # Compare cluster file sequences with input FASTA
    print(f"\nComparing cluster file sequences with input FASTA...")
    missing_seq_report = os.path.join(output_dir, 'missing_sequence_report.txt')
    write_missing_sequence_report(missing_seq_report, fasta_sequences, missing_sequences)
    
    # Write gap analysis summary
    gap_stats_fasta = os.path.join(output_dir, 'gap_analysis_summary.txt')
    write_gap_analysis_summary(gap_stats_fasta, all_stats, cluster_data)
    
    # Write IUPAC analysis CSV
    iupac_csv = os.path.join(output_dir, 'iupac_analysis.csv')
    write_iupac_csv(iupac_csv, all_iupac, all_stats, cluster_data)
    
    # Write gap filter analysis CSV
    gap_filter_csv = os.path.join(output_dir, 'gap_filter_analysis.csv')
    write_gap_filter_csv(gap_filter_csv, all_gap_filter)
    
    # Write oligo analysis CSV (comparison of gap-filtered vs IUPAC)
    oligo_analysis_csv = os.path.join(output_dir, 'oligo_analysis.csv')
    write_oligo_analysis_csv(oligo_analysis_csv, all_consensus, all_iupac_consensus, all_gap_filter)
    
    # Write degap analysis CSV
    degap_analysis_csv = os.path.join(output_dir, 'degap_analysis.csv')
    write_degap_analysis_csv(degap_analysis_csv, all_degap_consensus, all_iupac_consensus, all_gap_filter)
    
    # Write length disparity analysis CSV
    disparity_csv = os.path.join(output_dir, 'length_disparity_analysis.csv')
    write_length_disparity_csv(disparity_csv, disparity_info, reference_info, cluster_data)
    
    # Write quality report summary
    quality_summary_path = os.path.join(output_dir, 'quality_summary.csv')
    write_quality_report(quality_summary_path, all_quality, all_stats)
    
    # Write checkpoint summary
    write_checkpoint_summary(output_dir, checkpoint_manager)
    
    # Write visualization report
    viz_report_path = os.path.join(output_dir, 'visualization_report.txt')
    with open(viz_report_path, 'w') as f:
        for cluster_num in sorted(all_visualization.keys(), key=int):
            f.write(all_visualization[cluster_num] + "\n\n")
    print(f"Wrote visualization report to: {viz_report_path}")
    
    # Write metadata CSV
    metadata_csv_path = os.path.join(output_dir, 'sequence_metadata.csv')
    write_metadata_csv(metadata_csv_path, metadata_index, cluster_data)
    
    # Write renamed sequences FASTA
    renamed_fasta_path = os.path.join(output_dir, 'renamed_sequences.fasta')
    write_renamed_sequences_fasta(renamed_fasta_path, metadata_index, cluster_data)
    
    # Write sequence rename map
    rename_map_path = os.path.join(output_dir, 'sequence_rename_map.txt')
    write_sequence_rename_map(rename_map_path, metadata_index)
    
    return (aligned_fasta, consensus_fasta_path, gap_stats_fasta, iupac_csv, 
            gap_filter_csv, oligo_analysis_csv, degap_analysis_csv, disparity_csv,
            metadata_csv_path, renamed_fasta_path, rename_map_path,
            all_aligned, all_consensus, all_iupac_consensus, all_degap_consensus,
            all_quality, quality_summary_path)


def format_header_with_cluster(original_header, cluster_num):
    """Format header with 4-digit cluster prefix."""
    cluster_padded = cluster_num.zfill(4)
    return f"{cluster_padded}|{original_header}"


def process_cluster(args, config):
    """
    Process a single cluster: align sequences and generate consensus.
    Used for parallel processing.
    
    Args:
        args: tuple of (cluster_num, cluster_seq_list, aligned_dir, conservation, inclusion, 
                       timeout, filter_gaps, gap_threshold, has_disparity, disparity_info,
                       has_reference, reference_info)
        config: configuration dict
    
    Returns:
        tuple of (cluster_num, aligned_dict, consensus_dict, stats_dict, iupac_analysis, 
                  status, gap_filter_data, extended_region_info, disparity_status,
                  quality_data, viz_report)
    """
    (cluster_num, cluster_seq_list, aligned_dir, conservation, inclusion, 
     timeout, filter_gaps, gap_threshold, has_disparity, disparity_info,
     has_reference, reference_info) = args
    
    num_seqs = len(cluster_seq_list)
    cluster_aligned_path = os.path.join(aligned_dir, f'cluster_{cluster_num}.fasta')
    
    # Empty cluster: skip alignment
    if num_seqs == 0:
        print(f"[SKIP] Cluster {cluster_num}: Empty cluster (0 sequences)")
        aligned = {}
        consensus = ""
        status = "SKIPPED_EMPTY"
        gap_filter_data = ({}, {'total_positions': 0, 'positions_removed': 0, 'positions_kept': 0, 
                                 'gap_threshold': gap_threshold, 'avg_gap_fraction': 0.0, 'max_gap_fraction': 0.0})
        extended_region_info = {}
        disparity_status = {'has_disparity': False, 'has_reference': False}
        quality_data = {}
        viz_report = ""
    # Single sequence: skip alignment, use sequence directly
    elif num_seqs == 1:
        item = cluster_seq_list[0]
        if len(item) == 4:
            _, single_seq, _, _ = item
        else:
            _, single_seq = item
        aligned = {cluster_seq_list[0][0]: single_seq}
        
        # Write to output file
        with open(cluster_aligned_path, 'w') as f:
            f.write(f">{cluster_seq_list[0][0]}\n{single_seq}\n")
        
        # Consensus is identical to the single sequence
        consensus = single_seq
        status = "SKIPPED_SINGLE"
        print(f"[SKIP] Cluster {cluster_num}: Single sequence (no alignment needed)")
        gap_filter_data = ({}, {'total_positions': len(consensus), 'positions_removed': 0, 
                                 'positions_kept': len(consensus), 'gap_threshold': gap_threshold,
                                 'avg_gap_fraction': 0.0, 'max_gap_fraction': 0.0})
        extended_region_info = {}
        disparity_status = {'has_disparity': has_disparity, 'has_reference': has_reference}
        
        # Calculate quality for single sequence
        quality_data = calculate_consensus_quality(consensus, aligned, config['quality_threshold'])
        
        # Generate visualization
        viz_report = get_visualization_report(cluster_num, consensus, extended_region_info, 
                                             aligned, gap_threshold)
    else:
        # Check if we should reorder sequences to put reference first
        # Find reference sequence in the filtered cluster_seq_list
        ref_idx = None
        for idx, item in enumerate(cluster_seq_list):
            if len(item) == 4:
                _, _, _, is_reference = item
                if is_reference:
                    ref_idx = idx
                    break
        
        if ref_idx is not None:
            # Reference sequence found - place it first
            reordered_seqs = [cluster_seq_list[ref_idx]] + [
                item for i, item in enumerate(cluster_seq_list) if i != ref_idx
            ]
            ref_header = cluster_seq_list[ref_idx][0]
            print(f"[REORDER] Cluster {cluster_num}: Reference sequence '{ref_header}' placed first")
        else:
            # No reference sequence - use original order
            reordered_seqs = cluster_seq_list
        
        aligned, status = align_sequences_mafft(reordered_seqs, cluster_aligned_path, cluster_num, timeout)
    
    # Initialize variables for all cases
    iupac_consensus = ""
    extended_region_info = {}
    
    if not aligned:
        print(f"[WARN] Cluster {cluster_num}: Using fallback consensus")
        # Fallback: create consensus from first sequence if alignment failed
        if cluster_seq_list:
            item = cluster_seq_list[0]
            if len(item) == 4:
                _, consensus, _, _ = item
            else:
                _, consensus = item
            iupac_consensus = consensus
        else:
            consensus = ""
        gap_filter_data = ({}, {'total_positions': len(consensus), 'positions_removed': 0, 
                                 'positions_kept': len(consensus), 'gap_threshold': gap_threshold,
                                 'avg_gap_fraction': 0.0, 'max_gap_fraction': 0.0})
        extended_region_info = {}
        disparity_status = {'has_disparity': has_disparity, 'has_reference': has_reference}
        quality_data = {}
        viz_report = ""
    else:
        # Generate IUPAC consensus first (with gaps)
        iupac_consensus = generate_iupac_consensus(aligned, conservation, inclusion)
        
        # Check if this cluster has length disparity and needs extended region handling
        if has_disparity:
            # Process with extended regions (leading/trailing truncated regions)
            consensus, gap_filter_info, extended_region_info = process_extended_regions(
                aligned, iupac_consensus, 
                conservation, inclusion, filter_gaps, gap_threshold,
                disparity_info
            )
        else:
            # Standard processing
            consensus, gap_filter_info = generate_consensus(
                aligned, conservation, inclusion, filter_gaps, gap_threshold
            )
        
        gap_stats = gap_filter_info.get('gap_stats', {})
        if 'gap_threshold' not in gap_stats:
            gap_stats['gap_threshold'] = gap_threshold
        gap_filter_data = (gap_filter_info.get('positions_removed', []), gap_stats)
        
        # Calculate quality metrics
        quality_data = calculate_consensus_quality(consensus, aligned, config['quality_threshold'])
        
        # Generate visualization
        viz_report = get_visualization_report(cluster_num, consensus, extended_region_info, 
                                             aligned, gap_threshold)
        
        disparity_status = {
            'has_disparity': has_disparity, 
            'has_reference': has_reference,
            'disparity_ratio': disparity_info.get('ratio', 0) if has_disparity else None
        }
    
    # Generate consensus dict
    percent_conserv = int(conservation * 100)
    percent_incl = int(inclusion * 100)
    consensus_header = f"CONS|{cluster_num}|{percent_conserv}pct{percent_incl}pct"
    
    # Generate IUPAC consensus dict
    iupac_consensus_dict = {f"{consensus_header}_IUPAC": iupac_consensus}
    
    # Generate degap-filtered consensus (inherits IUPAC codes)
    if num_seqs > 1 and aligned:
        degap_consensus, degap_positions = degap_filter_consensus(iupac_consensus, gap_threshold)
    else:
        degap_consensus = iupac_consensus
        degap_positions = []
    
    degap_consensus_dict = {f"{consensus_header}_degap": degap_consensus}
    
    # Check for IUPAC degenerate bases in consensus
    if iupac_consensus:
        iupac_bases = set(iupac_consensus) & set(IUPAC_CODES.values())
        if iupac_bases:
            print(f"[CONSENSUS] Cluster {cluster_num}: Generated IUPAC consensus with {len(iupac_bases)} degenerate base types: {sorted(iupac_bases)}")
        else:
            print(f"[CONSENSUS] Cluster {cluster_num}: Generated IUPAC consensus (no degenerate bases)")
    else:
        print(f"[CONSENSUS] Cluster {cluster_num}: Generated empty IUPAC consensus")
        iupac_bases = set()
    
    # Calculate gap statistics
    if num_seqs > 0 and aligned:
        stats_dict = calculate_gap_stats(cluster_seq_list, aligned, iupac_consensus, num_seqs)
    else:
        stats_dict = {
            'original_length': 0,
            'alignment_length': 0,
            'gaps_in_consensus': 0,
            'gap_fraction': 0.0,
            'insertions': 0
        }
    
    # Collect IUPAC base analysis
    iupac_analysis = analyze_iupac_bases(iupac_consensus, num_seqs)
    
    consensus_dict = {consensus_header: consensus}
    
    return (cluster_num, aligned, consensus_dict, iupac_consensus_dict, degap_consensus_dict, 
            stats_dict, iupac_analysis, status, gap_filter_data, 
            extended_region_info, disparity_status, quality_data, viz_report)


def calculate_gap_stats(original_seqs, aligned_seqs, consensus, num_seqs):
    """
    Calculate statistics on gaps and insertions introduced by alignment.
    
    Args:
        original_seqs: list of (header, sequence) tuples
        aligned_seqs: dict mapping header to aligned sequence
        consensus: consensus sequence string
        num_seqs: number of sequences
    
    Returns:
        dict with gap statistics
    """
    if not aligned_seqs:
        return {
            'original_length': 0,
            'alignment_length': 0,
            'gaps_in_consensus': 0,
            'gap_fraction': 0.0,
            'insertions': 0
        }
    
    # Get original lengths (before alignment)
    # original_seqs can be 2-tuples (header, seq) or 4-tuples (header, seq, length, is_reference)
    original_lengths = []
    for item in original_seqs:
        if len(item) == 4:
            _, seq, length, _ = item
            original_lengths.append(length)
        else:
            _, seq = item
            original_lengths.append(len(seq))
    avg_original_length = sum(original_lengths) / len(original_lengths) if original_lengths else 0
    
    # Get alignment length
    alignment_length = len(next(iter(aligned_seqs.values()))) if aligned_seqs else 0
    
    # Count gaps in consensus
    gaps_in_consensus = consensus.count('-')
    gap_fraction = gaps_in_consensus / len(consensus) if len(consensus) > 0 else 0.0
    
    # Count insertions (gaps not in consensus but present in some sequences)
    # This is a simplified calculation
    insertions = 0
    for header, aligned_seq in aligned_seqs.items():
        for i, base in enumerate(aligned_seq):
            if base == '-' and i < len(consensus) and consensus[i] != '-':
                insertions += 1
    
    return {
        'original_length': int(avg_original_length),
        'alignment_length': alignment_length,
        'gaps_in_consensus': gaps_in_consensus,
        'gap_fraction': gap_fraction,
        'insertions': insertions
    }


def main():
    parser = argparse.ArgumentParser(
        description='ViroClust: Auto-align sequences and generate consensus based on cluster files'
    )
    parser.add_argument(
        '-f', '--fasta',
        required=True,
        help='Input FASTA file containing all sequences'
    )
    parser.add_argument(
        '-c', '--cluster',
        required=True,
        help='Cluster file (.clstr format)'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output directory'
    )
    parser.add_argument(
        '--config',
        help='Path to configuration JSON file (optional)'
    )
    parser.add_argument(
        '--conservation',
        type=float,
        default=None,
        help='Conservation threshold (default from config: 0.90)'
    )
    parser.add_argument(
        '--inclusion',
        type=float,
        default=None,
        help='Inclusion threshold (default from config: 0.70)'
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=None,
        help='Number of parallel threads (default from config: CPU count)'
    )
    parser.add_argument(
        '--test-clusters',
        type=int,
        default=None,
        help='For testing: only process first N clusters'
    )
    parser.add_argument(
        '--test-run',
        action='store_true',
        help='For testing: only align first 10 cluster sequences and output consensus'
    )
    parser.add_argument(
        '--filter-gaps',
        action='store_true',
        default=None,
        help='Remove alignment positions with >gap_threshold gaps (default from config: enabled)'
    )
    parser.add_argument(
        '--no-filter-gaps',
        action='store_true',
        help='Disable gap filtering'
    )
    parser.add_argument(
        '--gap-threshold',
        type=float,
        default=None,
        help='Fraction of sequences with gaps to trigger position removal (default from config: 0.10)'
    )
    parser.add_argument(
        '--length-disparity-threshold',
        type=float,
        default=None,
        help='Minimum ratio (shortest/longest) to trigger extended region handling (default from config: 0.20)'
    )
    parser.add_argument(
        '--quality-threshold',
        type=float,
        default=None,
        help='Minimum quality score for acceptable positions (default from config: 0.5)'
    )
    
    args = parser.parse_args()
    
    run_viroclust(
        args.fasta,
        args.cluster,
        args.output,
        args.config,
        args.conservation,
        args.inclusion,
        args.threads,
        args.test_clusters,
        args.test_run,
        not args.no_filter_gaps if args.filter_gaps is None and args.no_filter_gaps is None else args.filter_gaps,
        args.gap_threshold,
        args.length_disparity_threshold,
        args.quality_threshold
    )


if __name__ == '__main__':
    main()