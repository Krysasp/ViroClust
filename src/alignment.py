#!/usr/bin/env python3
"""
Alignment module for ViroClust.

Handles MAFFT alignment with optional reference sequence prioritization.
"""

import os
import subprocess
import tempfile


def align_sequences_mafft(seqs_with_headers, output_path, cluster_num, timeout=120):
    """
    Align sequences using MAFFT with reference sequence prioritization.
    
    If the first sequence in seqs_with_headers is marked as reference (has '*' marker),
    MAFFT will use it as the guide tree anchor for better alignment of divergent sequences.
    
    Args:
        seqs_with_headers: list of (header, sequence) tuples, with reference sequence first if available
        output_path: path to write aligned sequences
        cluster_num: cluster number for logging
        timeout: alignment timeout in seconds
    
    Returns:
        tuple of (aligned_dict, status_msg)
        - aligned_dict: dict mapping header to aligned sequence
        - status_msg: status message describing alignment outcome
    """
    num_seqs = len(seqs_with_headers)
    # Handle both 2-tuples (header, seq) and 4-tuples (header, seq, length, is_reference)
    seq_size = 0
    for item in seqs_with_headers:
        if len(item) == 4:
            _, seq, length, _ = item
            seq_size += length
        else:
            _, seq = item
            seq_size += len(seq)
    
    # Create temp input file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp_in:
        for item in seqs_with_headers:
            if len(item) == 4:
                header, seq, length, is_ref = item
            else:
                header, seq = item
            tmp_in.write(f">{header}\n{seq}\n")
        tmp_in_path = tmp_in.name
    
    # Show queued status
    cluster_info = f"Cluster {cluster_num}" if cluster_num else "Cluster"
    
    try:
        # Run MAFFT with --quiet to suppress verbose output
        # Using --thread 1 for predictable behavior in parallel processing
        result = subprocess.run(
            ['mafft', '--thread', '1', '--auto', '--maxiterate', '1000', '--localpair', '--quiet', tmp_in_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=timeout
        )
        
        # Check for errors
        if result.returncode != 0:
            status_msg = f"Alignment failed with return code {result.returncode}"
            print(f"[ALIGN] {cluster_info}: {status_msg}")
            os.remove(tmp_in_path)
            return {}, status_msg
        
        # Parse aligned output
        aligned = {}
        current_header = None
        current_seq = []
        
        for line in result.stdout.split('\n'):
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_header:
                    aligned[current_header] = ''.join(current_seq)
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line.upper())
        
        if current_header:
            aligned[current_header] = ''.join(current_seq)
        
        # Check if alignment produced output
        if not aligned:
            status_msg = "Alignment produced empty result"
            print(f"[ALIGN] {cluster_info}: {status_msg}")
            os.remove(tmp_in_path)
            return {}, status_msg
        
        # Write to output file
        # Preserve the order from seqs_with_headers (reference sequence first)
        with open(output_path, 'w') as f:
            for item in seqs_with_headers:
                if len(item) == 4:
                    header, seq, length, is_ref = item
                else:
                    header, seq = item
                if header in aligned:
                    f.write(f">{header}\n{aligned[header]}\n")
        
        status_msg = f"Completed ({len(aligned)} sequences, {sum(len(s) for s in aligned.values())} bases)"
        print(f"[ALIGN] {cluster_info}: {status_msg}")
        os.remove(tmp_in_path)
        return aligned, status_msg
        
    except subprocess.TimeoutExpired as e:
        status_msg = f"Alignment timeout after {timeout}s"
        print(f"[ALIGN] {cluster_info}: {status_msg}")
        # Try to read any partial output
        try:
            with open(tmp_in_path, 'r') as f:
                content = f.read()
                if content:
                    # Parse partial output if available
                    aligned = {}
                    current_header = None
                    current_seq = []
                    for line in content.split('\n'):
                        line = line.strip()
                        if not line:
                            continue
                        if line.startswith('>'):
                            if current_header:
                                aligned[current_header] = ''.join(current_seq)
                            current_header = line[1:]
                            current_seq = []
                        else:
                            current_seq.append(line.upper())
                    if current_header:
                        aligned[current_header] = ''.join(current_seq)
                    os.remove(tmp_in_path)
                    return aligned, f"Partial alignment ({len(aligned)} sequences)"
        except:
            pass
        os.remove(tmp_in_path)
        return {}, status_msg
    except Exception as e:
        status_msg = f"Alignment error: {str(e)}"
        print(f"[ALIGN] {cluster_info}: {status_msg}")
        os.remove(tmp_in_path)
        return {}, status_msg
