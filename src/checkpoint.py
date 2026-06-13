#!/usr/bin/env python3
"""
Checkpoint module for ViroClust.

Handles saving and restoring progress for interrupted runs.
"""

import json
import os
from typing import Dict, Any, List, Optional
from datetime import datetime


class CheckpointManager:
    """
    Manages checkpointing for ViroClust runs.
    
    Allows resuming from interruption by saving state periodically.
    """
    
    def __init__(self, output_dir: str, config: Dict[str, Any]):
        """
        Initialize checkpoint manager.
        
        Args:
            output_dir: output directory for checkpoints
            config: configuration dict
        """
        self.output_dir = output_dir
        self.config = config
        self.checkpoint_file = os.path.join(output_dir, 'checkpoint.json')
        self.enabled = config.get('enable_checkpointing', True)
        self.interval = config.get('checkpoint_interval', 10)
    
    def save(self, 
             cluster_results: Dict[str, Any],
             processed_clusters: List[str],
             start_time: datetime) -> None:
        """
        Save checkpoint state.
        
        Args:
            cluster_results: dict of completed cluster results
            processed_clusters: list of processed cluster numbers
            start_time: start time of the run
        """
        if not self.enabled:
            return
        
        checkpoint = {
            'config': self.config,
            'processed_clusters': processed_clusters,
            'cluster_results': cluster_results,
            'start_time': start_time.isoformat(),
            'last_update': datetime.now().isoformat(),
            'total_clusters': len(processed_clusters),
            'status': 'in_progress'
        }
        
        # Write to temp file first, then rename (atomic operation)
        temp_file = self.checkpoint_file + '.tmp'
        with open(temp_file, 'w') as f:
            json.dump(checkpoint, f, indent=2)
        os.rename(temp_file, self.checkpoint_file)
    
    def load(self) -> Optional[Dict[str, Any]]:
        """
        Load checkpoint if it exists.
        
        Returns:
            checkpoint dict or None if not found
        """
        if not os.path.exists(self.checkpoint_file):
            return None
        
        with open(self.checkpoint_file, 'r') as f:
            return json.load(f)
    
    def resume(self) -> tuple:
        """
        Resume from checkpoint if available.
        
        Returns:
            tuple of (cluster_results, processed_clusters, start_time, should_resume)
        """
        checkpoint = self.load()
        
        if checkpoint is None:
            return {}, [], datetime.now(), False
        
        return (
            checkpoint.get('cluster_results', {}),
            checkpoint.get('processed_clusters', []),
            datetime.fromisoformat(checkpoint.get('start_time', datetime.now().isoformat())),
            True
        )
    
    def complete(self, cluster_results: Dict[str, Any] = None, 
                 processed_clusters: List[str] = None) -> None:
        """
        Mark checkpoint as complete.
        
        Args:
            cluster_results: dict of completed cluster results (optional)
            processed_clusters: list of processed cluster numbers (optional)
        """
        if not self.enabled:
            return
        
        checkpoint = self.load()
        
        # If no existing checkpoint, create a new one
        if checkpoint is None and processed_clusters:
            checkpoint = {
                'config': self.config,
                'processed_clusters': processed_clusters,
                'cluster_results': cluster_results or {},
                'start_time': datetime.now().isoformat(),
                'last_update': datetime.now().isoformat(),
                'total_clusters': len(processed_clusters),
                'status': 'completed',
                'end_time': datetime.now().isoformat()
            }
        elif checkpoint:
            checkpoint['status'] = 'completed'
            checkpoint['end_time'] = datetime.now().isoformat()
            if cluster_results is not None:
                checkpoint['cluster_results'] = cluster_results
            if processed_clusters is not None:
                checkpoint['processed_clusters'] = processed_clusters
                checkpoint['total_clusters'] = len(processed_clusters)
            checkpoint['last_update'] = datetime.now().isoformat()
        else:
            # Nothing to checkpoint
            return
        
        with open(self.checkpoint_file, 'w') as f:
            json.dump(checkpoint, f, indent=2)
    
    def should_process_cluster(self, cluster_num: str) -> bool:
        """
        Check if cluster should be processed (not already completed).
        
        Args:
            cluster_num: cluster number to check
    
        Returns:
            True if cluster should be processed
        """
        if not self.enabled:
            return True
        
        checkpoint = self.load()
        if checkpoint is None:
            return True
        
        return cluster_num not in checkpoint.get('processed_clusters', [])


def write_checkpoint_summary(output_dir: str, checkpoint_manager: CheckpointManager) -> None:
    """
    Write checkpoint summary report.
    
    Args:
        output_dir: output directory
        checkpoint_manager: checkpoint manager instance
    """
    checkpoint = checkpoint_manager.load()
    
    if checkpoint is None:
        return
    
    summary_path = os.path.join(output_dir, 'checkpoint_summary.txt')
    
    with open(summary_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("CHECKPOINT SUMMARY\n")
        f.write("=" * 80 + "\n\n")
        
        f.write(f"Status: {checkpoint.get('status', 'unknown')}\n")
        f.write(f"Start time: {checkpoint.get('start_time', 'N/A')}\n")
        f.write(f"Last update: {checkpoint.get('last_update', 'N/A')}\n")
        f.write(f"Total clusters: {checkpoint.get('total_clusters', 0)}\n")
        f.write(f"Checkpointing enabled: {checkpoint_manager.enabled}\n")
        f.write(f"Checkpoint interval: {checkpoint_manager.interval}\n\n")
        
        if checkpoint.get('config'):
            f.write("Configuration:\n")
            for key, value in checkpoint['config'].items():
                if key not in ['threads', 'enable_progress_bars', 'enable_checkpointing']:
                    f.write(f"  {key}: {value}\n")
        f.write("\n")
        
        if checkpoint.get('processed_clusters'):
            f.write("Processed clusters:\n")
            for cluster_num in checkpoint['processed_clusters']:
                f.write(f"  - Cluster {cluster_num}\n")
        
        f.write("\n" + "=" * 80 + "\n")
    
    print(f"Wrote checkpoint summary to: {summary_path}")
