#!/usr/bin/env python3
"""
Progress module for ViroClust.

Handles progress bar display using tqdm.
"""

try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False


class ProgressManager:
    """
    Manages progress bar display for ViroClust.
    """
    
    def __init__(self, config: dict):
        """
        Initialize progress manager.
        
        Args:
            config: configuration dict with 'enable_progress_bars' key
        """
        self.enabled = config.get('enable_progress_bars', True)
    
    def create_cluster_progress(self, total: int, desc: str = "Clusters") -> 'ClusterProgress':
        """
        Create a progress manager for cluster processing.
        
        Args:
            total: total number of clusters
            desc: description for the progress bar
        
        Returns:
            ClusterProgress instance
        """
        return ClusterProgress(total, desc, self.enabled)
    
    def create_alignment_progress(self, total: int, desc: str = "Aligning") -> 'AlignmentProgress':
        """
        Create a progress manager for alignment tasks.
        
        Args:
            total: total number of alignments
            desc: description for the progress bar
        
        Returns:
            AlignmentProgress instance
        """
        return AlignmentProgress(total, desc, self.enabled)


class ClusterProgress:
    """
    Progress tracking for cluster processing.
    """
    
    def __init__(self, total: int, desc: str = "Clusters", enabled: bool = True):
        """
        Initialize cluster progress.
        
        Args:
            total: total number of clusters
            desc: description
            enabled: whether to show progress bar
        """
        self.total = total
        self.desc = desc
        self.enabled = enabled and TQDM_AVAILABLE
        self.processed = 0
        self.failed = 0
        self.skipped = 0
        self.progress_bar = None
    
    def __enter__(self):
        """Enter context manager."""
        if self.enabled:
            self.progress_bar = tqdm(
                total=self.total,
                desc=self.desc,
                unit='cluster',
                colour='blue'
            )
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Exit context manager."""
        if self.progress_bar:
            self.progress_bar.close()
    
    def update(self, cluster_num: str = None, status: str = "DONE"):
        """
        Update progress.
        
        Args:
            cluster_num: cluster number (optional)
            status: status string (DONE, FAILED, SKIPPED)
        """
        if self.enabled and self.progress_bar:
            self.progress_bar.update(1)
            if cluster_num:
                suffix = f"Cluster {cluster_num}: {status}"
                self.progress_bar.set_postfix_str(suffix)
        
        self.processed += 1
        if 'fail' in status.lower():
            self.failed += 1
        elif 'skip' in status.lower():
            self.skipped += 1
    
    def set_description(self, desc: str):
        """
        Set progress bar description.
        
        Args:
            desc: new description
        """
        if self.progress_bar:
            self.progress_bar.set_description(desc)
    
    def print_status(self, message: str):
        """
        Print a status message.
        
        Args:
            message: message to print
        """
        if self.enabled and self.progress_bar:
            self.progress_bar.write(message)


class AlignmentProgress:
    """
    Progress tracking for alignment tasks within a cluster.
    """
    
    def __init__(self, total: int, desc: str = "Aligning", enabled: bool = True):
        """
        Initialize alignment progress.
        
        Args:
            total: total number of sequences
            desc: description
            enabled: whether to show progress bar
        """
        self.total = total
        self.desc = desc
        self.enabled = enabled and TQDM_AVAILABLE
        self.processed = 0
        self.progress_bar = None
    
    def __enter__(self):
        """Enter context manager."""
        if self.enabled:
            self.progress_bar = tqdm(
                total=self.total,
                desc=self.desc,
                unit='seq',
                colour='green',
                leave=False
            )
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Exit context manager."""
        if self.progress_bar:
            self.progress_bar.close()
    
    def update(self, status: str = "OK"):
        """
        Update progress.
        
        Args:
            status: status string
        """
        if self.enabled and self.progress_bar:
            self.progress_bar.update(1)
        
        self.processed += 1
    
    def set_description(self, desc: str):
        """
        Set progress bar description.
        
        Args:
            desc: new description
        """
        if self.progress_bar:
            self.progress_bar.set_description(desc)


def format_progress_status(cluster_num: str, status: str, num_seqs: int, 
                          iupac_bases: int = 0) -> str:
    """
    Format a progress status message.
    
    Args:
        cluster_num: cluster number
        status: status string
        num_seqs: number of sequences
        iupac_bases: number of IUPAC bases (optional)
    
    Returns:
        formatted status string
    """
    base = f"[{cluster_num}] {num_seqs} seqs"
    if iupac_bases > 0:
        base += f" ({iupac_bases} IUPAC)"
    base += f" - {status}"
    return base
