"""
ViroClust package for sequence alignment and consensus generation.
"""

from .config import load_config, save_config, validate_config, calculate_timeout
from .parsers import (
    parse_fasta,
    parse_cluster_file,
    extract_fasta_metadata,
    build_sequence_metadata_index
)
from .alignment import align_sequences_mafft
from .consensus import generate_consensus, generate_iupac_consensus, IUPAC_CODES
from .disparity import detect_length_disparity, identify_reference_sequence
from .output import write_gap_analysis_summary, write_gap_filter_csv, write_iupac_csv
from .checkpoint import CheckpointManager, write_checkpoint_summary
from .progress import ProgressManager
from .quality import calculate_consensus_quality, get_quality_report
from .visualizer import get_visualization_report

__version__ = '1.0.0'
__author__ = 'ViroClust Team'

__all__ = [
    'load_config',
    'save_config',
    'validate_config',
    'calculate_timeout',
    'parse_fasta',
    'parse_cluster_file',
    'extract_fasta_metadata',
    'build_sequence_metadata_index',
    'align_sequences_mafft',
    'generate_consensus',
    'generate_iupac_consensus',
    'IUPAC_CODES',
    'detect_length_disparity',
    'identify_reference_sequence',
    'write_gap_analysis_summary',
    'write_gap_filter_csv',
    'write_iupac_csv',
    'CheckpointManager',
    'write_checkpoint_summary',
    'ProgressManager',
    'calculate_consensus_quality',
    'get_quality_report',
    'get_visualization_report',
]