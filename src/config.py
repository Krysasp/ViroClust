"""
Configuration module for ViroClust.

Handles configuration loading, saving, and validation.
"""

import json
import os
from typing import Dict, Any, Optional


DEFAULT_CONFIG = {
    'conservation': 0.90,
    'inclusion': 0.70,
    'threads': None,  # Auto-detect CPU count
    'filter_gaps': True,
    'gap_threshold': 0.10,
    'length_disparity_threshold': 0.20,
    'timeout_base': 120,
    'timeout_per_base': 5000,
    'max_timeout': 1800,
    'quality_threshold': 0.5,  # Default quality threshold
    'enable_progress_bars': True,
    'enable_checkpointing': True,
    'checkpoint_interval': 10,
}


def load_config(config_path: Optional[str] = None) -> Dict[str, Any]:
    """
    Load configuration from JSON file or use defaults.
    
    Args:
        config_path: path to config JSON file (optional)
    
    Returns:
        dict with configuration parameters
    """
    if config_path and os.path.exists(config_path):
        with open(config_path, 'r') as f:
            config = json.load(f)
        # Merge with defaults
        config = {**DEFAULT_CONFIG, **config}
    else:
        config = DEFAULT_CONFIG.copy()
    
    return config


def save_config(config: Dict[str, Any], config_path: str) -> None:
    """
    Save configuration to JSON file.
    
    Args:
        config: configuration dict
        config_path: output path for config file
    """
    with open(config_path, 'w') as f:
        json.dump(config, f, indent=2)


def validate_config(config: Dict[str, Any]) -> bool:
    """
    Validate configuration parameters.
    
    Args:
        config: configuration dict
    
    Returns:
        True if valid, False otherwise
    """
    errors = []
    
    if config['conservation'] < 0 or config['conservation'] > 1:
        errors.append(f"conservation must be between 0 and 1, got {config['conservation']}")
    
    if config['inclusion'] < 0 or config['inclusion'] > 1:
        errors.append(f"inclusion must be between 0 and 1, got {config['inclusion']}")
    
    if config['gap_threshold'] < 0 or config['gap_threshold'] > 1:
        errors.append(f"gap_threshold must be between 0 and 1, got {config['gap_threshold']}")
    
    if config['length_disparity_threshold'] < 0 or config['length_disparity_threshold'] > 1:
        errors.append(f"length_disparity_threshold must be between 0 and 1")
    
    if config['quality_threshold'] < 0 or config['quality_threshold'] > 1:
        errors.append(f"quality_threshold must be between 0 and 1, got {config['quality_threshold']}")
    
    if config['threads'] is not None and config['threads'] < 1:
        errors.append(f"threads must be positive, got {config['threads']}")
    
    if errors:
        print("Configuration validation errors:")
        for error in errors:
            print(f"  - {error}")
        return False
    
    return True


def calculate_timeout(total_bases: int, config: Dict[str, Any]) -> int:
    """
    Calculate timeout based on total sequence size.
    
    Args:
        total_bases: total number of bases across all sequences
        config: configuration dict
    
    Returns:
        timeout in seconds
    """
    timeout = max(
        config['timeout_base'],
        min(
            config['max_timeout'],
            total_bases // config['timeout_per_base']
        )
    )
    return timeout