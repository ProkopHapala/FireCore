"""
FireCore Configuration Utility

Loads configuration from config.json and environment variables.
Environment variables override config.json settings.

Usage:
    from pyBall.config_utils import get_config, get_path
    config = get_config()
    sk_path = get_path('dftb_sk_path')
"""

import os
import json
from pathlib import Path

# Default config file location (relative to this file)
DEFAULT_CONFIG_PATH = Path(__file__).parent.parent / 'firecore_config.json'

def get_config(config_path=None):
    """
    Load FireCore configuration from config.json.
    
    Priority (highest to lowest):
    1. Environment variables
    2. config.json file
    3. Default values
    
    Args:
        config_path: Optional path to config.json file
    
    Returns:
        dict with merged configuration
    """
    if config_path is None:
        config_path = DEFAULT_CONFIG_PATH
    
    # Load config file if it exists
    config = {}
    if Path(config_path).exists():
        with open(config_path, 'r') as f:
            config = json.load(f)
    
    # Override with environment variables
    env_overrides = {
        'dftb_sk_path': os.environ.get('DFTB_SK_PATH'),
        'dftb_exe': os.environ.get('DFTB_EXE'),
        'dftb_basis_path': os.environ.get('DFTB_BASIS_PATH'),
        'firecore_root': os.environ.get('FIRECORE_ROOT'),
        'fdata_dir': os.environ.get('FDATA_DIR'),
    }
    
    # Apply environment variable overrides to paths section
    if 'paths' not in config:
        config['paths'] = {}
    
    for key, env_value in env_overrides.items():
        if env_value is not None:
            config['paths'][key] = env_value
    
    # Set defaults for missing values
    if 'firecore_root' not in config['paths'] or config['paths']['firecore_root'] is None:
        config['paths']['firecore_root'] = str(Path(__file__).parent.parent)
    
    if 'dftb_basis_path' not in config['paths'] or config['paths']['dftb_basis_path'] is None:
        config['paths']['dftb_basis_path'] = str(Path(config['paths']['firecore_root']) / 'pyBall' / 'DFTB' / 'data')
    
    if 'fdata_dir' not in config['paths'] or config['paths']['fdata_dir'] is None:
        config['paths']['fdata_dir'] = str(Path(config['paths']['firecore_root']) / 'tests' / 'pyFireball' / 'Fdata')
    
    return config

def get_path(key, config_path=None, default=None):
    """
    Get a specific path from configuration.
    
    Args:
        key: Path key (e.g., 'dftb_sk_path', 'dftb_basis_path')
        config_path: Optional path to config.json file
        default: Optional default value if not found
    
    Returns:
        str path or None if not found
    """
    config = get_config(config_path)
    return config.get('paths', {}).get(key, default)

def get_dftb_basis_path(basis_name, config_path=None):
    """
    Get the path to a basis HSD file for a specific basis set.
    
    Args:
        basis_name: Name of basis set (e.g., 'mio-1-1', '3ob-3-1')
        config_path: Optional path to config.json file
    
    Returns:
        str path to basis HSD file or None if not found
    """
    config = get_config(config_path)
    
    # Check dftb.basis_sets section
    if 'dftb' in config and 'basis_sets' in config['dftb']:
        basis_info = config['dftb']['basis_sets'].get(basis_name)
        if basis_info and 'wfc_path' in basis_info:
            return basis_info['wfc_path']
    
    # Fall back to constructing from dftb_basis_path
    basis_path = get_path('dftb_basis_path', config_path)
    if basis_path:
        return str(Path(basis_path) / f'wfc.{basis_name}.hsd')
    
    return None

def get_dftb_sk_path(basis_name, config_path=None):
    """
    Get the path to Slater-Koster files for a specific basis set.
    
    Args:
        basis_name: Name of basis set (e.g., 'mio-1-1', '3ob-3-1')
        config_path: Optional path to config.json file
    
    Returns:
        str path to SK directory or None if not found
    """
    config = get_config(config_path)
    
    # Check dftb.basis_sets section
    if 'dftb' in config and 'basis_sets' in config['dftb']:
        basis_info = config['dftb']['basis_sets'].get(basis_name)
        if basis_info and 'sk_path' in basis_info:
            return basis_info['sk_path']
    
    # Fall back to constructing from dftb_sk_path
    sk_path = get_path('dftb_sk_path', config_path)
    if sk_path:
        # Try common subdirectories
        for subdir in [basis_name, f'mio/{basis_name}', f'library/{basis_name}']:
            candidate = Path(sk_path) / subdir
            if candidate.exists():
                return str(candidate)
        # Return the base sk_path and let caller handle subdirectory
        return sk_path
    
    return None

def print_config(config_path=None):
    """
    Print current configuration (useful for debugging).
    
    Args:
        config_path: Optional path to config.json file
    """
    config = get_config(config_path)
    print("[FireCore Configuration]")
    print(f"  Config file: {config_path or DEFAULT_CONFIG_PATH}")
    print(f"  Config exists: {Path(config_path or DEFAULT_CONFIG_PATH).exists()}")
    print("\n[Paths]")
    for key, value in config.get('paths', {}).items():
        print(f"  {key}: {value}")
    if 'dftb' in config:
        print("\n[DFTB]")
        if 'default_sk_set' in config['dftb']:
            print(f"  default_sk_set: {config['dftb']['default_sk_set']}")
        if 'available_sk_sets' in config['dftb']:
            print(f"  available_sk_sets: {config['dftb']['available_sk_sets']}")
