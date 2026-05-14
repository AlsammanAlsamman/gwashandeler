#!/usr/bin/env python3
"""
bioconfigme.py - Configuration loader for GWAS handler pipeline
Provides utilities to read and access analysis.yml and software.yml configuration files.
"""

import os
import yaml
from pathlib import Path


def _load_yaml(filepath):
    """Load a YAML file and return as dict."""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Config file not found: {filepath}")
    
    with open(filepath, 'r') as f:
        data = yaml.safe_load(f)
    
    return data if data is not None else {}


def get_results_dir():
    """
    Read results_dir from configs/analysis.yml.
    
    Returns:
        str: Absolute path to results directory
    """
    config_path = "configs/analysis.yml"
    config = _load_yaml(config_path)
    results_dir = config.get('results_dir')
    
    if not results_dir:
        raise ValueError("results_dir not found in configs/analysis.yml")
    
    return results_dir


def get_software_module(tool):
    """
    Get the module name for a given tool from configs/software.yml.
    
    Args:
        tool (str): Tool name (e.g., 'plink2', 'r', 'python')
    
    Returns:
        str: Module name (e.g., 'plink2/2.00a3.3lm')
    """
    config_path = "configs/software.yml"
    config = _load_yaml(config_path)
    
    if tool in config and 'module' in config[tool]:
        return config[tool]['module']
    
    raise ValueError(f"Tool '{tool}' or module not found in configs/software.yml")


def get_analysis_value(path_str):
    """
    Get a value from configs/analysis.yml using dot notation.
    
    Args:
        path_str (str): Path to value using dot notation (e.g., 'default_resources.mem_mb')
    
    Returns:
        Any: The value at the specified path
    """
    config_path = "configs/analysis.yml"
    config = _load_yaml(config_path)
    
    keys = path_str.split('.')
    value = config
    
    for key in keys:
        if isinstance(value, dict):
            value = value.get(key)
        else:
            raise ValueError(f"Cannot access key '{key}' in path '{path_str}'")
    
    return value


def get_gwastables():
    """
    Get all GWAS table definitions from configs/analysis.yml.
    
    Returns:
        list: List of gwastables entries
    """
    config_path = "configs/analysis.yml"
    config = _load_yaml(config_path)
    
    return config.get('gwastables', [])


def get_ref_panels():
    """
    Get reference panel paths from configs/analysis.yml.
    
    Returns:
        dict: Dictionary of reference panel names to paths
    """
    config_path = "configs/analysis.yml"
    config = _load_yaml(config_path)
    
    return config.get('ref_panels', {})


if __name__ == "__main__":
    # Test functions
    try:
        print("Results directory:", get_results_dir())
        print("GWAS tables:", len(get_gwastables()), "datasets")
        print("Reference panels:", list(get_ref_panels().keys()))
    except Exception as e:
        print(f"Error: {e}")
