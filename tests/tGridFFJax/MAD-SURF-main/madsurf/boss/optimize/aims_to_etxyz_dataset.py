#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created by edmanft
"""
import ase.io
import os
import argparse
# Set up argument parser
parser = argparse.ArgumentParser(description="Append trajectory data from a FHI-aims output file to a dataset file in extxyz format.")
parser.add_argument('aims_output_path', type=str, help='Path to the FHI-aims output file')
parser.add_argument('dataset_path', type=str, help='Path to the dataset file')
# Parse arguments
args = parser.parse_args()
aims_output_path = args.aims_output_path
dataset_path = args.dataset_path
# Check if the FHI-aims output file exists
if not os.path.exists(aims_output_path):
    raise FileNotFoundError(f"FHI-aims output file not found: {aims_output_path}")
# Check if the dataset file exists, if not create an empty file
if not os.path.exists(dataset_path):
    open(dataset_path, 'w').close()
# Read the FHI-aims output file
frames = ase.io.read(aims_output_path, index=':')
# Open the dataset file in append mode
with open(dataset_path, 'a') as dataset_file:
    for frame in frames:
        ase.io.write(dataset_file, frame, format='extxyz', append=True)