#!/usr/bin/env python3
import argparse
import sys
import os

sys.path.append("../../")
import pyBall.barrier_utils as bu

def main():
    parser = argparse.ArgumentParser(description="Cluster ensemble of DOFs using UMAP")
    parser.add_argument("--data", type=str, default="grid_search_out/ensemble_data.npz")
    parser.add_argument("--out_prefix", type=str, default="grid_search_out/ensemble_umap")
    parser.add_argument("--pn", action="store_true", help="Use FitREQ_PN instead of FitREQ")
    args = parser.parse_args()

    bu.init(use_pn=args.pn)
    bu.run_cluster_ensemble(args.data, args.out_prefix)

if __name__ == "__main__":
    main()
