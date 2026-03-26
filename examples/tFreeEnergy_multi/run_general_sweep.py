#!/usr/bin/env python3

import argparse
import json
import os
import subprocess
import shutil
import time

def run_cmd(cmd, cwd, log_path):
    print(f"Running: {' '.join(cmd)}")
    with open(log_path, "w") as f:
        p = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in p.stdout:
            f.write(line)
            # print(line, end="") # Optional: print to console
        return p.wait()

def main():
    parser = argparse.ArgumentParser(description="General Free Energy Sweep")
    parser.add_argument("--config", default="default_free_energy_config.json", help="JSON config file")
    args = parser.parse_args()

    with open(args.config, "r") as f:
        cfg = json.load(f)

    output_root = cfg.get("output_root", "results_sweep")
    common = cfg.get("common", {})
    param_sets = cfg.get("parameter_sets", [{}])

    os.makedirs(output_root, exist_ok=True)
    
    # Copy config to output for record
    shutil.copy2(args.config, os.path.join(output_root, "config.json"))

    summary_file = os.path.join(output_root, "summary.txt")
    with open(summary_file, "w") as sf:
        sf.write("Run Name | Status | Wall Time (s)\n")
        sf.write("-" * 40 + "\n")

    for pset in param_sets:
        # Merge common and specific parameters
        params = common.copy()
        params.update(pset)
        
        run_name = params.get("name", "unnamed_run")
        run_dir = os.path.join(output_root, run_name)
        os.makedirs(run_dir, exist_ok=True)
        
        print(f"\n>>> Starting run: {run_name}")
        
        cmd = ["python3", "run_ES.py"]
        # Add arguments from params
        for key, val in params.items():
            if key == "name": continue
            if isinstance(val, bool):
                if val: cmd.append(f"--{key}")
            elif key == "temperature":
                cmd.extend(["-T", str(val)])
            else:
                cmd.extend([f"--{key}", str(val)])

        log_path = os.path.join(run_dir, "run.log")
        t0 = time.time()
        rc = run_cmd(cmd, cwd=os.getcwd(), log_path=log_path)
        dt = time.time() - t0
        
        status = "OK" if rc == 0 else "FAILED"
        print(f"<<< Finished {run_name} with status {status} in {dt:.2f}s")
        
        with open(summary_file, "a") as sf:
            sf.write(f"{run_name} | {status} | {dt:.2f}\n")

        # Find the generated .dat file
        xyz_base = os.path.splitext(os.path.basename(params.get("xyz_name", "")))[0]
        dat_file = None
        for f in os.listdir("."):
            if f.endswith(".dat") and (xyz_base in f or "free_energy" in f):
                dat_file = f
                break

        if dat_file and rc == 0:
            print(f"Generating plot for {dat_file}...")
            plot_prefix = os.path.join(run_dir, f"{run_name}_plot")
            plot_cmd = [
                "python3", "plot_F_interactive.py",
                "--input", dat_file,
                "--output", plot_prefix,
                "--T", str(params.get("temperature", 300.0))
            ]
            # If we have N atoms information, we could pass it here as well
            run_cmd(plot_cmd, cwd=os.getcwd(), log_path=os.path.join(run_dir, "plot.log"))

        # Move all relevant files to run_dir
        for f in os.listdir("."):
            if f.endswith(".dat") or f.endswith(".html"):
                if xyz_base in f or "free_energy" in f or run_name in f:
                    target = os.path.join(run_dir, f)
                    if os.path.exists(target): os.remove(target)
                    shutil.move(f, target)

if __name__ == "__main__":
    main()
