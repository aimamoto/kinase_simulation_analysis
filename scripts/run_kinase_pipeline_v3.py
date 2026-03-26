#!/usr/bin/env python3
import subprocess
import sys
import os
import argparse
import shutil
from datetime import datetime

# Path Configuration
MODULE_DIR = "modules"
ARCHIVE_DIR = "archives"
OUTPUT_DIRS = ["cx_pocket_visualization_hmm"]

DEFAULT_OUTPUT_FILES =[
    "generated_matrix.csv", "proteins.yaml", "sequences.fasta", 
    "hmm_landmarks.json", "hmm_aloop_analysis_results.csv"
]

def archive_old_results(files_to_archive):
    """Archives old results, respecting files the user explicitly wants to use as inputs."""
    if not any(os.path.exists(f) for f in files_to_archive + OUTPUT_DIRS):
        return
        
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_archive_path = os.path.join(ARCHIVE_DIR, f"run_{timestamp}")
    os.makedirs(run_archive_path, exist_ok=True)
    
    print(f"\n>>> [ARCHIVE] Moving previous results to {run_archive_path}...")
    for item in files_to_archive + OUTPUT_DIRS:
        if os.path.exists(item):
            shutil.move(item, os.path.join(run_archive_path, item))

def run_step(script_name, args, description):
    script_path = os.path.join(MODULE_DIR, script_name)
    print(f"\n>>> [STEP] {description}...")
    if not os.path.exists(script_path):
        print(f"ERROR: {script_name} not found in {MODULE_DIR}/")
        sys.exit(1)
        
    cmd = ["python3", script_path] + args
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: {description} failed. Exit code: {e.returncode}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Kinase Structural Pipeline Wrapper")
    parser.add_argument("--resume", action="store_true", help="Skip discovery and resume from CSV audit")
    parser.add_argument("--use-yaml", action="store_true", help="Start from sequence extraction using existing proteins.yaml")
    parser.add_argument("--use-fasta", action="store_true", help="Start from HMM alignment using existing sequences.fasta")
    parser.add_argument("--no-archive", action="store_true", help="Skip archiving of old results entirely")
    parser.add_argument("-c", "--cores", type=int, default=8, help="Number of CPU cores to use (default: 8)")
    args = parser.parse_args()

    print("======================================================")
    print("   Kinase Structural Bioinformatic Pipeline")
    print("======================================================")

    files_to_archive = list(DEFAULT_OUTPUT_FILES)
    
    # Don't archive the matrix if we are actively trying to resume from it
    if args.resume and "generated_matrix.csv" in files_to_archive:
        files_to_archive.remove("generated_matrix.csv")
        
    if args.use_yaml and "proteins.yaml" in files_to_archive:
        files_to_archive.remove("proteins.yaml")
    if args.use_fasta:
        if "sequences.fasta" in files_to_archive: files_to_archive.remove("sequences.fasta")
        if "proteins.yaml" in files_to_archive: files_to_archive.remove("proteins.yaml")
        if "generated_matrix.csv" in files_to_archive: files_to_archive.remove("generated_matrix.csv")

    if not args.no_archive:
        archive_old_results(files_to_archive)

    run_config = not (args.use_yaml or args.use_fasta)
    run_fasta = not args.use_fasta

    if run_config:
        if not args.resume:
            run_step("generate_config.py", [], "Scanning directories for simulations")
            print("\n[?] Discovery Complete:\n    1. [Ready Mode] Proceed to full analysis.\n    2. [Audit Mode] Stop to edit 'generated_matrix.csv'.")
            choice = input("\nSelect an option (1 or 2): ").strip()
            if choice == "2":
                print("\n[*] PAUSED: Review 'generated_matrix.csv'.\n[*] When ready, run: python3 run_kinase_pipeline.py --resume")
                sys.exit(0)
        else:
            print("[*] Resuming pipeline from existing 'generated_matrix.csv'...")
            
        run_step("generate_config.py",["-m", "generated_matrix.csv"], "Generating proteins.yaml")

    if run_fasta:
        print("[*] Using 'proteins.yaml' for sequence extraction...")
        run_step("extract_fasta.py",["--out", "sequences.fasta"], "Extracting sequences to FASTA")

    if args.use_fasta:
        print("[*] Skipping extraction. Using user-provided 'sequences.fasta'...")

    run_step("extract_landmarks.py", ["--fasta", "sequences.fasta"], "HMM Aligning and Landmark extraction")
    
    # PASS CORES FLAG DYNAMICALLY HERE
    run_step("1_run_parallel_chimerax_hmm_v3.py", ["-c", str(args.cores)], f"Parallel ChimeraX analysis ({args.cores} workers)")

    print("\n======================================================\n   Pipeline Completed!\n======================================================")

if __name__ == "__main__":
    main()
