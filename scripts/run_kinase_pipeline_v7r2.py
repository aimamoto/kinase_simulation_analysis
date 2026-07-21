#!/usr/bin/env python3
import subprocess
import sys
import os
import re
import argparse
import shutil
import csv as pycsv
from datetime import datetime

# --- Path Configuration ---
MODULE_DIR = "modules"
ARCHIVE_DIR = "archives"
OUTPUT_DIRS = ["cx_viz_core", "cx_viz_allosteric"]

# --- Script Configuration ---
SCRIPT_DISCOVERY = "generate_config.py"
SCRIPT_FASTA     = "extract_fasta.py"
SCRIPT_LANDMARKS = "extract_landmarks.py"
SCRIPT_CHIMERAX  = "1_run_parallel_chimerax_hmm_v7r2.py"
SCRIPT_VIS       = "kinome_VISalign.py"
SCRIPT_AF3_METRICS = "extract_af3_metrics.py"

DEFAULT_OUTPUT_FILES = [
    "generated_matrix.csv", "proteins.yaml", "sequences.fasta",
    "hmm_landmarks.json", "hmm_kinase_analysis_results_v7r2.csv",
    "master_kinase_analysis_results_v7r2.csv"
]

def archive_old_results(files_to_archive):
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
        
def merge_csv_results(geom_csv, af3_csv, final_out_csv):
    if not os.path.exists(geom_csv):
        print(f"⚠️  [MERGE] Could not find geometry file: {geom_csv}")
        return
    if not os.path.exists(af3_csv):
        print(f"⚠️  [MERGE] Could not find AF3 file: {af3_csv}")
        return
        
    print("\n>>> [STEP] Merging geometric and AF3 confidence datasets...")
    
    # Read dynamically generated AF3 fields
    with open(af3_csv, 'r') as f:
        af3_reader = pycsv.DictReader(f)
        af3_fields = [fld for fld in af3_reader.fieldnames if fld != 'Directory']
        af3_data = {os.path.normpath(row['Directory']): row for row in af3_reader}
    
    with open(geom_csv, 'r') as f:
        geom_reader = pycsv.DictReader(f)
        geom_fields = geom_reader.fieldnames
        
        with open(final_out_csv, 'w', newline='') as out_f:
            writer = pycsv.DictWriter(out_f, fieldnames=geom_fields + af3_fields)
            writer.writeheader()
            
            for row in geom_reader:
                d = os.path.normpath(row['Directory'])
                if d in af3_data:
                    for fld in af3_fields:
                        row[fld] = af3_data[d].get(fld, "N/A")
                else:
                    for fld in af3_fields:
                        row[fld] = "N/A"
                writer.writerow(row)
                
    if os.path.exists(af3_csv): os.remove(af3_csv)
    if os.path.exists(geom_csv): os.remove(geom_csv)
    print(f"✅ Final comprehensive dataset compiled: {final_out_csv}")

def main():
    parser = argparse.ArgumentParser(description="Kinase Structural Pipeline Wrapper v7")
    parser.add_argument("--resume", action="store_true", help="Skip discovery and resume from CSV audit")
    parser.add_argument("--use-yaml", action="store_true", help="Start from sequence extraction using existing proteins.yaml")
    parser.add_argument("--use-fasta", action="store_true", help="Start from HMM alignment using existing sequences.fasta")
    parser.add_argument("--no-archive", action="store_true", help="Skip archiving of old results entirely")
    parser.add_argument("-c", "--cores", type=int, default=8, help="Number of CPU cores to use (default: 8)")
    parser.add_argument("-n", "--max-chains", type=int, default=None, help="Total number of chains per structure (Kinases + Co-Factors)")
    args = parser.parse_args()

    print("======================================================")
    print("   Kinase Structural Bioinformatic Pipeline v7r2")
    print("======================================================")

    files_to_archive = list(DEFAULT_OUTPUT_FILES)
    
    if args.resume and "generated_matrix.csv" in files_to_archive: files_to_archive.remove("generated_matrix.csv")
    if args.use_yaml and "proteins.yaml" in files_to_archive: files_to_archive.remove("proteins.yaml")
    if args.use_fasta:
        for f in ["sequences.fasta", "proteins.yaml", "generated_matrix.csv"]:
            if f in files_to_archive: files_to_archive.remove(f)

    if not args.no_archive:
        archive_old_results(files_to_archive)

    run_config = not (args.use_yaml or args.use_fasta)
    run_fasta = not args.use_fasta

    # --- STOICHIOMETRY PROMPT & SAFETY CATCH ---
    # Auto-detect maximum chains based on nested directory structures (a-, b-, c-, etc.)
    auto_detected_chains = 2
    for d in os.listdir("."):
        if os.path.isdir(d) and re.search(r'(?:^|_)[a-z]-', d, re.IGNORECASE):
            chain_count = len(re.findall(r'(?:^|_)[a-z]-', d, re.IGNORECASE))
            if chain_count > auto_detected_chains:
                auto_detected_chains = chain_count

    if run_fasta:
        if args.max_chains is None:
            print("\n[?] STOICHIOMETRY CONFIGURATION:")
            print(f"    Auto-detected maximum chains per simulation: {auto_detected_chains}")
            print("    How many total chains (Kinases + Target Co-factors) should be extracted per simulation?")
            print(f"    (e.g., enter '{auto_detected_chains}' based on your directory structures)")
            while True:
                ans = input(f"    Number of chains [default: {auto_detected_chains}]: ").strip()
                if not ans:
                    args.max_chains = auto_detected_chains
                    break
                try:
                    args.max_chains = int(ans)
                    if args.max_chains > 0:
                        break
                    print("    [!] Please enter a positive integer.")
                except ValueError:
                    print("    [!] Invalid input. Please enter a valid number.")
    else:
        # Fallback if running with --use-fasta but forgot to specify -n
        if args.max_chains is None:
            args.max_chains = auto_detected_chains

    if run_config:
        if not args.resume:
            run_step(SCRIPT_DISCOVERY, [], "Scanning directories for simulations")
            if not os.path.exists("generated_matrix.csv"):
                print("\n❌ CRITICAL: 'generated_matrix.csv' was not created.")
                sys.exit(1)
            
            print("\n[?] Discovery Complete:\n    1. [Ready Mode] Proceed to full analysis.\n    2. [Audit Mode] Stop to edit 'generated_matrix.csv'.")
            choice = input("\nSelect an option (1 or 2): ").strip()
            if choice == "2":
                print("\n[*] PAUSED: Review 'generated_matrix.csv'.\n[*] When ready, run: python3 run_kinase_pipeline_v7r2.py --resume")
                sys.exit(0)
        else:
            print("[*] Resuming pipeline from existing 'generated_matrix.csv'...")
            
        run_step(SCRIPT_DISCOVERY, ["-m", "generated_matrix.csv"], "Generating proteins.yaml")

    if run_fasta:
        print("[*] Using 'proteins.yaml' for sequence extraction...")
        run_step(SCRIPT_FASTA, [
            "--out", "sequences.fasta", 
            "--cores", str(args.cores),
            "--max-chains", str(args.max_chains)
        ], f"Parallel extraction of sequences ({args.cores} cores, max {args.max_chains} chains)")

    if args.use_fasta:
        print("[*] Skipping extraction. Using user-provided 'sequences.fasta'...")

    run_step(SCRIPT_LANDMARKS, ["--fasta", "sequences.fasta"], "HMM Aligning and Landmark extraction")
    run_step(SCRIPT_VIS, ["-i", "sequences.fasta", "-l", "hmm_landmarks.json"], "Generating 1D structural schematics")
    run_step(SCRIPT_CHIMERAX, ["-c", str(args.cores)], f"Parallel ChimeraX analysis ({args.cores} workers)")
    
    # AF3 Metrics
    run_step(SCRIPT_AF3_METRICS, ["--dir", ".", "--max-chains", str(args.max_chains), "--out", "temp_af3_metrics.csv"], "Extracting AF3 ipTM and PAE metrics")

    # Final Merge Execution
    merge_csv_results(
        geom_csv="hmm_kinase_analysis_results_v7r2.csv", 
        af3_csv="temp_af3_metrics.csv", 
        final_out_csv="master_kinase_analysis_results_v7r2.csv"
    )

    print("\n======================================================")
    print("   Pipeline Completed!")
    print("   Output: master_kinase_analysis_results_v7r2.csv")
    print("======================================================")

if __name__ == "__main__":
    main()
