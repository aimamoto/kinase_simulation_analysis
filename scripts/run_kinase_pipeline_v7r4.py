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
SCRIPT_CHIMERAX  = "1_run_parallel_chimerax_hmm_v7r4.py"
SCRIPT_VIS       = "kinome_VISalign.py"
SCRIPT_AF3_METRICS = "extract_af3_metrics.py"

DEFAULT_OUTPUT_FILES = [
    "generated_matrix.csv", "proteins.yaml", "sequences.fasta",
    "hmm_landmarks.json", "hmm_kinase_analysis_results_v7r4.csv",
    "master_kinase_analysis_results_v7r4.csv"
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
    # v7r4 (BLOCKERS.md B12.6 item 1): these were warn-and-return in v7r3, so a run that
    # produced no geometry still printed "Pipeline Completed!" and exited 0 with no master
    # CSV written at all. Same failure class as B12.5. Stop instead.
    if not os.path.exists(geom_csv):
        print(f"ERROR: [MERGE] Could not find geometry file: {geom_csv}", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(af3_csv):
        print(f"ERROR: [MERGE] Could not find AF3 file: {af3_csv}", file=sys.stderr)
        sys.exit(1)
        
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
            
            n_geom = 0
            for row in geom_reader:
                d = os.path.normpath(row['Directory'])
                if d in af3_data:
                    for fld in af3_fields:
                        row[fld] = af3_data[d].get(fld, "N/A")
                else:
                    for fld in af3_fields:
                        row[fld] = "N/A"
                writer.writerow(row)
                n_geom += 1

    # v7r4: this merge is a LEFT join FROM the geometry CSV, so it is 1:1 by construction --
    # which is exactly why it propagated B12.5's gap into the master without comment. Coverage
    # is checked upstream in 1_run_parallel_chimerax_hmm_v7r4.py; this only confirms the join
    # itself lost nothing.
    with open(final_out_csv, 'r') as f:
        n_master = sum(1 for _ in pycsv.DictReader(f))
    if n_master != n_geom:
        print(f"ERROR: [MERGE] master has {n_master} rows but geometry had {n_geom}.",
              file=sys.stderr)
        sys.exit(1)
    print(f"[*] Merge check: {n_master} rows carried through from the geometry CSV.")

    if os.path.exists(af3_csv): os.remove(af3_csv)
    if os.path.exists(geom_csv): os.remove(geom_csv)
    print(f"✅ Final comprehensive dataset compiled: {final_out_csv}")

def find_af3_confidence_files(root="."):
    # AF3 confidence side files under any naming: local runs write summary_confidences.json /
    # <job>_summary_confidences.json / confidences.json, the AF3 server writes
    # fold_<job>_summary_confidences_0.json and fold_<job>_full_data_0.json.
    found = []
    for dirpath, dirs, files in os.walk(root):
        dirs[:] = [d for d in dirs if not d.startswith('.') and not d.lower().startswith('archive')
                   and d not in (MODULE_DIR, "temp_chimerax_chunks", "temp_fasta_chunks")]
        for f in files:
            fl = f.lower()
            if fl.endswith(".json") and ("confidences" in fl or "full_data" in fl):
                found.append(os.path.join(dirpath, f))
    return sorted(found)

def resolve_missing_af3(af3_csv, max_chains):
    """v7r4 (2026-10-01). Tell non-AF3 input apart from AF3 output under other names.

    extract_af3_metrics.py writes nothing when it finds no file literally named `model.cif`
    (B12.1), and since B12.6 item 1 the merge stops on a missing AF3 file. Together they failed
    every MD-trajectory and experimental-structure run at the last step, with the geometry
    complete (seen on 3D7T and the AURKA-TPX2 MD run; v7r3 ended at the geometry CSV instead).
    With no AF3 confidence files anywhere, the input is not AF3: write a header-only AF3 table so
    the merge fills ipTM/pTM/PAE with N/A and the run completes. With confidence files present
    under other names, the AF3 metrics really are being lost, which is the case the hard stop is
    for: stop and say how to fix it.
    """
    if os.path.exists(af3_csv): return
    conf = find_af3_confidence_files(".")
    if conf:
        print("\nERROR: AF3 confidence files found but no 'model.cif' -- non-standard AF3 file names (B12.1).",
              file=sys.stderr)
        for p in conf[:5]:
            print(f"         {p}", file=sys.stderr)
        if len(conf) > 5:
            print(f"         ... and {len(conf) - 5} more", file=sys.stderr)
        print("       ipTM/pTM/PAE would be lost. Rename each model to model.cif with summary_confidences.json and",
              file=sys.stderr)
        print("       confidences.json beside it, or remove the JSON files if AF3 metrics are not wanted.",
              file=sys.stderr)
        print("       The geometry CSV hmm_kinase_analysis_results_v7r4.csv is complete and has been kept.",
              file=sys.stderr)
        sys.exit(1)

    # Same column names as extract_af3_metrics.py builds for this chain count.
    labels = ['A', 'B', 'C', 'D'][:max_chains]
    fields = ["Directory", "ipTM", "pTM"]
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            fields += [f"PAE_{labels[i]}_to_{labels[j]}", f"PAE_{labels[j]}_to_{labels[i]}",
                       f"PAE_Mean_{labels[i]}{labels[j]}"]
    with open(af3_csv, 'w', newline='') as f:
        pycsv.writer(f).writerow(fields)
    print("[i] No AF3 models ('model.cif') found and no AF3 confidence files present -- treating input as")
    print("    non-AF3 (MD / experimental). ipTM, pTM and PAE columns will be N/A.")

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
    print("   Kinase Structural Bioinformatic Pipeline v7r4")
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
                print("\n[*] PAUSED: Review 'generated_matrix.csv'.\n[*] When ready, run: python3 run_kinase_pipeline_v7r4.py --resume")
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
    resolve_missing_af3("temp_af3_metrics.csv", args.max_chains)

    # Final Merge Execution
    merge_csv_results(
        geom_csv="hmm_kinase_analysis_results_v7r4.csv", 
        af3_csv="temp_af3_metrics.csv", 
        final_out_csv="master_kinase_analysis_results_v7r4.csv"
    )

    print("\n======================================================")
    print("   Pipeline Completed!")
    print("   Output: master_kinase_analysis_results_v7r4.csv")
    print("======================================================")

if __name__ == "__main__":
    main()
