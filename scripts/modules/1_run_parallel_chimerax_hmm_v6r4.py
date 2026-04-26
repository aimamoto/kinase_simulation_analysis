#!/usr/bin/env python3
"""
Master script to run ChimeraX structural analysis in parallel.
Usage: python3 1_run_parallel_chimerax_hmm_v6r3.py -c 24
"""

import os
import glob
import subprocess
import sys
import multiprocessing
import argparse
import math
import shutil  
import pandas as pd

# --- CONFIGURATION ---
WORKER_SCRIPT = os.path.join("modules", "chimerax_hmm_worker_v6r4.py")
FINAL_CSV_NAME = "hmm_kinase_analysis_results_v6.csv"
CHUNK_DIR = "temp_chimerax_chunks"

def ensure_hmm_landmarks():
    """Checks for the landmarks JSON and triggers the 0_ prep script if missing."""
    if os.path.exists("hmm_landmarks.json"):
        print("Found hmm_landmarks.json. Proceeding to structural analysis...")
        return

    print("hmm_landmarks.json not found. Launching 0_extract_hmm_landmarks.py...")
    try:
        subprocess.run(["python3", "0_extract_hmm_landmarks.py"], check=True)
    except subprocess.CalledProcessError:
        print("ERROR: Landmark generation failed. Halting pipeline.")
        sys.exit(1)
    except FileNotFoundError:
        print("ERROR: 0_extract_hmm_landmarks.py not found in this directory.")
        sys.exit(1)

def filter_cif_files(cif_files):
    """
    Smart-filters CIF files. If a directory contains AF3 nested seeds 
    (e.g., 'seed-m_sample-n/model.cif'), we ignore the redundant top-level 
    '*_model.cif' in that same directory. 
    For flat directories without nested seeds, we keep everything.
    """
    dirs_with_nested_seeds = set()
    
    # 1. Identify directories containing nested AF3 models
    for cif in cif_files:
        cif_norm = cif.replace("\\", "/") # Windows/Linux compatibility
        parts = cif_norm.split("/")
        
        if len(parts) >= 2:
            filename = parts[-1]
            parent = parts[-2]
            
            # Check if it looks like an AF3 nested model
            if filename == "model.cif" and ("seed" in parent.lower() or "sample" in parent.lower()):
                # The grandparent is the main job directory
                grandparent = "/".join(parts[:-2]) if len(parts) > 2 else "."
                dirs_with_nested_seeds.add(grandparent)
                
    # 2. Filter out redundant top-ranked models
    filtered_files =[]
    for cif in cif_files:
        cif_norm = cif.replace("\\", "/")
        parts = cif_norm.split("/")
        filename = parts[-1]
        parent_dir = "/".join(parts[:-1]) if len(parts) > 1 else "."
        
        # If this is a summary *_model.cif AND its parent directory has nested seeds, drop it
        is_redundant = False
        if filename.endswith("_model.cif") and parent_dir in dirs_with_nested_seeds:
            is_redundant = True
            
        if is_redundant:
            print(f"[*] Orchestrator: Excluded redundant AF3 top-ranked model -> {cif}")
        else:
            filtered_files.append(cif)
            
    return filtered_files

def create_chunks(cif_files, num_cores):
    """Splits the list of CIF files into smaller text files for each worker."""
    if not os.path.exists(CHUNK_DIR):
        os.makedirs(CHUNK_DIR)
        
    print(f"[*] Orchestrator: Partitioning {len(cif_files)} CIF files across {num_cores} workers...")
    chunk_size = math.ceil(len(cif_files) / num_cores)
    chunk_files =[]
    
    for i in range(num_cores):
        start_idx = i * chunk_size
        end_idx = start_idx + chunk_size
        chunk = cif_files[start_idx:end_idx]
        
        if not chunk:
            continue
            
        chunk_path = os.path.join(CHUNK_DIR, f"chunk_{i}.txt")
        with open(chunk_path, "w") as f:
            for cif in chunk:
                f.write(f"{cif}\n")
        chunk_files.append(chunk_path)
        
    print(f"[*] Orchestrator: Successfully generated {len(chunk_files)} chunk files in '{CHUNK_DIR}/'.")
    return chunk_files

def run_chimerax_chunk(chunk_file):
    """Spawns a headless ChimeraX instance to process a specific chunk."""
    chunk_name = os.path.basename(chunk_file)
    print(f"  -> [Worker START] Booting headless ChimeraX for {chunk_name}...")
    
    env = os.environ.copy()
    env["CHIMERAX_CHUNK"] = chunk_file

    cmd =[
        "chimerax", 
        "--nogui", 
        WORKER_SCRIPT
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, env=env)
        
        if result.returncode != 0:
            print(f"  [!] [Worker ERROR] ChimeraX failed on {chunk_name}")
            print(f"Details:\n{result.stderr}\n{result.stdout}")
            return False
            
        print(f"  <- [Worker DONE] Successfully processed {chunk_name}.")
        return True
        
    except FileNotFoundError:
        print("  [!] [Worker ERROR] Could not find 'chimerax' in PATH.")
        return False
    except Exception as e:
        print(f"  [!] [Worker ERROR] Unexpected failure on {chunk_name}. Reason: {e}")
        return False

def merge_csvs():
    """Combines all the individual worker CSVs into one final dataset and sorts by Simulation_ID."""
    print("[*] Orchestrator: Scanning for worker CSV outputs...")
    
    # Check both the main folder AND the temp folder specifically for the v6 suffix
    csv_files = glob.glob("*_results_v6.csv") + glob.glob(os.path.join(CHUNK_DIR, "*_results_v6.csv"))
    
    if FINAL_CSV_NAME in csv_files:
        csv_files.remove(FINAL_CSV_NAME)
        
    if not csv_files:
        print("[!] Orchestrator WARNING: No worker output CSVs found to merge.")
        return

    print(f"[*] Orchestrator: Merging {len(csv_files)} worker outputs into '{FINAL_CSV_NAME}'...")
    
    df_list =[]
    for csv_file in csv_files:
        try:
            df = pd.read_csv(csv_file)
            df_list.append(df)
        except pd.errors.EmptyDataError:
            print(f"  [!] Skipping empty file: {csv_file}")
            
    if not df_list:
        print("[!] Orchestrator WARNING: No valid data found in worker CSVs.")
        return

    # Combine all worker dataframes
    final_df = pd.concat(df_list, ignore_index=True)
    
    # Sort by Simulation_ID
    if "Simulation_ID" in final_df.columns:
        print("[*] Orchestrator: Sorting dataset by 'Simulation_ID'...")
        final_df = final_df.sort_values(by="Simulation_ID")
    else:
        print("  [!] Warning: 'Simulation_ID' column not found. Skipping sorting step.")
        
    # Export to final CSV
    final_df.to_csv(FINAL_CSV_NAME, index=False)
    print(f"[*] Orchestrator: Final dataset successfully compiled and saved to '{FINAL_CSV_NAME}'.")
    
    # Cleanup individual worker CSVs
    print("[*] Orchestrator: Removing intermediate worker CSV files...")
    for csv_file in csv_files:
        os.remove(csv_file)

def main():
    parser = argparse.ArgumentParser(description="Run ChimeraX HMM Pipeline in Parallel")
    parser.add_argument("-c", "--cores", type=int, default=8, help="Number of CPU cores to use (default: 8)")
    args = parser.parse_args()
    
    system_cores = multiprocessing.cpu_count()
    if args.cores > system_cores:
        print(f"[*] Orchestrator WARNING: Requested {args.cores} cores, but system only has {system_cores}.")
        print(f"[*] Orchestrator: Capping at {system_cores} cores to prevent system crash.")
        args.cores = system_cores

    print("\n=======================================================")
    print("   Starting Parallel ChimeraX Structural Pipeline")
    print("=======================================================\n")

    # 1. Guarantee the HMM landmarks exist before any workers wake up
    ensure_hmm_landmarks()

    # 2. Find the AlphaFold3 outputs recursively
    cif_files = glob.glob("**/*.cif", recursive=True)
    if not cif_files:
        print("ERROR: No .cif files found in the current directory or subdirectories.")
        sys.exit(1)
        
    print(f"[*] Orchestrator: Discovered {len(cif_files)} total CIF files.")
    
    cif_files = filter_cif_files(cif_files)
    print(f"[*] Orchestrator: {len(cif_files)} CIF files remaining after AF3 redundancy filter.")

    # 3. Split the work
    chunk_files = create_chunks(cif_files, args.cores)

    # 4. Run the workers in parallel
    print(f"[*] Orchestrator: Dispatching tasks to multiprocessing pool ({args.cores} cores)...")
    with multiprocessing.Pool(processes=args.cores) as pool:
        pool.map(run_chimerax_chunk, chunk_files)
    print("[*] Orchestrator: All workers have reported back.")

    # 5. Clean up and merge
    merge_csvs()
    
    # 6. Cleanup temp folder
    print(f"[*] Orchestrator: Cleaning up temporary directory '{CHUNK_DIR}/'...")
    if os.path.exists(CHUNK_DIR):
        shutil.rmtree(CHUNK_DIR)
        
    print("\n=======================================================")
    print("   ChimeraX Execution Fully Completed!")
    print("=======================================================\n")

if __name__ == "__main__":
    main()
