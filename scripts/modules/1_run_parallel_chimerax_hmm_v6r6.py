#!/usr/bin/env python3
"""
Master script to run ChimeraX structural analysis in parallel.
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
WORKER_SCRIPT = os.path.join("modules", "chimerax_hmm_worker_v6r6.py")
FINAL_CSV_NAME = "hmm_kinase_analysis_results_v6r6.csv"
CHUNK_DIR = "temp_chimerax_chunks"

def ensure_hmm_landmarks():
    if os.path.exists("hmm_landmarks.json"): return
    try: subprocess.run(["python3", "0_extract_hmm_landmarks.py"], check=True)
    except Exception: sys.exit(1)

def filter_cif_files(cif_files):
    dirs_with_nested_seeds = set()
    for cif in cif_files:
        cif_norm = cif.replace("\\", "/")
        parts = cif_norm.split("/")
        if len(parts) >= 2 and parts[-1] == "model.cif" and ("seed" in parts[-2].lower() or "sample" in parts[-2].lower()):
            dirs_with_nested_seeds.add("/".join(parts[:-2]) if len(parts) > 2 else ".")
                
    filtered_files = []
    for cif in cif_files:
        cif_norm = cif.replace("\\", "/")
        parts = cif_norm.split("/")
        parent_dir = "/".join(parts[:-1]) if len(parts) > 1 else "."
        
        is_redundant = parts[-1].endswith("_model.cif") and parent_dir in dirs_with_nested_seeds
        if not is_redundant: filtered_files.append(cif)
            
    return filtered_files

def create_chunks(cif_files, num_cores):
    if not os.path.exists(CHUNK_DIR): os.makedirs(CHUNK_DIR)
    chunk_size = math.ceil(len(cif_files) / num_cores)
    chunk_files = []
    
    for i in range(num_cores):
        chunk = cif_files[i * chunk_size : (i + 1) * chunk_size]
        if not chunk: continue
        chunk_path = os.path.join(CHUNK_DIR, f"chunk_{i}.txt")
        with open(chunk_path, "w") as f:
            for cif in chunk: f.write(f"{cif}\n")
        chunk_files.append(chunk_path)
    return chunk_files

def run_chimerax_chunk(chunk_file):
    env = os.environ.copy()
    env["CHIMERAX_CHUNK"] = chunk_file
    cmd = ["chimerax", "--nogui", WORKER_SCRIPT]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, env=env)
        if result.returncode != 0: return False
        return True
    except Exception: return False

def merge_csvs():
    csv_files = glob.glob("*_results_v6r6.csv") + glob.glob(os.path.join(CHUNK_DIR, "*_results_v6r6.csv"))
    if FINAL_CSV_NAME in csv_files: csv_files.remove(FINAL_CSV_NAME)
    if not csv_files: return

    df_list = []
    for csv_file in csv_files:
        try: df_list.append(pd.read_csv(csv_file))
        except pd.errors.EmptyDataError: pass
            
    if not df_list: return

    final_df = pd.concat(df_list, ignore_index=True)
    if "Simulation_ID" in final_df.columns: final_df = final_df.sort_values(by="Simulation_ID")
        
    final_df.to_csv(FINAL_CSV_NAME, index=False)
    for csv_file in csv_files: os.remove(csv_file)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--cores", type=int, default=8)
    args = parser.parse_args()
    
    system_cores = multiprocessing.cpu_count()
    if args.cores > system_cores: args.cores = system_cores

    ensure_hmm_landmarks()
    cif_files = filter_cif_files(glob.glob("**/*.cif", recursive=True))
    if not cif_files: sys.exit(1)

    chunk_files = create_chunks(cif_files, args.cores)
    with multiprocessing.Pool(processes=args.cores) as pool: pool.map(run_chimerax_chunk, chunk_files)

    merge_csvs()
    if os.path.exists(CHUNK_DIR): shutil.rmtree(CHUNK_DIR)

if __name__ == "__main__":
    main()
