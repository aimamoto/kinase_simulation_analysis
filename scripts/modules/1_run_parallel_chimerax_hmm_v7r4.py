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
WORKER_SCRIPT = os.path.join("modules", "chimerax_hmm_worker_v7r4.py")
FINAL_CSV_NAME = "hmm_kinase_analysis_results_v7r4.csv"
CHUNK_DIR = "temp_chimerax_chunks"

def ensure_hmm_landmarks():
    if os.path.exists("hmm_landmarks.json"): return
    try: subprocess.run(["python3", "0_extract_hmm_landmarks.py"], check=True)
    except Exception: sys.exit(1)

def filter_cif_files(cif_files):
    # Skip any path under an 'archive*' subdirectory (case-insensitive)
    cif_files = [c for c in cif_files
                 if not any(p.lower().startswith('archive') for p in c.replace('\\', '/').split('/'))]
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
    # v7r4 (BLOCKERS.md B12.6 item 2): return the chunk name and the failure
    # reason alongside the bool. v7r3 returned a bare bool that pool.map then
    # discarded, so a worker that died took its chunk's models with it silently.
    env = os.environ.copy()
    env["CHIMERAX_CHUNK"] = chunk_file
    cmd = ["chimerax", "--nogui", WORKER_SCRIPT]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, env=env)
        if result.returncode != 0:
            tail = (result.stderr or result.stdout or "").strip().splitlines()
            why = tail[-1] if tail else "no output"
            return (chunk_file, False, f"exit {result.returncode}: {why}")
        return (chunk_file, True, "")
    except Exception as e:
        return (chunk_file, False, f"{type(e).__name__}: {e}")

def merge_csvs():
    # Only merge per-chunk worker outputs (chunk_<i>_results_v7r3.csv). A broad
    # "*_results_v7*.csv" also matches the aggregate FINAL/master files, which
    # would concatenate a stale full dataset into the merge and duplicate rows.
    csv_files = glob.glob("chunk_*_results_v7r4*.csv") + glob.glob(os.path.join(CHUNK_DIR, "chunk_*_results_v7*.csv"))
    if FINAL_CSV_NAME in csv_files: csv_files.remove(FINAL_CSV_NAME)
    if not csv_files: return None

    df_list = []
    for csv_file in csv_files:
        try: df_list.append(pd.read_csv(csv_file))
        except pd.errors.EmptyDataError: pass
            
    if not df_list: return None

    final_df = pd.concat(df_list, ignore_index=True)
    if "Simulation_ID" in final_df.columns: final_df = final_df.sort_values(by="Simulation_ID")
        
    final_df.to_csv(FINAL_CSV_NAME, index=False)
    for csv_file in csv_files: os.remove(csv_file)
    return final_df

def verify_coverage(final_df, cif_files):
    """v7r4 (BLOCKERS.md B12.6 item 1). Every model handed to the workers must
    appear in the merged geometry CSV.

    On 260918 a whole 240-model arm produced zero rows and the run still exited
    0 with 'Pipeline Completed!' (B12.5). Nothing downstream could see the loss:
    the wrapper's merge into the master is a LEFT join FROM this CSV, so the
    master simply inherits the gap. Rows are per chain, so the exact count is
    not known ahead of time -- but model COVERAGE is, and that is the invariant
    that actually broke. Fail loudly and name the missing directories.
    """
    if final_df is None or final_df.empty:
        print(f"\n[!] FATAL: no rows were produced for any of {len(cif_files)} models.", file=sys.stderr)
        return False
    if not {"Directory", "File"}.issubset(final_df.columns):
        print("\n[!] Cannot verify coverage: merged CSV lacks Directory/File.", file=sys.stderr)
        return False

    # The worker records Directory/File from the same path string it was handed,
    # so reconstruct the expected keys the same way.
    expected = {(os.path.dirname(c), os.path.basename(c)) for c in cif_files}
    got = set(zip(final_df["Directory"].astype(str), final_df["File"].astype(str)))
    missing = sorted(expected - got)
    if not missing:
        print(f"[*] Coverage check: all {len(expected)} models produced rows "
              f"({len(final_df)} chain rows).")
        return True

    # Roll up to the arm (the parent of the per-seed directory) -- a lost arm is 240
    # directories on a run like 260918, and the arm is what the reader needs to see.
    by_arm = {}
    for d, f in missing:
        arm = os.path.dirname(d) or d
        by_arm[arm] = by_arm.get(arm, 0) + 1
    print("\n" + "=" * 75, file=sys.stderr)
    print(f" [!] FATAL: {len(missing)} of {len(expected)} models produced NO rows.", file=sys.stderr)
    print("", file=sys.stderr)
    print("     Missing models, by arm:", file=sys.stderr)
    for arm in sorted(by_arm):
        print(f"       {by_arm[arm]:6d}  {arm}", file=sys.stderr)
    print("", file=sys.stderr)
    print("     First missing directories:", file=sys.stderr)
    for d, f in missing[:10]:
        print(f"       {d}", file=sys.stderr)
    if len(missing) > 10:
        print(f"       ... and {len(missing) - 10} more", file=sys.stderr)
    print("", file=sys.stderr)
    print("     The commonest cause is B12.5: the landmark gate found no candidate", file=sys.stderr)
    print("     for these models, every chain fell through to the co-factor branch,", file=sys.stderr)
    print("     and the model emitted nothing. v7r4 falls back to the full landmark", file=sys.stderr)
    print("     set when the name gate matches nothing, so if you are seeing this,", file=sys.stderr)
    print("     check sequences.fasta and hmm_landmarks.json cover these chains at", file=sys.stderr)
    print("     all, and read the per-chunk CSVs left in " + CHUNK_DIR + ".", file=sys.stderr)
    print("=" * 75, file=sys.stderr)
    return False


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
    with multiprocessing.Pool(processes=args.cores) as pool:
        results = pool.map(run_chimerax_chunk, chunk_files)

    # v7r4 (B12.6 item 2): act on the worker return codes instead of dropping them.
    failed = [(c, why) for c, ok, why in results if not ok]
    if failed:
        print("\n" + "=" * 75, file=sys.stderr)
        print(f" [!] {len(failed)} of {len(chunk_files)} ChimeraX workers failed:", file=sys.stderr)
        for c, why in failed:
            print(f"       {c}  --  {why}", file=sys.stderr)
        print("=" * 75, file=sys.stderr)

    final_df = merge_csvs()
    covered = verify_coverage(final_df, cif_files)

    # v7r4: the chunk directory is the only evidence of what was lost, so keep it
    # when anything went wrong. v7r3 deleted it unconditionally.
    if failed or not covered:
        print(f"[!] Keeping {CHUNK_DIR} for diagnosis; not deleting it.", file=sys.stderr)
        sys.exit(1)
    if os.path.exists(CHUNK_DIR): shutil.rmtree(CHUNK_DIR)

if __name__ == "__main__":
    main()
