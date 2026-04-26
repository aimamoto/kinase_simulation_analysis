#!/usr/bin/env python3
import subprocess
import sys
import os
import csv
import itertools
import glob
import re
import pandas as pd

def get_unique_values(csv_path, column_name):
    values = set()
    try:
        with open(csv_path, 'r', encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            for row in reader:
                if column_name in row and row[column_name].strip():
                    values.add(row[column_name].strip())
        return sorted(list(values), key=lambda x: int(x.split()[-1]) if x.split()[-1].isdigit() else x)
    except Exception as e:
        return []

def main():
    print("="*75)
    print(" 📐 Kinase Allostery Discovery Pipeline (Dual-Engine)")
    print("="*75)

    if not os.path.isdir("modules"):
        print("\n[!] Error: 'modules' directory not found.")
        sys.exit(1)

    print("\n[?] Choose your starting point:")
    print("    1. Run the Full Core Pipeline (Phases 1-7)")
    print("    2. Skip directly to Post-Hoc Differential Analysis (Phase 8)")
    
    while True:
        choice = input("    -> Enter 1 or 2: ").strip()
        if choice in ['1', '2']: break
    
    run_posthoc = False
    base_dir = "." 
    engine_choice = "1"

    if choice == '1':
        raw_candidates = glob.glob("**/*_results_v6.csv", recursive=True) + glob.glob("**/*_results.csv", recursive=True)
        ignore_pattern = re.compile(r"temp_chimerax_chunks|archive|old", re.IGNORECASE)
        csv_candidates = [f for f in raw_candidates if not ignore_pattern.search(f.replace("\\", "/"))]
        
        if not csv_candidates:
            print("\n[!] Error: Could not find any valid HMM output CSVs.")
            sys.exit(1)
            
        temp_list_path = ".temp_active_csv_list.txt"
        with open(temp_list_path, "w") as f:
            for csv_file in csv_candidates: f.write(f"{os.path.abspath(csv_file)}\n")
        
        engine_choice = input("\n[?] Engine: (1) Standard or (2) ERBB Asymmetric Dimer? -> ").strip()
        engine_script = "multimer_core_engine.R" if engine_choice == '1' else "erbb_asymmetric_engine.R"
        
        target_protein = input("\n[?] Enter the target protein (e.g., BRAF, SRC): ").strip().upper()
        cluster_method = input("\n[?] Clustering: (1) GMM or (2) K-Means? [Default: 1] -> ").strip()
        cluster_method = "kmeans" if cluster_method == "2" else "gmm"
        
        print("\n" + "-"*75)
        print(f" [*] EXECUTING: Engine | Target: {target_protein} | Clustering: {cluster_method.upper()}")
        print("-" * 75)
        
        core_run = subprocess.run(["Rscript", os.path.join("modules", engine_script), target_protein, cluster_method, temp_list_path])
        if os.path.exists(temp_list_path): os.remove(temp_list_path)
        
        if core_run.returncode != 0: sys.exit(1)
            
        engine_suffix = "" if engine_choice == '1' else "_ERBB"
        base_dir = f"plots_and_stats_{target_protein}_{cluster_method.upper()}{engine_suffix}"
        
        run_posthoc = input(f"\n[?] Run Phase 8 (Post-Hoc Analysis) on '{base_dir}'? (y/n): ").strip().lower() in ['y', 'yes']

    elif choice == '2':
        base_dir = input("\n[?] Enter directory containing Phase 7 CSV [Default: .]: ").strip() or "."
        run_posthoc = True

    # --- PHASE 8: POST-HOC BINGO MODE ---
    if run_posthoc:
        csv_path = os.path.join(base_dir, "Phase7_Complete_Structural_Metadata.csv")
        if not os.path.exists(csv_path): sys.exit(1)
        
        with open(csv_path, 'r', encoding='utf-8-sig') as f:
            headers = next(csv.reader(f))
            
        comp_mode = input("\n[?] Compare: (1) Meta-States or (2) Sub-Populations (Condition)? -> ").strip()
        group_col = "Macro_State" if comp_mode == '1' else ("Group_Name" if "Group_Name" in headers else "Condition")
        available_targets = get_unique_values(csv_path, group_col)
        
        mode = input("\n[?] Mode: (1) Lazy (2) Sequential (3) Manual (4) Bingo Matrix -> ").strip()
        pairs_to_run = []
        
        if mode == '4':
            exp_csv_path = "experiment.csv"
            df_exp = pd.read_csv(exp_csv_path)

            df_exp.columns = [str(c).lower().strip() for c in df_exp.columns]
            if 'protein_a' in df_exp.columns: df_exp.rename(columns={'protein_a': 'chain_a'}, inplace=True)
            if 'protein_b' in df_exp.columns: df_exp.rename(columns={'protein_b': 'chain_b'}, inplace=True)
            if 'condition' in df_exp.columns:
                if 'condition_a' not in df_exp.columns: df_exp['condition_a'] = df_exp['condition']
                if 'condition_b' not in df_exp.columns: df_exp['condition_b'] = df_exp['condition']

            for col in ['chain_a', 'mutation_a', 'ptm_a', 'condition_a', 'chain_b', 'mutation_b', 'ptm_b', 'condition_b']:
                if col not in df_exp.columns: df_exp[col] = ""
                else: df_exp[col] = df_exp[col].fillna("").astype(str).str.strip()

            # Translates Bingo Mode logic to match Phase 1's physical Apo/Holo conversions
            def get_translated_states(row):
                base_a = "-".join([p for p in [row['chain_a'], row['mutation_a'], row['ptm_a']] if p and p != 'nan'])
                base_b = "-".join([p for p in [row['chain_b'], row['mutation_b'], row['ptm_b']] if p and p != 'nan'])
                cond = str(row.get('condition_a', '')).lower()
                
                if '0atp' in cond: return [f"{base_a}-apo\n{base_b}-apo"]
                if '2atp' in cond: return [f"{base_a}-holo\n{base_b}-holo"]
                if '1atp' in cond: 
                    # Add both possible structural outcomes for Bingo detection
                    return [f"{base_a}-holo\n{base_b}-apo", f"{base_a}-apo\n{base_b}-holo"]
                
                return [f"{base_a}-{cond}\n{base_b}-{cond}"]

            df_exp['Expected_Conditions'] = df_exp.apply(get_translated_states, axis=1)

            seen_pairs = set()
            records = df_exp.to_dict('records')

            # We mathematically check distances between logical groups
            for i in range(len(records)):
                for j in range(i + 1, len(records)):
                    diffs = 0
                    for col in ['chain_a', 'mutation_a', 'ptm_a', 'chain_b', 'mutation_b', 'ptm_b']:
                        if records[i][col] != records[j][col]: diffs += 1
                        
                    # If logic is 1 step away mathematically...
                    if diffs <= 1:
                        for target1 in records[i]['Expected_Conditions']:
                            for target2 in records[j]['Expected_Conditions']:
                                if target1 != target2 and target1 in available_targets and target2 in available_targets:
                                    pair_hash = tuple(sorted([target1, target2]))
                                    if pair_hash not in seen_pairs:
                                        seen_pairs.add(pair_hash)
                                        pairs_to_run.append((target1, target2))
            
            print(f"\n[*] Bingo Mode matched {len(pairs_to_run)} valid physical 1-variable transitions.")
            if not pairs_to_run: sys.exit(1)

        elif mode == '1': pairs_to_run = list(itertools.combinations(available_targets, 2))
        
        success_count = 0
        for state_a, state_b in pairs_to_run:
            print(f"    -> Running Phase 8 comparison...")
            subprocess.run(["Rscript", os.path.join("modules", "posthoc_differentiate.R"), base_dir, state_a.replace('\n', '___'), state_b.replace('\n', '___'), group_col], stdout=subprocess.DEVNULL)
            success_count += 1
        print(f"\n[✓] Phase 8 Complete! Successfully generated {success_count} comparisons.")

if __name__ == "__main__":
    try: main()
    except KeyboardInterrupt: sys.exit(0)
