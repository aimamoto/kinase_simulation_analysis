
#!/usr/bin/env python3
import subprocess
import sys
import os
import csv
import itertools
import glob
import re
import pandas as pd
from collections import Counter

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
        return[]

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
        # Strictly search the current directory to prevent duplicate nested parsing.
        raw_candidates = list(set(glob.glob("*_results_v7*.csv") + glob.glob("*_results.csv")))
        ignore_pattern = re.compile(r"temp_chimerax_chunks|archive|old", re.IGNORECASE)
        csv_candidates = [f for f in raw_candidates if not ignore_pattern.search(f.replace("\\", "/"))]
        
        if not csv_candidates:
            print("\n[!] Error: Could not find any valid HMM output CSVs in the working directory.")
            sys.exit(1)
            
        print(f"\n[*] Discovered {len(csv_candidates)} root result CSV files.")
        
        # Quick Data Previewer
        print("\n[*] Available Data Preview (Unique Naming Formats):")
        sample_ids = set()
        for c in csv_candidates:
            with open(c, 'r', encoding='utf-8-sig') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    sim_id = str(row.get('Simulation_ID', ''))
                    if sim_id:
                        base = sim_id.split('_')[0] 
                        sample_ids.add(base)
        print("    " + ", ".join(sorted(list(sample_ids))[:15]))

        print("\n[Optional] Filter the dataset to analyze a specific biological cohort.")
        print("    Syntax: Use commas for 'AND', pipes '|' for 'OR', minus '-' to 'EXCLUDE'.")
        print("    Prefix with '=' for WT Anchor (e.g., =EGFR finds 'egfrcat' but ignores 'egfrcat-t790m').")
        print("    Prefix with '~' for Regex (e.g., ~EGFR(?![_-])).")
        print("    Use ':' for Dimer Pairs (e.g., =EGFR:=EGFR | =EGFR:HER2).")
        filter_str = input("    -> Filter by keyword(s) [Leave blank to use all]: ").strip()
        
        filter_suffix = ""
        if filter_str:
            raw_conditions = [c.strip() for c in filter_str.split(',')]
            inc_blocks = []
            exc_keys = []
            
            for c in raw_conditions:
                if not c: continue
                if c.startswith('-'):
                    exc_keys.append(c.lstrip('-').strip())
                else:
                    or_terms = [t.strip() for t in c.split('|') if t.strip()]
                    if or_terms: inc_blocks.append(or_terms)
            
            filtered_rows = []
            for c in csv_candidates:
                with open(c, 'r', encoding='utf-8-sig') as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        row_str_orig = str(row.get('Directory', '')) + " " + str(row.get('Simulation_ID', ''))
                        row_str_lower = row_str_orig.lower()
                        
                        has_inc = True
                        for block in inc_blocks:
                            block_match = False
                            for term in block:
                                sub_terms = [s.strip() for s in term.split(':')]
                                sub_counts = Counter(sub_terms)
                                term_match = True
                                for sub, req_count in sub_counts.items():
                                    if sub.startswith('~'):
                                        try:
                                            if len(re.findall(sub[1:], row_str_orig, re.IGNORECASE)) < req_count:
                                                term_match = False; break
                                        except re.error as e:
                                            print(f"\n[!] Invalid regex '{sub[1:]}': {e}")
                                            sys.exit(1)
                                    elif sub.startswith('='):
                                        strict_term = re.escape(sub[1:])
                                        pattern = rf"(?<![a-zA-Z]){strict_term}(?![a-zA-Z0-9]*-)"
                                        if len(re.findall(pattern, row_str_orig, re.IGNORECASE)) < req_count:
                                            term_match = False; break
                                    else:
                                        if row_str_lower.count(sub.lower()) < req_count:
                                            term_match = False; break
                                if term_match:
                                    block_match = True
                                    break
                            if not block_match:
                                has_inc = False
                                break
                        
                        has_exc = False
                        for k in exc_keys:
                            sub_terms = [s.strip() for s in k.split(':')]
                            sub_counts = Counter(sub_terms)
                            k_match = True
                            for sub, req_count in sub_counts.items():
                                if sub.startswith('~'):
                                    try:
                                        if len(re.findall(sub[1:], row_str_orig, re.IGNORECASE)) < req_count:
                                            k_match = False; break
                                    except re.error: pass
                                elif sub.startswith('='):
                                    strict_sub = re.escape(sub[1:])
                                    pattern = rf"(?<![a-zA-Z]){strict_sub}(?![a-zA-Z0-9]*-)"
                                    if len(re.findall(pattern, row_str_orig, re.IGNORECASE)) < req_count:
                                        k_match = False; break
                                else:
                                    if row_str_lower.count(sub.lower()) < req_count:
                                        k_match = False; break
                            if k_match:
                                has_exc = True
                                break
                        
                        if has_inc and not has_exc:
                            filtered_rows.append(row)
            
            if not filtered_rows:
                print("\n[!] Filter too strict. No structural models matched your criteria. Exiting.")
                sys.exit(1)
                
            print(f"    -> [✓] Filter applied. Kept {len(filtered_rows)} structural models.")
            
            master_header = []
            for row in filtered_rows:
                for k in row.keys():
                    if k not in master_header and k is not None:
                        master_header.append(k)
            
            temp_csv = ".temp_filtered_results.csv"
            with open(temp_csv, 'w', newline='', encoding='utf-8') as f:
                writer = csv.DictWriter(f, fieldnames=master_header, extrasaction='ignore')
                writer.writeheader()
                writer.writerows(filtered_rows)
                
            csv_candidates = [temp_csv]
            
            clean_filename = lambda s: re.sub(r'[^A-Za-z0-9_]', '', s).replace(':', 'AND')
            suffix_parts = ["_OR_".join(clean_filename(t.upper()) for t in block) for block in inc_blocks]
            suffix_parts += [f"NO_{clean_filename(k.upper())}" for k in exc_keys]
            if suffix_parts: 
                filter_suffix = "_" + "_".join(suffix_parts).replace(" ", "")
                
        temp_list_path = ".temp_active_csv_list.txt"
        with open(temp_list_path, "w") as f:
            for csv_file in csv_candidates: f.write(f"{os.path.abspath(csv_file)}\n")
        
        print("\n[?] Choose the Biological Engine Type:")
        print("    1. Standard Kinase Engine (e.g., BRAF, SRC, Symmetric Homodimers)")
        print("    2. Asymmetric Dimer Engine (e.g., Massive ERBB Network Discovery)")
        while True:
            engine_choice = input("    -> Enter 1 or 2: ").strip()
            if engine_choice in ['1', '2']: break
            
        if engine_choice == '2':
            print("\n" + "="*75)
            print(" [!] Asymmetric Dimer (ERBB) Engine: UNDER CONSTRUCTION")
            print("     Module 2 receptor-family (ERBB) analysis is not yet released.")
            print("="*75)
            if os.path.exists(temp_list_path): os.remove(temp_list_path)
            if os.path.exists(".temp_filtered_results.csv"): os.remove(".temp_filtered_results.csv")
            sys.exit(0)

        engine_script = "multimer_core_engine.R"
        
        if engine_choice == '1':
            while True:
                print("\n[?] Enter the Target Protein Name (e.g., BRAF, SRC):")
                target_protein = input("    -> Target Protein: ").strip().upper()
                if target_protein and not target_protein.isdigit(): break
        else:
            print("\n[?] Enter a Project/Network Name for Output Labeling.")
            print("    (The script will automatically discover & parse all proteins inside the data)")
            target_protein = input("    -> Project Name [Default: ERBB_NETWORK]: ").strip().upper()
            if not target_protein: target_protein = "ERBB_NETWORK"
            
        print("\n[?] Choose your Meta-Stable State Discovery logic:")
        print("    1. Gaussian Mixture Models (Allows non-spherical state basins) [Default]")
        print("    2. K-Means + Gap Statistic (Strict spherical basins)")
        while True:
            cluster_method = input("    -> Clustering (1 or 2): ").strip()
            if not cluster_method: cluster_method = '1'
            if cluster_method in ['1', '2']: break
            
        cluster_method = "kmeans" if cluster_method == "2" else "gmm"
        
        engine_suffix = "" if engine_choice == '1' else "_ERBB"
        base_dir = f"plots_and_stats_{target_protein}_{cluster_method.upper()}{engine_suffix}{filter_suffix}"
        
        print("\n" + "-"*75)
        print(f" [*] EXECUTING: {'Asymmetric Engine' if engine_choice=='2' else 'Standard Engine'}")
        print(f" [*] OUTPUT DIR: {base_dir}")
        print("-" * 75)
        
        os.environ["CUSTOM_OUT_DIR"] = base_dir
        core_run = subprocess.run(["Rscript", os.path.join("modules", engine_script), target_protein, cluster_method, temp_list_path])
        
        if os.path.exists(temp_list_path): os.remove(temp_list_path)
        if os.path.exists(".temp_filtered_results.csv"): os.remove(".temp_filtered_results.csv")
            
        if core_run.returncode != 0: sys.exit(1)
            
        run_posthoc = input(f"\n[?] Run Phase 8 (Post-Hoc Analysis) on '{base_dir}'? (y/n): ").strip().lower() in ['y', 'yes']

    elif choice == '2':
        while True:
            user_dir = input("\n[?] Enter Phase 7 directory [Default: .]: ").strip()
            base_dir = user_dir if user_dir else "."
            if os.path.exists(os.path.join(base_dir, "Phase7_Complete_Structural_Metadata.csv")): break
            print(f"    [!] Could not find metadata CSV in {base_dir}.")
        run_posthoc = True

    if run_posthoc:
        family_dirs = sorted(glob.glob(os.path.join(base_dir, "Family_*")))
        if family_dirs:
            print("\n[?] Which dataset would you like to analyze for Phase 8?")
            print("    0. ALL FAMILIES (Entire Network - Allows Cross-Family Comparisons)")
            for idx, fdir in enumerate(family_dirs):
                print(f"    {idx + 1}. {os.path.basename(fdir)} (Isolated)")
            while True:
                try:
                    f_idx = int(input(f"    -> Select a folder (0-{len(family_dirs)}): ").strip())
                    if f_idx == 0: break
                    elif 1 <= f_idx <= len(family_dirs):
                        base_dir = family_dirs[f_idx - 1]
                        break
                except ValueError: pass
        
        csv_path = os.path.join(base_dir, "Phase7_Complete_Structural_Metadata.csv")
        if not os.path.exists(csv_path): sys.exit(1)
        
        with open(csv_path, 'r', encoding='utf-8-sig') as f: headers = next(csv.reader(f))
            
        print("\n[?] What would you like to compare in Phase 8?")
        print("    1. Meta-Stable States (e.g., State 1 vs State 2)")
        print("    2. Specific Sub-Populations (e.g., WT Homodimer vs L858R Heterodimer)")
        while True:
            comp_mode = input("    -> Enter 1 or 2: ").strip()
            if comp_mode in ['1', '2']: break
            
        group_col = "Macro_State" if comp_mode == '1' else ("Group_Name" if "Group_Name" in headers else "Condition")
        available_targets = get_unique_values(csv_path, group_col)
        
        print("\n[?] Select Effect Size Metric for Phase 8 Volcano Plots:")
        print("    1. Cohen's d (Parametric, standard but assumes normality)")
        print("    2. Rank-Biserial Correlation (Non-parametric, robust)")
        while True:
            eff_choice = input("    -> Enter 1 or 2 [Default: 2]: ").strip()
            if not eff_choice: eff_choice = '2'
            if eff_choice in ['1', '2']: break
        eff_metric_arg = "cohens_d" if eff_choice == '1' else "wilcox"

        mode = input("\n[?] Mode: (1) Lazy (Auto-run all) (2) Sequential (3) Manual -> ").strip()
        
        pairs_to_run =[]
        if mode == '1': 
            pairs_to_run = list(itertools.combinations(available_targets, 2))
            if group_col == "Group_Name":
                pairs_to_run = [p for p in pairs_to_run if "@" in p[0] and "@" in p[1] and p[0].split("@")[-1] == p[1].split("@")[-1]]
        elif mode == '2':
            pairs_to_run = [(available_targets[i], available_targets[i+1]) for i in range(len(available_targets)-1)]
        elif mode == '3':
            print("\n[Available Targets:]")
            for idx, tgt in enumerate(available_targets): print(f"  {idx + 1}. {tgt.replace('@', ' | ')}")
            while True:
                try:
                    idx1 = int(input(f"\n    -> Select Target A (1-{len(available_targets)}): ")) - 1
                    idx2 = int(input(f"    -> Select Target B (1-{len(available_targets)}): ")) - 1
                    pairs_to_run = [(available_targets[idx1], available_targets[idx2])]
                    break
                except ValueError: pass

        success_count = 0
        print("-" * 75)
        for state_a, state_b in pairs_to_run:
            print(f"    -> Running Phase 8 comparison: {state_a.replace('@', ' | ')} vs {state_b.replace('@', ' | ')}")
            result = subprocess.run(["Rscript", os.path.join("modules", "posthoc_differentiate.R"), 
                                     base_dir, state_a.replace('\n', '___'), state_b.replace('\n', '___'), group_col, eff_metric_arg],
                                     stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            if result.returncode == 0: success_count += 1
            else: print(f"       [!] Warning: Data might be too sparse for this plot.")
                
        print(f"\n[✓] Phase 8 Complete! Successfully generated {success_count}/{len(pairs_to_run)} comparisons.")

if __name__ == "__main__":
    try: main()
    except KeyboardInterrupt: sys.exit(0)

