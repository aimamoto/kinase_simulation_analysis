#!/usr/bin/env python3
import subprocess
import sys
import os
import csv
import itertools
import glob
import re

def get_unique_values(csv_path, column_name):
    """Parses the Phase 7 CSV to find all unique values in a specified column."""
    values = set()
    try:
        with open(csv_path, 'r', encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            for row in reader:
                if column_name in row and row[column_name].strip():
                    values.add(row[column_name].strip())
        return sorted(list(values), key=lambda x: int(x.split()[-1]) if x.split()[-1].isdigit() else x)
    except Exception as e:
        print(f"    [!] Could not parse {column_name} from CSV: {e}")
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
        print("    [!] Invalid choice.")
    
    run_posthoc = False
    base_dir = "." 
    engine_choice = "1"

    if choice == '1':
        # Single Source of Truth: Discover all CSVs in Python
        raw_candidates = glob.glob("**/*_results_v6.csv", recursive=True) + glob.glob("**/*_results.csv", recursive=True)
        
        # Rigorous filtering logic
        ignore_pattern = re.compile(r"temp_chimerax_chunks|archive|old", re.IGNORECASE)
        csv_candidates = [f for f in raw_candidates if not ignore_pattern.search(f.replace("\\", "/"))]
        
        if not csv_candidates:
            print("\n[!] Error: Could not find any valid HMM analysis results CSVs in the active directories.")
            sys.exit(1)
            
        print(f"\n[*] Discovered {len(csv_candidates)} valid HMM output CSV(s) across active directories.")
        
        # Write exactly these files to a temp list to pass to R
        temp_list_path = ".temp_active_csv_list.txt"
        with open(temp_list_path, "w") as f:
            for csv_file in csv_candidates:
                f.write(f"{os.path.abspath(csv_file)}\n")
        
        print("\n[?] Choose the Biological Engine Type:")
        print("    1. Standard Kinase Engine (e.g., CSK, SRC, Monomers)")
        print("    2. Asymmetric Dimer Engine (e.g., EGFR, ERBB, competitive binding)")
        
        while True:
            engine_choice = input("    -> Enter 1 or 2: ").strip()
            if engine_choice in ['1', '2']: break
            print("    [!] Invalid choice.")
            
        engine_script = "multimer_core_engine.R" if engine_choice == '1' else "erbb_asymmetric_engine.R"
        engine_label = "Standard Engine" if engine_choice == '1' else "ERBB Asymmetric Engine"
        
        if engine_choice == '1':
            while True:
                target_protein = input("\n[?] Enter the primary target protein (e.g., CSK, SRC): ").strip().upper()
                if not target_protein: print("    [!] You must enter a target protein.")
                elif target_protein.isdigit(): print("    [!] Protein name cannot be just a number.")
                else: break
        else:
            print("\n[?] Enter the name of the base protein that contains the mutations/variants.")
            print("    (e.g., If comparing ERBB2 vs EGFR-L858R and ERBB3 vs EGFR-T790M, enter 'EGFR')")
            while True:
                target_protein = input("    -> Target Protein: ").strip().upper()
                if not target_protein: print("    [!] You must enter a target protein.")
                elif target_protein.isdigit(): print("    [!] Protein name cannot be just a number.")
                else: break
        
        print("\n[?] Choose your Meta-Stable State Discovery logic:")
        print("    1. Gaussian Mixture Models (Allows non-spherical basins)[Default]")
        print("    2. K-Means + Gap Statistic (Strict spherical basins)")
        
        while True:
            logic_choice = input("    -> Enter 1 or 2[Default: 1]: ").strip()
            if not logic_choice:
                logic_choice = '1'
                break
            elif logic_choice in ['1', '2']: break
        
        cluster_method = "kmeans" if logic_choice == "2" else "gmm"
        
        print("\n" + "-"*75)
        print(f" [*] EXECUTING: {engine_label} | Target: {target_protein} | Clustering: {cluster_method.upper()}")
        print("-" * 75)
        input("[Press Enter to acknowledge and begin execution] ")
        
        core_script = os.path.join("modules", engine_script)
        
        # Pass the exact temp list file instead of a directory path
        core_run = subprocess.run(["Rscript", core_script, target_protein, cluster_method, temp_list_path])
        
        # Clean up the temp list
        if os.path.exists(temp_list_path):
            os.remove(temp_list_path)
        
        if core_run.returncode != 0:
            print("\n[!] Core pipeline failed. Exiting.")
            sys.exit(1)
            
        engine_suffix = "" if engine_choice == '1' else "_ERBB"
        base_dir = f"plots_and_stats_{target_protein}_{cluster_method.upper()}{engine_suffix}"
        
        while True:
            proceed = input(f"\n[?] Core pipeline complete. Output saved to '{base_dir}'. \n    Run Phase 8 (Post-Hoc Analysis)? (y/n): ").strip().lower()
            if proceed in['y', 'yes', 'n', 'no']:
                run_posthoc = (proceed in['y', 'yes'])
                break
            
        if not run_posthoc:
            print("\n[✓] Pipeline finished. Exiting.")

    elif choice == '2':
        while True:
            print("\n[?] Enter the directory containing the Phase 7 CSV")
            user_dir = input("    -> Directory [Default: current folder]: ").strip()
            base_dir = user_dir if user_dir else "."
            csv_check = os.path.join(base_dir, "Phase7_Complete_Structural_Metadata.csv")
            if os.path.exists(csv_check): break
            print(f"    [!] Could not find '{csv_check}'.")
        run_posthoc = True

    # --- PHASE 8: DUAL-MODE POST-HOC ANALYSIS ---
    if run_posthoc:
        
        if engine_choice == '2' and choice == '1':
            family_dirs = sorted(glob.glob(os.path.join(base_dir, "Family_*")))
            if family_dirs:
                print("\n[?] Which Dimer Family would you like to analyze for Phase 8?")
                for idx, fdir in enumerate(family_dirs):
                    print(f"    {idx + 1}. {os.path.basename(fdir)}")
                while True:
                    try:
                        f_idx = int(input(f"    -> Select a folder (1-{len(family_dirs)}): ").strip()) - 1
                        if 0 <= f_idx < len(family_dirs):
                            base_dir = family_dirs[f_idx]
                            break
                    except ValueError: pass
        
        csv_path = os.path.join(base_dir, "Phase7_Complete_Structural_Metadata.csv")
        
        if not os.path.exists(csv_path):
            print(f"\n[!] Error: Cannot find '{csv_path}'. Exiting.")
            sys.exit(1)
        
        print("\n[?] What would you like to compare in Phase 8?")
        print("    1. Meta-Stable States (e.g., State 1 vs State 2)")
        print("    2. Specific Sub-Populations (Requires a 'Group_Name' or 'Condition' column)")
        
        while True:
            comp_mode = input("    -> Enter 1 or 2: ").strip()
            if comp_mode in['1', '2']: break
            print("    [!] Invalid choice.")
            
        with open(csv_path, 'r', encoding='utf-8-sig') as f:
            headers = next(csv.reader(f))
            
        if comp_mode == '1': group_col = "Macro_State"
        else: group_col = "Group_Name" if "Group_Name" in headers else "Condition"
            
        available_targets = get_unique_values(csv_path, group_col)
        
        if not available_targets:
            print(f"\n[!] No valid entries found for {group_col} in the CSV. Exiting.")
            sys.exit(1)
            
        print(f"\n[*] Discovered {len(available_targets)} valid {group_col}s in {base_dir}")
        
        print("\n[?] How would you like to run the differential analysis?")
        print("    1. Lazy Mode (Auto-run ALL unique pairwise comparisons)")
        print("    2. Sequential Mode (1 vs 2, 2 vs 3, etc.)")
        print("    3. Manual Mode (Pick from a numbered list)")
        
        while True:
            mode = input("    -> Enter 1, 2, or 3: ").strip()
            if mode in['1', '2', '3']: break
            
        pairs_to_run =[]
        
        if mode == '1':
            pairs_to_run = list(itertools.combinations(available_targets, 2))
        elif mode == '2':
            pairs_to_run =[(available_targets[i], available_targets[i+1]) for i in range(len(available_targets)-1)]
        elif mode == '3':
            print(f"\n    [Available {group_col}s:]")
            for idx, val in enumerate(available_targets):
                print(f"      {idx + 1}. {val.replace(chr(10), ' + ')}")
                
            while True:
                try:
                    s1_idx = int(input(f"\n    -> Enter the NUMBER of the FIRST target (1-{len(available_targets)}): ").strip()) - 1
                    s2_idx = int(input(f"    -> Enter the NUMBER of the SECOND target (1-{len(available_targets)}): ").strip()) - 1
                    
                    if 0 <= s1_idx < len(available_targets) and 0 <= s2_idx < len(available_targets):
                        pairs_to_run =[(available_targets[s1_idx], available_targets[s2_idx])]
                        break
                    else:
                        print("    [!] Invalid numbers. Try again.")
                except ValueError:
                    print("    [!] Please enter valid numbers.")

        print(f"\n[*] Executing {len(pairs_to_run)} pairwise comparisons...")
        print("-" * 75)
        
        posthoc_script = os.path.join("modules", "posthoc_differentiate.R")
        success_count = 0
        
        for state_a, state_b in pairs_to_run:
            display_a = state_a.replace('\n', ' + ')
            display_b = state_b.replace('\n', ' + ')
            print(f"    -> Running: {display_a}  VS  {display_b}...")
            
            safe_a = state_a.replace('\n', '___')
            safe_b = state_b.replace('\n', '___')
            
            posthoc_run = subprocess.run(["Rscript", posthoc_script, base_dir, safe_a, safe_b, group_col],
                stdout=subprocess.DEVNULL, 
                stderr=subprocess.STDOUT
            )
            if posthoc_run.returncode == 0: success_count += 1
            else: print(f"       [!] Warning: Error generating plot for this pair. Data might be too sparse.")

        print("-" * 75)
        print(f"\n[✓] Phase 8 Complete! Successfully generated {success_count}/{len(pairs_to_run)} comparisons in '{base_dir}/Phase8_Volcanos/'.")

if __name__ == "__main__":
    try: main()
    except KeyboardInterrupt:
        print("\n\n[!] Pipeline aborted by user.")
        sys.exit(0)
