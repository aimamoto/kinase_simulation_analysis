#!/usr/bin/env python3
import os
import re
import csv
import argparse
import sys
from datetime import datetime

try:
    import yaml
except ImportError:
    print("\n❌ ERROR: The 'yaml' module is missing.")
    sys.exit(1)

DEFAULT_MATRIX = "generated_matrix.csv"
DEFAULT_YAML = "proteins.yaml"

def clean_protein_name(raw_name):
    """Universally cleans structural states, PTMs, and construct tags, preserving mutations."""
    name = re.sub(r'^\d+[_\-|]', '', raw_name)
    
    # 1. Broaden to catch apo, holo, 0atp, 1atp, gtp, etc., ANYWHERE in the string
    name = re.sub(r'[-_](apo|holo|\d*atp|\d*adp|\d*amp|\d*gtp|\d*gdp|\d*anp)\b', '', name, flags=re.IGNORECASE)
    
    # 2. Smartly remove construct tags (wt, cattail, cat, tail) ONLY at the end of a block/word
    # This matches the pristine regex used in the R Engine.
    name = re.sub(r'(?i)(cattail|cat|tail|wt)(?=_|-|$)', '', name)
    
    # 3. Catch all phosphorylation types (pY159, pS494, pT491) ANYWHERE in the string
    name = re.sub(r'[-_]p[sty]?\d+\b', '', name, flags=re.IGNORECASE)
    
    # 4. Clean up dangling dashes
    name = re.sub(r'[-_]+', '-', name).strip('-')
    
    return name if name else raw_name

def scan_directories():
    combos, singles = set(), set()
    IGNORE_DIRS = {'modules', 'archives', 'temp_chimerax_chunks', 'old', 'cx_viz_core', 'cx_viz_allosteric'}

    print("[*] Scanning directories recursively (Massive Mode Enabled)...")
    for root, dirs, files in os.walk("."):
        # Filter ignored directories in-place to speed up the recursive walk
        dirs[:] = [d for d in dirs if not d.startswith('.') and d not in IGNORE_DIRS and not d.startswith('plots_and_stats')]

        # Check if this specific folder contains structure files
        has_structure = any(f.lower().endswith(('.cif', '.pdb')) for f in files)
        if not has_structure:
            continue
            
        # Regex search the entire path for the interaction pattern
        # e.g., ./egfr-egfr/a-egfrcat_b-egfrcat_1atp/model.cif
        multi_match = re.search(r"a-([^_/\\]+)_b-([^_/\\]+)", root, re.IGNORECASE)
        
        if multi_match:
            combos.add((multi_match.group(1), multi_match.group(2)))
            continue
        
        # Logic for single structure fallback
        path_parts = root.replace('\\', '/').split('/')
        base_name = path_parts[-1]
        
        # If we are deep inside a seed or model folder, step up to the descriptive name
        if re.search(r'(seed|model|fold|sample|unrelaxed|relaxed)', base_name, re.IGNORECASE) and len(path_parts) > 1:
            base_name = path_parts[-2]
        
        singles.add(base_name.split('_')[0])

    if not combos and not singles: 
        return None, None, None
        
    chain_a = sorted(list(set(c[0] for c in combos)))
    chain_b = sorted(list(set(c[1] for c in combos)))
    
    matrix = {a: {b: '' for b in chain_b} for a in chain_a}
    for a, b in combos: matrix[a][b] = 'x'
    return matrix, chain_b, sorted(list(singles))

def save_matrix_csv(matrix, chain_b, singles, filename=DEFAULT_MATRIX):
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["# MULTIMER MATRIX (Chain A rows, Chain B columns)"])
        writer.writerow([""] + chain_b)
        for a, row_data in matrix.items():
            writer.writerow([a] + [row_data[b] for b in chain_b])
            
        writer.writerow([])
        writer.writerow(["# SINGLE PROTEIN RUNS"])
        for s in singles: writer.writerow([s, "x"])
    print(f"✅ Scanned directories and saved matrix to: {filename}")

def load_matrix_and_make_yaml(csv_path):
    patterns = {}
    with open(csv_path, 'r') as f:
        reader = csv.reader(f)
        current_section = None
        chain_b_vars =[]
        for row in reader:
            if not row or not row[0]: continue
            if "# MULTIMER" in row[0]:
                current_section = "MULTI"
                chain_b_vars = next(reader)[1:]
                continue
            elif "# SINGLE" in row[0]:
                current_section = "SINGLE"
                continue
            
            if current_section == "MULTI":
                a_var = row[0]
                for i, val in enumerate(row[1:]):
                    if val.lower() == 'x':
                        b_var = chain_b_vars[i]
                        pattern = f"a-{a_var}_b-{b_var}"
                        clean_a = clean_protein_name(a_var)
                        clean_b = clean_protein_name(b_var)
                        patterns[pattern] = {"proteins":[{"name": clean_a.lower()}, {"name": clean_b.lower()}]}
            elif current_section == "SINGLE":
                if row[1].lower() == 'x':
                    patterns[row[0]] = {"proteins":[{"name": clean_protein_name(row[0]).lower()}]}
    
    config = {"metadata": {"generated": datetime.now().isoformat(), "tool": "generate_config.py"}, "pattern_matches": patterns}
    with open(DEFAULT_YAML, "w") as f: yaml.dump(config, f, sort_keys=False)
    print(f"✅ Generated {DEFAULT_YAML} from {csv_path} ({len(patterns)} total patterns).")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--matrix", help="Path to an existing CSV matrix")
    args = parser.parse_args()

    if args.matrix: load_matrix_and_make_yaml(args.matrix)
    else:
        m, b_vars, s_vars = scan_directories()
        if m or s_vars:
            save_matrix_csv(m, b_vars, s_vars)
            print(f"👉 Review {DEFAULT_MATRIX}, then run pipeline with --resume")
        else: print("❌ No valid simulation directories found containing .cif/.pdb files.")
