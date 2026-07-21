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
    name = re.sub(r'^[a-zA-Z]-', '', name)  # Strip chain tags like 'a-', 'b-', 'c-'
    
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
        dirs[:] = [d for d in dirs if not d.startswith('.') and d not in IGNORE_DIRS and not d.startswith('plots_and_stats') and not d.lower().startswith('archive')]

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
        
        # Strip cosmetic single-letter prefixes (like a-, b-) from the raw extracted name
        raw_single = base_name.split('_')[0]
        clean_single = re.sub(r'^[a-zA-Z]-', '', raw_single)
        
        singles.add(clean_single)

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
    approved_multi = set()
    approved_single = set()
    
    # 1. Parse the audited 2D Matrix
    with open(csv_path, 'r') as f:
        reader = csv.reader(f)
        current_section = None
        chain_b_vars = []
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
                        approved_multi.add((a_var, chain_b_vars[i]))
            elif current_section == "SINGLE":
                if row[1].lower() == 'x':
                    approved_single.add(row[0])
                    
    # 2. Rescan directories to build the fully-resolved n-dimensional patterns (a, b, c...)
    patterns = {}
    IGNORE_DIRS = {'modules', 'archives', 'temp_chimerax_chunks', 'old', 'cx_viz_core', 'cx_viz_allosteric'}
    
    for root, dirs, files in os.walk("."):
        dirs[:] = [d for d in dirs if not d.startswith('.') and d not in IGNORE_DIRS and not d.startswith('plots_and_stats') and not d.lower().startswith('archive')]
        if not any(f.lower().endswith(('.cif', '.pdb')) for f in files): continue
            
        path_parts = root.replace('\\', '/').split('/')
        base_name = path_parts[-1]
        if re.search(r'(seed|model|fold|sample|unrelaxed|relaxed)', base_name, re.IGNORECASE) and len(path_parts) > 1:
            base_name = path_parts[-2]
            
        tokens = base_name.split("_")
        chain_tokens = [p for p in tokens if re.match(r"^[a-z]-", p, re.IGNORECASE)]
        
        if len(chain_tokens) >= 2:
            a_name = chain_tokens[0].split("-", 1)[1]
            b_name = chain_tokens[1].split("-", 1)[1]
            # If A and B were approved in the matrix, add the ENTIRE chain list to YAML
            if (a_name, b_name) in approved_multi:
                pattern_key = "_".join(chain_tokens)
                proteins = [{"name": clean_protein_name(t.split("-", 1)[1]).lower()} for t in chain_tokens]
                patterns[pattern_key] = {"proteins": proteins}
        elif len(chain_tokens) == 1:
            single_name = chain_tokens[0].split("-", 1)[1]
            if single_name in approved_single:
                patterns[chain_tokens[0]] = {"proteins": [{"name": clean_protein_name(single_name).lower()}]}
        else:
            clean_single = re.sub(r'^[a-zA-Z]-', '', base_name.split('_')[0])
            if clean_single in approved_single:
                patterns[base_name] = {"proteins": [{"name": clean_protein_name(clean_single).lower()}]}

    config = {"metadata": {"generated": datetime.now().isoformat(), "tool": "generate_config.py"}, "pattern_matches": patterns}
    with open(DEFAULT_YAML, "w") as f: yaml.dump(config, f, sort_keys=False)
    print(f"✅ Generated {DEFAULT_YAML} from {csv_path} ({len(patterns)} fully resolved multimeric patterns).")

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
