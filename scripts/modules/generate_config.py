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
    """Universally cleans structural states and construct tags, preserving mutations."""
    name = re.sub(r'^\d+[_\-|]', '', raw_name)
    name = re.sub(r'[-_](apo|holo|py\d+|\d*atp)$', '', name, flags=re.IGNORECASE)
    name = re.sub(r'(wt|cattail|cat|tail)', '', name, flags=re.IGNORECASE)
    name = re.sub(r'[-_]+', '-', name).strip('-')
    return name if name else raw_name

def scan_directories():
    combos, singles = set(), set()
    IGNORE_DIRS = {'modules', 'archives', 'temp_chimerax_chunks', 'old'}

    for d in os.listdir("."):
        if not os.path.isdir(d) or d.startswith('.') or d in IGNORE_DIRS: continue
        has_structure = any(f.lower().endswith(('.cif', '.pdb')) for root, dirs, files in os.walk(d) for f in files)
        if not has_structure: continue
            
        multi_match = re.search(r"a-([^_]+)_b-([^_]+)", d, re.IGNORECASE)
        if multi_match:
            combos.add((multi_match.group(1), multi_match.group(2)))
            continue
        
        parts = d.split("_")
        if parts: singles.add(parts[0])
    
    if not combos and not singles: return None, None, None
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
