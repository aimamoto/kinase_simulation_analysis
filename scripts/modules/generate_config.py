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
    print("👉 Please install it by running: pip install pyyaml\n")
    sys.exit(1)

# Default filenames
DEFAULT_MATRIX = "generated_matrix.csv"
DEFAULT_YAML = "proteins.yaml"

def scan_directories():
    """Scans folders to find unique Chain A and Chain B combinations."""
    combos = set()
    singles = set()
    
    # Folders that the pipeline creates or we know we want to ignore
    IGNORE_DIRS = {'modules', 'archives', 'cx_pocket_visualization_hmm', 'temp_chimerax_chunks', 'old'}

    for d in os.listdir("."):
        # 1. Skip if it's not a directory, is hidden, or is in the ignore list
        if not os.path.isdir(d) or d.startswith('.') or d in IGNORE_DIRS:
            continue
            
        # 2. Verify it actually contains structural files (.cif or .pdb)
        has_structure = False
        for root, dirs, files in os.walk(d):
            if any(f.lower().endswith('.cif') or f.lower().endswith('.pdb') for f in files):
                has_structure = True
                break
                
        if not has_structure:
            continue
            
        # 3. Multimer detection: a-NAME_b-NAME
        # The [^_]+ safely stops capturing at the first underscore, 
        # ignoring _1atp_20260302_040645 etc.
        multi_match = re.search(r"a-([^_]+)_b-([^_]+)", d, re.IGNORECASE)
        if multi_match:
            combos.add((multi_match.group(1), multi_match.group(2)))
            continue
        
        # 4. Single protein fallback (only reached if structural files exist)
        parts = d.split("_")
        if parts:
            singles.add(parts[0])
    
    if not combos and not singles:
        return None, None, None
    
    chain_a = sorted(list(set(c[0] for c in combos)))
    chain_b = sorted(list(set(c[1] for c in combos)))
    
    # Build matrix for multimers
    matrix = {a: {b: '' for b in chain_b} for a in chain_a}
    for a, b in combos:
        matrix[a][b] = 'x'
        
    return matrix, chain_b, sorted(list(singles))

def save_matrix_csv(matrix, chain_b, singles, filename=DEFAULT_MATRIX):
    """Writes the discovered combinations to a CSV file."""
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Section 1: Multimer Matrix
        writer.writerow(["# MULTIMER MATRIX (Chain A rows, Chain B columns)"])
        writer.writerow([""] + chain_b)
        for a, row_data in matrix.items():
            writer.writerow([a] + [row_data[b] for b in chain_b])
            
        # Section 2: Single Protein List
        writer.writerow([])
        writer.writerow(["# SINGLE PROTEIN RUNS"])
        for s in singles:
            writer.writerow([s, "x"])
            
    print(f"✅ Scanned directories and saved matrix to: {filename}")

def load_matrix_and_make_yaml(csv_path):
    """Reads the CSV and generates the proteins.yaml pattern file strictly."""
    patterns = {}
    
    with open(csv_path, 'r') as f:
        reader = csv.reader(f)
        current_section = None
        chain_b_vars = []
        
        for row in reader:
            if not row or not row[0]: continue
            
            # Identify Sections
            if "# MULTIMER" in row[0]:
                current_section = "MULTI"
                headers = next(reader)
                chain_b_vars = headers[1:]
                continue
            elif "# SINGLE" in row[0]:
                current_section = "SINGLE"
                continue
            
            # Process Data STRICTLY (No hyphenation logic)
            if current_section == "MULTI":
                a_var = row[0]
                for i, val in enumerate(row[1:]):
                    if val.lower() == 'x':
                        b_var = chain_b_vars[i]
                        pattern = f"a-{a_var}_b-{b_var}"
                        
                        # --- CLEAN NAMES ---
                        # re.sub replaces ALL occurrences, so 'src-wtcat-py159-holo' becomes 'src-wtcat'
                        clean_a = re.sub(r'-(apo|holo|py\d+)', '', a_var, flags=re.IGNORECASE)
                        clean_b = re.sub(r'-(apo|holo|py\d+)', '', b_var, flags=re.IGNORECASE)
                        
                        patterns[pattern] = {
                            "proteins":[{"name": clean_a.lower()}, {"name": clean_b.lower()}]
                        }
            
            elif current_section == "SINGLE":
                name = row[0]
                if row[1].lower() == 'x':
                    # --- CLEAN NAME ---
                    clean_name = re.sub(r'-(apo|holo|py\d+|\d*atp)', '', name, flags=re.IGNORECASE)
                    
                    patterns[name] = {
                        "proteins":[{"name": clean_name.lower()}]
                    }
    
    config = {
        "metadata": {"generated": datetime.now().isoformat(), "tool": "generate_config.py"},
        "pattern_matches": patterns
    }
    
    with open(DEFAULT_YAML, "w") as f:
        yaml.dump(config, f, sort_keys=False)
    print(f"✅ Generated {DEFAULT_YAML} from {csv_path} ({len(patterns)} total patterns).")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--matrix", help="Path to an existing CSV matrix")
    args = parser.parse_args()

    if args.matrix:
        load_matrix_and_make_yaml(args.matrix)
    else:
        m, b_vars, s_vars = scan_directories()
        if m or s_vars:
            save_matrix_csv(m, b_vars, s_vars)
            print(f"👉 Review {DEFAULT_MATRIX}, then run pipeline with --resume")
        else:
            print("❌ No valid simulation directories found containing .cif/.pdb files.")
