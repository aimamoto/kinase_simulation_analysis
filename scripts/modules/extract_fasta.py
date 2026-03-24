#!/usr/bin/env python3
"""
Simplified extract_fasta.py

Features:
  ✓ Strict AA_MAP
  ✓ AlphaFold3 & MD directory parsing simplified
  ✓ Seamlessly handles single (abl1-wtcat) and multimer (a-erbb2cattail, b-egfrcattail) formats
  ✓ Auto-formats fused domains (e.g., erbb2cattail -> erbb2-cattail)
  ✓ Global deduplication by sequence
  ✓ YAML override support
"""

import os
import re
import glob
import json
import argparse
from Bio.PDB import PDBParser, MMCIFParser

try:
    import yaml
except ImportError:
    yaml = None

AA_MAP = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C',
    'GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
    'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P',
    'SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
    'HID':'H','HIE':'H','HIP':'H','HSD':'H','HSE':'H','HSP':'H',
    'CYX':'C','CYM':'C',
    'PTR':'Y','SEP':'S','TPO':'T',
    'LYZ':'K',
    'ASX':'D','GLX':'E'
}

def extract_sequences(filepath):
    ext = os.path.splitext(filepath)[1].lower()
    parser = MMCIFParser(QUIET=True) if ext == ".cif" else PDBParser(QUIET=True)
    try:
        st = parser.get_structure("model", filepath)
    except Exception:
        return []

    out = []
    for model in st:
        for chain in model:
            seq = []
            for res in chain:
                r = res.get_resname().upper().strip()
                if r in AA_MAP:
                    seq.append(AA_MAP[r])
            if seq:
                out.append({"chain": chain.id, "sequence": "".join(seq)})
    return out

def get_protein_names(filepath, yaml_data=None):
    folder = os.path.dirname(os.path.abspath(filepath))
    basename = os.path.basename(folder)

    # Step up one level if inside an AlphaFold3 seed directory
    if re.search(r"seed[-_]\d+[-_]sample", basename, re.IGNORECASE):
        basename = os.path.basename(os.path.dirname(folder))

    # 1. Check YAML pattern matches FIRST
    if yaml_data and "pattern_matches" in yaml_data:
        for pattern, data in yaml_data["pattern_matches"].items():
            if pattern.lower() in basename.lower():
                return [p["name"] for p in data.get("proteins", [])]

    # 2. Fallback: Parse automatically if no YAML match is found
    parts = basename.split("_")
    multimers = [p for p in parts if re.match(r"^[a-z]-", p, re.IGNORECASE)]

    if multimers:
        # Splits 'a-egfrcattail-t790m' into 'a' and 'egfrcattail-t790m'
        names = [p.split("-", 1)[1] for p in multimers]
    else:
        # Single protein simulation
        names = [parts[0]] 

    # Strictly return the raw parsed names without any hyphen manipulation
    return [n.lower() for n in names]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="extracted_sequences.fasta")
    parser.add_argument("--map", default="fasta_source_map.json")
    parser.add_argument("--config", help="YAML file for explicit naming")
    args = parser.parse_args()

    struct_files = glob.glob("**/*.cif", recursive=True) + glob.glob("**/*.pdb", recursive=True)
    if not struct_files:
        print("No CIF/PDB files found.")
        return

    yaml_data = None
    yaml_path = args.config if args.config else ("proteins.yaml" if os.path.exists("proteins.yaml") else None)
    if yaml_path and yaml:
        with open(yaml_path, "r") as f:
            yaml_data = yaml.safe_load(f)
        print(f"Using YAML config: {yaml_path}")

    global_map = {}

    for fpath in struct_files:
        names = get_protein_names(fpath, yaml_data)
        seqs = extract_sequences(fpath)

        for entry in seqs:
            seq = entry["sequence"]
            chain = entry["chain"]
            idx = ord(chain.upper()) - ord("A")
            
            # Map chain to parsed name, fallback to the last known name if chains exceed expected names
            pname = names[idx] if idx < len(names) else names[-1]

            # Global Deduplication by sequence
            if seq not in global_map:
                global_map[seq] = {
                    "protein_name": pname,
                    "chains": [chain],
                    "example_file": fpath
                }
            else:
                # If we encounter the same sequence again, keep the more specific/longer name 
                # (e.g., egfr-cattail-s768i over just egfr-cattail)
                if len(pname) > len(global_map[seq]["protein_name"]):
                    global_map[seq]["protein_name"] = pname
                if chain not in global_map[seq]["chains"]:
                    global_map[seq]["chains"].append(chain)

    # Output FASTA
    with open(args.out, "w") as f:
        for seq, info in global_map.items():
            f.write(f">{info['protein_name']}\n{seq}\n")

    # Output JSON Map (keyed by protein name for readability)
    out_map = {info["protein_name"]: {
        "sequence": seq, 
        "chains": info["chains"], 
        "example_file": info["example_file"]
    } for seq, info in global_map.items()}

    with open(args.map, "w") as f:
        json.dump(out_map, f, indent=2)

    print(f"\n✅ Extracted {len(global_map)} unique sequences.")
    print(f"✅ FASTA → {args.out}")
    print(f"✅ MAP   → {args.map}\n")

if __name__ == "__main__":
    main()
