#!/usr/bin/env python3
import os
import re
import glob
import json
import argparse
import subprocess
import sys
import math
import shutil
from Bio.PDB import PDBParser, MMCIFParser

try: import yaml
except ImportError: yaml = None

AA_MAP = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
    'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
    'HID':'H','HIE':'H','HIP':'H','HSD':'H','HSE':'H','HSP':'H','CYX':'C','CYM':'C','PTR':'Y','SEP':'S','TPO':'T',
    'LYZ':'K','ASX':'D','GLX':'E'
}

def clean_protein_name(raw_name):
    name = re.sub(r'^\d+[_\-|]', '', raw_name)
    name = re.sub(r'[-_](apo|holo|\d*atp|\d*adp|\d*amp|\d*gtp|\d*gdp|\d*anp)$', '', name, flags=re.IGNORECASE)
    name = re.sub(r'(wt|cattail|cat|tail)', '', name, flags=re.IGNORECASE)
    name = re.sub(r'[-_]p[sty]?\d+', '', name, flags=re.IGNORECASE) 
    name = re.sub(r'[-_]+', '-', name).strip('-')
    return name if name else raw_name

def extract_sequences(filepath):
    ext = os.path.splitext(filepath)[1].lower()
    parser = MMCIFParser(QUIET=True) if ext == ".cif" else PDBParser(QUIET=True)
    try: st = parser.get_structure("model", filepath)
    except Exception: return []
    out = []
    for model in st:
        for chain in model:
            seq = [AA_MAP[res.get_resname().upper().strip()] for res in chain if res.get_resname().upper().strip() in AA_MAP]
            if seq: out.append({"chain": chain.id, "sequence": "".join(seq)})
    return out

def get_protein_names(filepath, yaml_data=None):
    folder = os.path.dirname(os.path.abspath(filepath))
    basename = os.path.basename(folder)
    if re.search(r"seed[-_]\d+[-_]sample", basename, re.IGNORECASE): 
        basename = os.path.basename(os.path.dirname(folder))
    
    if yaml_data and "pattern_matches" in yaml_data:
        for pattern in sorted(yaml_data["pattern_matches"].keys(), key=len, reverse=True):
            if pattern.lower() in basename.lower(): 
                names = [p["name"] for p in yaml_data["pattern_matches"][pattern].get("proteins",[])]
                return {"is_multimer": len(names) > 1, "names": names}
                
    parts = basename.split("_")
    multimers = [p for p in parts if re.match(r"^[a-z]-", p, re.IGNORECASE)]
    names = [p.split("-", 1)[1] for p in multimers] if multimers else [parts[0]]
    
    return {
        "is_multimer": bool(multimers), 
        "names": [clean_protein_name(n).lower() for n in names]
    }

def run_worker(chunk_file, config_file, max_chains):
    worker_id = re.search(r'chunk_(\d+)', chunk_file).group(1) if re.search(r'chunk_(\d+)', chunk_file) else "?"
    
    with open(chunk_file, 'r') as f:
        files = [line.strip() for line in f if line.strip()]
        
    print(f"  -> [Worker {worker_id}] Booted up. Parsing {len(files)} files...", flush=True)

    yaml_data = None
    if config_file and os.path.exists(config_file) and yaml:
        with open(config_file, "r") as f: yaml_data = yaml.safe_load(f)

    results = []
    for fpath in files:
        name_data = get_protein_names(fpath, yaml_data)
        names = name_data["names"]
        
        seqs = extract_sequences(fpath)[:max_chains]
        
        for entry in seqs:
            seq, chain = entry["sequence"], entry["chain"]
            idx = ord(chain.upper()) - ord("A")
            pname = names[idx] if idx < len(names) else names[-1]
                
            results.append({"seq": seq, "chain": chain, "pname": pname, "file": fpath})
            
    out_json = chunk_file.replace(".txt", ".json")
    with open(out_json, "w") as f: json.dump(results, f)
    print(f"  <- [Worker {worker_id}] Finished! Extracted {len(results)} sequence chains.", flush=True)

def run_orchestrator(args):
    print("\n[*] Orchestrator: Scanning for structural files (.cif / .pdb)...")
    struct_files = glob.glob("**/*.cif", recursive=True) + glob.glob("**/*.pdb", recursive=True)
    if not struct_files: 
        print("[!] Orchestrator: No CIF/PDB files found.")
        return

    chunk_dir = "temp_fasta_chunks"
    os.makedirs(chunk_dir, exist_ok=True)

    num_cores = min(args.cores, len(struct_files))
    chunk_size = math.ceil(len(struct_files) / num_cores)
    print(f"[*] Orchestrator: Discovered {len(struct_files)} files. Partitioning across {num_cores} parallel workers...", flush=True)

    processes = []
    for i in range(num_cores):
        chunk_files = struct_files[i*chunk_size : (i+1)*chunk_size]
        if not chunk_files: continue
        
        chunk_txt = os.path.join(chunk_dir, f"chunk_{i}.txt")
        with open(chunk_txt, "w") as f: f.write("\n".join(chunk_files))
            
        cmd = [sys.executable, __file__, "--chunk", chunk_txt, "--max-chains", str(args.max_chains)]
        if args.config: cmd.extend(["--config", args.config])
        elif os.path.exists("proteins.yaml"): cmd.extend(["--config", "proteins.yaml"])
        
        p = subprocess.Popen(cmd)
        processes.append(p)

    for p in processes: p.wait()

    print(f"\n[*] Orchestrator: Workers completed. Merging and deduplicating sequence data...", flush=True)
    global_map = {}
    json_chunks = glob.glob(os.path.join(chunk_dir, "*.json"))
    
    for j_file in json_chunks:
        with open(j_file, 'r') as f:
            data = json.load(f)
            for item in data:
                seq, pname, chain, fpath = item["seq"], item["pname"], item["chain"], item["file"]
                if seq not in global_map: global_map[seq] = {"protein_name": pname, "chains": [chain], "example_file": fpath}
                else:
                    if len(pname) > len(global_map[seq]["protein_name"]): global_map[seq]["protein_name"] = pname
                    if chain not in global_map[seq]["chains"]: global_map[seq]["chains"].append(chain)

    with open(args.out, "w") as f:
        for seq, info in global_map.items(): f.write(f">{info['protein_name'].upper()}\n{seq}\n")

    out_map = {info["protein_name"]: {"sequence": seq, "chains": info["chains"], "example_file": info["example_file"]} for seq, info in global_map.items()}
    with open(args.map, "w") as f: json.dump(out_map, f, indent=2)

    shutil.rmtree(chunk_dir)
    print(f"✅ Extracted {len(global_map)} globally unique sequences to {args.out}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="extracted_sequences.fasta")
    parser.add_argument("--map", default="fasta_source_map.json")
    parser.add_argument("--config", help="YAML file for explicit naming")
    parser.add_argument("-c", "--cores", type=int, default=8, help="Number of CPU cores")
    parser.add_argument("--max-chains", type=int, default=2, help="Limit extraction to the first N chains")
    parser.add_argument("--chunk", type=str, help=argparse.SUPPRESS)
    args = parser.parse_args()

    if args.chunk: run_worker(args.chunk, args.config, args.max_chains)
    else: run_orchestrator(args)

if __name__ == "__main__":
    main()
