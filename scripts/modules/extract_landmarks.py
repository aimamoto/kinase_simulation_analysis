#!/usr/bin/env python3
import os
import re
import json
import glob
import argparse
import subprocess
import sys

DEFAULT_HMM = os.path.expanduser("~/pfam/Pkinase.hmm")
DEFAULT_OUTFILE = "hmm_landmarks.json"
TEMP_A2M = "_temp_aligned.a2m"

TARGET_NODES = {30: "k", 48: "c", 52: "rs1", 56: "n99", 61: "v104", 62: "k105", 63: "rs2", 64: "e107", 75: "m118", 77: "m120", 78: "e121", 113: "i150", 119: "y156", 121: "hrd", 142: "f", 167: "ape", 179: "d220"}

def clean_protein_name(raw_name):
    name = raw_name.split('/')[0] 
    name = re.sub(r'^\d+[_\-|]', '', name)
    name = re.sub(r'^[a-zA-Z]-', '', name)  # Strip chain tags like 'a-', 'b-', 'c-'
    name = re.sub(r'[-_](apo|holo|py\d+|\d*atp)$', '', name, flags=re.IGNORECASE)
    name = re.sub(r'(wt|cattail|cat|tail)', '', name, flags=re.IGNORECASE)
    name = re.sub(r'[-_]+', '-', name).strip('-')
    return name.upper() if name else raw_name.upper()

def extract_nodes_from_aligned_seq(a2m_seq: str, fasta_name: str):
    seq_index, hmm_index = 0, 1  
    mapping = {v: None for v in TARGET_NODES.values()}
    for ch in a2m_seq:
        if ch.isalpha():
            if ch.isupper():
                if hmm_index in TARGET_NODES: mapping[TARGET_NODES[hmm_index]] = seq_index
                seq_index += 1; hmm_index += 1
            else: seq_index += 1
        elif ch == "-": hmm_index += 1
    mapping["type"] = clean_protein_name(fasta_name)
    return mapping

def run_hmmalign(hmm_path: str, fasta_path: str) -> str:
    # Filter out short peptides before running HMMer to keep JSON clean
    filtered_fasta = "_temp_filtered.fasta"
    with open(fasta_path, "r") as f_in, open(filtered_fasta, "w") as f_out:
        name, seq_lines = None, []
        for line in f_in:
            line = line.strip()
            if line.startswith(">"):
                if name and sum(len(x) for x in seq_lines) >= 100:
                    f_out.write(f">{name}\n{''.join(seq_lines)}\n")
                name, seq_lines = line[1:], []
            elif line:
                seq_lines.append(line)
        if name and sum(len(x) for x in seq_lines) >= 100:
            f_out.write(f">{name}\n{''.join(seq_lines)}\n")
            
    cmd = ["hmmalign", "--outformat", "A2M", "-o", TEMP_A2M, hmm_path, filtered_fasta]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if os.path.exists(filtered_fasta): os.remove(filtered_fasta)
    return TEMP_A2M

def parse_a2m(a2m_file: str) -> dict:
    landmarks, current_name, current_seq_chunks = {}, None, []
    with open(a2m_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                if current_name is not None: landmarks[current_name] = extract_nodes_from_aligned_seq("".join(current_seq_chunks), current_name)
                current_name = line[1:].strip().split()[0]
                current_seq_chunks = []
            else: current_seq_chunks.append(line.strip())
        if current_name is not None: landmarks[current_name] = extract_nodes_from_aligned_seq("".join(current_seq_chunks), current_name)
    return landmarks

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", type=str)
    parser.add_argument("--hmm", type=str, default=DEFAULT_HMM)
    parser.add_argument("--out", type=str, default=DEFAULT_OUTFILE)
    args = parser.parse_args()
    
    fastas = glob.glob("*.fasta") + glob.glob("*.fa") + glob.glob("*.faa")
    fasta_path = args.fasta if args.fasta else fastas[0]
    
    a2m_file = run_hmmalign(args.hmm, fasta_path)
    with open(args.out, "w") as f: json.dump(parse_a2m(a2m_file), f, indent=2)
    if os.path.exists(a2m_file): os.remove(a2m_file)
    print(f"Successfully extracted landmarks to {args.out}")

if __name__ == "__main__":
    main()
