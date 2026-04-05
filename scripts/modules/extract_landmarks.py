#!/usr/bin/env python3
"""
extract_landmarks.py
Generate HMM-based kinase landmark positions from protein FASTA sequences.
"""
import os
import re
import json
import glob
import argparse
import subprocess
import sys

# Version 4 Upgrade

DEFAULT_HMM = os.path.expanduser("~/pfam/Pkinase.hmm")
DEFAULT_OUTFILE = "hmm_landmarks.json"
TEMP_A2M = "_temp_aligned.a2m"

# ================================================================
#  CANONICAL & ALLOSTERIC KINASE LANDMARKS (1-based HMM indexing)
# ================================================================
TARGET_NODES = {
    30:  "k",      # VAIK Lysine
    48:  "c",      # alphaC-helix Glutamate
    52:  "rs1",    # R-spine 1 (PKA: L95)
    56:  "n99",    # aC-b4 Loop Anchor (PKA: N99)
    61:  "v104",   # Shell / R-C spine bridge (PKA: V104)
    62:  "k105",   # aC-b4 Toggle Switch (PKA: K105)
    63:  "rs2",    # R-spine 2 (PKA: L106)
    64:  "e107",   # K105 Toggle Target 1 (PKA: E107)
    75:  "m118",   # Shell (PKA: M118)
    77:  "m120",   # Shell (PKA: M120)
    78:  "e121",   # K105 Toggle Target 2 (PKA: E121)
    113: "i150",   # Bridge to alphaE (PKA: I150)
    119: "y156",   # alphaE Anchor to N99 (PKA: Y156)
    121: "hrd",    # HRD Histidine
    142: "f",      # DFG Phenylalanine
    167: "ape",    # APE Glutamate
    179: "d220"    # alphaF Scaffold / Cat-loop Anchor (PKA: D220)
}

def extract_nodes_from_aligned_seq(a2m_seq: str, fasta_name: str):
    seq_index = 0  
    hmm_index = 1  
    mapping = {v: None for v in TARGET_NODES.values()}

    for ch in a2m_seq:
        if ch.isalpha():
            if ch.isupper():
                if hmm_index in TARGET_NODES:
                    mapping[TARGET_NODES[hmm_index]] = seq_index
                seq_index += 1
                hmm_index += 1
            else:
                seq_index += 1
        elif ch == "-":
            hmm_index += 1

    base = re.split(r"[-_]", fasta_name)[0].upper()
    mapping["type"] = base
    return mapping

def run_hmmalign(hmm_path: str, fasta_path: str) -> str:
    if not os.path.exists(hmm_path):
        print(f"ERROR: HMM profile not found at: {hmm_path}"); sys.exit(1)
    cmd = ["hmmalign", "--outformat", "A2M", "-o", TEMP_A2M, hmm_path, fasta_path]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return TEMP_A2M

def parse_a2m(a2m_file: str) -> dict:
    landmarks = {}
    current_name = None
    current_seq_chunks = []
    with open(a2m_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                if current_name is not None:
                    landmarks[current_name] = extract_nodes_from_aligned_seq("".join(current_seq_chunks), current_name)
                current_name = line[1:].strip().split()[0]
                current_seq_chunks = []
            else:
                current_seq_chunks.append(line.strip())
        if current_name is not None:
            landmarks[current_name] = extract_nodes_from_aligned_seq("".join(current_seq_chunks), current_name)
    return landmarks

def find_fasta_file(user_specified: str = None) -> str:
    if user_specified: return user_specified
    fastas = glob.glob("*.fasta") + glob.glob("*.fa") + glob.glob("*.faa")
    if len(fastas) == 1: return fastas[0]
    return fastas[int(input("Select the FASTA file number to use: ")) - 1]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", type=str)
    parser.add_argument("--hmm", type=str, default=DEFAULT_HMM)
    parser.add_argument("--out", type=str, default=DEFAULT_OUTFILE)
    args = parser.parse_args()

    fasta_path = find_fasta_file(args.fasta)
    a2m_file = run_hmmalign(args.hmm, fasta_path)
    landmarks = parse_a2m(a2m_file)

    with open(args.out, "w") as f:
        json.dump(landmarks, f, indent=2)
    if os.path.exists(a2m_file): os.remove(a2m_file)
    print(f"Successfully extracted v4 landmarks to {args.out}")

if __name__ == "__main__":
    main()
