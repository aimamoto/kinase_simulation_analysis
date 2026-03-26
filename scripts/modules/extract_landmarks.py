#!/usr/bin/env python3
"""
extract_landmarks.py

Generate HMM-based kinase landmark positions from protein FASTA sequences.

This script:
  1. Detects (or accepts) a FASTA file
  2. Runs hmmalign using the user's Pkinase.hmm profile
  3. Parses the A2M alignment output
  4. Extracts the canonical kinase landmarks + Hydrophobic Shell/Bridges
  5. Outputs: hmm_landmarks.json

Usage:
    python extract_landmarks.py
    python extract_landmarks.py --fasta my_sequences.fasta
    python extract_landmarks.py --hmm /path/to/Pkinase.hmm
"""

import os
import re
import json
import glob
import argparse
import subprocess
import sys


# ================================================================
#  DEFAULT LOCATIONS
# ================================================================

DEFAULT_HMM = os.path.expanduser("~/pfam/Pkinase.hmm")
DEFAULT_OUTFILE = "hmm_landmarks.json"
TEMP_A2M = "_temp_aligned.a2m"


# ================================================================
#  CANONICAL KINASE LANDMARK POSITIONS (1-based HMM indexing)
#  Updated to include the Hydrophobic Shell and Bridging residues
# ================================================================

TARGET_NODES = {
    30:  "k",      # VAIK Lysine
    48:  "c",      # alphaC-helix Glutamate
    52:  "rs1",    # R-spine 1 (PKA eq: L95)
    61:  "v104",   # Shell / R-C spine bridge (PKA eq: V104)
    63:  "rs2",    # R-spine 2 (PKA eq: L106)
    75:  "m118",   # Shell (PKA eq: M118)
    77:  "m120",   # Shell (PKA eq: M120)
    113: "i150",   # Bridge to alphaE (PKA eq: I150)
    121: "hrd",    # HRD Histidine
    142: "f",      # DFG Phenylalanine
    167: "ape"     # APE Glutamate
}


# ================================================================
#  FUNCTIONS
# ================================================================

def extract_nodes_from_aligned_seq(a2m_seq: str, fasta_name: str):
    """
    Given an aligned A2M sequence string and the FASTA header name,
    map the sequence positions to the canonical kinase landmark indices.
    """

    seq_index = 0  # position in actual protein sequence
    hmm_index = 1  # alignment index (1-based)
    mapping = {v: None for v in TARGET_NODES.values()}

    for ch in a2m_seq:
        if ch.isalpha():
            # Uppercase => aligned AA in HMM
            if ch.isupper():
                if hmm_index in TARGET_NODES:
                    key = TARGET_NODES[hmm_index]
                    mapping[key] = seq_index
                seq_index += 1
                hmm_index += 1
            else:
                # lowercase => insertion
                seq_index += 1
        elif ch == "-":
            hmm_index += 1
        else:
            # unexpected character (should not happen)
            pass

    # Infer protein "type" from FASTA name (strip after _ or - etc.)
    base = re.split(r"[-_]", fasta_name)[0].upper()
    mapping["type"] = base

    return mapping


def run_hmmalign(hmm_path: str, fasta_path: str) -> str:
    """
    Run hmmalign with:
       hmmalign --outformat A2M -o _temp_aligned.a2m HMM FASTA

    Returns the output A2M filename.
    """

    if not os.path.exists(hmm_path):
        print(f"ERROR: HMM profile not found at: {hmm_path}")
        sys.exit(1)

    cmd = [
        "hmmalign",
        "--outformat", "A2M",
        "-o", TEMP_A2M,
        hmm_path,
        fasta_path
    ]

    try:
        subprocess.run(cmd, check=True,
                       stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        print("\nERROR: hmmalign failed — please check your HMMER installation or input FASTA.")
        sys.exit(1)

    return TEMP_A2M


def parse_a2m(a2m_file: str) -> dict:
    """
    Parse the A2M alignment file produced by hmmalign.
    """
    landmarks = {}
    current_name = None
    current_seq_chunks = []

    with open(a2m_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                # Save previous record
                if current_name is not None:
                    full_aln = "".join(current_seq_chunks)
                    landmarks[current_name] = extract_nodes_from_aligned_seq(full_aln, current_name)

                # Start a new one
                current_name = line[1:].strip().split()[0]
                current_seq_chunks = []
            else:
                current_seq_chunks.append(line.strip())

        # Save last entry
        if current_name is not None:
            full_aln = "".join(current_seq_chunks)
            landmarks[current_name] = extract_nodes_from_aligned_seq(full_aln, current_name)

    return landmarks


def find_fasta_file(user_specified: str = None) -> str:
    """
    Return the FASTA file path.
    """

    if user_specified:
        if not os.path.exists(user_specified):
            print(f"ERROR: Provided FASTA file '{user_specified}' does not exist.")
            sys.exit(1)
        return user_specified

    fastas = glob.glob("*.fasta") + glob.glob("*.fa") + glob.glob("*.faa")

    if not fastas:
        print("ERROR: No FASTA file found in working directory.")
        sys.exit(1)

    if len(fastas) == 1:
        print(f"Using FASTA file: {fastas[0]}")
        return fastas[0]

    # Multiple FASTA files found — ask user which one
    print("\nMultiple FASTA files found:")
    for i, f in enumerate(fastas, 1):
        print(f"  {i}. {f}")
    choice = input("Select the FASTA file number to use: ")

    try:
        idx = int(choice) - 1
        assert 0 <= idx < len(fastas)
        return fastas[idx]
    except:
        print("Invalid selection.")
        sys.exit(1)


# ================================================================
#  MAIN
# ================================================================

def main():
    parser = argparse.ArgumentParser(description="Extract kinase HMM landmarks from FASTA.")
    parser.add_argument("--fasta", type=str,
                        help="Path to FASTA file. If omitted, auto-detect.")
    parser.add_argument("--hmm", type=str,
                        default=DEFAULT_HMM,
                        help=f"Path to Pkinase.hmm (default: {DEFAULT_HMM})")
    parser.add_argument("--out", type=str,
                        default=DEFAULT_OUTFILE,
                        help=f"Output JSON file (default: {DEFAULT_OUTFILE})")

    args = parser.parse_args()

    fasta_path = find_fasta_file(args.fasta)

    print("\nRunning hmmalign...")
    a2m_file = run_hmmalign(args.hmm, fasta_path)

    print("Parsing A2M alignment...")
    landmarks = parse_a2m(a2m_file)

    # Write JSON
    with open(args.out, "w") as f:
        json.dump(landmarks, f, indent=2)

    # Cleanup
    if os.path.exists(a2m_file):
        os.remove(a2m_file)

    print(f"\nSuccessfully extracted landmarks for {len(landmarks)} sequences.")
    print(f"Saved to: {args.out}\n")


if __name__ == "__main__":
    main()
