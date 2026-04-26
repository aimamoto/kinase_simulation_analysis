#!/usr/bin/env python3
import os
import json
import numpy as np
import csv
import glob
import argparse

def get_chain_lengths_from_fasta(fasta_path):
    """Dynamically parses the FASTA file to get the exact lengths of the first two chains."""
    lengths =[]
    current_len = 0
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_len > 0:
                    lengths.append(current_len)
                    current_len = 0
            elif line:
                current_len += len(line)
        if current_len > 0:
            lengths.append(current_len)
    
    if len(lengths) < 2:
        print("⚠️ Warning: Less than 2 chains found in FASTA. Defaulting PAE slice to 300 aa.")
        return 300, 300
    return lengths[0], lengths[1]

def get_iptm(summary_json_path):
    try:
        with open(summary_json_path, 'r') as f:
            data = json.load(f)
        return round(float(data.get('iptm', 0.0)), 3), round(float(data.get('ptm', 0.0)), 3)
    except Exception:
        return "N/A", "N/A"

def extract_interchain_pae(confidences_json_path, lenA, lenB):
    try:
        with open(confidences_json_path, 'r') as f:
            data = json.load(f)
        if 'pae' not in data:
            return "N/A", "N/A"
            
        pae_matrix = np.array(data['pae'])
        pae_AB = pae_matrix[0:lenA, lenA : lenA + lenB]
        pae_BA = pae_matrix[lenA : lenA + lenB, 0:lenA]
        return round(float(np.mean(pae_AB)), 2), round(float(np.mean(pae_BA)), 2)
    except Exception:
        return "N/A", "N/A"

def main():
    parser = argparse.ArgumentParser(description="Extract AF3 Confidence Metrics")
    parser.add_argument("-d", "--dir", default=".", help="Root directory")
    parser.add_argument("-f", "--fasta", default="sequences.fasta", help="FASTA to determine chain lengths")
    parser.add_argument("-o", "--out", default="temp_af3_metrics.csv", help="Temp Output CSV file")
    args = parser.parse_args()

    lenA, lenB = get_chain_lengths_from_fasta(args.fasta)
    print(f"[*] Detected Chain A length: {lenA}, Chain B length: {lenB}")

    cif_files = glob.glob(os.path.join(args.dir, "**", "model.cif"), recursive=True)
    if not cif_files:
        print("❌ No 'model.cif' files found.")
        return

    csv_data =[]
    for cif_path in cif_files:
        sim_dir = os.path.dirname(cif_path)
        summary_json = os.path.join(sim_dir, "summary_confidences.json")
        confidences_json = os.path.join(sim_dir, "confidences.json")
        
        iptm, ptm = get_iptm(summary_json) if os.path.exists(summary_json) else ("N/A", "N/A")
        
        if os.path.exists(confidences_json):
            pae_ab, pae_ba = extract_interchain_pae(confidences_json, lenA, lenB)
            pae_mean = round((pae_ab + pae_ba) / 2.0, 2) if pae_ab != "N/A" else "N/A"
        else:
            pae_ab, pae_ba, pae_mean = "N/A", "N/A", "N/A"
            
        csv_data.append({
            "Directory": os.path.normpath(sim_dir),
            "ipTM": iptm, "pTM": ptm,
            "PAE_ChainA_to_ChainB": pae_ab, "PAE_ChainB_to_ChainA": pae_ba, "PAE_Mean_Interface": pae_mean
        })

    with open(args.out, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["Directory", "ipTM", "pTM", "PAE_ChainA_to_ChainB", "PAE_ChainB_to_ChainA", "PAE_Mean_Interface"])
        writer.writeheader()
        writer.writerows(csv_data)

if __name__ == "__main__":
    main()
