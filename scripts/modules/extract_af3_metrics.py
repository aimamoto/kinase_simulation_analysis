#!/usr/bin/env python3
"""
=============================================================================
EVALUATION GUIDE: ALPHAFOLD 3 INTERACTION METRICS
=============================================================================
1. ipTM (Interface Predicted Template Modeling score):
   - Scale: 0 to 1 (Higher is better).
   - > 0.8: High confidence in the predicted protein-protein interface.
   - 0.6 - 0.8: Moderate confidence; interaction is likely, but structural 
     details/rotamers may be flexible or require orthogonal validation.
   - < 0.6: Low confidence; potential transient interaction or artifact.

2. Inter-chain PAE (Predicted Aligned Error):
   - Scale: 0 to ~31 Å (Lower is better).
   - Evaluates the predicted distance error between residues of Chain X and Y.
   - < 5 Å: Highly confident, rigid interface.
   - 5 - 10 Å: Probable interaction, possible domain flexibility.
   - > 15 Å: Chains are likely non-interacting or randomly oriented in space.
=============================================================================
"""

import os
import json
import numpy as np
import csv
import glob
import argparse
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed

from Bio.PDB import MMCIFParser
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

def get_chain_lengths_from_cif(cif_path, max_chains=4):
    """Extracts the exact number of residues per protein chain for this specific model."""
    parser = MMCIFParser(QUIET=True)
    try:
        structure = parser.get_structure("model", cif_path)
    except Exception:
        return []
        
    lengths = []
    for model in structure:
        for chain in model:
            # Count only standard amino acids (hetero flag is ' ')
            res_count = sum(1 for res in chain if res.id[0] == ' ')
            if res_count > 0:
                lengths.append(res_count)
        break  # only process the first model
    return lengths[:max_chains]

def get_iptm(summary_json_path):
    # Parse ipTM and pTM independently: a single-entity model (e.g. an apo
    # monomer) has no interface, so AF3 writes iptm=null. Coercing that shared
    # with ptm would discard the valid ptm too, so handle each field on its own.
    def _num(v):
        try:
            return round(float(v), 3)
        except (TypeError, ValueError):
            return "N/A"
    try:
        with open(summary_json_path, 'r') as f:
            data = json.load(f)
    except Exception:
        return "N/A", "N/A"
    return _num(data.get('iptm')), _num(data.get('ptm'))

def extract_interchain_pae(confidences_json_path, chain_lengths):
    try:
        with open(confidences_json_path, 'r') as f:
            data = json.load(f)
        if 'pae' not in data:
            return {}
            
        pae_matrix = np.array(data['pae'])
        total_len = sum(chain_lengths)
        
        # STRICT SAFETY CHECK: Prevent empty slices on Monomers or matrix bounds
        if len(chain_lengths) < 2 or pae_matrix.shape[0] < total_len:
            return {}
            
        labels = ['A', 'B', 'C', 'D'][:len(chain_lengths)]
        starts = [0]
        for l in chain_lengths[:-1]:
            starts.append(starts[-1] + l)
            
        pae_results = {}
        for i in range(len(chain_lengths)):
            for j in range(i + 1, len(chain_lengths)):
                start_i, end_i = starts[i], starts[i] + chain_lengths[i]
                start_j, end_j = starts[j], starts[j] + chain_lengths[j]
                
                pae_ij = pae_matrix[start_i:end_i, start_j:end_j]
                pae_ji = pae_matrix[start_j:end_j, start_i:end_i]
                
                if pae_ij.size == 0 or pae_ji.size == 0:
                    continue
                    
                mean_ij = round(float(np.mean(pae_ij)), 2)
                mean_ji = round(float(np.mean(pae_ji)), 2)
                mean_interface = round((mean_ij + mean_ji) / 2.0, 2)
                
                pair_lbl = f"{labels[i]}{labels[j]}"
                pae_results[f"PAE_{labels[i]}_to_{labels[j]}"] = mean_ij
                pae_results[f"PAE_{labels[j]}_to_{labels[i]}"] = mean_ji
                pae_results[f"PAE_Mean_{pair_lbl}"] = mean_interface
                
        return pae_results
    except Exception:
        return {}

def process_single_cif(cif_path, max_chains, expected_pae_keys):
    """Worker function to process a single simulation directory."""
    sim_dir = os.path.dirname(cif_path)
    summary_json = os.path.join(sim_dir, "summary_confidences.json")
    confidences_json = os.path.join(sim_dir, "confidences.json")
    
    iptm, ptm = get_iptm(summary_json) if os.path.exists(summary_json) else ("N/A", "N/A")
    res = {"Directory": os.path.normpath(sim_dir), "ipTM": iptm, "pTM": ptm}
    
    # Initialize expected metric keys to N/A
    for k in expected_pae_keys:
        res[k] = "N/A"
        
    chain_lengths = get_chain_lengths_from_cif(cif_path, max_chains)
    
    if os.path.exists(confidences_json) and len(chain_lengths) >= 2:
        extracted_pae = extract_interchain_pae(confidences_json, chain_lengths)
        for k, v in extracted_pae.items():
            if k in expected_pae_keys:
                res[k] = v
            
    return res

def main():
    parser = argparse.ArgumentParser(description="Extract AF3 Confidence Metrics")
    parser.add_argument("-d", "--dir", default=".", help="Root directory")
    parser.add_argument("-n", "--max-chains", type=int, default=2, help="Max chains to evaluate")
    parser.add_argument("-o", "--out", default="temp_af3_metrics.csv", help="Temp Output CSV file")
    parser.add_argument("-c", "--cores", type=int, default=None, help="Number of cores to use")
    args = parser.parse_args()

    # Generate expected PAE column names dynamically based on the global max chains passed
    labels = ['A', 'B', 'C', 'D'][:args.max_chains]
    expected_pae_keys = []
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            expected_pae_keys.extend([
                f"PAE_{labels[i]}_to_{labels[j]}", 
                f"PAE_{labels[j]}_to_{labels[i]}", 
                f"PAE_Mean_{labels[i]}{labels[j]}"
            ])

    cif_files = glob.glob(os.path.join(args.dir, "**", "model.cif"), recursive=True)
    # Skip any path under an 'archive*' subdirectory (case-insensitive)
    cif_files = [f for f in cif_files
                 if not any(p.lower().startswith('archive') for p in f.replace('\\', '/').split('/'))]
    if not cif_files:
        print("❌ No 'model.cif' files found.")
        return

    max_workers = args.cores if args.cores else min(8, multiprocessing.cpu_count())
    print(f"[*] Extracting metrics from {len(cif_files)} files using {max_workers} parallel workers...")

    csv_data = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_single_cif, cif, args.max_chains, expected_pae_keys): cif for cif in cif_files}
        for future in as_completed(futures):
            try:
                res = future.result()
                csv_data.append(res)
            except Exception as e:
                print(f"Error processing a file: {e}")

    fieldnames = ["Directory", "ipTM", "pTM"] + expected_pae_keys
    with open(args.out, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(csv_data)

if __name__ == "__main__":
    main()
