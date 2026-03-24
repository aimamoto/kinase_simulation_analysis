import json
import sys

## USAGE
# python3 debug_landmarks.py SRC_and_related_kinases_cat.fasta

def debug_landmarks(fasta_file, json_file="hmm_landmarks.json"):
    # 1. Load the FASTA sequences
    sequences = {}
    current_name = None
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                current_name = line[1:].strip().split()[0]
                sequences[current_name] = ""
            elif current_name:
                sequences[current_name] += line.strip()

    # 2. Load the HMM landmarks
    try:
        with open(json_file, 'r') as f:
            landmarks = json.load(f)
    except FileNotFoundError:
        print("Error: hmm_landmarks.json not found.")
        return

    # 3. Print the mapped residues
    print(f"{'Sequence Name':<30} | {'VAIK (K)':<10} | {'aC (E)':<10} | {'HRD (H)':<10} | {'DFG (F)':<10} | {'APE (E)':<10}")
    print("-" * 90)

    for name, lm in landmarks.items():
        seq = sequences.get(name)
        if not seq: continue
        
        try:
            # Grab the actual amino acid letter at the mapped index
            k_res = f"{seq[lm['k']]}{lm['k']}" if lm.get('k') is not None else "None"
            c_res = f"{seq[lm['c']]}{lm['c']}" if lm.get('c') is not None else "None"
            hrd_res = f"{seq[lm['hrd']]}{lm['hrd']}" if lm.get('hrd') is not None else "None"
            dfg_res = f"{seq[lm['f']]}{lm['f']}" if lm.get('f') is not None else "None"
            ape_res = f"{seq[lm['ape']]}{lm['ape']}" if lm.get('ape') is not None else "None"
            
            print(f"{name:<30} | {k_res:<10} | {c_res:<10} | {hrd_res:<10} | {dfg_res:<10} | {ape_res:<10}")
        except IndexError:
            print(f"{name:<30} | ERROR: Index out of bounds")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 debug_landmarks.py <input.fasta>")
        sys.exit(1)
    debug_landmarks(sys.argv[1])
