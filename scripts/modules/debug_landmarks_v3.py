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
    # Expanded header to include the hydrophobic shell and bridge residues
    header = f"{'Sequence Name':<30} | {'VAIK (K)':<8} | {'aC (E)':<8} | {'V104':<6} | {'M118':<6} | {'M120':<6} | {'I150':<6} | {'HRD (H)':<8} | {'DFG (F)':<8} | {'APE (E)':<8}"
    print(header)
    print("-" * len(header))

    for name, lm in landmarks.items():
        seq = sequences.get(name)
        if not seq: continue
        
        try:
            # Grab the actual amino acid letter at the mapped index
            k_res = f"{seq[lm['k']]}{lm['k']}" if lm.get('k') is not None else "None"
            c_res = f"{seq[lm['c']]}{lm['c']}" if lm.get('c') is not None else "None"
            
            # New Shell and Bridge residues
            v104_res = f"{seq[lm['v104']]}{lm['v104']}" if lm.get('v104') is not None else "None"
            m118_res = f"{seq[lm['m118']]}{lm['m118']}" if lm.get('m118') is not None else "None"
            m120_res = f"{seq[lm['m120']]}{lm['m120']}" if lm.get('m120') is not None else "None"
            i150_res = f"{seq[lm['i150']]}{lm['i150']}" if lm.get('i150') is not None else "None"
            
            hrd_res = f"{seq[lm['hrd']]}{lm['hrd']}" if lm.get('hrd') is not None else "None"
            dfg_res = f"{seq[lm['f']]}{lm['f']}" if lm.get('f') is not None else "None"
            ape_res = f"{seq[lm['ape']]}{lm['ape']}" if lm.get('ape') is not None else "None"
            
            print(f"{name:<30} | {k_res:<8} | {c_res:<8} | {v104_res:<6} | {m118_res:<6} | {m120_res:<6} | {i150_res:<6} | {hrd_res:<8} | {dfg_res:<8} | {ape_res:<8}")
        except IndexError:
            print(f"{name:<30} | ERROR: Index out of bounds")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 debug_landmarks.py <input.fasta>")
        sys.exit(1)
    debug_landmarks(sys.argv[1])
