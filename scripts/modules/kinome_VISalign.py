"""
=============================================================================
Kinome Structural Pipeline: Schematics & MSA Panels
=============================================================================
Description:
1. Generates abstract 1D topological schematics for each input kinase.
2. Generates a wrapped, annotated Multiple Sequence Alignment (MSA) panel 
   showing all sequences with structural annotations mapped to the columns.

Execution (In Pipeline):
    python3 kinome_VISalign.py -i sequences.fasta -l hmm_landmarks.json
=============================================================================
"""

import re
import os
import sys
import json
import argparse
import subprocess
import math
import matplotlib
matplotlib.use('Agg') # Headless mode
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path

# =============================================================================
# CONFIG: Macro-Structures & Dimensions
# =============================================================================
DEFAULT_DOMAIN_LENGTH = 310
LOBE_GRAY = "#eeeeee"

LOBES = [
    {"name": "N-Lobe", "start": 1, "end": 116, "color": LOBE_GRAY},
    {"name": "C-Lobe", "start": 126, "end": DEFAULT_DOMAIN_LENGTH, "color": LOBE_GRAY}
]

# =============================================================================
# CONFIG: Topology & Motifs
# =============================================================================
HELIX_COLOR = "#ff4d4d"
SHEET_COLOR = "#66b3ff"

TOPOLOGY = [
    {"name": "αA", "start": 15, "end": 35, "color": HELIX_COLOR, "type": "helix"},
    {"name": "β1", "start": 43, "end": 49, "color": SHEET_COLOR, "type": "sheet"},
    {"name": "β2", "start": 52, "end": 57, "color": SHEET_COLOR, "type": "sheet"},
    {"name": "β3", "start": 65, "end": 74, "color": SHEET_COLOR, "type": "sheet"},
    {"name": "αB", "start": 76, "end": 80, "color": HELIX_COLOR, "type": "helix"},
    {"name": "αC", "start": 84, "end": 97, "color": HELIX_COLOR, "type": "helix"},
    {"name": "β4", "start": 104, "end": 108, "color": SHEET_COLOR, "type": "sheet"},
    {"name": "β5", "start": 111, "end": 115, "color": SHEET_COLOR, "type": "sheet"},
    {"name": "αD", "start": 127, "end": 135, "color": HELIX_COLOR, "type": "helix"},
    {"name": "αE", "start": 138, "end": 160, "color": HELIX_COLOR, "type": "helix"},
    {"name": "β6", "start": 170, "end": 174, "color": SHEET_COLOR, "type": "sheet"},
    {"name": "β7", "start": 179, "end": 183, "color": SHEET_COLOR, "type": "sheet"},
    {"name": "β8", "start": 188, "end": 193, "color": SHEET_COLOR, "type": "sheet"},
    {"name": "β9", "start": 200, "end": 205, "color": SHEET_COLOR, "type": "sheet"},
    {"name": "αF", "start": 216, "end": 236, "color": HELIX_COLOR, "type": "helix"},
    {"name": "αG", "start": 245, "end": 255, "color": HELIX_COLOR, "type": "helix"},
    {"name": "αH", "start": 261, "end": 273, "color": HELIX_COLOR, "type": "helix"},
    {"name": "αI", "start": 282, "end": 286, "color": HELIX_COLOR, "type": "helix"},
    {"name": "αJ", "start": 288, "end": 294, "color": HELIX_COLOR, "type": "helix"}
]

MOTIFS = [
    {"name": "Gly-rich loop", "start": 47, "end": 55, "color": "#3498db"},
    {"name": "αC-Helix (Reg)", "start": 84, "end": 97, "color": "#9b59b6"},
    {"name": "Hinge", "start": 118, "end": 125, "color": "#f1c40f"},
    {"name": "Catalytic Loop", "start": 161, "end": 170, "color": "#e74c3c"},
    {"name": "Activation Loop", "start": 184, "end": 208, "color": "#2ecc71"},
    {"name": "αF-Helix (Core)", "start": 216, "end": 236, "color": "#34495e"}
]

# Canonical Nodes from HMM extraction to map visibly
KEY_LANDMARKS = {
    "k": {"label": "K(VAIK)", "color": "#8e44ad"},
    "c": {"label": "E(αC)", "color": "#e67e22"},
    "hrd": {"label": "HRD", "color": "#d35400"},
    "f": {"label": "DFG", "color": "#c0392b"},
    "ape": {"label": "APE", "color": "#27ae60"}
}

# =============================================================================
# Subprocess: HMM Alignment
# =============================================================================
def run_hmmalign(fasta_path, hmm_path):
    if not os.path.isfile(fasta_path):
        sys.exit(f"Error: Input file '{fasta_path}' not found.")
    
    cmd = ["hmmalign", "--trim", hmm_path, fasta_path]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except FileNotFoundError:
        sys.exit("Error: 'hmmalign' not found in PATH.")
    except subprocess.CalledProcessError as e:
        sys.exit(f"Error: hmmalign failed.\n{e.stderr}")

    seqs = {}
    for line in result.stdout.splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line == "//": continue
        parts = line.split()
        if len(parts) == 2:
            name, seq = parts
            if name not in seqs: seqs[name] = ""
            seqs[name] += seq
    return seqs

# =============================================================================
# Hybrid Coordinate Mapper (HMM JSON -> MSA Columns)
# =============================================================================
def map_coordinates_to_alignment(aligned_seqs, fasta_path, json_path=None, reference_id="PKA"):
    """
    Translates structural coordinates to MSA columns.
    Calculates the N-terminal trimming offset to perfectly sync JSON mapping to the trimmed MSA.
    """
    # 1. Read Raw FASTA to determine how much hmmalign trimmed
    full_seqs = {}
    with open(fasta_path, 'r') as f:
        name = ""
        for line in f:
            if line.startswith(">"): 
                name = line[1:].strip().split()[0]
            else: 
                full_seqs[name] = full_seqs.get(name, "") + line.strip()

    # 2. Identify the Reference Sequence
    ref_id = next((k for k in aligned_seqs.keys() if reference_id.lower() in k.lower()), None)
    if not ref_id:
        ref_id = list(aligned_seqs.keys())[0]
        print(f"Warning: Reference '{reference_id}' missing! Defaulting to '{ref_id}'.")
    
    ref_seq = aligned_seqs[ref_id]
    
    # Map gapped columns to true ungapped sequence strings
    ungapped_to_gapped = []
    ungapped_seq = ""
    for col_idx, char in enumerate(ref_seq):
        if char.isalpha():  
            ungapped_to_gapped.append(col_idx)
            ungapped_seq += char.upper()

    # 3. Calculate Trimming Offset (Crucial for JSON alignment)
    offset = 0
    full_ref = next((seq for name, seq in full_seqs.items() if ref_id in name), None)
    if full_ref:
        offset = full_ref.upper().find(ungapped_seq)
        if offset == -1: 
            offset = 0 # Fallback if mismatch occurs

    mapped_landmarks = {}
    canonical_offset = 0
    dfg_ungapped_idx = None
    ape_ungapped_idx = None

    # --- Mode A: HMM Landmarks JSON ---
    if json_path and os.path.isfile(json_path):
        print(f"Loading HMM landmarks from {json_path} (Correcting for --trim offset: {offset} residues)...")
        with open(json_path, 'r') as f:
            lm_data = json.load(f)
        
        json_key = next((k for k in lm_data.keys() if k in ref_id or ref_id in k), None)
        if json_key:
            for node, ungapped_idx in lm_data[json_key].items():
                if isinstance(ungapped_idx, int):
                    adjusted_idx = ungapped_idx - offset
                    if 0 <= adjusted_idx < len(ungapped_to_gapped):
                        mapped_landmarks[node] = ungapped_to_gapped[adjusted_idx]
                        if node == "f": dfg_ungapped_idx = adjusted_idx
                        if node == "ape": ape_ungapped_idx = adjusted_idx

            if dfg_ungapped_idx is not None:
                canonical_offset = 185 - dfg_ungapped_idx  

    # --- Mode B: Regex Fallback ---
    else:
        print("Warning: hmm_landmarks.json not found. Falling back to Regex mapping.")
        dfg_match = re.search(r'D[FWL]G', ungapped_seq)
        ape_match = re.search(r'.APE', ungapped_seq)
        
        if dfg_match:
            dfg_ungapped_idx = dfg_match.start() + 1
            mapped_landmarks["f"] = ungapped_to_gapped[dfg_ungapped_idx]
            canonical_offset = 185 - dfg_ungapped_idx
            
        if ape_match:
            ape_ungapped_idx = ape_match.end() - 1
            mapped_landmarks["ape"] = ungapped_to_gapped[ape_ungapped_idx]

    # Apply Canonical Topology Mapping
    residue_to_column = { (i + canonical_offset): col_idx for i, col_idx in enumerate(ungapped_to_gapped) }

    def get_col(canonical_res):
        if canonical_res in residue_to_column: return residue_to_column[canonical_res]
        closest_res = min(residue_to_column.keys(), key=lambda k: abs(k - canonical_res))
        return residue_to_column[closest_res]

    for f in TOPOLOGY:
        f["start"], f["end"] = get_col(f["start"]), get_col(f["end"])

    for m in MOTIFS:
        if "Activation Loop" in m["name"] and dfg_ungapped_idx and ape_ungapped_idx:
            m["start"] = ungapped_to_gapped[max(0, dfg_ungapped_idx - 1)]
            m["end"] = ungapped_to_gapped[ape_ungapped_idx]
        else:
            m["start"], m["end"] = get_col(m["start"]), get_col(m["end"])
            
    print("Successfully mapped alignment coordinates.")
    return mapped_landmarks

# =============================================================================
# Visualizer 1: Abstract Schematic 
# =============================================================================
def draw_helix(ax, x, y, w, h, c):
    ax.add_patch(patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.15", ec='#cc0000', fc=c, zorder=2))

def draw_sheet(ax, x, y, w, h, c):
    hw = min(w * 0.4, 3.0) 
    p = [(Path.MOVETO, (x, y+h*.2)), (Path.LINETO, (x+w-hw, y+h*.2)), (Path.LINETO, (x+w-hw, y)), 
         (Path.LINETO, (x+w, y+h/2)), (Path.LINETO, (x+w-hw, y+h)), (Path.LINETO, (x+w-hw, y+h*.8)), 
         (Path.LINETO, (x, y+h*.8)), (Path.CLOSEPOLY, (x, y+h*.2))]
    cd, vt = zip(*p)
    ax.add_patch(patches.PathPatch(Path(vt, cd), fc=c, ec='#0059b3', zorder=2))

def plot_abstract_schematic(domain_len, seq_id, out_file, landmarks_dict):
    fig, ax = plt.subplots(figsize=(16, 5)) 
    
    # Lobes
    for lobe in LOBES:
        end_pos = domain_len if lobe["name"] == "C-Lobe" and lobe["end"] == DEFAULT_DOMAIN_LENGTH else lobe["end"]
        ax.add_patch(patches.Rectangle((lobe["start"], -1.3), end_pos - lobe["start"], 2.6, lw=0, fc=lobe["color"], zorder=0, alpha=0.8))
        ax.text((lobe["start"] + end_pos)/2, 1.25, lobe["name"], ha='center', va='bottom', fontsize=14, fontweight='bold', color='gray')

    ax.hlines(0, 0, domain_len, color='black', lw=1.5, zorder=1)
    ax.text(-2, 0, seq_id, ha='right', va='center', fontsize=12, fontweight='bold', clip_on=False)

    # Topology
    for f in TOPOLOGY:
        if f["start"] > domain_len: continue
        w = f["end"] - f["start"]
        draw_helix(ax, f["start"], 0.1, w, 0.5, f["color"]) if f["type"] == "helix" else draw_sheet(ax, f["start"], 0.1, w, 0.5, f["color"])
        ax.text((f["start"] + f["end"])/2, 0.7, f["name"], ha='center', va='bottom', fontsize=10, fontweight='bold')

    # Motifs
    stagger = [-0.6, -0.9] * 3
    for i, m in enumerate(MOTIFS):
        if m["start"] > domain_len: continue
        w = m["end"] - m["start"]
        ax.add_patch(patches.Rectangle((m["start"], -0.4), w, 0.3, ec='black', fc=m["color"], zorder=2, alpha=0.9))
        ax.plot([(m["start"]+m["end"])/2]*2, [-0.4, stagger[i]], color='gray', ls=':', zorder=1)
        ax.text((m["start"]+m["end"])/2, stagger[i]-0.05, m["name"], ha='center', va='top', fontsize=10)

    # Plot Canonical Landmarks
    for node, col in landmarks_dict.items():
        if node in KEY_LANDMARKS and col <= domain_len:
            info = KEY_LANDMARKS[node]
            ax.axvline(x=col, ymin=0.25, ymax=0.55, color=info["color"], ls="--", lw=1.5, zorder=1)
            ax.text(col, 0.80, info["label"], color=info["color"], fontsize=10, fontweight='bold', ha='center', va='bottom', rotation=35)

    ax.set(xlim=(-5, domain_len + 10), ylim=(-1.5, 1.6), yticks=[])
    ax.spines[['top','right','left']].set_visible(False)
    plt.tight_layout(rect=[0.1, 0, 1, 1]) 
    plt.savefig(os.path.abspath(out_file), format="pdf", bbox_inches='tight')
    plt.close(fig)
    
# =============================================================================
# Visualizer 2: Annotated MSA Text Panel
# =============================================================================
def plot_msa_panel(aligned_seqs, domain_len, out_file, landmarks_dict):
    cols_per_row = 80
    num_seqs = len(aligned_seqs)
    num_chunks = math.ceil(domain_len / cols_per_row)
    
    lines_per_chunk = num_seqs + 5  
    fig_height = max(5, (num_chunks * lines_per_chunk) * 0.3)
    
    fig, ax = plt.subplots(figsize=(18, fig_height))
    ax.set_xlim(-15, cols_per_row + 2) 
    ax.set_ylim(-(num_chunks * lines_per_chunk), 2)
    
    y_cursor = 0
    
    for chunk in range(num_chunks):
        start_col = chunk * cols_per_row
        end_col = min(start_col + cols_per_row, domain_len)
        
        # --- 1. Draw Landmarks Above Block ---
        for node, col in landmarks_dict.items():
            if node in KEY_LANDMARKS and start_col <= col < end_col:
                local_col = col - start_col
                info = KEY_LANDMARKS[node]
                ax.text(local_col, y_cursor + 0.8, f"▼\n{info['label']}", 
                        color=info["color"], ha='center', va='bottom', fontsize=8, fontweight='bold', linespacing=0.8)

        # --- 2. Draw Topology Shapes Above Alignment ---
        for f in TOPOLOGY:
            if f["end"] >= start_col and f["start"] < end_col:
                local_start_idx = max(0, f["start"] - start_col)
                local_end_idx = min(cols_per_row - 1, f["end"] - start_col)
                w = local_end_idx - local_start_idx + 1
                
                ax.add_patch(patches.Rectangle((local_start_idx - 0.5, y_cursor), w, 0.6, ec='black', fc=f["color"], zorder=2))
                ax.text(local_start_idx + w/2 - 0.5, y_cursor + 0.3, f["name"], ha='center', va='center', fontsize=9, fontweight='bold')
        
        y_cursor -= 1 
        
        # --- 3. Print Monospaced Sequence Alignment ---
        for seq_id, seq in aligned_seqs.items():
            chunk_str = seq[start_col:end_col]
            ax.text(-2, y_cursor, seq_id[:12], ha='right', va='center', fontfamily='monospace', fontsize=11, fontweight='bold')
            
            for x_idx, char in enumerate(chunk_str):
                color = "lightgray" if char == "-" else "black"
                ax.text(x_idx, y_cursor, char, ha='center', va='center', fontfamily='monospace', fontsize=11, color=color)
            y_cursor -= 1
            
        # --- 4. Draw Functional Motifs Below Alignment ---
        for m in MOTIFS:
            if m["end"] >= start_col and m["start"] < end_col:
                local_start_idx = max(0, m["start"] - start_col)
                local_end_idx = min(cols_per_row - 1, m["end"] - start_col)
                w = local_end_idx - local_start_idx + 1 
                
                ax.add_patch(patches.Rectangle((local_start_idx - 0.5, y_cursor - 0.1), w, 0.4, ec='black', fc=m["color"], zorder=2, alpha=0.8))
                ax.text(local_start_idx + w/2 - 0.5, y_cursor - 0.6, m["name"], ha='center', va='top', fontsize=8)
        
        y_cursor -= 3.0 # Gap before the next chunk block begins

    ax.axis('off') 
    plt.title("Annotated Kinome Alignment Panel", fontsize=16, fontweight='bold', pad=20)
    plt.tight_layout()
    plt.savefig(os.path.abspath(out_file), format="pdf", bbox_inches='tight')
    plt.close(fig)
    print(f"Saved MSA Panel: {os.path.abspath(out_file)}")

# =============================================================================
# Execution Block
# =============================================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True, help="Path to input FASTA")
    parser.add_argument("-p", "--hmm", type=str, default=os.path.expanduser("~/pfam/Pkinase.hmm"), help="Path to Pkinase.hmm")
    parser.add_argument("-l", "--landmarks", type=str, default="hmm_landmarks.json", help="Path to hmm_landmarks.json (Output from extract_landmarks.py)")
    args = parser.parse_args()
    
    aligned_seqs = run_hmmalign(args.input, args.hmm)
    alignment_length = len(list(aligned_seqs.values())[0])
    
    mapped_landmarks = map_coordinates_to_alignment(
        aligned_seqs=aligned_seqs, 
        fasta_path=args.input, 
        json_path=args.landmarks, 
        reference_id="PKA"
    )
    
    print(f"\nAligned {len(aligned_seqs)} sequences ({alignment_length} cols).")
    
    # 1. Generate individual abstract schematics (now PDF)
    for seq_id in aligned_seqs.keys():
        safe_name = "".join([c if c.isalnum() else "_" for c in seq_id])
        plot_abstract_schematic(alignment_length, seq_id, f"{safe_name}_schematic.pdf", mapped_landmarks)
        print(f"Saved Schematic: {safe_name}_schematic.pdf")
        
    # 2. Generate the unified MSA text panel (now PDF)
    plot_msa_panel(aligned_seqs, alignment_length, "MSA_Annotated_Panel.pdf", mapped_landmarks)
    print("="*70 + "\nPipeline complete!\n")
