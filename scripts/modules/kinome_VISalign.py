import re
import os
import sys
import json
import argparse
import subprocess
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path

DEFAULT_DOMAIN_LENGTH = 310
LOBES = [{"name": "N-Lobe", "start": 1, "end": 116, "color": "#eeeeee"}, {"name": "C-Lobe", "start": 126, "end": DEFAULT_DOMAIN_LENGTH, "color": "#eeeeee"}]
TOPOLOGY = [
    {"name": "αA", "start": 15, "end": 35, "color": "#ff4d4d", "type": "helix"},
    {"name": "β1", "start": 43, "end": 49, "color": "#66b3ff", "type": "sheet"},
    {"name": "β2", "start": 52, "end": 57, "color": "#66b3ff", "type": "sheet"},
    {"name": "β3", "start": 65, "end": 74, "color": "#66b3ff", "type": "sheet"},
    {"name": "αB", "start": 76, "end": 80, "color": "#ff4d4d", "type": "helix"},
    {"name": "αC", "start": 84, "end": 97, "color": "#ff4d4d", "type": "helix"},
    {"name": "β4", "start": 104, "end": 108, "color": "#66b3ff", "type": "sheet"},
    {"name": "β5", "start": 111, "end": 115, "color": "#66b3ff", "type": "sheet"},
    {"name": "αD", "start": 127, "end": 135, "color": "#ff4d4d", "type": "helix"},
    {"name": "αE", "start": 138, "end": 160, "color": "#ff4d4d", "type": "helix"},
    {"name": "β6", "start": 170, "end": 174, "color": "#66b3ff", "type": "sheet"},
    {"name": "β7", "start": 179, "end": 183, "color": "#66b3ff", "type": "sheet"},
    {"name": "β8", "start": 188, "end": 193, "color": "#66b3ff", "type": "sheet"},
    {"name": "β9", "start": 200, "end": 205, "color": "#66b3ff", "type": "sheet"},
    {"name": "αF", "start": 216, "end": 236, "color": "#ff4d4d", "type": "helix"},
    {"name": "αG", "start": 245, "end": 255, "color": "#ff4d4d", "type": "helix"},
    {"name": "αH", "start": 261, "end": 273, "color": "#ff4d4d", "type": "helix"},
    {"name": "αI", "start": 282, "end": 286, "color": "#ff4d4d", "type": "helix"},
    {"name": "αJ", "start": 288, "end": 294, "color": "#ff4d4d", "type": "helix"}
]
MOTIFS = [
    {"name": "Gly-rich loop", "start": 47, "end": 55, "color": "#3498db"},
    {"name": "αC-Helix (Reg)", "start": 84, "end": 97, "color": "#9b59b6"},
    {"name": "Hinge", "start": 118, "end": 125, "color": "#f1c40f"},
    {"name": "Catalytic Loop", "start": 161, "end": 170, "color": "#e74c3c"},
    {"name": "Activation Loop", "start": 184, "end": 208, "color": "#2ecc71"},
    {"name": "αF-Helix (Core)", "start": 216, "end": 236, "color": "#34495e"}
]
KEY_LANDMARKS = {"k": {"label": "K(VAIK)", "color": "#8e44ad"}, "c": {"label": "E(αC)", "color": "#e67e22"}, "hrd": {"label": "HRD", "color": "#d35400"}, "f": {"label": "DFG", "color": "#c0392b"}, "ape": {"label": "APE", "color": "#27ae60"}}

def clean_display_name(raw_name):
    name = raw_name.split('/')[0]
    name = re.sub(r'^\d+[_\-|]', '', name)
    name = re.sub(r'[-_](apo|holo|py\d+|\d*atp)$', '', name, flags=re.IGNORECASE)
    name = re.sub(r'(wt|cattail|cat|tail)', '', name, flags=re.IGNORECASE)
    name = re.sub(r'[-_]+', '-', name).strip('-')
    return name.upper() if name else raw_name.upper()

def run_hmmalign(fasta_path, hmm_path):
    cmd = ["hmmalign", "--trim", hmm_path, fasta_path]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    seqs = {}
    for line in result.stdout.splitlines():
        if not line.strip() or line.startswith("#") or line == "//": continue
        parts = line.split()
        if len(parts) == 2:
            name, seq = parts
            seqs[name] = seqs.get(name, "") + seq
    return seqs

def map_coordinates_to_alignment(aligned_seqs, fasta_path, json_path=None, reference_id="PKA"):
    full_seqs = {}
    with open(fasta_path, 'r') as f:
        name = ""
        for line in f:
            if line.startswith(">"): name = line[1:].strip().split()[0]
            else: full_seqs[name] = full_seqs.get(name, "") + line.strip()

    ref_id = next((k for k in aligned_seqs.keys() if reference_id.lower() in k.lower()), list(aligned_seqs.keys())[0])
    ref_seq = aligned_seqs[ref_id]
    
    ungapped_to_gapped = []
    ungapped_seq = ""
    for col_idx, char in enumerate(ref_seq):
        if char.isalpha():  
            ungapped_to_gapped.append(col_idx)
            ungapped_seq += char.upper()

    offset = 0
    full_ref = next((seq for name, seq in full_seqs.items() if ref_id in name), None)
    if full_ref: offset = max(0, full_ref.upper().find(ungapped_seq))

    mapped_landmarks, canonical_offset, dfg_ungapped_idx, ape_ungapped_idx = {}, 0, None, None

    if json_path and os.path.isfile(json_path):
        with open(json_path, 'r') as f: lm_data = json.load(f)
        json_key = next((k for k in lm_data.keys() if k in ref_id or ref_id in k), None)
        if json_key:
            for node, ungapped_idx in lm_data[json_key].items():
                if isinstance(ungapped_idx, int):
                    adjusted_idx = ungapped_idx - offset
                    if 0 <= adjusted_idx < len(ungapped_to_gapped):
                        mapped_landmarks[node] = ungapped_to_gapped[adjusted_idx]
                        if node == "f": dfg_ungapped_idx = adjusted_idx
                        if node == "ape": ape_ungapped_idx = adjusted_idx
            if dfg_ungapped_idx is not None: canonical_offset = 185 - dfg_ungapped_idx  
    else:
        dfg_match, ape_match = re.search(r'D[FWL]G', ungapped_seq), re.search(r'.APE', ungapped_seq)
        if dfg_match:
            dfg_ungapped_idx = dfg_match.start() + 1
            mapped_landmarks["f"] = ungapped_to_gapped[dfg_ungapped_idx]
            canonical_offset = 185 - dfg_ungapped_idx
        if ape_match: mapped_landmarks["ape"] = ungapped_to_gapped[ape_match.end() - 1]

    residue_to_column = { (i + canonical_offset): col_idx for i, col_idx in enumerate(ungapped_to_gapped) }
    def get_col(c_res): return residue_to_column.get(c_res, residue_to_column[min(residue_to_column.keys(), key=lambda k: abs(k - c_res))] if residue_to_column else 0)

    for f in TOPOLOGY: f["start"], f["end"] = get_col(f["start"]), get_col(f["end"])
    for m in MOTIFS:
        if "Activation Loop" in m["name"] and dfg_ungapped_idx and ape_ungapped_idx:
            m["start"], m["end"] = ungapped_to_gapped[max(0, dfg_ungapped_idx - 1)], ungapped_to_gapped[ape_ungapped_idx]
        else: m["start"], m["end"] = get_col(m["start"]), get_col(m["end"])
    return mapped_landmarks

def draw_helix(ax, x, y, w, h, c): ax.add_patch(patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.15", ec='#cc0000', fc=c, zorder=2))
def draw_sheet(ax, x, y, w, h, c):
    hw = min(w * 0.4, 3.0) 
    p = [(Path.MOVETO, (x, y+h*.2)), (Path.LINETO, (x+w-hw, y+h*.2)), (Path.LINETO, (x+w-hw, y)), (Path.LINETO, (x+w, y+h/2)), (Path.LINETO, (x+w-hw, y+h)), (Path.LINETO, (x+w-hw, y+h*.8)), (Path.LINETO, (x, y+h*.8)), (Path.CLOSEPOLY, (x, y+h*.2))]
    cd, vt = zip(*p)
    ax.add_patch(patches.PathPatch(Path(vt, cd), fc=c, ec='#0059b3', zorder=2))

def plot_abstract_schematic(domain_len, display_name, out_file, landmarks_dict):
    fig, ax = plt.subplots(figsize=(16, 5)) 
    for lobe in LOBES:
        ax.add_patch(patches.Rectangle((lobe["start"], -1.3), (domain_len if lobe["name"] == "C-Lobe" else lobe["end"]) - lobe["start"], 2.6, lw=0, fc=lobe["color"], zorder=0, alpha=0.8))
        ax.text((lobe["start"] + (domain_len if lobe["name"] == "C-Lobe" else lobe["end"]))/2, 1.25, lobe["name"], ha='center', va='bottom', fontsize=14, fontweight='bold', color='gray')
    ax.hlines(0, 0, domain_len, color='black', lw=1.5, zorder=1)
    ax.text(-2, 0, display_name, ha='right', va='center', fontsize=14, fontweight='bold', clip_on=False)
    for f in TOPOLOGY:
        if f["start"] > domain_len: continue
        draw_helix(ax, f["start"], 0.1, f["end"]-f["start"], 0.5, f["color"]) if f["type"] == "helix" else draw_sheet(ax, f["start"], 0.1, f["end"]-f["start"], 0.5, f["color"])
        ax.text((f["start"] + f["end"])/2, 0.7, f["name"], ha='center', va='bottom', fontsize=10, fontweight='bold')
    stagger = [-0.6, -0.9] * 3
    for i, m in enumerate(MOTIFS):
        if m["start"] > domain_len: continue
        ax.add_patch(patches.Rectangle((m["start"], -0.4), m["end"]-m["start"], 0.3, ec='black', fc=m["color"], zorder=2, alpha=0.9))
        ax.plot([(m["start"]+m["end"])/2]*2, [-0.4, stagger[i]], color='gray', ls=':', zorder=1)
        ax.text((m["start"]+m["end"])/2, stagger[i]-0.05, m["name"], ha='center', va='top', fontsize=10)
    for node, col in landmarks_dict.items():
        if node in KEY_LANDMARKS and col <= domain_len:
            info = KEY_LANDMARKS[node]
            ax.axvline(x=col, ymin=0.25, ymax=0.55, color=info["color"], ls="--", lw=1.5, zorder=1)
            ax.text(col, 0.80, info["label"], color=info["color"], fontsize=10, fontweight='bold', ha='center', va='bottom', rotation=35)
    ax.set(xlim=(-5, domain_len + 10), ylim=(-1.5, 1.6), yticks=[])
    ax.spines[['top','right','left']].set_visible(False)
    plt.tight_layout(rect=[0.1, 0, 1, 1]); plt.savefig(os.path.abspath(out_file), format="pdf", bbox_inches='tight'); plt.close(fig)

def plot_msa_panel(aligned_seqs, domain_len, out_file, landmarks_dict):
    cols_per_row = 80
    num_chunks = math.ceil(domain_len / cols_per_row)
    fig, ax = plt.subplots(figsize=(18, max(5, (num_chunks * (len(aligned_seqs) + 5)) * 0.3)))
    ax.set_xlim(-18, cols_per_row + 2); ax.set_ylim(-(num_chunks * (len(aligned_seqs) + 5)), 2)
    y_cursor = 0
    for chunk in range(num_chunks):
        start_col, end_col = chunk * cols_per_row, min((chunk + 1) * cols_per_row, domain_len)
        for node, col in landmarks_dict.items():
            if node in KEY_LANDMARKS and start_col <= col < end_col:
                ax.text(col - start_col, y_cursor + 0.8, f"▼\n{KEY_LANDMARKS[node]['label']}", color=KEY_LANDMARKS[node]["color"], ha='center', va='bottom', fontsize=8, fontweight='bold', linespacing=0.8)
        for f in TOPOLOGY:
            if f["end"] >= start_col and f["start"] < end_col:
                ls, le = max(0, f["start"] - start_col), min(cols_per_row - 1, f["end"] - start_col)
                ax.add_patch(patches.Rectangle((ls - 0.5, y_cursor), le - ls + 1, 0.6, ec='black', fc=f["color"], zorder=2))
                ax.text(ls + (le - ls + 1)/2 - 0.5, y_cursor + 0.3, f["name"], ha='center', va='center', fontsize=9, fontweight='bold')
        y_cursor -= 1 
        for seq_id, seq in aligned_seqs.items():
            ax.text(-2, y_cursor, clean_display_name(seq_id)[:18], ha='right', va='center', fontfamily='monospace', fontsize=11, fontweight='bold')
            for x_idx, char in enumerate(seq[start_col:end_col]): ax.text(x_idx, y_cursor, char, ha='center', va='center', fontfamily='monospace', fontsize=11, color="lightgray" if char == "-" else "black")
            y_cursor -= 1
        for m in MOTIFS:
            if m["end"] >= start_col and m["start"] < end_col:
                ls, le = max(0, m["start"] - start_col), min(cols_per_row - 1, m["end"] - start_col)
                ax.add_patch(patches.Rectangle((ls - 0.5, y_cursor - 0.1), le - ls + 1, 0.4, ec='black', fc=m["color"], zorder=2, alpha=0.8))
                ax.text(ls + (le - ls + 1)/2 - 0.5, y_cursor - 0.6, m["name"], ha='center', va='top', fontsize=8)
        y_cursor -= 3.0 
    ax.axis('off'); plt.title("Annotated Kinome Alignment Panel", fontsize=16, fontweight='bold', pad=20); plt.tight_layout(); plt.savefig(os.path.abspath(out_file), format="pdf", bbox_inches='tight'); plt.close(fig)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True)
    parser.add_argument("-p", "--hmm", type=str, default=os.path.expanduser("~/pfam/Pkinase.hmm"))
    parser.add_argument("-l", "--landmarks", type=str, default="hmm_landmarks.json")
    args = parser.parse_args()
    
    aligned_seqs = run_hmmalign(args.input, args.hmm)

    # --- FILTER OUT NON-KINASES (CO-FACTORS) ---
    if os.path.exists(args.landmarks):
        with open(args.landmarks, 'r') as f:
            lm_data = json.load(f)
        
        filtered_seqs = {}
        for seq_id, seq in aligned_seqs.items():
            # Only keep sequences that have the conserved DFG ('f') landmark
            if seq_id in lm_data and lm_data[seq_id].get("f") is not None:
                filtered_seqs[seq_id] = seq
            else:
                print(f"[*] Excluding non-kinase sequence from visualizations: {seq_id}")
        
        aligned_seqs = filtered_seqs

    if not aligned_seqs:
        print("[!] No valid kinases found for visualization. Exiting module.")
        sys.exit(0)

    alignment_length = len(list(aligned_seqs.values())[0])
    mapped_landmarks = map_coordinates_to_alignment(aligned_seqs, args.input, args.landmarks, "PKA")
    
    for seq_id in aligned_seqs.keys():
        safe_filename = "".join([c if c.isalnum() else "_" for c in seq_id])
        plot_abstract_schematic(alignment_length, clean_display_name(seq_id), f"{safe_filename}_schematic.pdf", mapped_landmarks)
        
    plot_msa_panel(aligned_seqs, alignment_length, "MSA_Annotated_Panel.pdf", mapped_landmarks)
