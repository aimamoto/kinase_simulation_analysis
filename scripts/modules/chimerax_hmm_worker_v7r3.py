import os
import csv
import json
import numpy as np
import re
import difflib
from typing import Dict, List, Tuple, Optional, Any
from chimerax.core.commands import run

# ==============================================================================
# UPDATE LOG (v6r5 -> v6r6)
# Date: May 18, 2026
#
# 1. Universal Steric Co-Factor Integration (CDK-Cyclin Support):
#    - Non-kinase chains are no longer discarded. They are preserved as `cofactors`.
#    - Distances from the kinase's αC-helix and Activation Loop are tracked.
#    - CoFactor_Name, CoFactor_aC_Dist, CoFactor_ActLoop_Dist added to CSV.
#
# 2. Co-Factor Role Assignment Fix:
#    - Kinases successfully bound (<8.0 Å) to a non-kinase chain will correctly 
#      display "Co-factor Bound" as their Role, preventing misleading "Unpaired" statuses.
#
# 3. Cosmetic CXC Fix:
#    - Strips 'A-' or 'B-' prefixes natively from the ChimeraX visualization text.
#
# UPDATE LOG (v6r7 -> v7)
# Date: June 4, 2026
#
# 1. Added Cleft Gaping (Roof-to-Floor Distance) metric (Roof CA to HRD-Asp CA)
# 2. Added Mg2+ Hijacking metric (Inhibitory Phosphate O/P to Mg2+ ion)
# 3. Added Substrate Clearance Angle (Roof CA - ATP_PG - HRD-Asp CA)
#
# UPDATE LOG (v7r2 -> v7r3)   *** CHANGES REPORTED VALUES -- NOT A DROP-IN FOR v7r2 ***
# Date: July 29, 2026
# Prompted by a concordance check against Kincore-standalone2 (Dunbrack lab) on PDB 3D7T.
# v7r2 outputs are NOT reproduced by v7r3; re-run any dataset before comparing across versions.
#
# 1. D1_Dist anchor corrected (AFFECTS Spatial + State).
#    Modi & Dunbrack, PNAS 2019 (116:6818) define D1 = dist(aC-Glu(+4)-CA, DFG-Phe-CZ).
#    v7r2 measured from the aC-Glu ITSELF, while still applying the published 11/14 A
#    cutoffs, which are calibrated on the (+4) anchor. That compressed D1's dynamic range,
#    put the 11 A cut mid-distribution, inflated the fall-through "Outlier" class and could
#    report a true DFGout as DFGinter. D2_Dist was already correct and is unchanged.
#    Verified on 3D7T: D1 now 4.85 A (CSK) / 12.81 A (SRC), matching Kincore exactly.
#
# 2. ActLoop_CT anchor corrected. v7r2 used HRD+6, which is not a conserved position
#    (Arg in CSK, Ala in SRC). The C-terminal contact partner is the HRD arginine, HRD+1,
#    as in Modi & Dunbrack's APE9-Arg contact. The 5.5 A all-atom cutoff is retained as
#    AlloQuant's own sensitivity choice (Kincore uses 6.0 A; 8.0 A for the TYR group).
#
# 3. Missing-density guard added for ActLoop_NT / ActLoop_CT (see residues_contiguous).
#    Landmark indices address the RESOLVED-residue list, so a disordered activation loop
#    was invisible to index arithmetic and v7r2 silently measured across the gap. Offsets
#    are now validated against deposited residue numbering and yield N/A when broken.
#    No effect on AlphaFold3 models (complete chains); matters for experimental structures.
# ==============================================================================

# --- CONFIGURATION ---
SEARCH_DIR = "."  
VIZ_OUT_DIR_CORE = "cx_viz_core"
VIZ_OUT_DIR_ALLO = "cx_viz_allosteric"
FILE_EXTENSION = ".cif"
LANDMARKS_JSON = "hmm_landmarks.json"

AA_MAP = {
    'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLN':'Q', 'GLU':'E', 
    'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 
    'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V',
    'PTR':'Y', 'SEP':'S', 'TPO':'T', 'HIE':'H', 'HID':'H', 'HIP':'H', 'LYZ':'K'
}
LIGAND_NAMES = ["ATP", "ADP", "ANP", "ACP", "AGS", "AMP", "GTP", "GDP", "STU"]

# --- MATH & ATOM HELPERS ---
def normalize_sim_name(name: str) -> str:
    clean = re.sub(r'_?(model|seed|rank|pred|unrelaxed|relaxed)[\-_]?\d+', '', name, flags=re.IGNORECASE)
    clean = re.sub(r'_?(unrelaxed|relaxed)', '', clean, flags=re.IGNORECASE)
    clean = re.sub(r'^(fold|job|run)_', '', clean, flags=re.IGNORECASE)
    clean = re.sub(r'[-_](apo|holo|py\d+|\d*atp)$', '', clean, flags=re.IGNORECASE)
    clean = re.sub(r'(wt|cattail|cat|tail)', '', clean, flags=re.IGNORECASE)
    clean = re.sub(r'[-_]+', '-', clean).strip('-')
    return clean.strip('_')

def get_atom(residue: Any, atom_name: str) -> Optional[np.ndarray]:
    if residue is None: return None
    atoms = residue.atoms[residue.atoms.names == atom_name]
    if len(atoms) > 0: return atoms[0].scene_coord
    return None

def get_sidechain_atoms(residue: Any) -> Optional[Any]:
    if residue is None: return None
    bb_names = {'N', 'CA', 'C', 'O', 'OXT'}
    mask = ~np.isin(residue.atoms.names, list(bb_names))
    sc_atoms = residue.atoms[mask]
    if len(sc_atoms) == 0: return residue.atoms[residue.atoms.names == 'CA']
    return sc_atoms

def get_center_of_mass(residues: Any) -> Optional[np.ndarray]:
    if hasattr(residues, 'scene_coords'): coords = residues.scene_coords
    elif hasattr(residues, 'atoms'): coords = residues.atoms.scene_coords
    else: return None
    if len(coords) == 0: return None
    return np.mean(coords, axis=0)

def calculate_dihedral(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray, p4: np.ndarray) -> Optional[float]:
    if any(x is None for x in [p1, p2, p3, p4]): return None
    b0 = -1.0 * (p2 - p1); b1 = p3 - p2; b2 = p4 - p3
    b1_norm = np.linalg.norm(b1)
    if b1_norm == 0: return None
    b1 /= b1_norm
    v = b0 - np.dot(b0, b1) * b1; w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w); y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))

def calculate_angle(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> Optional[float]:
    if any(x is None for x in [p1, p2, p3]): return None
    v1 = p1 - p2
    v2 = p3 - p2
    n1 = np.linalg.norm(v1)
    n2 = np.linalg.norm(v2)
    if n1 == 0 or n2 == 0: return None
    cosine_angle = np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0)
    return float(np.degrees(np.arccos(cosine_angle)))

def get_min_sc_dist(res1: Any, res2: Any) -> float:
    sc1 = get_sidechain_atoms(res1); sc2 = get_sidechain_atoms(res2)
    if sc1 is None or sc2 is None or len(sc1) == 0 or len(sc2) == 0: return 999.9
    diff = sc1.scene_coords[:, np.newaxis, :] - sc2.scene_coords[np.newaxis, :, :]
    return np.min(np.linalg.norm(diff, axis=2))

def get_res_min_dist(res1: Any, res2: Any) -> float:
    if res1 is None or res2 is None: return 999.9
    c1 = res1.atoms.scene_coords
    c2 = res2.atoms.scene_coords
    if len(c1) == 0 or len(c2) == 0: return 999.9
    diff = c1[:, np.newaxis, :] - c2[np.newaxis, :, :]
    return float(np.min(np.linalg.norm(diff, axis=2)))

def get_min_dist(res_group_1: Any, res_group_2: Any) -> float:
    c1 = res_group_1.atoms.scene_coords
    c2 = res_group_2.atoms.scene_coords
    if len(c1) == 0 or len(c2) == 0: return 999.9
    diff = c1[:, np.newaxis, :] - c2[np.newaxis, :, :]
    return float(np.min(np.linalg.norm(diff, axis=2)))

def clean_protein_type(raw_type: str) -> str:
    name = re.sub(r'[_|-]?chain[_|-]?[A-Za-z0-9]+$', '', raw_type, flags=re.IGNORECASE)
    name = re.sub(r'^\d+[_\-|]', '', name)
    name = re.sub(r'^[a-zA-Z]-', '', name) 
    tokens = re.split(r'[-_]', name)
    ignore_tokens = {'wt', 'cat', 'cattail', 'tail', 'wtcat', 'catwt', 'kd', 'apo', 'holo'}
    clean_tokens = []
    for t in tokens:
        if not t: continue
        t_lower = t.lower()
        if t_lower in ignore_tokens: continue
        if re.match(r'^p[sty]?\d+$', t_lower): continue
        if re.match(r'^\d*(atp|adp|amp|gdp|gtp|anp|stu)$', t_lower): continue
        clean_tokens.append(t.upper())
    return "-".join(clean_tokens) if clean_tokens else raw_type.upper()

def is_erbb_system(sim_id: str, chain_types: List[str]) -> bool:
    search_string = (sim_id + " " + " ".join(chain_types)).lower()
    return any(x in search_string for x in ['erbb', 'egfr', 'her2', 'her3', 'her4'])

def analyze_dimer_interface(chain_a_data: Dict, chain_b_data: Dict, use_erbb_terminology: bool) -> Tuple[str, str, float, float]:
    res_a, lm_a = chain_a_data['residues'], chain_a_data['landmarks']
    res_b, lm_b = chain_b_data['residues'], chain_b_data['landmarks']
    
    if not lm_a or not lm_b or lm_a.get('f') is None or lm_b.get('f') is None:
        return "Indeterminate", "Indeterminate", 999.9, 999.9

    def get_lobes(res, lm):
        if lm.get('k') is not None and lm.get('hrd') is not None:
            split_idx = lm['hrd'] - 20
            start_idx = max(0, lm['k'] - 35)
            end_idx = min(len(res), lm['ape'] + 40 if lm.get('ape') else lm['f'] + 60)
        else:
            split_idx = lm['f']
            start_idx = max(0, lm['f'] - 180)
            end_idx = min(len(res), lm['f'] + 180)
        return res[start_idx : split_idx], res[split_idx : end_idx]

    core_n_a, core_c_a = get_lobes(res_a, lm_a)
    core_n_b, core_c_b = get_lobes(res_b, lm_b)

    dist_AC_BN = get_min_dist(core_c_a, core_n_b)
    dist_BC_AN = get_min_dist(core_c_b, core_n_a)

    role_a, role_b = "Unpaired", "Unpaired"
    CONTACT_THRESHOLD = 8.0 
    lbl_sym = "Symmetric"
    lbl_act, lbl_rec = ("Activator", "Receiver") if use_erbb_terminology else ("C-lobe Donor", "N-lobe Receiver")

    if dist_AC_BN < CONTACT_THRESHOLD and dist_BC_AN < CONTACT_THRESHOLD: role_a = role_b = lbl_sym
    elif dist_AC_BN < CONTACT_THRESHOLD: role_a, role_b = lbl_act, lbl_rec
    elif dist_BC_AN < CONTACT_THRESHOLD: role_a, role_b = lbl_rec, lbl_act
        
    return role_a, role_b, round(dist_AC_BN, 2), round(dist_BC_AN, 2)

def get_phi_psi(res_prev: Any, res: Any, res_next: Any) -> Tuple[Optional[float], Optional[float]]:
    if any(r is None for r in [res_prev, res, res_next]): return None, None
    phi = calculate_dihedral(get_atom(res_prev, "C"), get_atom(res, "N"), get_atom(res, "CA"), get_atom(res, "C"))
    psi = calculate_dihedral(get_atom(res, "N"), get_atom(res, "CA"), get_atom(res, "C"), get_atom(res_next, "N"))
    return phi, psi

def assign_ramachandran_region(phi: float, psi: float) -> str:
    if phi is None or psi is None: return "X"
    phi = (phi + 180) % 360 - 180; psi = (psi + 180) % 360 - 180
    if phi < 0: return "A" if -100 <= psi <= 50 else "B"
    else: return "L" if -50 <= psi <= 100 else "E"

def get_rotamer(chi1: float) -> str:
    if chi1 is None: return "Unknown"
    if chi1 < 0: chi1 += 360
    if 0 <= chi1 < 120: return "plus"
    elif 120 <= chi1 < 240: return "trans"
    else: return "minus"

def get_sequence_and_residues(model: Any, chain_id: str) -> Tuple[str, Any]:
    chain_res = model.residues[model.residues.chain_ids == chain_id]
    seq_str = "".join([AA_MAP.get(rname, 'X') for rname in chain_res.names])
    return seq_str, chain_res

def match_cofactor_to_fasta(chain_seq: str, fasta_seqs: Dict, expected_name: str = None) -> str:
    best_name = "Unknown_CoFactor"
    best_ratio = 0.0
    for fasta_name, fseq in fasta_seqs.items():
        sm = difflib.SequenceMatcher(None, chain_seq, fseq.upper(), autojunk=False)
        total_identical = sum(block.size for block in sm.get_matching_blocks())
        ratio = total_identical / len(chain_seq) if len(chain_seq) > 0 else 0.0
        
        if expected_name:
            norm_exp = normalize_sim_name(expected_name).lower()
            norm_fas = normalize_sim_name(fasta_name).lower()
            if norm_exp == norm_fas: ratio += 1.0 
            elif expected_name.lower() in fasta_name.lower() or norm_exp in norm_fas: ratio += 0.5
                
        if ratio > best_ratio:
            best_ratio = ratio
            best_name = fasta_name
            
    if best_ratio >= 0.30: 
        raw_name = clean_protein_type(best_name)
        # Strip trailing unique identifiers added by extract_fasta.py (e.g., _1, _2)
        return re.sub(r'_\d+$', '', raw_name)
    return "Unknown_CoFactor"

def residues_contiguous(residues: Any, i: int, j: int) -> bool:
    """True only if residues i..j are consecutively numbered in the deposited coordinates.

    Landmark indices address the list of RESOLVED residues, so a disordered stretch (common
    in crystallographic activation loops) is invisible to indexing and a naive scan silently
    measures ACROSS the gap. Kincore reports 999 -> 'None' in that situation; we return N/A.
    """
    lo, hi = (i, j) if i <= j else (j, i)
    if lo < 0 or hi >= len(residues): return False
    try:
        return all(residues[k + 1].number == residues[k].number + 1 for k in range(lo, hi))
    except Exception:
        return False


def analyze_activation_loop_dynamic(residues: Any, lm: Dict) -> Tuple[str, str]:
    nt_label, ct_label = "NT_Unk", "CT_Unk"
    idx_f, idx_hrd, idx_ape = lm.get('f'), lm.get('hrd'), lm.get('ape')
    if idx_f is None or idx_hrd is None: return "N/A", "N/A"

    # NB the two endpoints of each contact sit on opposite sides of the activation loop, so they
    # must NOT be required to be contiguous with each other. What has to be verified is that each
    # landmark-plus-offset index still lands on the residue it is supposed to name.
    idx_xhrd = idx_hrd - 1
    if 0 <= idx_xhrd < len(residues) and residues_contiguous(residues, idx_xhrd, idx_hrd):
        nt_dists = []
        for offset in range(3, 7):
            idx_nt_scan = idx_f + offset
            if idx_nt_scan < len(residues) and (idx_ape is None or idx_nt_scan < idx_ape) \
               and residues_contiguous(residues, idx_f, idx_nt_scan):
                nt_dists.append(get_res_min_dist(residues[idx_nt_scan], residues[idx_xhrd]))
        nt_label = ("NTin" if min(nt_dists) <= 5.5 else "NTout") if nt_dists else "N/A"
    else:
        nt_label = "N/A"

    if idx_ape is not None:
        # The C-terminal anchor is the HRD arginine (HRD+1), as in Modi & Dunbrack's APE9-Arg
        # contact -- not HRD+6, which is not a conserved position (it is an Arg in CSK but an
        # Ala in SRC). NB the 5.5 A all-atom cutoff remains AlloQuant's own sensitivity choice
        # and is deliberately looser than Kincore's 6.0 A (8.0 A for the TYR group).
        idx_arg = idx_hrd + 1
        if 0 <= idx_arg < len(residues) and residues_contiguous(residues, idx_hrd, idx_arg):
            ct_dists = []
            scan_limit = max(idx_f + 2, idx_ape - 9)
            for offset in range(idx_ape - 6, scan_limit - 1, -1):
                if 0 <= offset < len(residues) and residues_contiguous(residues, offset, idx_ape):
                    ct_dists.append(get_res_min_dist(residues[offset], residues[idx_arg]))
            ct_label = ("CTin" if min(ct_dists) <= 5.5 else "CTout") if ct_dists else "N/A"
        else:
            ct_label = "N/A"
    return nt_label, ct_label

def analyze_spines(residues: Any, lm: Dict, ligand_atoms: Any) -> Tuple[str, str]:
    r_spine = "Missing"; c_spine = "No Ligand"
    try:
        if all(lm.get(x) is not None for x in ['hrd', 'f', 'rs1', 'rs2']):
            d43 = get_min_sc_dist(residues[lm['hrd']], residues[lm['f']])
            d31 = get_min_sc_dist(residues[lm['f']], residues[lm['rs1']])
            d12 = get_min_sc_dist(residues[lm['rs1']], residues[lm['rs2']])
            r_spine = "Intact" if (d43 < 4.5 and d31 < 4.5 and d12 < 4.5) else "Broken"
    except IndexError: pass

    if ligand_atoms is not None and len(ligand_atoms) > 0 and lm.get('k') is not None and lm['k'] >= 3:
        try:
            vaik_coords = []
            for offset in range(4):
                res_idx = lm['k'] - offset
                if res_idx >= 0: vaik_coords.extend(residues[res_idx].atoms.scene_coords)
            if vaik_coords:
                diff = np.array(vaik_coords)[:, np.newaxis, :] - ligand_atoms.scene_coords[np.newaxis, :, :]
                min_dist = np.min(np.linalg.norm(diff, axis=2))
                c_spine = "Intact" if min_dist < 6.0 else "Ligand Distant"
        except (IndexError, AttributeError): pass
    return r_spine, c_spine

def analyze_core_bridges(residues: Any, lm: Dict) -> Tuple[Any, Any, Any, str]:
    v104_dist, i150_dist, shell_dist, shell_state = "N/A", "N/A", "N/A", "Unknown"
    try:
        if all(lm.get(x) is not None for x in ['v104', 'rs2']): v104_dist = round(get_min_sc_dist(residues[lm['v104']], residues[lm['rs2']]), 2)
        if all(lm.get(x) is not None for x in ['i150', 'hrd']): i150_dist = round(get_min_sc_dist(residues[lm['i150']], residues[lm['hrd']]), 2)
        if all(lm.get(x) is not None for x in ['m118', 'm120']):
            shell_dist = round(get_min_sc_dist(residues[lm['m118']], residues[lm['m120']]), 2)
            shell_state = "Packed" if shell_dist < 5.0 else "Loose"
    except IndexError: pass
    return v104_dist, i150_dist, shell_dist, shell_state

def analyze_alphaC_beta4_loop(residues: Any, lm: Dict) -> Tuple[Any, Any, Any, Any, Any]:
    y156_n99_d, k105_e107_d, k105_e121_d, k105_n99_d, d220_hrd_d = ["N/A"] * 5
    try:
        if all(lm.get(x) is not None for x in ['y156', 'n99']): y156_n99_d = round(get_res_min_dist(residues[lm['y156']], residues[lm['n99']]), 2)
        if all(lm.get(x) is not None for x in ['k105', 'e107']): k105_e107_d = round(get_min_sc_dist(residues[lm['k105']], residues[lm['e107']]), 2)
        if all(lm.get(x) is not None for x in ['k105', 'e121']): k105_e121_d = round(get_min_sc_dist(residues[lm['k105']], residues[lm['e121']]), 2)
        if all(lm.get(x) is not None for x in ['k105', 'n99']): k105_n99_d = round(get_min_sc_dist(residues[lm['k105']], residues[lm['n99']]), 2)
        if all(lm.get(x) is not None for x in ['d220', 'hrd']): d220_hrd_d = round(get_min_sc_dist(residues[lm['d220']], residues[lm['hrd']]), 2)
    except IndexError: pass
    return y156_n99_d, k105_e107_d, k105_e121_d, k105_n99_d, d220_hrd_d

def get_core_name(rel_path: str) -> str:
    base = os.path.splitext(rel_path)[0]
    parts = base.replace('\\', '/').split('/')
    seen_tokens = []
    for part in parts:
        for t in part.split('_'):
            if t and t not in seen_tokens: seen_tokens.append(t)
    return "_".join(seen_tokens)

def get_best_landmark_for_chain(chain_seq: str, cid: str, candidate_lms: List, fasta_seqs: Dict, sim_id: str = "", expected_name: str = None) -> Tuple[str, Optional[Dict]]:
    best_lm = None; best_name = "Unknown"; best_score = -1.0; used_ratio = False
    best_ratio = 0.0  
    
    for fasta_name, coords in candidate_lms:
        struct_score = 0
        if coords.get('f') is not None and coords['f'] < len(chain_seq):
            if chain_seq[coords['f']] in ['F', 'W', 'L', 'Y']: struct_score += 2
            if coords['f'] >= 1 and chain_seq[coords['f']-1] == 'D': struct_score += 2
            if coords['f'] + 1 < len(chain_seq) and chain_seq[coords['f']+1] == 'G': struct_score += 2
        if coords.get('hrd') is not None and coords['hrd'] + 2 < len(chain_seq):
            if chain_seq[coords['hrd']] in ['H', 'Y', 'F']: struct_score += 1
            if chain_seq[coords['hrd']+1] in ['R', 'C', 'H', 'K']: struct_score += 1
            if chain_seq[coords['hrd']+2] == 'D': struct_score += 2

        ratio = 0.0; ref_seq = None
        for header, fseq in fasta_seqs.items():
            clean_header = re.sub(r'^\d+[_\-|]', '', header)
            if fasta_name.lower() in clean_header.lower() or clean_header.lower() in fasta_name.lower():
                ref_seq = fseq; break
        
        if ref_seq:
            sm = difflib.SequenceMatcher(None, chain_seq, ref_seq.upper(), autojunk=False)
            total_identical = sum(block.size for block in sm.get_matching_blocks())
            ratio = total_identical / len(chain_seq) if len(chain_seq) > 0 else 0.0
            used_ratio = True

        # [CRITICAL SCORING FIX] Make sequence identity the dominant factor. 
        # A 100% match gets 1000 pts. A point mutant gets ~990 pts. Mismatches get 0.
        composite = (ratio * 1000) + struct_score
        
        if fasta_name.endswith(f"_{cid}") or fasta_name.endswith(f"-{cid}") or f"chain_{cid}" in fasta_name.lower(): 
            composite += 50
            
        if sim_id:
            sim_tokens = set(re.split(r'[-_]', sim_id.lower()))
            fasta_tokens = set(re.split(r'[-_]', fasta_name.lower()))
            overlap = len(sim_tokens.intersection(fasta_tokens))
            composite += overlap * 2  
            for ft in fasta_tokens:
                if ft in {'apo', 'holo'} or re.match(r'^\d*(atp|adp|amp|gdp|gtp|anp)$', ft) or re.match(r'^p[sty]?\d+$', ft):
                    if ft not in sim_tokens: composite -= 15  

        # Only apply YAML folder-name assumptions if the sequence identity is remotely plausible (>30%)
        if expected_name and ratio > 0.30:
            norm_exp = normalize_sim_name(expected_name).lower()
            norm_fas = normalize_sim_name(fasta_name).lower()
            if norm_exp == norm_fas: composite += 500
            elif expected_name.lower() in fasta_name.lower() or norm_exp in norm_fas: composite += 250

        if composite > best_score:
            best_score = composite; best_lm = coords; best_name = fasta_name; best_ratio = ratio

    if used_ratio:
        if best_ratio >= 0.30: return best_name, best_lm
    else:
        if best_score >= 3.0: return best_name, best_lm
    return "Unknown", None

# --- MASTER EXECUTION ---
def process_model(session, full_cif_path: str, base_dir: str, out_dir_core: str, out_dir_allo: str, all_landmarks: Dict, fasta_seqs: Dict, yaml_patterns: Dict = None) -> List[Dict]:
    if yaml_patterns is None: yaml_patterns = {}
    csv_rows = []
    rel_to_base = os.path.relpath(full_cif_path, base_dir)
    sim_id = get_core_name(rel_to_base)
    
    expected_proteins = []
    for pattern in sorted(yaml_patterns.keys(), key=len, reverse=True):
        if sim_id.startswith(pattern) or pattern in sim_id:
            expected_proteins = [p.get('name') for p in yaml_patterns[pattern].get('proteins', [])]
            break

    # --- NEW HIGHER-ORDER COMPLEX FIX ---
    # If the YAML config failed to capture tertiary/quaternary targets (e.g. c-shc-ch1), 
    # dynamically extract the stoichiometry and identities from the sim_id namespace.
    dynamic_expected = []
    for token in sim_id.split('_'):
        if re.match(r'^[a-z]-', token, re.IGNORECASE):
            dynamic_expected.append(token[2:])
            
    if len(dynamic_expected) > len(expected_proteins):
        expected_proteins = dynamic_expected
    # ------------------------------------
    
    sim_base = normalize_sim_name(sim_id)
    candidate_lms = []
    for fasta_name, coords in all_landmarks.items():
        fasta_base = normalize_sim_name(fasta_name)
        if fasta_base.lower() in sim_base.lower() or sim_base.lower() in fasta_base.lower(): candidate_lms.append((fasta_name, coords))
        elif fasta_name.lower() in sim_id.lower() or fasta_name.replace("_", "-").lower() in sim_id.lower(): candidate_lms.append((fasta_name, coords))
            
    candidate_lms = list({name: coords for name, coords in candidate_lms}.items())
    if not candidate_lms and not fasta_seqs: return csv_rows
    
    models = run(session, f"open '{full_cif_path}'", log=False)
    if not models: return csv_rows
    model = models[0]
    
    cids = sorted(list(set(model.residues.chain_ids)))
    ligand_mask = np.isin(np.array(model.residues.names), LIGAND_NAMES)
    ligands = model.residues[ligand_mask]
    
    chain_data = {}
    cofactors = {}
    has_kinase = False
    
    valid_cids = [cid for cid in cids if len(get_sequence_and_residues(model, cid)[0]) >= 8]
    
    for cid in valid_cids:
        seq, res = get_sequence_and_residues(model, cid)
        expected_type = expected_proteins[valid_cids.index(cid)] if expected_proteins and valid_cids.index(cid) < len(expected_proteins) else None
        lm_name, lm = get_best_landmark_for_chain(seq, cid, candidate_lms, fasta_seqs, sim_id, expected_name=expected_type)
        
        # --- UNIVERSAL CO-FACTOR DETECTION ---
        if lm is None or lm.get('f') is None or lm['f'] >= len(res):
            cf_name = match_cofactor_to_fasta(seq, fasta_seqs, expected_type)
            if expected_type and cf_name == "Unknown_CoFactor": cf_name = clean_protein_type(expected_type)
            cofactors[cid] = {'name': cf_name, 'residues': res}
            continue
            
        has_kinase = True
        spec = f"/{cid}"
        chain_ligand_res, chain_ligand_atoms = None, None
        
        if len(ligands) > 0:
            anchor_pos = get_atom(res[lm['f']], "CA")
            if anchor_pos is None and lm.get('k') is not None: anchor_pos = get_atom(res[lm['k']], "CA")
            if anchor_pos is not None:
                for lig_res in ligands:
                    if np.linalg.norm(anchor_pos - get_center_of_mass(lig_res.atoms)) < 15.0:
                        chain_ligand_res = lig_res; chain_ligand_atoms = lig_res.atoms; break

        r_spine, c_spine = analyze_spines(res, lm, chain_ligand_atoms)
        nt_loop, ct_loop = analyze_activation_loop_dynamic(res, lm)
        v104_d, i150_d, shell_d, shell_state = analyze_core_bridges(res, lm)
        y156_n99_d, k105_e107_d, k105_e121_d, k105_n99_d, d220_hrd_d = analyze_alphaC_beta4_loop(res, lm)
        
        spatial_label, dihedral_label, chelix_label = "Unknown", "Unknown", "Unknown"
        raw_phi_d, raw_psi_d = "N/A", "N/A"
        d1_val, d2_val, sb_dist_val, hrd_atp_dist_val = "N/A", "N/A", "N/A", "N/A"
        dfg_mg_dist_val, dfg_atp_dist_val, ploop_atp_dist_val = "N/A", "N/A", "N/A"
        aCb4_aE_dist_val, spine_bridge_dist_val = "N/A", "N/A"
        cleft_gape_val, mg_hijack_val, clearance_angle_val = "N/A", "N/A", "N/A"
        
        closest_d_o_name, atp_spec = None, None
        closest_dfg_mg_atom, closest_mg_spec = None, None
        aCb4_aE_cmd = None
        ploop_res_nums, ac_b4_res_nums, aE_res_nums, sb_oe_names = [], [], [], []

        # --- Substrate Angle & Cleft Gaping (Roof vs Floor) ---
        roof_res = None
        if lm.get('k') is not None and lm['k'] > 10:
            start_search = max(0, lm['k'] - 35)
            matches = list(re.finditer(r'G.G', seq[start_search:lm['k'] - 5]))
            if matches:
                p_idx = start_search + matches[-1].start()
                # Expand search to reliably hit T14/Y15 equivalents
                ploop_residues = res[p_idx : p_idx + 8] 
            else:
                ploop_residues = res[max(0, lm['k'] - 25):max(0, lm['k'] - 15)]
                
            ploop_res_nums = [r.number for r in ploop_residues]

            # Prioritize Inhibitory Phosphates -> T/Y/S -> Center of P-Loop
            for r in ploop_residues:
                if r.name in ['PTR', 'SEP', 'TPO']:
                    roof_res = r
                    break
            if roof_res is None:
                for r in ploop_residues:
                    if r.name in ['TYR', 'Y', 'THR', 'T', 'SER', 'S']:
                        roof_res = r
                        break
            if roof_res is None and len(ploop_residues) > 0:
                roof_res = ploop_residues[len(ploop_residues) // 2]
                
            # Metric 1: Cleft Gaping (Roof-to-Floor Dist)
            floor_res = res[lm['hrd'] + 2] if (lm.get('hrd') is not None and lm['hrd'] + 2 < len(res)) else None
            roof_ca = get_atom(roof_res, "CA") if roof_res else None
            floor_ca = get_atom(floor_res, "CA") if floor_res else None
            
            if roof_ca is not None and floor_ca is not None:
                cleft_gape_val = round(float(np.linalg.norm(roof_ca - floor_ca)), 2)

            # Metric 2: Mg Hijacking (Inhibitory Phosphate to Mg2+)
            mg_atoms = model.atoms[model.atoms.names == 'MG']
            if len(mg_atoms) > 0 and roof_res is not None:
                target_atoms = roof_res.atoms[np.isin(roof_res.atoms.names, ['P', 'O1P', 'O2P', 'O3P', 'OH', 'OG', 'OG1', 'O'])]
                if len(target_atoms) == 0: target_atoms = roof_res.atoms
                dists = np.linalg.norm(target_atoms.scene_coords[:, np.newaxis, :] - mg_atoms.scene_coords[np.newaxis, :, :], axis=2)
                mg_hijack_val = round(float(np.min(dists)), 2)

            # Metric 3: Substrate Clearance Angle (Roof_CA -- ATP_PG -- HRD_CA)
            if roof_ca is not None and floor_ca is not None and chain_ligand_atoms is not None:
                target_phos = []
                for phos_group in [['PG', 'O1G', 'O2G', 'O3G'],['PB', 'O1B', 'O2B', 'O3B'],['PA', 'O1A', 'O2A', 'O3A']]:
                    if any(np.isin(chain_ligand_atoms.names, phos_group)): target_phos = phos_group; break
                if target_phos:
                    lig_phosphates = chain_ligand_atoms[np.isin(chain_ligand_atoms.names, target_phos)]
                    if len(lig_phosphates) > 0:
                        atp_pg = lig_phosphates.scene_coords[0] 
                        angle = calculate_angle(roof_ca, atp_pg, floor_ca)
                        if angle is not None:
                            clearance_angle_val = round(float(angle), 2)

            if len(ploop_residues) > 0 and chain_ligand_atoms is not None:
                lig_phosphates = chain_ligand_atoms[np.isin(chain_ligand_atoms.names, ['PG', 'O1G', 'O2G', 'O3G', 'PB', 'O1B', 'O2B', 'O3B', 'PA', 'O1A', 'O2A', 'O3A'])]
                ploop_coords = []
                for pr in ploop_residues: ploop_coords.extend(pr.atoms[np.isin(pr.atoms.names, ['N', 'CA', 'C', 'CB'])].scene_coords)
                if len(lig_phosphates) > 0 and len(ploop_coords) > 0: ploop_atp_dist_val = round(float(np.min(np.linalg.norm(np.array(ploop_coords)[:, np.newaxis, :] - lig_phosphates.scene_coords[np.newaxis, :, :], axis=2))), 2)

        # Standard Core Distances
        if lm.get('hrd') is not None and lm['hrd'] + 2 < len(res) and chain_ligand_atoms is not None:
            hrd_sc = get_sidechain_atoms(res[lm['hrd'] + 2])
            d_oxygens = hrd_sc[np.isin(hrd_sc.names, [n for n in hrd_sc.names if n.startswith('O') or n.startswith('N')])]
            target_phos = []
            for phos_group in [['PG', 'O1G', 'O2G', 'O3G'],['PB', 'O1B', 'O2B', 'O3B'],['PA', 'O1A', 'O2A', 'O3A']]:
                if any(np.isin(chain_ligand_atoms.names, phos_group)): target_phos = phos_group; break
            if target_phos:
                lig_phosphates = chain_ligand_atoms[np.isin(chain_ligand_atoms.names, target_phos)]
                if len(d_oxygens) > 0 and len(lig_phosphates) > 0:
                    dists = np.linalg.norm(d_oxygens.scene_coords[:, np.newaxis, :] - lig_phosphates.scene_coords[np.newaxis, :, :], axis=2)
                    min_idx = np.unravel_index(np.argmin(dists), dists.shape)
                    hrd_atp_dist_val = round(float(dists[min_idx]), 2)
                    closest_d_o_name = d_oxygens.names[min_idx[0]]
                    atp_spec = f"/{chain_ligand_res.chain_id}:{chain_ligand_res.number}@{lig_phosphates.names[min_idx[1]]}"

        dfg_d_res = res[lm['f'] - 1] if lm.get('f') is not None and lm['f'] >= 1 else None
        if dfg_d_res is not None:
            dfg_sc = get_sidechain_atoms(dfg_d_res)
            d_oxygens = dfg_sc[np.isin(dfg_sc.names, [n for n in dfg_sc.names if n.startswith('O') or n.startswith('N')])]
            mg_atoms = model.atoms[model.atoms.names == 'MG']
            
            if len(d_oxygens) > 0 and len(mg_atoms) > 0:
                dists = np.linalg.norm(d_oxygens.scene_coords[:, np.newaxis, :] - mg_atoms.scene_coords[np.newaxis, :, :], axis=2)
                min_idx = np.unravel_index(np.argmin(dists), dists.shape)
                if float(dists[min_idx]) <= 15.0:
                    dfg_mg_dist_val = round(float(dists[min_idx]), 2)
                    closest_dfg_mg_atom = d_oxygens.names[min_idx[0]]
                    closest_mg_spec = f"/{mg_atoms[min_idx[1]].residue.chain_id}:{mg_atoms[min_idx[1]].residue.number}@{mg_atoms[min_idx[1]].name}"
                
            if len(d_oxygens) > 0 and chain_ligand_atoms is not None:
                lig_phosphates = chain_ligand_atoms[np.isin(chain_ligand_atoms.names, ['PG', 'O1G', 'O2G', 'O3G', 'PB', 'O1B', 'O2B', 'O3B', 'PA', 'O1A', 'O2A', 'O3A'])]
                if len(lig_phosphates) > 0: dfg_atp_dist_val = round(float(np.min(np.linalg.norm(d_oxygens.scene_coords[:, np.newaxis, :] - lig_phosphates.scene_coords[np.newaxis, :, :], axis=2))), 2)

        if lm.get('c') is not None and lm.get('hrd') is not None:
            ac_b4_res, aE_res = res[lm['c'] + 8 : lm['c'] + 15], res[max(0, lm['hrd'] - 25) : max(0, lm['hrd'] - 10)]
            if len(ac_b4_res) > 0 and len(aE_res) > 0:
                ac_b4_res_nums, aE_res_nums = [r.number for r in ac_b4_res], [r.number for r in aE_res]
                best_dist = 999.9
                for r_ac in ac_b4_res:
                    for r_ae in aE_res:
                        ac_polar = r_ac.atoms[np.isin(r_ac.atoms.names, [n for n in r_ac.atoms.names if n.startswith('N') or n.startswith('O')])]
                        ae_polar = r_ae.atoms[np.isin(r_ae.atoms.names, [n for n in r_ae.atoms.names if n.startswith('N') or n.startswith('O')])]
                        if len(ac_polar) > 0 and len(ae_polar) > 0:
                            dists = np.linalg.norm(ac_polar.scene_coords[:, np.newaxis, :] - ae_polar.scene_coords[np.newaxis, :, :], axis=2)
                            min_idx = np.unravel_index(np.argmin(dists), dists.shape)
                            if float(dists[min_idx]) < best_dist:
                                best_dist, aCb4_aE_dist_val = float(dists[min_idx]), round(float(dists[min_idx]), 2)
                                aCb4_aE_cmd = f"distance {spec}:{r_ac.number}@{ac_polar.names[min_idx[0]]} {spec}:{r_ae.number}@{ae_polar.names[min_idx[1]]} color white radius 0.05"

                if chain_ligand_atoms is not None and len(chain_ligand_atoms) > 0:
                    c1, c2 = ac_b4_res.atoms.scene_coords, chain_ligand_atoms.scene_coords
                    if len(c1) > 0 and len(c2) > 0: spine_bridge_dist_val = round(float(np.min(np.linalg.norm(c1[:, np.newaxis, :] - c2[np.newaxis, :, :], axis=2))), 2)
                elif lm.get('k') is not None and lm['k'] >= 3:
                    spine_bridge_dist_val = round(get_res_min_dist(ac_b4_res, res[lm['k'] - 3]), 2)

        k_ca = get_atom(res[lm['k']], "CA") if lm.get('k') is not None else None
        # D1 anchor is the FOURTH residue after the aC-Glu -- dist(aC-Glu(+4)-CA, DFG-Phe-CZ)
        # per Modi & Dunbrack, PNAS 2019 (116:6818). lm['c'] is the aC-Glu itself, which is
        # correct for the salt bridge below but NOT for D1; the 11/14 A cutoffs are calibrated
        # on the (+4) anchor, so using lm['c'] here compresses D1 and mislabels DFGout/inter.
        c4_idx = lm['c'] + 4 if lm.get('c') is not None else None
        c_ca = get_atom(res[c4_idx], "CA") if (c4_idx is not None and 0 <= c4_idx < len(res)) else None
        f_res = res[lm['f']] if lm.get('f') is not None else None
        f_cz = get_atom(f_res, "CZ")
        if f_cz is None and f_res is not None:
            f_sc, f_ca = get_sidechain_atoms(f_res), get_atom(f_res, "CA")
            if f_sc is not None and f_ca is not None and len(f_sc) > 0: f_cz = f_sc.scene_coords[np.argmax(np.linalg.norm(f_sc.scene_coords - f_ca, axis=1))] 

        if all(x is not None for x in [k_ca, c_ca, f_cz]):
            d1, d2 = float(np.linalg.norm(c_ca - f_cz)), float(np.linalg.norm(k_ca - f_cz))
            d1_val, d2_val = round(d1, 2), round(d2, 2)
            if d1 <= 11.0 and d2 >= 11.0: spatial_label = "DFGin"
            elif d1 > 11.0 and d2 <= 14.0: spatial_label = "DFGout"
            elif d1 <= 11.0 and d2 <= 11.0: spatial_label = "DFGinter"
            else: spatial_label = "Outlier"

        k_nz = get_atom(res[lm['k']], "NZ") if lm.get('k') is not None else None
        c_res = res[lm['c']] if lm.get('c') is not None else None
        if k_nz is not None and c_res is not None:
            c_sc = get_sidechain_atoms(c_res)
            c_polar_atoms = c_sc[np.isin(c_sc.names, [n for n in c_sc.names if n.startswith('O') or n.startswith('N')])]
            if len(c_polar_atoms) > 0:
                dists = np.linalg.norm(c_polar_atoms.scene_coords - k_nz, axis=1)
                sb_dist_val, sb_oe_names = round(float(dists[np.argmin(dists)]), 2), c_polar_atoms.names[np.argsort(dists)[:2]].tolist()

        k_cb, c_cb = get_atom(res[lm['k']], "CB") if lm.get('k') is not None else None, get_atom(res[lm['c']], "CB") if lm.get('c') is not None else None
        if k_cb is not None and c_cb is not None: chelix_label = "In" if np.linalg.norm(k_cb - c_cb) <= 10.0 else "Out"

        idx_f = lm['f']; idx_d = idx_f - 1; idx_x = idx_f - 2  
        if idx_x - 1 >= 0 and idx_f + 1 < len(res):
            phi_x, psi_x = get_phi_psi(res[idx_x-1], res[idx_x], res[idx_d])
            phi_d, psi_d = get_phi_psi(res[idx_x], res[idx_d], res[idx_f])
            phi_f, psi_f = get_phi_psi(res[idx_d], res[idx_f], res[idx_f+1])
            if phi_d is not None and psi_d is not None: raw_phi_d, raw_psi_d = round(phi_d, 2), round(psi_d, 2)
            rama_x, rama_d, rama_f = assign_ramachandran_region(phi_x, psi_x), assign_ramachandran_region(phi_d, psi_d), assign_ramachandran_region(phi_f, psi_f)
            rot_f = get_rotamer(calculate_dihedral(get_atom(res[idx_f], "N"), get_atom(res[idx_f], "CA"), get_atom(res[idx_f], "CB"), get_atom(res[idx_f], "CG")))
            if "X" not in [rama_x, rama_d, rama_f] and rot_f != "Unknown": dihedral_label = f"{rama_x}{rama_d}{rama_f}{rot_f}"

        state = "Unknown"
        if spatial_label == "DFGin":
            if dihedral_label == "BLAminus": state = "Active (BLAminus)"
            elif dihedral_label == "ABAminus": state = "Active-Like (ABAminus)"
            elif dihedral_label == "BLBplus": state = "Inactive (BLBplus)"
            elif dihedral_label == "BLBtrans": state = "Inactive (BLBtrans)"
            elif dihedral_label == "BLAplus": state = "Inactive (BLAplus)"
            elif dihedral_label == "BLBminus": state = "Inactive (BLBminus)"
            else: state = f"DFGin ({dihedral_label})"
        elif spatial_label == "DFGout" and dihedral_label == "BBAminus": state = "Inactive (BBAminus)"
        elif spatial_label == "DFGinter" and dihedral_label == "BABtrans": state = "Inactive (BABtrans)"
        else: state = f"{spatial_label} ({dihedral_label})"

        chain_data[cid] = {
            "residues": res, "landmarks": lm, 
            "meta": {"Type": lm_name, "State": state, "CHelix": chelix_label,
                     "RSpine": r_spine, "CSpine": c_spine, "Spatial": spatial_label, "Dihedral": dihedral_label,
                     "ActLoop_NT": nt_loop, "ActLoop_CT": ct_loop, "Phi_D": raw_phi_d, "Psi_D": raw_psi_d,
                     "Cleft_Gape_Dist": cleft_gape_val, "Mg_Hijack_Dist": mg_hijack_val, "Substrate_Clearance_Angle": clearance_angle_val,
                     "D1_Dist": d1_val, "D2_Dist": d2_val, "SB_Dist": sb_dist_val, 
                     "HRD_ATP_Dist": hrd_atp_dist_val, "DFG_Mg_Dist": dfg_mg_dist_val,
                     "DFG_ATP_Dist": dfg_atp_dist_val, "PLoop_ATP_Dist": ploop_atp_dist_val,
                     "V104_RS2_Dist": v104_d, "I150_HRD_Dist": i150_d, "Shell_M118_M120_Dist": shell_d,
                     "Shell_State": shell_state,
                     "Y156_N99_Dist": y156_n99_d, "K105_E107_Dist": k105_e107_d, 
                     "K105_E121_Dist": k105_e121_d, "K105_N99_Dist": k105_n99_d, "D220_HRD_Dist": d220_hrd_d,
                     "aCb4_aE_Dist": aCb4_aE_dist_val, "Spine_Bridge_Dist": spine_bridge_dist_val, "aCb4_aE_cmd": aCb4_aE_cmd,
                     "atp_spec": atp_spec, "closest_d_o_name": closest_d_o_name, "dfg_d_res": dfg_d_res,
                     "closest_dfg_mg_atom": closest_dfg_mg_atom, "closest_mg_spec": closest_mg_spec,
                     "ploop_res_nums": ploop_res_nums, "sb_oe_names": sb_oe_names,
                     "ac_b4_res_nums": ac_b4_res_nums, "aE_res_nums": aE_res_nums}
        }

    # =========================================================================
    # --- UNIVERSAL CO-FACTOR DISTANCE MAPPING ---
    # =========================================================================
    for cid, data in chain_data.items():
        res, lm = data['residues'], data['landmarks']
        cf_name_val, cf_aC_dist_val, cf_actloop_dist_val, cf_global_dist = "None", "N/A", "N/A", 999.9
        
        if cofactors:
            best_cf_dist = 999.9
            best_cf = None
            for cf_id, cf_data in cofactors.items():
                d = get_min_dist(res, cf_data['residues'])
                if d < best_cf_dist:
                    best_cf_dist, best_cf = d, cf_data
                    
            if best_cf and best_cf_dist <= 15.0:
                cf_name_val = best_cf['name']
                cf_global_dist = round(best_cf_dist, 2)
                if lm.get('c') is not None:
                    aC_res = res[max(0, lm['c'] - 4) : min(len(res), lm['c'] + 5)]
                    cf_aC_dist_val = round(get_min_dist(aC_res, best_cf['residues']), 2)
                if lm.get('f') is not None and lm.get('ape') is not None and lm['ape'] < len(res):
                    act_res = res[lm['f'] : lm['ape'] + 1]
                    cf_actloop_dist_val = round(get_min_dist(act_res, best_cf['residues']), 2)
                    
        data['meta']['CoFactor_Name'] = cf_name_val
        data['meta']['CoFactor_Global_Dist'] = cf_global_dist
        data['meta']['CoFactor_aC_Dist'] = cf_aC_dist_val
        data['meta']['CoFactor_ActLoop_Dist'] = cf_actloop_dist_val

    # =========================================================================
    # --- DIMER INTERFACE DETECTION & ROLE ASSIGNMENTS ---
    # =========================================================================
    chain_types = [data['meta']['Type'] for data in chain_data.values()]
    use_erbb_terms = is_erbb_system(sim_id, chain_types)

    chain_roles = {cid: "Unpaired" for cid in chain_data}
    partner_map = {cid: "None" for cid in chain_data}
    int_dist_AC_BN = {cid: "N/A" for cid in chain_data}
    int_dist_BC_AN = {cid: "N/A" for cid in chain_data}
    
    cids_valid = list(chain_data.keys())
    for i in range(len(cids_valid)):
        for j in range(i+1, len(cids_valid)):
            cA, cB = cids_valid[i], cids_valid[j]
            if chain_data[cA]['landmarks'] and chain_data[cB]['landmarks']:
                role_a, role_b, dAC_BN, dBC_AN = analyze_dimer_interface(chain_data[cA], chain_data[cB], use_erbb_terms)

                if role_a != "Unpaired":
                    chain_roles[cA], chain_roles[cB] = role_a, role_b
                    partner_map[cA], partner_map[cB] = cB, cA
                    int_dist_AC_BN[cA], int_dist_AC_BN[cB] = dAC_BN, dBC_AN 
                    int_dist_BC_AN[cA], int_dist_BC_AN[cB] = dBC_AN, dAC_BN

    # --- Post-Processing: Co-Factor Partner Override ---
    for cid in cids_valid:
        if chain_roles[cid] == "Unpaired":
            cf_name = chain_data[cid]['meta'].get('CoFactor_Name', 'None')
            cf_dist = chain_data[cid]['meta'].get('CoFactor_Global_Dist', 999.9)
            if cf_name != "None" and cf_dist <= 8.0:
                chain_roles[cid] = "Co-factor Bound"
                partner_map[cid] = cf_name

    # =========================================================================
    # --- SPLIT VISUALIZATION GENERATORS ---
    # =========================================================================
    display_names = [re.sub(r'^[A-Za-z]-', '', data['meta']['Type']) for data in chain_data.values()]
    protein_name = " / ".join(sorted(list(set(display_names))))
    
    rel_to_out_core = os.path.relpath(full_cif_path, out_dir_core).replace('\\', '/')
    rel_to_out_allo = os.path.relpath(full_cif_path, out_dir_allo).replace('\\', '/')

    cxc_core = [
        f"# KinCore Core Visualization Macro for {sim_id}", f"open {rel_to_out_core}",
        "graphics silhouettes true", "graphics silhouettes width 1.5", "lighting soft", "color white", 
        "hide atoms", "show cartoons", "transparency 50 cartoons", "style stick", "2dlabels delete all",
        "\n# --- On-Screen Titles & Legend ---",
        f"2dlabels create title_prot text '{protein_name}' color white size 24 bold true xpos 0.05 ypos 0.94",
        "2dlabels create subtitle_prot text '(Catalytic Core)' color white size 20 xpos 0.05 ypos 0.90",
        "2dlabels create leg_title text 'Color Legend:' color white size 20 bold true xpos 0.80 ypos 0.86",
        "2dlabels create leg_aloop text 'A-Loop (DFG/Dynamic)' color coral size 16 xpos 0.80 ypos 0.83",
        "2dlabels create leg_ploop text 'P-Loop (Gly-Rich)' color deep sky blue size 16 xpos 0.80 ypos 0.80",
        "2dlabels create leg_rspine text 'R-Spine' color medium purple size 16 xpos 0.80 ypos 0.77",
        "2dlabels create leg_sb text 'Salt Bridge (K-C)' color spring green size 16 xpos 0.80 ypos 0.74",
        "2dlabels create leg_cat text 'Catalytic HRD-Asp' color red size 16 xpos 0.80 ypos 0.71",
        "2dlabels create leg_dfg text 'DFG-Asp Coordination' color dodger blue size 16 xpos 0.80 ypos 0.68",
        "2dlabels create leg_lig text 'ATP / Magnesium' color gold size 16 xpos 0.80 ypos 0.65"
    ]

    cxc_allo = [
        f"# KinCore Allosteric Macro for {sim_id}", f"open {rel_to_out_allo}",
        "graphics silhouettes true", "graphics silhouettes width 1.5", "lighting soft", "color white", 
        "hide atoms", "show cartoons", "transparency 50 cartoons", "style stick", "2dlabels delete all",
        "\n# --- On-Screen Titles & Legend ---",
        f"2dlabels create title_prot text '{protein_name}' color white size 24 bold true xpos 0.05 ypos 0.94",
        "2dlabels create subtitle_prot text '(aC-b4 Allosteric Network)' color white size 20 xpos 0.05 ypos 0.90",
        "2dlabels create leg_title text 'Color Legend:' color white size 20 bold true xpos 0.80 ypos 0.86",
        "2dlabels create leg_acb4 text 'aC-b4 Loop (Allosteric Bridge)' color hot pink size 16 xpos 0.80 ypos 0.83",
        "2dlabels create leg_ae text 'aE Helix (Core Anchor)' color yellow size 16 xpos 0.80 ypos 0.80",
        "2dlabels create leg_cat text 'Catalytic HRD-Asp' color red size 16 xpos 0.80 ypos 0.77",
        "2dlabels create leg_lig text 'ATP / Magnesium' color gold size 16 xpos 0.80 ypos 0.74"
    ]

    if cofactors:
        cxc_core.append("2dlabels create leg_cofactor text 'Steric Co-Factor' color dark cyan size 16 xpos 0.80 ypos 0.62")
        cxc_allo.append("2dlabels create leg_cofactor text 'Steric Co-Factor' color dark cyan size 16 xpos 0.80 ypos 0.71")
        for cf_id, cf_data in cofactors.items():
            spec = f"/{cf_id}"
            cxc_core.extend([f"color {spec} dark cyan", f"transparency {spec} 30 cartoons"])
            cxc_allo.extend([f"color {spec} dark cyan", f"transparency {spec} 30 cartoons"])

    if len(ligands) > 0:
        lig_cmds = [
            "\n# --- Ligands & Interactions ---", "show :ATP,ADP,ANP,ACP,AGS,AMP,GTP,GDP,STU",
            "color :ATP,ADP,ANP,ACP,AGS,AMP,GTP,GDP,STU byhetero", "color :ATP,ADP,ANP,ACP,AGS,AMP,GTP,GDP,STU@C* gold",
            "label :ATP,ADP,ANP,ACP,AGS,AMP,GTP,GDP,STU residues color gold height 1.5",
            "hbonds :ATP,ADP,ANP,ACP,AGS,AMP,GTP,GDP,STU restrict protein color cyan radius 0.05",
            "show :MG", "color :MG green", "style :MG sphere", "size :MG atomRadius 1.0"
        ]
        cxc_core.extend(lig_cmds); cxc_allo.extend(lig_cmds)

    y_offset = 0.84

    for cid, data in chain_data.items():
        meta = data['meta']
        csv_rows.append({
            "Simulation_ID": sim_id, "Directory": os.path.dirname(full_cif_path), "File": os.path.basename(full_cif_path), 
            "Chain": cid, "Type": meta['Type'], "State": meta['State'],
            "Role": chain_roles[cid], "Partner": partner_map[cid],
            "Interface_C_Lobe_Donor_Dist": int_dist_AC_BN[cid], "Interface_N_Lobe_Rec_Dist": int_dist_BC_AN[cid],
            "CoFactor_Name": meta.get('CoFactor_Name', 'None'),
            "CoFactor_aC_Dist": meta.get('CoFactor_aC_Dist', 'N/A'),
            "CoFactor_ActLoop_Dist": meta.get('CoFactor_ActLoop_Dist', 'N/A'),
            "R_Spine": meta['RSpine'], "C_Spine": meta['CSpine'], "C_Helix": meta['CHelix'],
            "Shell_State": meta['Shell_State'], "Spatial": meta['Spatial'], "Dihedral": meta['Dihedral'],
            "ActLoop_NT": meta.get('ActLoop_NT', 'N/A'), "ActLoop_CT": meta.get('ActLoop_CT', 'N/A'), 
            "Phi_D": meta.get('Phi_D', 'N/A'), "Psi_D": meta.get('Psi_D', 'N/A'),
            "Cleft_Gape_Dist": meta.get('Cleft_Gape_Dist', 'N/A'),
            "Mg_Hijack_Dist": meta.get('Mg_Hijack_Dist', 'N/A'),
            "Substrate_Clearance_Angle": meta.get('Substrate_Clearance_Angle', 'N/A'),
            "D1_Dist": meta['D1_Dist'], "D2_Dist": meta['D2_Dist'], "SB_Dist": meta['SB_Dist'],
            "HRD_ATP_Dist": meta['HRD_ATP_Dist'], "DFG_Mg_Dist": meta.get('DFG_Mg_Dist', 'N/A'), 
            "DFG_ATP_Dist": meta.get('DFG_ATP_Dist', 'N/A'), "PLoop_ATP_Dist": meta.get('PLoop_ATP_Dist', 'N/A'),
            "aCb4_aE_Dist": meta.get('aCb4_aE_Dist', 'N/A'), "Spine_Bridge_Dist": meta.get('Spine_Bridge_Dist', 'N/A'),
            "V104_RS2_Dist": meta['V104_RS2_Dist'], "I150_HRD_Dist": meta['I150_HRD_Dist'], "Shell_M118_M120_Dist": meta['Shell_M118_M120_Dist'],
            "Y156_N99_Dist": meta['Y156_N99_Dist'], "K105_E107_Dist": meta['K105_E107_Dist'], 
            "K105_E121_Dist": meta['K105_E121_Dist'], "K105_N99_Dist": meta['K105_N99_Dist'], "D220_HRD_Dist": meta['D220_HRD_Dist']
        })

        lm, res = data['landmarks'], data['residues']
        f_num, d_num = res[lm['f']].number, res[lm['f']].number - 1
        hrd_d_num = res[lm['hrd'] + 2].number if (lm.get('hrd') is not None and lm['hrd'] + 2 < len(res)) else None
        spec = f"/{cid}"
        
        clean_type = re.sub(r'^[A-Za-z]-', '', meta['Type'])
        role_str = f" | {chain_roles[cid]}" if chain_roles[cid] not in ["N/A", "Unpaired"] else ""
        state_label = f"2dlabels create State_{cid} text 'Chain {cid} ({clean_type}{role_str}): {meta['State']}' color white size 20 xpos 0.05 ypos {y_offset}"
        cxc_core.append(state_label); cxc_allo.append(state_label)

        cxc_core.extend([f"color {spec} #d3d3d3", f"show {spec}:{d_num}-{f_num+1}", f"color {spec}:{d_num}-{f_num+1} dark orange"])
        if meta.get('ploop_res_nums'): cxc_core.extend([f"color {spec}:{meta['ploop_res_nums'][0]}-{meta['ploop_res_nums'][-1]} deep sky blue"])
        if lm.get('ape') is not None: cxc_core.extend([f"color {spec}:{f_num+1}-{res[lm['ape']].number} coral", f"show {spec}:{res[lm['ape']].number} & sidechain"])
            
        if all(lm.get(x) is not None for x in ['hrd', 'f', 'rs1', 'rs2']):
            r1, r2, r3, r4 = res[lm['rs2']].number, res[lm['rs1']].number, res[lm['f']].number, res[lm['hrd']].number
            cxc_core.extend([f"color {spec}:{r1},{r2},{r3},{r4} medium purple", f"show {spec}:{r1},{r2},{r3},{r4} & sidechain",
                             f"distance {spec}:{r1}@CA {spec}:{r2}@CA color medium purple radius 0.05",
                             f"distance {spec}:{r2}@CA {spec}:{r3}@CA color medium purple radius 0.05",
                             f"distance {spec}:{r3}@CA {spec}:{r4}@CA color medium purple radius 0.05"])
            
        if lm.get('k') is not None and lm.get('c') is not None:
            k_num, c_num = res[lm['k']].number, res[lm['c']].number
            cxc_core.extend([f"color {spec}:{k_num},{c_num} spring green", f"show {spec}:{k_num},{c_num} & sidechain"])
            for oe_name in meta.get('sb_oe_names', []): cxc_core.append(f"distance {spec}:{k_num}@NZ {spec}:{c_num}@{oe_name} color spring green radius 0.05")

        if all(lm.get(x) is not None for x in ['v104', 'm118', 'm120', 'rs2']):
            v104, m118, m120, rs2 = res[lm['v104']].number, res[lm['m118']].number, res[lm['m120']].number, res[lm['rs2']].number
            cxc_core.extend([f"color {spec}:{v104},{m118},{m120} teal", f"show {spec}:{v104},{m118},{m120} & sidechain",
                             f"distance {spec}:{m118}@CA {spec}:{m120}@CA color teal radius 0.05", f"distance {spec}:{v104}@CA {spec}:{rs2}@CA color teal radius 0.05"])
            
        if lm.get('i150') is not None and lm.get('hrd') is not None:
            cxc_core.extend([f"color {spec}:{res[lm['i150']].number} tan", f"show {spec}:{res[lm['i150']].number} & sidechain",
                             f"distance {spec}:{res[lm['i150']].number}@CA {spec}:{res[lm['hrd']].number}@CA color tan radius 0.05"])

        cxc_allo.extend([f"color {spec} #d3d3d3"])
        if meta.get('ac_b4_res_nums') and meta.get('aE_res_nums'):
            ac_b4_n1, ac_b4_n2 = meta['ac_b4_res_nums'][0], meta['ac_b4_res_nums'][-1]
            aE_n1, aE_n2 = meta['aE_res_nums'][0], meta['aE_res_nums'][-1]
            cxc_allo.extend([f"\n# Allosteric Network: aC-b4 Loop & aE Helix", f"color {spec}:{ac_b4_n1}-{ac_b4_n2} hot pink",
                             f"color {spec}:{aE_n1}-{aE_n2} yellow", f"show {spec}:{ac_b4_n1}-{ac_b4_n2},{aE_n1}-{aE_n2} & sidechain"])
            if meta.get('aCb4_aE_cmd'): cxc_allo.append(meta['aCb4_aE_cmd'])

        if lm.get('y156') is not None and lm.get('n99') is not None:
            y156, n99 = res[lm['y156']].number, res[lm['n99']].number
            cxc_allo.extend([f"\n# aE Anchor", f"color {spec}:{y156} lime green", f"show {spec}:{y156} & sidechain", f"distance {spec}:{y156}@CA {spec}:{n99}@CA color lime green radius 0.05"])

        if hrd_d_num:
            cat_cmd = [f"show {spec}:{hrd_d_num} & sidechain", f"color {spec}:{hrd_d_num} & sidechain red"]
            cxc_core.extend(cat_cmd); cxc_allo.extend(cat_cmd)
            if meta.get('atp_spec') and meta.get('closest_d_o_name') and meta.get('HRD_ATP_Dist') != "N/A" and meta['HRD_ATP_Dist'] < 10.0:
                dist_cmd = f"distance {spec}:{hrd_d_num}@{meta['closest_d_o_name']} {meta['atp_spec']} color magenta radius 0.05"
                cxc_core.append(dist_cmd); cxc_allo.append(dist_cmd)

        if meta.get('dfg_d_res'):
            dfg_d_num = meta['dfg_d_res'].number
            dfg_cmd = [f"show {spec}:{dfg_d_num} & sidechain", f"color {spec}:{dfg_d_num} & sidechain dodger blue"]
            cxc_core.extend(dfg_cmd); cxc_allo.extend(dfg_cmd)
            if meta.get('DFG_Mg_Dist') != "N/A" and meta.get('closest_dfg_mg_atom') and meta.get('closest_mg_spec') and meta['DFG_Mg_Dist'] < 10.0:
                dist_cmd = f"distance {spec}:{dfg_d_num}@{meta['closest_dfg_mg_atom']} {meta['closest_mg_spec']} color green radius 0.05"
                cxc_core.append(dist_cmd); cxc_allo.append(dist_cmd)
        
        y_offset -= 0.04

    if has_kinase:
        cxc_core.extend(["\n# Final Polish", "hide H"]); cxc_allo.extend(["\n# Final Polish", "hide H"])
        view_targets = " ".join([f"/{cid}:{data['residues'][data['landmarks']['f']].number}" 
                                 for cid, data in chain_data.items() if data['landmarks'] and data['landmarks'].get('f') is not None])
        if view_targets: 
            cxc_core.extend([f"view {view_targets}", "zoom 1.2", "clip off"])
            cxc_allo.extend([f"view {view_targets}", "zoom 1.2", "clip off"])
        
        with open(os.path.join(out_dir_core, f"viz_core_{sim_id}.cxc"), 'w') as out_f: out_f.write("\n".join(cxc_core))
        with open(os.path.join(out_dir_allo, f"viz_allo_{sim_id}.cxc"), 'w') as out_f: out_f.write("\n".join(cxc_allo))

    run(session, "close session", log=False)
    return csv_rows

def parse_yaml_fallback(filepath: str) -> Dict:
    patterns = {}
    if not os.path.exists(filepath): return patterns
    with open(filepath, 'r') as yf:
        current_pattern, in_pm = None, False
        for line in yf:
            stripped = line.strip()
            if not stripped: continue
            if line.startswith('pattern_matches:'): in_pm = True
            elif in_pm and line.startswith('  ') and not line.startswith('    '):
                current_pattern = stripped.rstrip(':')
                patterns[current_pattern] = {'proteins': []}
            elif in_pm and line.startswith('    - name:'):
                if current_pattern: patterns[current_pattern]['proteins'].append({'name': stripped.split('name:')[1].strip()})
    return patterns

def main(session):
    base_dir = os.path.abspath(SEARCH_DIR)
    out_dir_core, out_dir_allo = os.path.join(base_dir, VIZ_OUT_DIR_CORE), os.path.join(base_dir, VIZ_OUT_DIR_ALLO)
    os.makedirs(out_dir_core, exist_ok=True); os.makedirs(out_dir_allo, exist_ok=True)
    
    chunk_list_file = os.environ.get("CHIMERAX_CHUNK")
    if not chunk_list_file or not os.path.exists(chunk_list_file): run(session, "quit"); return

    out_csv_name = f"{os.path.splitext(os.path.basename(chunk_list_file))[0]}_results_v7r3.csv"
    landmarks_file = os.path.join(base_dir, LANDMARKS_JSON)
    if not os.path.exists(landmarks_file): run(session, "quit"); return
        
    with open(landmarks_file, 'r') as f: all_landmarks = json.load(f)
    with open(chunk_list_file, 'r') as f: files_to_process = list(dict.fromkeys([line.strip() for line in f if line.strip()]))

    fasta_path, fasta_seqs = os.path.join(base_dir, "sequences.fasta"), {}
    if os.path.exists(fasta_path):
        with open(fasta_path, 'r') as f:
            name, seq_lines = None, []
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if name: fasta_seqs[name] = "".join(seq_lines)
                    name, seq_lines = line[1:].strip().split()[0], []
                elif line: seq_lines.append(line)
            if name: fasta_seqs[name] = "".join(seq_lines)

    yaml_patterns = {}
    yaml_path = os.path.join(base_dir, "proteins.yaml")
    if os.path.exists(yaml_path):
        try:
            import yaml
            with open(yaml_path, 'r') as yf:
                ydata = yaml.safe_load(yf)
                if ydata and 'pattern_matches' in ydata: yaml_patterns = ydata['pattern_matches']
        except Exception: yaml_patterns = parse_yaml_fallback(yaml_path)

    cols = ["Simulation_ID", "Directory", "File", "Chain", "Type", "State", "Role", "Partner",
            "CoFactor_Name", "CoFactor_aC_Dist", "CoFactor_ActLoop_Dist",
            "Interface_C_Lobe_Donor_Dist", "Interface_N_Lobe_Rec_Dist",
            "R_Spine", "C_Spine", "C_Helix", "Shell_State", "Spatial", "Dihedral", "ActLoop_NT", "ActLoop_CT", 
            "Phi_D", "Psi_D", 
            "Cleft_Gape_Dist", "Mg_Hijack_Dist", "Substrate_Clearance_Angle", # NEW METRICS
            "D1_Dist", "D2_Dist", "SB_Dist", 
            "HRD_ATP_Dist", "DFG_Mg_Dist", "DFG_ATP_Dist", "PLoop_ATP_Dist",
            "aCb4_aE_Dist", "Spine_Bridge_Dist",
            "V104_RS2_Dist", "I150_HRD_Dist", "Shell_M118_M120_Dist",
            "Y156_N99_Dist", "K105_E107_Dist", "K105_E121_Dist", "K105_N99_Dist", "D220_HRD_Dist"]
           
    all_rows = list() 
    for full_cif_path in files_to_process:
        try: all_rows.extend(process_model(session, full_cif_path, base_dir, out_dir_core, out_dir_allo, all_landmarks, fasta_seqs, yaml_patterns))
        except Exception as e: print(f"Error processing {full_cif_path}: {e}")

    with open(os.path.join(base_dir, out_csv_name), 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=cols)
        writer.writeheader()
        writer.writerows(all_rows)

    run(session, "quit")

main(session)
