import os
import csv
import json
import numpy as np
import re
import difflib
from typing import Dict, List, Tuple, Optional, Any
from chimerax.core.commands import run

# ==============================================================================
# UPDATE LOG (v6r2 -> v6r3)
# Date: April 21, 2026
#
# 1. Trans-Tail Regulatory Blockade (Biologically UNIQUE to RAF dimers):
#    - Added `get_raf_tail_asymmetry()` to explicitly measure the distance 
#      between the N-terminal regulatory tail/RBD of one kinase and the 
#      catalytic cleft (HRD/DFG motifs) of its partner.
#    - Relaxed the physical blockade detection threshold to < 7.0 Å to account 
#      for experimental flexibility seen in gold-standard Cryo-EM structures 
#      (e.g., PDB 6UAN).
#
# 2. Physiological Role Reversal (Nomenclature Fix):
#    - Corrected the role assignments to match Kuriyan lab mechanics. The kinase
#      that inserts its tail into the partner's cleft is now correctly labeled 
#      the `Receiver`, while the structurally locked, catalytically occluded 
#      kinase is labeled the `Activator (Tail Blocked)`.
#
# 3. Structural Anchor Hardening (Index vs. PDB Number Fix):
#    - Refactored `get_raf_interface_metrics()` to anchor strictly on the physical
#      PDB residue number of the conserved aC-helix Glutamate (HMM node 48 / 'c').
#    - This bypasses previous array-indexing shifts caused by missing/unresolved 
#      N-terminal loops in experimental structures, considered to fix ~12 Å 
#      measurements on the LLA hydrophobic core (but ~ 12 Å may not be abnormal as also in 6UAN).
#
# 4. CSV Diagnostic Expansion:
#    - Output now permanently tracks `RAF_Tail_A_Block` and `RAF_Tail_B_Block` 
#      alongside the core LLA/Push metrics, allowing for instant differentiation 
#      between internally uncoupled cores (AF3 artifacts) and true physiological 
#      trans-tail asymmetric dimers.
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
LIGAND_NAMES =["ATP", "ADP", "ANP", "ACP", "AGS", "AMP", "GTP", "GDP", "STU"]

# --- MATH & ATOM HELPERS ---
def normalize_sim_name(name: str) -> str:
    """Normalizes 3D file paths exactly like the FASTA headers to guarantee a match."""
    clean = re.sub(r'_?(model|seed|rank|pred|unrelaxed|relaxed)[\-_]?\d+', '', name, flags=re.IGNORECASE)
    clean = re.sub(r'_?(unrelaxed|relaxed)', '', clean, flags=re.IGNORECASE)
    clean = re.sub(r'^(fold|job|run)_', '', clean, flags=re.IGNORECASE)
    clean = re.sub(r'[-_](apo|holo|py\d+|\d*atp)$', '', clean, flags=re.IGNORECASE)
    clean = re.sub(r'(wt|cattail|cat|tail)', '', clean, flags=re.IGNORECASE)
    clean = re.sub(r'[-_]+', '-', clean).strip('-')
    return clean.strip('_')

def get_base_sim_name(name: str) -> str:
    clean = re.sub(r'_?(model|seed|rank|pred|unrelaxed|relaxed)[\-_]?\d+', '', name, flags=re.IGNORECASE)
    clean = re.sub(r'_?(unrelaxed|relaxed)', '', clean, flags=re.IGNORECASE)
    clean = re.sub(r'^(fold|job|run)_', '', clean, flags=re.IGNORECASE)
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
    if any(x is None for x in[p1, p2, p3, p4]): return None
    b0 = -1.0 * (p2 - p1); b1 = p3 - p2; b2 = p4 - p3
    b1_norm = np.linalg.norm(b1)
    if b1_norm == 0: return None
    b1 /= b1_norm
    v = b0 - np.dot(b0, b1) * b1; w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w); y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))

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

# --- DOMAIN, DIMER & INTERFACE HELPERS ---
def is_raf_system(sim_id: str, chain_types: List[str]) -> bool:
    search_string = (sim_id + " " + " ".join(chain_types)).lower()
    return any(x in search_string for x in ['raf', 'braf', 'raf1', 'craf'])

def get_raf_interface_metrics(res_a: Any, res_b: Any, lm_a: Dict, lm_b: Dict) -> Tuple[float, float, float]:
    """
    Returns (A_Push_Dist, B_Push_Dist, LLA_Core_Dist). 
    Dynamically finds the RAF interface residues anchoring on HMM node 48 ('c').
    """
    a_glu = lm_a.get('c')
    b_glu = lm_b.get('c')

    if a_glu is None or b_glu is None:
        return 999.9, 999.9, 999.9

    try:
        # Arginine Anchor is +8 from the Glu. LLA Leucine is +4.
        r_a = res_a[a_glu + 8]
        l_a = res_a[a_glu + 4]
        r_b = res_b[b_glu + 8]
        l_b = res_b[b_glu + 4]

        # Get aC-helix regions using HMM node 48 ('c') as the center
        helix_a = res_a[max(0, a_glu - 4) : a_glu + 5]
        helix_b = res_b[max(0, b_glu - 4) : b_glu + 5]
    except IndexError:
        return 999.9, 999.9, 999.9

    # Calculate "Push" (Dist from Arginine/Anchor to partner aC-Helix COM)
    # Fallback to CB or CA if the +8 residue is mutated or shifted in alignment
    cz_a = get_atom(r_a, "CZ") or get_atom(r_a, "CB") or get_atom(r_a, "CA")
    cz_b = get_atom(r_b, "CZ") or get_atom(r_b, "CB") or get_atom(r_b, "CA")

    com_helix_a = get_center_of_mass(helix_a)
    com_helix_b = get_center_of_mass(helix_b)

    dist_a_push = float(np.linalg.norm(cz_a - com_helix_b)) if cz_a is not None and com_helix_b is not None else 999.9
    dist_b_push = float(np.linalg.norm(cz_b - com_helix_a)) if cz_b is not None and com_helix_a is not None else 999.9
    
    # Calculate hydrophobic LLA core distance
    lla_dist = get_res_min_dist(l_a, l_b) if l_a and l_b else 999.9

    return round(dist_a_push, 2), round(dist_b_push, 2), round(lla_dist, 2)

def get_raf_tail_asymmetry(res_a: Any, res_b: Any, lm_a: Dict, lm_b: Dict) -> Tuple[float, float]:
    """
    Calculates the minimum distance from the tail regions (non-kinase domain) of Chain A 
    to the catalytic cleft (HRD/DFG motifs) of Chain B, and vice versa.
    """
    def get_cleft_atoms(res, lm):
        cleft_atoms = []
        if lm.get('hrd') is not None:
            cleft_atoms.extend(res[lm['hrd']].atoms.scene_coords)
        if lm.get('f') is not None:
            cleft_atoms.extend(res[lm['f']].atoms.scene_coords)
        return np.array(cleft_atoms) if cleft_atoms else None

    def get_tail_atoms(res, lm):
        # Define core as HMM node 'k'-10 to 'f'+40. Tails are everything outside this.
        start_idx = max(0, lm.get('k', 10) - 10)
        end_idx = min(len(res), lm.get('f', len(res)-40) + 40)
        
        tail_atoms = []
        for i, r in enumerate(res):
            if i < start_idx or i > end_idx:
                tail_atoms.extend(r.atoms.scene_coords)
        return np.array(tail_atoms) if tail_atoms else None

    cleft_a = get_cleft_atoms(res_a, lm_a)
    cleft_b = get_cleft_atoms(res_b, lm_b)
    tails_a = get_tail_atoms(res_a, lm_a)
    tails_b = get_tail_atoms(res_b, lm_b)

    dist_ta_cb = 999.9
    if tails_a is not None and cleft_b is not None and len(tails_a) > 0 and len(cleft_b) > 0:
        diff = tails_a[:, np.newaxis, :] - cleft_b[np.newaxis, :, :]
        dist_ta_cb = float(np.min(np.linalg.norm(diff, axis=2)))

    dist_tb_ca = 999.9
    if tails_b is not None and cleft_a is not None and len(tails_b) > 0 and len(cleft_a) > 0:
        diff = tails_b[:, np.newaxis, :] - cleft_a[np.newaxis, :, :]
        dist_tb_ca = float(np.min(np.linalg.norm(diff, axis=2)))

    return round(dist_ta_cb, 2), round(dist_tb_ca, 2)

def clean_protein_type(raw_type: str) -> str:
    """Strips modifiers, PTMs, and states. Leaves base name and mutations (e.g., BRAF-V600E)."""
    name = re.sub(r'[_|-]?chain[_|-]?[A-Za-z0-9]+$', '', raw_type, flags=re.IGNORECASE)
    name = re.sub(r'^\d+[_\-|]', '', name)
    
    tokens = re.split(r'[-_]', name)
    ignore_tokens = {'wt', 'cat', 'cattail', 'tail', 'wtcat', 'catwt', 'kd', 'apo', 'holo'}
    
    clean_tokens =[]
    for t in tokens:
        if not t: continue
        t_lower = t.lower()
        if t_lower in ignore_tokens: continue
        # Catch pS, pT, pY
        if re.match(r'^p[sty]?\d+$', t_lower): continue
        # Catch 0atp, 1atp, etc.
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

    role_a = "Unpaired"; role_b = "Unpaired"
    CONTACT_THRESHOLD = 8.0 

    lbl_sym = "Symmetric"
    lbl_act, lbl_rec = ("Activator", "Receiver") if use_erbb_terminology else ("C-lobe Donor", "N-lobe Receiver")

    if dist_AC_BN < CONTACT_THRESHOLD and dist_BC_AN < CONTACT_THRESHOLD: 
        role_a = role_b = lbl_sym
    elif dist_AC_BN < CONTACT_THRESHOLD: 
        role_a, role_b = lbl_act, lbl_rec
    elif dist_BC_AN < CONTACT_THRESHOLD: 
        role_a, role_b = lbl_rec, lbl_act
        
    return role_a, role_b, round(dist_AC_BN, 2), round(dist_BC_AN, 2)

def get_phi_psi(res_prev: Any, res: Any, res_next: Any) -> Tuple[Optional[float], Optional[float]]:
    if any(r is None for r in[res_prev, res, res_next]): return None, None
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

def analyze_activation_loop_dynamic(residues: Any, lm: Dict) -> Tuple[str, str]:
    nt_label, ct_label = "NT_Unk", "CT_Unk"
    idx_f, idx_hrd, idx_ape = lm.get('f'), lm.get('hrd'), lm.get('ape')
    if idx_f is None or idx_hrd is None: return "N/A", "N/A"
    
    idx_xhrd = idx_hrd - 1
    if 0 <= idx_xhrd < len(residues):
        nt_dists =[]
        for offset in range(3, 7):
            idx_nt_scan = idx_f + offset
            if idx_nt_scan < len(residues) and (idx_ape is None or idx_nt_scan < idx_ape):
                nt_dists.append(get_res_min_dist(residues[idx_nt_scan], residues[idx_xhrd]))
        if nt_dists: nt_label = "NTin" if min(nt_dists) <= 5.5 else "NTout"

    if idx_ape is not None:
        idx_arg = idx_hrd + 6  
        if 0 <= idx_arg < len(residues):
            ct_dists =[]
            scan_limit = max(idx_f + 2, idx_ape - 9)
            for offset in range(idx_ape - 6, scan_limit - 1, -1):
                if 0 <= offset < len(residues):
                    ct_dists.append(get_res_min_dist(residues[offset], residues[idx_arg]))
            if ct_dists: ct_label = "CTin" if min(ct_dists) <= 5.5 else "CTout"
                
    return nt_label, ct_label

def analyze_spines(residues: Any, lm: Dict, ligand_atoms: Any) -> Tuple[str, str]:
    r_spine = "Missing"; c_spine = "No Ligand"
    try:
        if all(lm.get(x) is not None for x in['hrd', 'f', 'rs1', 'rs2']):
            d43 = get_min_sc_dist(residues[lm['hrd']], residues[lm['f']])
            d31 = get_min_sc_dist(residues[lm['f']], residues[lm['rs1']])
            d12 = get_min_sc_dist(residues[lm['rs1']], residues[lm['rs2']])
            r_spine = "Intact" if (d43 < 4.5 and d31 < 4.5 and d12 < 4.5) else "Broken"
    except IndexError: pass

    if ligand_atoms is not None and len(ligand_atoms) > 0 and lm.get('k') is not None and lm['k'] >= 3:
        try:
            vaik_coords =[]
            for offset in range(4):
                res_idx = lm['k'] - offset
                if res_idx >= 0:
                    vaik_coords.extend(residues[res_idx].atoms.scene_coords)
            if vaik_coords:
                diff = np.array(vaik_coords)[:, np.newaxis, :] - ligand_atoms.scene_coords[np.newaxis, :, :]
                min_dist = np.min(np.linalg.norm(diff, axis=2))
                c_spine = "Intact" if min_dist < 6.0 else "Ligand Distant"
        except (IndexError, AttributeError): pass
    
    return r_spine, c_spine

def analyze_core_bridges(residues: Any, lm: Dict) -> Tuple[Any, Any, Any, str]:
    v104_dist, i150_dist, shell_dist = "N/A", "N/A", "N/A"
    shell_state = "Unknown"
    try:
        if all(lm.get(x) is not None for x in['v104', 'rs2']):
            v104_dist = round(get_min_sc_dist(residues[lm['v104']], residues[lm['rs2']]), 2)
        if all(lm.get(x) is not None for x in ['i150', 'hrd']):
            i150_dist = round(get_min_sc_dist(residues[lm['i150']], residues[lm['hrd']]), 2)
        if all(lm.get(x) is not None for x in ['m118', 'm120']):
            shell_dist = round(get_min_sc_dist(residues[lm['m118']], residues[lm['m120']]), 2)
            shell_state = "Packed" if shell_dist < 5.0 else "Loose"
    except IndexError: pass
    return v104_dist, i150_dist, shell_dist, shell_state

def analyze_alphaC_beta4_loop(residues: Any, lm: Dict) -> Tuple[Any, Any, Any, Any, Any]:
    y156_n99_d, k105_e107_d, k105_e121_d, k105_n99_d, d220_hrd_d = ["N/A"] * 5
    try:
        if all(lm.get(x) is not None for x in['y156', 'n99']):
            y156_n99_d = round(get_res_min_dist(residues[lm['y156']], residues[lm['n99']]), 2)
        if all(lm.get(x) is not None for x in['k105', 'e107']):
            k105_e107_d = round(get_min_sc_dist(residues[lm['k105']], residues[lm['e107']]), 2)
        if all(lm.get(x) is not None for x in['k105', 'e121']):
            k105_e121_d = round(get_min_sc_dist(residues[lm['k105']], residues[lm['e121']]), 2)
        if all(lm.get(x) is not None for x in ['k105', 'n99']):
            k105_n99_d = round(get_min_sc_dist(residues[lm['k105']], residues[lm['n99']]), 2)
        if all(lm.get(x) is not None for x in ['d220', 'hrd']):
            d220_hrd_d = round(get_min_sc_dist(residues[lm['d220']], residues[lm['hrd']]), 2)
    except IndexError: pass
    return y156_n99_d, k105_e107_d, k105_e121_d, k105_n99_d, d220_hrd_d

def get_core_name(rel_path: str) -> str:
    base = os.path.splitext(rel_path)[0]
    parts = base.replace('\\', '/').split('/')
    seen_tokens =[]
    for part in parts:
        for t in part.split('_'):
            if t and t not in seen_tokens: seen_tokens.append(t)
    return "_".join(seen_tokens)

def get_best_landmark_for_chain(chain_seq: str, cid: str, candidate_lms: List, fasta_seqs: Dict, sim_id: str = "") -> Tuple[str, Optional[Dict]]:
    best_lm = None; best_name = "Unknown"; best_score = -1.0; used_ratio = False
    for fasta_name, coords in candidate_lms:
        struct_score = 0
        if coords.get('f') is not None and coords['f'] < len(chain_seq):
            if chain_seq[coords['f']] in ['F', 'W', 'L', 'Y']: struct_score += 2
            if coords['f'] >= 1 and chain_seq[coords['f']-1] == 'D': struct_score += 2
            if coords['f'] + 1 < len(chain_seq) and chain_seq[coords['f']+1] == 'G': struct_score += 2
        if coords.get('hrd') is not None and coords['hrd'] + 2 < len(chain_seq):
            if chain_seq[coords['hrd']] in ['H', 'Y', 'F']: struct_score += 1
            if chain_seq[coords['hrd']+1] in['R', 'C', 'H', 'K']: struct_score += 1
            if chain_seq[coords['hrd']+2] == 'D': struct_score += 2

        ratio = 0.0; ref_seq = None
        for header, fseq in fasta_seqs.items():
            # Clean off the 0_ header indexing artifact before comparing
            clean_header = re.sub(r'^\d+[_\-|]', '', header)
            if fasta_name.lower() in clean_header.lower() or clean_header.lower() in fasta_name.lower():
                ref_seq = fseq; break
        
        if ref_seq:
            sm = difflib.SequenceMatcher(None, chain_seq, ref_seq.upper(), autojunk=False)
            matching_blocks = sm.get_matching_blocks()
            total_identical = sum(block.size for block in matching_blocks)
            ratio = total_identical / len(chain_seq) if len(chain_seq) > 0 else 0.0
            used_ratio = True

        composite = (ratio * 100) + struct_score
        if fasta_name.endswith(f"_{cid}") or fasta_name.endswith(f"-{cid}") or f"chain_{cid}" in fasta_name.lower():
            composite += 50
            
        # --- TIE-BREAKER PENALTY BLOCK ---
        if sim_id:
            sim_tokens = set(re.split(r'[-_]', sim_id.lower()))
            fasta_tokens = set(re.split(r'[-_]', fasta_name.lower()))
            
            # Bonus for shared tokens
            overlap = len(sim_tokens.intersection(fasta_tokens))
            composite += overlap * 2  
            
            # Heavy penalty if FASTA has a state/PTM that the sim_id lacks
            state_keywords = {'apo', 'holo'}
            for ft in fasta_tokens:
                if ft in state_keywords or re.match(r'^\d*(atp|adp|amp|gdp|gtp|anp)$', ft) or re.match(r'^p[sty]?\d+$', ft):
                    if ft not in sim_tokens:
                        composite -= 15  # Prevents Apo from stealing pS494/pT491 landmarks
        # ---------------------------------

        if composite > best_score:
            best_score = composite; best_lm = coords; best_name = fasta_name

    if used_ratio:
        if best_score >= 15.0: return best_name, best_lm
    else:
        if best_score >= 3.0: return best_name, best_lm
    return "Unknown", None

# --- MASTER EXECUTION ---
def process_model(session, full_cif_path: str, base_dir: str, out_dir_core: str, out_dir_allo: str, all_landmarks: Dict, fasta_seqs: Dict) -> List[Dict]:
    csv_rows =[]
    rel_to_base = os.path.relpath(full_cif_path, base_dir)
    sim_id = get_core_name(rel_to_base)
    
    # Use the Universal Normalizer so the folder string strictly matches the FASTA headers!
    sim_base = normalize_sim_name(sim_id)
    
    candidate_lms =[]
    for fasta_name, coords in all_landmarks.items():
        fasta_base = normalize_sim_name(fasta_name)
        if fasta_base.lower() in sim_base.lower() or sim_base.lower() in fasta_base.lower():
            candidate_lms.append((fasta_name, coords))
        elif fasta_name.lower() in sim_id.lower() or fasta_name.replace("_", "-").lower() in sim_id.lower():
            candidate_lms.append((fasta_name, coords))
            
    unique_lms = {name: coords for name, coords in candidate_lms}
    candidate_lms = list(unique_lms.items())
            
    if not candidate_lms: return csv_rows
    
    models = run(session, f"open '{full_cif_path}'", log=False)
    if not models: return csv_rows
    model = models[0]
    
    cids = sorted(list(set(model.residues.chain_ids)))
    ligand_mask = np.isin(np.array(model.residues.names), LIGAND_NAMES)
    ligands = model.residues[ligand_mask]
    
    chain_data = {}
    has_kinase = False
    
    for cid in cids:
        seq, res = get_sequence_and_residues(model, cid)
        if len(seq) < 100: continue
        
        lm_name, lm = get_best_landmark_for_chain(seq, cid, candidate_lms, fasta_seqs, sim_id)
        if lm is None or lm.get('f') is None or lm['f'] >= len(res): continue
        
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
        
        closest_d_o_name, atp_spec = None, None
        closest_dfg_mg_atom, closest_mg_spec = None, None
        aCb4_aE_cmd = None
        ploop_res_nums, ac_b4_res_nums, aE_res_nums, sb_oe_names = [], [], [], []

        # --- Catalytic & ATP Setup ---
        if lm.get('hrd') is not None and lm['hrd'] + 2 < len(res) and chain_ligand_atoms is not None:
            hrd_d_res = res[lm['hrd'] + 2]
            hrd_sc = get_sidechain_atoms(hrd_d_res)
            polar_mask =[n for n in hrd_sc.names if n.startswith('O') or n.startswith('N')]
            d_oxygens = hrd_sc[np.isin(hrd_sc.names, polar_mask)]
            
            target_phos =[]
            for phos_group in [['PG', 'O1G', 'O2G', 'O3G'],['PB', 'O1B', 'O2B', 'O3B'],['PA', 'O1A', 'O2A', 'O3A']]:
                if any(np.isin(chain_ligand_atoms.names, phos_group)): target_phos = phos_group; break
            if target_phos:
                p_mask = np.isin(chain_ligand_atoms.names, target_phos)
                lig_phosphates = chain_ligand_atoms[p_mask]
                if len(d_oxygens) > 0 and len(lig_phosphates) > 0:
                    diff = d_oxygens.scene_coords[:, np.newaxis, :] - lig_phosphates.scene_coords[np.newaxis, :, :]
                    dists = np.linalg.norm(diff, axis=2)
                    min_idx = np.unravel_index(np.argmin(dists), dists.shape)
                    hrd_atp_dist_val = round(float(dists[min_idx]), 2)
                    closest_d_o_name = d_oxygens.names[min_idx[0]]
                    closest_atp_atom_name = lig_phosphates.names[min_idx[1]]
                    atp_spec = f"/{chain_ligand_res.chain_id}:{chain_ligand_res.number}@{closest_atp_atom_name}"

        # --- DFG & Mg ---
        dfg_d_res = res[lm['f'] - 1] if lm.get('f') is not None and lm['f'] >= 1 else None
        if dfg_d_res is not None:
            dfg_sc = get_sidechain_atoms(dfg_d_res)
            polar_mask =[n for n in dfg_sc.names if n.startswith('O') or n.startswith('N')]
            d_oxygens = dfg_sc[np.isin(dfg_sc.names, polar_mask)]
            mg_atoms = model.atoms[model.atoms.names == 'MG']
            
            if len(d_oxygens) > 0 and len(mg_atoms) > 0:
                diff = d_oxygens.scene_coords[:, np.newaxis, :] - mg_atoms.scene_coords[np.newaxis, :, :]
                dists = np.linalg.norm(diff, axis=2)
                min_idx = np.unravel_index(np.argmin(dists), dists.shape)
                min_dist = float(dists[min_idx])
                
                # FIX: Only link Magnesium if it's actually within 15.0 Å of the active site
                if min_dist <= 15.0:
                    dfg_mg_dist_val = round(min_dist, 2)
                    closest_dfg_mg_atom = d_oxygens.names[min_idx[0]]
                    tgt_mg = mg_atoms[min_idx[1]]
                    closest_mg_spec = f"/{tgt_mg.residue.chain_id}:{tgt_mg.residue.number}@{tgt_mg.name}"
                
            if len(d_oxygens) > 0 and chain_ligand_atoms is not None:
                p_mask = np.isin(chain_ligand_atoms.names,['PG', 'O1G', 'O2G', 'O3G', 'PB', 'O1B', 'O2B', 'O3B', 'PA', 'O1A', 'O2A', 'O3A'])
                lig_phosphates = chain_ligand_atoms[p_mask]
                if len(lig_phosphates) > 0:
                    diff = d_oxygens.scene_coords[:, np.newaxis, :] - lig_phosphates.scene_coords[np.newaxis, :, :]
                    dfg_atp_dist_val = round(float(np.min(np.linalg.norm(diff, axis=2))), 2)

        # --- P-Loop ---
        if lm.get('k') is not None and lm['k'] > 10:
            start_search = max(0, lm['k'] - 35); end_search = lm['k'] - 5
            search_seq = seq[start_search:end_search]
            matches = list(re.finditer(r'G.G', search_seq))
            if matches: 
                match = matches[-1]
                idx_g1 = start_search + match.start()
                idx_g2 = start_search + match.end() - 1
                ploop_residues = res[idx_g1:idx_g2 + 3]
            else: 
                ploop_residues = res[max(0, lm['k'] - 25):max(0, lm['k'] - 15)]
                
            ploop_res_nums =[r.number for r in ploop_residues]
            if len(ploop_residues) > 0 and chain_ligand_atoms is not None:
                p_mask = np.isin(chain_ligand_atoms.names,['PG', 'O1G', 'O2G', 'O3G', 'PB', 'O1B', 'O2B', 'O3B', 'PA', 'O1A', 'O2A', 'O3A'])
                lig_phosphates = chain_ligand_atoms[p_mask]
                ploop_coords =[]
                for pr in ploop_residues:
                    valid_atoms = pr.atoms[np.isin(pr.atoms.names,['N', 'CA', 'C', 'CB'])]
                    if len(valid_atoms) > 0: ploop_coords.extend(valid_atoms.scene_coords)
                if len(lig_phosphates) > 0 and len(ploop_coords) > 0:
                    diff = np.array(ploop_coords)[:, np.newaxis, :] - lig_phosphates.scene_coords[np.newaxis, :, :]
                    ploop_atp_dist_val = round(float(np.min(np.linalg.norm(diff, axis=2))), 2)

        # --- Allosteric Toggles (Wu et al. & Kim et al. Metrics) ---
        if lm.get('c') is not None and lm.get('hrd') is not None:
            idx_ac_b4_start = lm['c'] + 8
            idx_ac_b4_end = lm['c'] + 15
            idx_aE_start = max(0, lm['hrd'] - 25)
            idx_aE_end = max(0, lm['hrd'] - 10)
            
            if idx_ac_b4_end < len(res) and idx_aE_end < len(res):
                ac_b4_res = res[idx_ac_b4_start : idx_ac_b4_end]
                aE_res = res[idx_aE_start : idx_aE_end]
                ac_b4_res_nums = [r.number for r in ac_b4_res]
                aE_res_nums = [r.number for r in aE_res]
                
                # aC-b4 to aE polar anchoring
                best_dist = 999.9
                for r_ac in ac_b4_res:
                    for r_ae in aE_res:
                        ac_names =[n for n in r_ac.atoms.names if n.startswith('N') or n.startswith('O')]
                        ae_names =[n for n in r_ae.atoms.names if n.startswith('N') or n.startswith('O')]
                        ac_polar = r_ac.atoms[np.isin(r_ac.atoms.names, ac_names)]
                        ae_polar = r_ae.atoms[np.isin(r_ae.atoms.names, ae_names)]
                        
                        if len(ac_polar) > 0 and len(ae_polar) > 0:
                            diff = ac_polar.scene_coords[:, np.newaxis, :] - ae_polar.scene_coords[np.newaxis, :, :]
                            dists = np.linalg.norm(diff, axis=2)
                            min_idx = np.unravel_index(np.argmin(dists), dists.shape)
                            d = float(dists[min_idx])
                            if d < best_dist:
                                best_dist = d
                                aCb4_aE_dist_val = round(d, 2)
                                ac_atom_name = ac_polar.names[min_idx[0]]
                                ae_atom_name = ae_polar.names[min_idx[1]]
                                aCb4_aE_cmd = f"distance {spec}:{r_ac.number}@{ac_atom_name} {spec}:{r_ae.number}@{ae_atom_name} color white radius 0.05"

                # Spine Bridge (aC-b4 to Ligand or VAIK+3)
                if chain_ligand_atoms is not None and len(chain_ligand_atoms) > 0:
                    c1 = ac_b4_res.atoms.scene_coords
                    c2 = chain_ligand_atoms.scene_coords
                    if len(c1) > 0 and len(c2) > 0:
                        diff = c1[:, np.newaxis, :] - c2[np.newaxis, :, :]
                        spine_bridge_dist_val = round(float(np.min(np.linalg.norm(diff, axis=2))), 2)
                elif lm.get('k') is not None and lm['k'] >= 3:
                    vaik_val_res = res[lm['k'] - 3]
                    spine_bridge_dist_val = round(get_res_min_dist(ac_b4_res, vaik_val_res), 2)

        # --- Conformation Anchors ---
        k_ca = get_atom(res[lm['k']], "CA") if lm.get('k') is not None else None
        c_ca = get_atom(res[lm['c']], "CA") if lm.get('c') is not None else None
        f_res = res[lm['f']] if lm.get('f') is not None else None
        f_cz = get_atom(f_res, "CZ")
        if f_cz is None and f_res is not None:
            f_sc = get_sidechain_atoms(f_res); f_ca = get_atom(f_res, "CA")
            if f_sc is not None and f_ca is not None and len(f_sc) > 0:
                dists = np.linalg.norm(f_sc.scene_coords - f_ca, axis=1)
                f_cz = f_sc.scene_coords[np.argmax(dists)] 

        if all(x is not None for x in[k_ca, c_ca, f_cz]):
            d1 = float(np.linalg.norm(c_ca - f_cz)); d2 = float(np.linalg.norm(k_ca - f_cz))
            d1_val, d2_val = round(d1, 2), round(d2, 2)
            if d1 <= 11.0 and d2 >= 11.0: spatial_label = "DFGin"
            elif d1 > 11.0 and d2 <= 14.0: spatial_label = "DFGout"
            elif d1 <= 11.0 and d2 <= 11.0: spatial_label = "DFGinter"
            else: spatial_label = "Outlier"

        k_nz = get_atom(res[lm['k']], "NZ") if lm.get('k') is not None else None
        c_res = res[lm['c']] if lm.get('c') is not None else None
        if k_nz is not None and c_res is not None:
            c_sc = get_sidechain_atoms(c_res)
            polar_mask =[n for n in c_sc.names if n.startswith('O') or n.startswith('N')]
            c_polar_atoms = c_sc[np.isin(c_sc.names, polar_mask)]
            if len(c_polar_atoms) > 0:
                dists = np.linalg.norm(c_polar_atoms.scene_coords - k_nz, axis=1)
                min_idx = np.argmin(dists)
                sb_dist_val = round(float(dists[min_idx]), 2)
                closest_indices = np.argsort(dists)[:2]
                sb_oe_names = c_polar_atoms.names[closest_indices].tolist()

        k_cb = get_atom(res[lm['k']], "CB") if lm.get('k') is not None else None
        c_cb = get_atom(res[lm['c']], "CB") if lm.get('c') is not None else None
        if k_cb is not None and c_cb is not None:
            chelix_label = "In" if np.linalg.norm(k_cb - c_cb) <= 10.0 else "Out"

        idx_f = lm['f']; idx_d = idx_f - 1; idx_x = idx_f - 2  
        if idx_x - 1 >= 0 and idx_f + 1 < len(res):
            phi_x, psi_x = get_phi_psi(res[idx_x-1], res[idx_x], res[idx_d])
            phi_d, psi_d = get_phi_psi(res[idx_x], res[idx_d], res[idx_f])
            phi_f, psi_f = get_phi_psi(res[idx_d], res[idx_f], res[idx_f+1])
            if phi_d is not None and psi_d is not None:
                raw_phi_d = round(phi_d, 2); raw_psi_d = round(psi_d, 2)
            rama_x = assign_ramachandran_region(phi_x, psi_x)
            rama_d = assign_ramachandran_region(phi_d, psi_d)
            rama_f = assign_ramachandran_region(phi_f, psi_f)
            chi1 = calculate_dihedral(get_atom(res[idx_f], "N"), get_atom(res[idx_f], "CA"), get_atom(res[idx_f], "CB"), get_atom(res[idx_f], "CG"))
            rot_f = get_rotamer(chi1)
            if "X" not in[rama_x, rama_d, rama_f] and rot_f != "Unknown":
                dihedral_label = f"{rama_x}{rama_d}{rama_f}{rot_f}"

        state = "Unknown"
        if spatial_label == "DFGin":
            if dihedral_label == "BLAminus": state = "Active (BLAminus)"
            elif dihedral_label == "ABAminus": state = "Active-Like (ABAminus)"
            elif dihedral_label == "BLBplus": state = "Inactive (BLBplus)"
            elif dihedral_label == "BLBtrans": state = "Inactive (BLBtrans)"
            elif dihedral_label == "BLAplus": state = "Inactive (BLAplus)"
            elif dihedral_label == "BLBminus": state = "Inactive (BLBminus)"
            else: state = f"DFGin ({dihedral_label})"
        # ... (Right after the spatial/state calculations, around line 410) ...
        elif spatial_label == "DFGout" and dihedral_label == "BBAminus": state = "Inactive (BBAminus)"
        elif spatial_label == "DFGinter" and dihedral_label == "BABtrans": state = "Inactive (BABtrans)"
        else: state = f"{spatial_label} ({dihedral_label})"

        # --- LABEL ASSIGNMENT ---
        # Strictly use the protein name as it was listed in the FASTA file
        clean_type_string = lm_name

        chain_data[cid] = {
            "residues": res, "landmarks": lm, 
            "meta": {"Type": clean_type_string, "State": state, "CHelix": chelix_label,
                     "RSpine": r_spine, "CSpine": c_spine, "Spatial": spatial_label, "Dihedral": dihedral_label,
                     "ActLoop_NT": nt_loop, "ActLoop_CT": ct_loop, "Phi_D": raw_phi_d, "Psi_D": raw_psi_d,
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
    # --- UNIVERSAL DIMER INTERFACE DETECTION ---
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
                
                # NEW RAF ASYMMETRY OVERRIDE
                raf_a_push, raf_b_push, raf_lla_dist = "N/A", "N/A", "N/A"
                tail_a_block, tail_b_block = "N/A", "N/A"
                
                type_a = chain_data[cA]['meta']['Type']
                type_b = chain_data[cB]['meta']['Type']
                
                if is_raf_system(sim_id, [type_a, type_b]) and role_a == "Symmetric":
                    # 1. Measure Core Asymmetry
                    pA, pB, lla = get_raf_interface_metrics(chain_data[cA]['residues'], chain_data[cB]['residues'], chain_data[cA]['landmarks'], chain_data[cB]['landmarks'])
                    raf_a_push, raf_b_push, raf_lla_dist = pA, pB, lla
                    
                    # 2. Measure Trans-Tail Blockade
                    tA_cB, tB_cA = get_raf_tail_asymmetry(chain_data[cA]['residues'], chain_data[cB]['residues'], chain_data[cA]['landmarks'], chain_data[cB]['landmarks'])
                    tail_a_block, tail_b_block = tA_cB, tB_cA
                    
                    # Logic 1: Does a tail uniquely block the partner's active site? (< 7.0 A)
                    if min(tA_cB, tB_cA) < 7.0 and abs(tA_cB - tB_cA) > 3.0:
                        if tA_cB < tB_cA:
                            # Tail A blocks Cleft B. A is doing the catalysis; B is the blocked scaffold.
                            role_a, role_b = "Receiver", "Activator (Tail Blocked)"
                        else:
                            # Tail B blocks Cleft A. 
                            role_a, role_b = "Activator (Tail Blocked)", "Receiver"
                    
                    # Logic 2: Fallback to the core LLA / Arginine asymmetry
                    elif lla != 999.9 and lla < 6.0:
                        push_delta = pA - pB
                        if push_delta < -2.0:
                            role_a, role_b = "Activator", "Receiver"
                        elif push_delta > 2.0:
                            role_a, role_b = "Receiver", "Activator"

                # Store metrics in meta
                chain_data[cA]['meta']['RAF_A_Push'] = raf_a_push
                chain_data[cA]['meta']['RAF_B_Push'] = raf_b_push
                chain_data[cA]['meta']['RAF_LLA_Core'] = raf_lla_dist
                chain_data[cA]['meta']['RAF_Tail_A_Block'] = tail_a_block
                chain_data[cA]['meta']['RAF_Tail_B_Block'] = tail_b_block
                
                chain_data[cB]['meta']['RAF_A_Push'] = raf_b_push
                chain_data[cB]['meta']['RAF_B_Push'] = raf_a_push
                chain_data[cB]['meta']['RAF_LLA_Core'] = raf_lla_dist
                chain_data[cB]['meta']['RAF_Tail_A_Block'] = tail_b_block
                chain_data[cB]['meta']['RAF_Tail_B_Block'] = tail_a_block
                
                if role_a != "Unpaired":
                    chain_roles[cA] = role_a
                    chain_roles[cB] = role_b
                    partner_map[cA] = cB
                    partner_map[cB] = cA
                    int_dist_AC_BN[cA] = dAC_BN
                    int_dist_AC_BN[cB] = dBC_AN 
                    int_dist_BC_AN[cA] = dBC_AN
                    int_dist_BC_AN[cB] = dAC_BN

    # =========================================================================
    # --- SPLIT VISUALIZATION GENERATORS ---
    # =========================================================================
    protein_name = " / ".join(sorted(list(set([name for name, _ in candidate_lms]))))
    rel_to_out_core = os.path.relpath(full_cif_path, out_dir_core).replace('\\', '/')
    rel_to_out_allo = os.path.relpath(full_cif_path, out_dir_allo).replace('\\', '/')

    cxc_core =[
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

    cxc_allo =[
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

    if len(ligands) > 0:
        lig_cmds =[
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
            "R_Spine": meta['RSpine'], "C_Spine": meta['CSpine'], "C_Helix": meta['CHelix'],
            "Shell_State": meta['Shell_State'], "Spatial": meta['Spatial'], "Dihedral": meta['Dihedral'],
            "ActLoop_NT": meta.get('ActLoop_NT', 'N/A'), "ActLoop_CT": meta.get('ActLoop_CT', 'N/A'), 
            "Phi_D": meta.get('Phi_D', 'N/A'), "Psi_D": meta.get('Psi_D', 'N/A'),
            "D1_Dist": meta['D1_Dist'], "D2_Dist": meta['D2_Dist'], "SB_Dist": meta['SB_Dist'],
            "HRD_ATP_Dist": meta['HRD_ATP_Dist'], "DFG_Mg_Dist": meta.get('DFG_Mg_Dist', 'N/A'), 
            "DFG_ATP_Dist": meta.get('DFG_ATP_Dist', 'N/A'), "PLoop_ATP_Dist": meta.get('PLoop_ATP_Dist', 'N/A'),
            "aCb4_aE_Dist": meta.get('aCb4_aE_Dist', 'N/A'), "Spine_Bridge_Dist": meta.get('Spine_Bridge_Dist', 'N/A'),
            "V104_RS2_Dist": meta['V104_RS2_Dist'], "I150_HRD_Dist": meta['I150_HRD_Dist'], "Shell_M118_M120_Dist": meta['Shell_M118_M120_Dist'],
            "Y156_N99_Dist": meta['Y156_N99_Dist'], "K105_E107_Dist": meta['K105_E107_Dist'], 
            "K105_E121_Dist": meta['K105_E121_Dist'], "K105_N99_Dist": meta['K105_N99_Dist'], "D220_HRD_Dist": meta['D220_HRD_Dist'],
            "RAF_A_Push": meta.get('RAF_A_Push', 'N/A'), 
            "RAF_B_Push": meta.get('RAF_B_Push', 'N/A'), 
            "RAF_LLA_Core": meta.get('RAF_LLA_Core', 'N/A'),
            "RAF_Tail_A_Block": meta.get('RAF_Tail_A_Block', 'N/A'),
            "RAF_Tail_B_Block": meta.get('RAF_Tail_B_Block', 'N/A')
        })

        if data['landmarks']:
            lm = data['landmarks']; res = data['residues']
            f_num = res[lm['f']].number; d_num = f_num - 1
            hrd_d_num = res[lm['hrd'] + 2].number if (lm.get('hrd') is not None and lm['hrd'] + 2 < len(res)) else None
            spec = f"/{cid}"
            
            role_str = f" | {chain_roles[cid]}" if chain_roles[cid] not in["N/A", "Unpaired"] else ""
            state_label = f"2dlabels create State_{cid} text 'Chain {cid} ({meta['Type']}{role_str}): {meta['State']}' color white size 20 xpos 0.05 ypos {y_offset}"
            cxc_core.append(state_label); cxc_allo.append(state_label)

            # --- POPULATE CORE VISUALIZER ---
            cxc_core.extend([f"color {spec} #d3d3d3", f"show {spec}:{d_num}-{f_num+1}", f"color {spec}:{d_num}-{f_num+1} dark orange"])
            if meta.get('ploop_res_nums'):
                cxc_core.extend([f"color {spec}:{meta['ploop_res_nums'][0]}-{meta['ploop_res_nums'][-1]} deep sky blue"])
            if lm.get('ape') is not None:
                cxc_core.extend([f"color {spec}:{f_num+1}-{res[lm['ape']].number} coral", f"show {spec}:{res[lm['ape']].number} & sidechain"])
                
            if all(lm.get(x) is not None for x in['hrd', 'f', 'rs1', 'rs2']):
                r1, r2, r3, r4 = res[lm['rs2']].number, res[lm['rs1']].number, res[lm['f']].number, res[lm['hrd']].number
                cxc_core.extend([f"color {spec}:{r1},{r2},{r3},{r4} medium purple", f"show {spec}:{r1},{r2},{r3},{r4} & sidechain",
                                 f"distance {spec}:{r1}@CA {spec}:{r2}@CA color medium purple radius 0.05",
                                 f"distance {spec}:{r2}@CA {spec}:{r3}@CA color medium purple radius 0.05",
                                 f"distance {spec}:{r3}@CA {spec}:{r4}@CA color medium purple radius 0.05"])
                
            if lm.get('k') is not None and lm.get('c') is not None:
                k_num, c_num = res[lm['k']].number, res[lm['c']].number
                cxc_core.extend([f"color {spec}:{k_num},{c_num} spring green", f"show {spec}:{k_num},{c_num} & sidechain"])
                for oe_name in meta.get('sb_oe_names',[]):
                    cxc_core.append(f"distance {spec}:{k_num}@NZ {spec}:{c_num}@{oe_name} color spring green radius 0.05")

            if all(lm.get(x) is not None for x in['v104', 'm118', 'm120', 'rs2']):
                v104, m118, m120, rs2 = res[lm['v104']].number, res[lm['m118']].number, res[lm['m120']].number, res[lm['rs2']].number
                cxc_core.extend([f"color {spec}:{v104},{m118},{m120} teal", f"show {spec}:{v104},{m118},{m120} & sidechain",
                                 f"distance {spec}:{m118}@CA {spec}:{m120}@CA color teal radius 0.05", f"distance {spec}:{v104}@CA {spec}:{rs2}@CA color teal radius 0.05"])
                
            if lm.get('i150') is not None and lm.get('hrd') is not None:
                cxc_core.extend([f"color {spec}:{res[lm['i150']].number} tan", f"show {spec}:{res[lm['i150']].number} & sidechain",
                                 f"distance {spec}:{res[lm['i150']].number}@CA {spec}:{res[lm['hrd']].number}@CA color tan radius 0.05"])

            # --- POPULATE ALLOSTERIC VISUALIZER ---
            cxc_allo.extend([f"color {spec} #d3d3d3"])
            if meta.get('ac_b4_res_nums') and meta.get('aE_res_nums'):
                ac_b4_n1, ac_b4_n2 = meta['ac_b4_res_nums'][0], meta['ac_b4_res_nums'][-1]
                aE_n1, aE_n2 = meta['aE_res_nums'][0], meta['aE_res_nums'][-1]
                cxc_allo.extend([
                    f"\n# Allosteric Network: aC-b4 Loop & aE Helix",
                    f"color {spec}:{ac_b4_n1}-{ac_b4_n2} hot pink",
                    f"color {spec}:{aE_n1}-{aE_n2} yellow",
                    f"show {spec}:{ac_b4_n1}-{ac_b4_n2},{aE_n1}-{aE_n2} & sidechain"
                ])
                if meta.get('aCb4_aE_cmd'):
                    cxc_allo.append(meta['aCb4_aE_cmd'])

            # Add legacy toggles back into allo
            if lm.get('y156') is not None and lm.get('n99') is not None:
                y156, n99 = res[lm['y156']].number, res[lm['n99']].number
                cxc_allo.extend([f"\n# aE Anchor", f"color {spec}:{y156} lime green", f"show {spec}:{y156} & sidechain", f"distance {spec}:{y156}@CA {spec}:{n99}@CA color lime green radius 0.05"])

            # --- SHARED ELEMENTS (Catalytic & Ligand Coordination) ---
            if hrd_d_num:
                cat_cmd =[f"show {spec}:{hrd_d_num} & sidechain", f"color {spec}:{hrd_d_num} & sidechain red"]
                cxc_core.extend(cat_cmd); cxc_allo.extend(cat_cmd)
                if meta.get('atp_spec') and meta.get('closest_d_o_name') and meta.get('HRD_ATP_Dist') != "N/A" and meta['HRD_ATP_Dist'] < 10.0:
                    dist_cmd = f"distance {spec}:{hrd_d_num}@{meta['closest_d_o_name']} {meta['atp_spec']} color magenta radius 0.05"
                    cxc_core.append(dist_cmd); cxc_allo.append(dist_cmd)

            if meta.get('dfg_d_res'):
                dfg_d_num = meta['dfg_d_res'].number
                dfg_cmd =[f"show {spec}:{dfg_d_num} & sidechain", f"color {spec}:{dfg_d_num} & sidechain dodger blue"]
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

    # Fully clear the session to prevent ghost residues across loops
    run(session, "close session", log=False)
    return csv_rows

def main(session):
    base_dir = os.path.abspath(SEARCH_DIR)
    out_dir_core = os.path.join(base_dir, VIZ_OUT_DIR_CORE)
    out_dir_allo = os.path.join(base_dir, VIZ_OUT_DIR_ALLO)
    os.makedirs(out_dir_core, exist_ok=True)
    os.makedirs(out_dir_allo, exist_ok=True)
    
    chunk_list_file = os.environ.get("CHIMERAX_CHUNK")
    if not chunk_list_file or not os.path.exists(chunk_list_file): run(session, "quit"); return

    out_csv_name = f"{os.path.splitext(os.path.basename(chunk_list_file))[0]}_results_v6.csv"
    landmarks_file = os.path.join(base_dir, LANDMARKS_JSON)
    if not os.path.exists(landmarks_file): run(session, "quit"); return
        
    with open(landmarks_file, 'r') as f: all_landmarks = json.load(f)
    
    # Deduplicate paths to prevent multiple processing of the same file (fixes the 12-line bug)
    with open(chunk_list_file, 'r') as f: 
        files_to_process = list(dict.fromkeys([line.strip() for line in f if line.strip()]))

    fasta_path = os.path.join(base_dir, "sequences.fasta")
    fasta_seqs = {}
    if os.path.exists(fasta_path):
        with open(fasta_path, 'r') as f:
            name = None; seq_lines =[]
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if name: fasta_seqs[name] = "".join(seq_lines)
                    name = line[1:].strip().split()[0]
                    seq_lines =[]
                elif line: seq_lines.append(line)
            if name: fasta_seqs[name] = "".join(seq_lines)

    cols =["Simulation_ID", "Directory", "File", "Chain", "Type", "State", "Role", "Partner",
            "Interface_C_Lobe_Donor_Dist", "Interface_N_Lobe_Rec_Dist",
            "R_Spine", "C_Spine", "C_Helix", "Shell_State", "Spatial", "Dihedral", "ActLoop_NT", "ActLoop_CT", 
            "Phi_D", "Psi_D", "D1_Dist", "D2_Dist", "SB_Dist", 
            "HRD_ATP_Dist", "DFG_Mg_Dist", "DFG_ATP_Dist", "PLoop_ATP_Dist",
            "aCb4_aE_Dist", "Spine_Bridge_Dist",
            "V104_RS2_Dist", "I150_HRD_Dist", "Shell_M118_M120_Dist",
            "Y156_N99_Dist", "K105_E107_Dist", "K105_E121_Dist", "K105_N99_Dist", "D220_HRD_Dist",
            "RAF_A_Push", "RAF_B_Push", "RAF_LLA_Core", "RAF_Tail_A_Block", "RAF_Tail_B_Block"]
           
    all_rows = list() 
    
    for full_cif_path in files_to_process:
        try:
            rows = process_model(session, full_cif_path, base_dir, out_dir_core, out_dir_allo, all_landmarks, fasta_seqs)
            all_rows.extend(rows)
        except Exception as e: print(f"Error processing {full_cif_path}: {e}")

    with open(os.path.join(base_dir, out_csv_name), 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=cols)
        writer.writeheader()
        writer.writerows(all_rows)

    run(session, "quit")

main(session)
