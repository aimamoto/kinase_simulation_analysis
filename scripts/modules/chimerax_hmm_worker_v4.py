import os
import csv
import json
import numpy as np
import re
import difflib  
from chimerax.core.commands import run

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
LIGAND_NAMES =["ATP", "ADP", "ANP", "ACP", "AGS", "AMP", "GTP", "GDP"]

# --- MATH & ATOM HELPERS ---
def get_base_sim_name(name):
    clean = re.sub(r'_?(model|seed|rank|pred|unrelaxed|relaxed)[\-_]?\d+', '', name, flags=re.IGNORECASE)
    clean = re.sub(r'_?(unrelaxed|relaxed)', '', clean, flags=re.IGNORECASE)
    clean = re.sub(r'^(fold|job|run)_', '', clean, flags=re.IGNORECASE)
    return clean.strip('_')

def get_atom(residue, atom_name):
    if residue is None: return None
    atoms = residue.atoms[residue.atoms.names == atom_name]
    if len(atoms) > 0: return atoms[0].scene_coord
    return None

def get_sidechain_atoms(residue):
    if residue is None: return None
    bb_names = {'N', 'CA', 'C', 'O', 'OXT'}
    mask = ~np.isin(residue.atoms.names, list(bb_names))
    sc_atoms = residue.atoms[mask]
    if len(sc_atoms) == 0: return residue.atoms[residue.atoms.names == 'CA']
    return sc_atoms

def get_center_of_mass(residues):
    if hasattr(residues, 'scene_coords'): coords = residues.scene_coords
    elif hasattr(residues, 'atoms'): coords = residues.atoms.scene_coords
    else: return None
    if len(coords) == 0: return None
    return np.mean(coords, axis=0)

def calculate_dihedral(p1, p2, p3, p4):
    if any(x is None for x in[p1, p2, p3, p4]): return None
    b0 = -1.0 * (p2 - p1); b1 = p3 - p2; b2 = p4 - p3
    b1_norm = np.linalg.norm(b1)
    if b1_norm == 0: return None
    b1 /= b1_norm
    v = b0 - np.dot(b0, b1) * b1; w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w); y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))

def get_min_sc_dist(res1, res2):
    sc1 = get_sidechain_atoms(res1); sc2 = get_sidechain_atoms(res2)
    if sc1 is None or sc2 is None or len(sc1) == 0 or len(sc2) == 0: return 999.9
    diff = sc1.scene_coords[:, np.newaxis, :] - sc2.scene_coords[np.newaxis, :, :]
    return np.min(np.linalg.norm(diff, axis=2))

def get_res_min_dist(res1, res2):
    if res1 is None or res2 is None: return 999.9
    c1 = res1.atoms.scene_coords
    c2 = res2.atoms.scene_coords
    if len(c1) == 0 or len(c2) == 0: return 999.9
    diff = c1[:, np.newaxis, :] - c2[np.newaxis, :, :]
    return float(np.min(np.linalg.norm(diff, axis=2)))

def get_phi_psi(res_prev, res, res_next):
    if any(r is None for r in [res_prev, res, res_next]): return None, None
    phi = calculate_dihedral(get_atom(res_prev, "C"), get_atom(res, "N"), get_atom(res, "CA"), get_atom(res, "C"))
    psi = calculate_dihedral(get_atom(res, "N"), get_atom(res, "CA"), get_atom(res, "C"), get_atom(res_next, "N"))
    return phi, psi

def assign_ramachandran_region(phi, psi):
    if phi is None or psi is None: return "X"
    phi = (phi + 180) % 360 - 180; psi = (psi + 180) % 360 - 180
    if phi < 0: return "A" if -100 <= psi <= 50 else "B"
    else: return "L" if -50 <= psi <= 100 else "E"

def get_rotamer(chi1):
    if chi1 is None: return "Unknown"
    if chi1 < 0: chi1 += 360
    if 0 <= chi1 < 120: return "plus"
    elif 120 <= chi1 < 240: return "trans"
    else: return "minus"

def get_sequence_and_residues(model, chain_id):
    chain_res = model.residues[model.residues.chain_ids == chain_id]
    seq_str = "".join([AA_MAP.get(rname, 'X') for rname in chain_res.names])
    return seq_str, chain_res

def analyze_activation_loop_dynamic(residues, lm):
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

def analyze_spines(residues, lm, ligand_atoms):
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

def analyze_core_bridges(residues, lm):
    v104_dist, i150_dist, shell_dist = "N/A", "N/A", "N/A"
    shell_state = "Unknown"
    try:
        if all(lm.get(x) is not None for x in ['v104', 'rs2']):
            v104_dist = round(get_min_sc_dist(residues[lm['v104']], residues[lm['rs2']]), 2)
        if all(lm.get(x) is not None for x in ['i150', 'hrd']):
            i150_dist = round(get_min_sc_dist(residues[lm['i150']], residues[lm['hrd']]), 2)
        if all(lm.get(x) is not None for x in ['m118', 'm120']):
            shell_dist = round(get_min_sc_dist(residues[lm['m118']], residues[lm['m120']]), 2)
            shell_state = "Packed" if shell_dist < 5.0 else "Loose"
    except IndexError: pass
    return v104_dist, i150_dist, shell_dist, shell_state

def analyze_alphaC_beta4_loop(residues, lm):
    y156_n99_d, k105_e107_d, k105_e121_d, k105_n99_d, d220_hrd_d = ["N/A"] * 5
    try:
        if all(lm.get(x) is not None for x in ['y156', 'n99']):
            y156_n99_d = round(get_res_min_dist(residues[lm['y156']], residues[lm['n99']]), 2)
        if all(lm.get(x) is not None for x in ['k105', 'e107']):
            k105_e107_d = round(get_min_sc_dist(residues[lm['k105']], residues[lm['e107']]), 2)
        if all(lm.get(x) is not None for x in ['k105', 'e121']):
            k105_e121_d = round(get_min_sc_dist(residues[lm['k105']], residues[lm['e121']]), 2)
        if all(lm.get(x) is not None for x in ['k105', 'n99']):
            k105_n99_d = round(get_min_sc_dist(residues[lm['k105']], residues[lm['n99']]), 2)
        if all(lm.get(x) is not None for x in ['d220', 'hrd']):
            d220_hrd_d = round(get_min_sc_dist(residues[lm['d220']], residues[lm['hrd']]), 2)
    except IndexError: pass
    return y156_n99_d, k105_e107_d, k105_e121_d, k105_n99_d, d220_hrd_d

def get_core_name(rel_path):
    base = os.path.splitext(rel_path)[0]
    parts = base.replace('\\', '/').split('/')
    seen_tokens =[]
    for part in parts:
        for t in part.split('_'):
            if t and t not in seen_tokens: seen_tokens.append(t)
    return "_".join(seen_tokens)

def find_wt_reference(sim_id, wt_dir):
    if not wt_dir or not os.path.exists(wt_dir): return None
    prefix_match = re.split(r'[-_]', sim_id)[0].upper()
    for f in os.listdir(wt_dir):
        if f.endswith(".cif") or f.endswith(".pdb"):
            f_upper = f.upper()
            if prefix_match in f_upper and "WT" in f_upper:
                return os.path.abspath(os.path.join(wt_dir, f))
    return None

def get_best_landmark_for_chain(chain_seq, candidate_lms, fasta_seqs):
    best_lm = None
    best_name = "Unknown"
    best_score = -1
    best_ratio = -1.0
    used_ratio = False
    
    for fasta_name, coords in candidate_lms:
        ref_seq = None
        for header, fseq in fasta_seqs.items():
            if fasta_name.lower() == header.lower() or fasta_name.lower() in header.lower():
                ref_seq = fseq
                break
        
        if ref_seq:
            ratio = difflib.SequenceMatcher(None, chain_seq, ref_seq).quick_ratio()
            if ratio > best_ratio:
                best_ratio = ratio
                best_lm = coords
                best_name = fasta_name
            used_ratio = True
        
        if not used_ratio:
            score = 0
            if coords.get('f') is not None and coords['f'] < len(chain_seq):
                if chain_seq[coords['f']] in['F', 'W', 'L', 'Y']: score += 2
                if coords['f'] >= 1 and chain_seq[coords['f']-1] == 'D': score += 2
                if coords['f'] + 1 < len(chain_seq) and chain_seq[coords['f']+1] == 'G': score += 2
            if coords.get('hrd') is not None and coords['hrd'] + 2 < len(chain_seq):
                if chain_seq[coords['hrd']] in ['H', 'Y', 'F']: score += 1
                if chain_seq[coords['hrd']+1] in ['R', 'C', 'H', 'K']: score += 1
                if chain_seq[coords['hrd']+2] == 'D': score += 2
                
            if score > best_score:
                best_score = score
                best_lm = coords
                best_name = fasta_name

    if used_ratio and best_ratio > 0.30: return best_name, best_lm
    if not used_ratio and best_score >= 3: return best_name, best_lm
    return "Unknown", None

# --- MASTER EXECUTION ---
def process_model(session, full_cif_path, base_dir, out_dir_core, out_dir_allo, wt_dir, all_landmarks, fasta_seqs):
    csv_rows =[]
    rel_to_base = os.path.relpath(full_cif_path, base_dir)
    sim_id = get_core_name(rel_to_base)
    sim_base = get_base_sim_name(sim_id)
    
    candidate_lms =[]
    for fasta_name, coords in all_landmarks.items():
        fasta_base = get_base_sim_name(fasta_name)
        if fasta_base.lower() in sim_base.lower() or sim_base.lower() in fasta_base.lower():
            candidate_lms.append((fasta_name, coords))
        elif fasta_name.lower() in sim_id.lower() or fasta_name.replace("_", "-").lower() in sim_id.lower():
            candidate_lms.append((fasta_name, coords))
            
    unique_lms = {}
    for name, coords in candidate_lms:
        unique_lms[name] = coords
    candidate_lms = list(unique_lms.items())
            
    if not candidate_lms: 
        return csv_rows
    
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
        
        lm_name, lm = get_best_landmark_for_chain(seq, candidate_lms, fasta_seqs)
        if lm is None or lm.get('f') is None or lm['f'] >= len(res): continue
            
        has_kinase = True
        chain_ligand_res, chain_ligand_atoms = None, None
        
        if len(ligands) > 0:
            anchor_pos = get_atom(res[lm['f']], "CA")
            if anchor_pos is None and lm.get('k') is not None: anchor_pos = get_atom(res[lm['k']], "CA")
            if anchor_pos is not None:
                for lig_res in ligands:
                    if np.linalg.norm(anchor_pos - get_center_of_mass(lig_res.atoms)) < 15.0:
                        chain_ligand_res = lig_res
                        chain_ligand_atoms = lig_res.atoms
                        break

        r_spine, c_spine = analyze_spines(res, lm, chain_ligand_atoms)
        nt_loop, ct_loop = analyze_activation_loop_dynamic(res, lm)
        v104_d, i150_d, shell_d, shell_state = analyze_core_bridges(res, lm)
        y156_n99_d, k105_e107_d, k105_e121_d, k105_n99_d, d220_hrd_d = analyze_alphaC_beta4_loop(res, lm)
        
        spatial_label, dihedral_label = "Unknown", "Unknown"
        raw_phi_d, raw_psi_d = "N/A", "N/A"
        d1_val, d2_val, sb_dist_val, hrd_atp_dist_val = "N/A", "N/A", "N/A", "N/A"
        closest_d_o_name, atp_spec = None, None
        
        if lm.get('hrd') is not None and lm['hrd'] + 2 < len(res) and chain_ligand_atoms is not None:
            hrd_d_res = res[lm['hrd'] + 2]
            hrd_sc = get_sidechain_atoms(hrd_d_res)
            polar_mask =[n for n in hrd_sc.names if n.startswith('O') or n.startswith('N')]
            d_oxygens = hrd_sc[np.isin(hrd_sc.names, polar_mask)]
            
            target_phos = []
            for phos_group in [['PG', 'O1G', 'O2G', 'O3G'],['PB', 'O1B', 'O2B', 'O3B'],['PA', 'O1A', 'O2A', 'O3A']]:
                if any(np.isin(chain_ligand_atoms.names, phos_group)):
                    target_phos = phos_group; break
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

        dfg_mg_dist_val, dfg_atp_dist_val = "N/A", "N/A"
        closest_dfg_mg_atom = None
        closest_mg_spec = None
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
                dfg_mg_dist_val = round(float(dists[min_idx]), 2)
                closest_dfg_mg_atom = d_oxygens.names[min_idx[0]]
                tgt_mg = mg_atoms[min_idx[1]]
                closest_mg_spec = f"/{tgt_mg.residue.chain_id}:{tgt_mg.residue.number}@{tgt_mg.name}"
                
            if len(d_oxygens) > 0 and chain_ligand_atoms is not None:
                p_mask = np.isin(chain_ligand_atoms.names,['PG', 'O1G', 'O2G', 'O3G', 'PB', 'O1B', 'O2B', 'O3B', 'PA', 'O1A', 'O2A', 'O3A'])
                lig_phosphates = chain_ligand_atoms[p_mask]
                if len(lig_phosphates) > 0:
                    diff = d_oxygens.scene_coords[:, np.newaxis, :] - lig_phosphates.scene_coords[np.newaxis, :, :]
                    dfg_atp_dist_val = round(float(np.min(np.linalg.norm(diff, axis=2))), 2)

        ploop_atp_dist_val = "N/A"
        ploop_res_nums =[]
        if lm.get('k') is not None and lm['k'] > 10:
            start_search = max(0, lm['k'] - 35)
            end_search = lm['k'] - 5
            search_seq = seq[start_search:end_search]
            matches = list(re.finditer(r'G.G', search_seq))
            if matches: 
                match = matches[-1]
                idx_g1 = start_search + match.start()
                idx_g2 = start_search + match.end() - 1
                ploop_residues = res[idx_g1:idx_g2 + 3]
            else: 
                ploop_residues = res[max(0, lm['k'] - 25):max(0, lm['k'] - 15)]
                
            ploop_res_nums = [r.number for r in ploop_residues]
            if len(ploop_residues) > 0 and chain_ligand_atoms is not None:
                p_mask = np.isin(chain_ligand_atoms.names,['PG', 'O1G', 'O2G', 'O3G', 'PB', 'O1B', 'O2B', 'O3B', 'PA', 'O1A', 'O2A', 'O3A'])
                lig_phosphates = chain_ligand_atoms[p_mask]
                ploop_coords =[]
                for pr in ploop_residues:
                    valid_atoms = pr.atoms[np.isin(pr.atoms.names, ['N', 'CA', 'C', 'CB'])]
                    if len(valid_atoms) > 0: ploop_coords.extend(valid_atoms.scene_coords)
                if len(lig_phosphates) > 0 and len(ploop_coords) > 0:
                    diff = np.array(ploop_coords)[:, np.newaxis, :] - lig_phosphates.scene_coords[np.newaxis, :, :]
                    ploop_atp_dist_val = round(float(np.min(np.linalg.norm(diff, axis=2))), 2)

        k_ca = get_atom(res[lm['k']], "CA") if lm.get('k') is not None else None
        c_ca = get_atom(res[lm['c']], "CA") if lm.get('c') is not None else None
        f_res = res[lm['f']] if lm.get('f') is not None else None
        f_cz = get_atom(f_res, "CZ")
        if f_cz is None and f_res is not None:
            f_sc = get_sidechain_atoms(f_res)
            f_ca = get_atom(f_res, "CA")
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
        sb_oe_names =[]
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
        chelix_label = "N/A"
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
        elif spatial_label == "DFGout" and dihedral_label == "BBAminus": state = "Inactive (BBAminus)"
        elif spatial_label == "DFGinter" and dihedral_label == "BABtrans": state = "Inactive (BABtrans)"
        else: state = f"{spatial_label} ({dihedral_label})"

        chain_data[cid] = {
            "residues": res, "landmarks": lm, 
            "meta": {"Type": lm['type'], "State": state, "CHelix": chelix_label,
                     "RSpine": r_spine, "CSpine": c_spine, "Spatial": spatial_label, "Dihedral": dihedral_label,
                     "ActLoop_NT": nt_loop, "ActLoop_CT": ct_loop, "Phi_D": raw_phi_d, "Psi_D": raw_psi_d,
                     "D1_Dist": d1_val, "D2_Dist": d2_val, "SB_Dist": sb_dist_val, 
                     "HRD_ATP_Dist": hrd_atp_dist_val, "DFG_Mg_Dist": dfg_mg_dist_val,
                     "DFG_ATP_Dist": dfg_atp_dist_val, "PLoop_ATP_Dist": ploop_atp_dist_val,
                     "V104_RS2_Dist": v104_d, "I150_HRD_Dist": i150_d, "Shell_M118_M120_Dist": shell_d,
                     "Shell_State": shell_state,
                     "Y156_N99_Dist": y156_n99_d, "K105_E107_Dist": k105_e107_d, 
                     "K105_E121_Dist": k105_e121_d, "K105_N99_Dist": k105_n99_d, "D220_HRD_Dist": d220_hrd_d,
                     "atp_spec": atp_spec, "closest_d_o_name": closest_d_o_name, "dfg_d_res": dfg_d_res,
                     "closest_dfg_mg_atom": closest_dfg_mg_atom, "closest_mg_spec": closest_mg_spec,
                     "ploop_res_nums": ploop_res_nums, "sb_oe_names": sb_oe_names}
        }

    # =========================================================================
    # --- SPLIT VISUALIZATION GENERATORS ---
    # =========================================================================
    
    protein_name = " / ".join(sorted(list(set([name for name, _ in candidate_lms]))))
    rel_to_out_core = os.path.relpath(full_cif_path, out_dir_core).replace('\\', '/')
    rel_to_out_allo = os.path.relpath(full_cif_path, out_dir_allo).replace('\\', '/')

    # 1. SETUP CORE MACRO
    cxc_core = [
        f"# KinCore Visualization Macro for {sim_id}", f"open {rel_to_out_core}",
        "graphics silhouettes true", "graphics silhouettes width 1.5", "lighting soft", "color white", 
        "hide atoms", "show cartoons", "transparency 50 cartoons", "style stick", "2dlabels delete all",
        "\n# --- On-Screen Titles & Legend ---",
        f"2dlabels create title_prot text '{protein_name} (Catalytic Core)' color white size 24 bold true xpos 0.05 ypos 0.92",
        "2dlabels create leg_title text 'Color Legend:' color white size 20 bold true xpos 0.80 ypos 0.90",
        "2dlabels create leg_aloop text 'A-Loop (DFG/Dynamic)' color coral size 16 xpos 0.80 ypos 0.86",
        "2dlabels create leg_ploop text 'P-Loop (Gly-Rich)' color deep sky blue size 16 xpos 0.80 ypos 0.83",
        "2dlabels create leg_rspine text 'R-Spine' color medium purple size 16 xpos 0.80 ypos 0.80",
        "2dlabels create leg_sb text 'Salt Bridge (K-C)' color spring green size 16 xpos 0.80 ypos 0.77",
        "2dlabels create leg_cat text 'Catalytic HRD-Asp' color red size 16 xpos 0.80 ypos 0.74",
        "2dlabels create leg_dfg text 'DFG-Asp Coordination' color dodger blue size 16 xpos 0.80 ypos 0.71",
        "2dlabels create leg_shell text 'Hydrophobic Shell' color teal size 16 xpos 0.80 ypos 0.68",
        "2dlabels create leg_bridge text 'Core Bridges (V104/I150)' color tan size 16 xpos 0.80 ypos 0.65",
        "2dlabels create leg_lig text 'ATP / Magnesium' color gold size 16 xpos 0.80 ypos 0.62"
    ]

    # 2. SETUP ALLOSTERIC MACRO
    cxc_allo = [
        f"# Allosteric Communication Macro for {sim_id}", f"open {rel_to_out_allo}",
        "graphics silhouettes true", "graphics silhouettes width 1.5", "lighting soft", "color white", 
        "hide atoms", "show cartoons", "transparency 50 cartoons", "style stick", "2dlabels delete all",
        "\n# --- On-Screen Titles & Legend ---",
        f"2dlabels create title_prot text '{protein_name} (aC-b4 Allosteric Network)' color white size 24 bold true xpos 0.05 ypos 0.92",
        "2dlabels create leg_title text 'Color Legend:' color white size 20 bold true xpos 0.80 ypos 0.90",
        "2dlabels create leg_k105 text 'aC-b4 Sensor (K105)' color hot pink size 16 xpos 0.80 ypos 0.86",
        "2dlabels create leg_toggles text 'Toggle Targets (N99/E107/E121)' color pink size 16 xpos 0.80 ypos 0.83",
        "2dlabels create leg_y156 text 'aE Anchor (Y156)' color lime green size 16 xpos 0.80 ypos 0.80",
        "2dlabels create leg_d220 text 'aF Scaffold (D220)' color cyan size 16 xpos 0.80 ypos 0.77",
        "2dlabels create leg_cat text 'Catalytic HRD-Asp' color red size 16 xpos 0.80 ypos 0.74",
        "2dlabels create leg_lig text 'ATP / Magnesium' color gold size 16 xpos 0.80 ypos 0.71"
    ]

    if len(ligands) > 0:
        lig_cmds = [
            "\n# --- Ligands & Interactions ---", "show :ATP,ADP,ANP,ACP,AGS,AMP,GTP,GDP",
            "color :ATP,ADP,ANP,ACP,AGS,AMP,GTP,GDP byhetero", "color :ATP,ADP,ANP,ACP,AGS,AMP,GTP,GDP@C* gold",
            "label :ATP,ADP,ANP,ACP,AGS,AMP,GTP,GDP residues color gold height 1.5",
            "hbonds :ATP,ADP,ANP,ACP,AGS,AMP,GTP,GDP restrict protein color cyan radius 0.05",
            "show :MG", "color :MG green", "style :MG sphere", "size :MG atomRadius 1.0"
        ]
        cxc_core.extend(lig_cmds); cxc_allo.extend(lig_cmds)

    y_offset = 0.88

    for cid, data in chain_data.items():
        meta = data['meta']
        csv_rows.append({
            "Simulation_ID": sim_id, "Directory": os.path.dirname(full_cif_path), "File": os.path.basename(full_cif_path), 
            "Chain": cid, "Type": meta['Type'], "State": meta['State'],
            "R_Spine": meta['RSpine'], "C_Spine": meta['CSpine'], "C_Helix": meta['CHelix'],
            "Shell_State": meta['Shell_State'], "Spatial": meta['Spatial'], "Dihedral": meta['Dihedral'],
            "ActLoop_NT": meta.get('ActLoop_NT', 'N/A'), "ActLoop_CT": meta.get('ActLoop_CT', 'N/A'), 
            "Phi_D": meta.get('Phi_D', 'N/A'), "Psi_D": meta.get('Psi_D', 'N/A'),
            "D1_Dist": meta['D1_Dist'], "D2_Dist": meta['D2_Dist'], "SB_Dist": meta['SB_Dist'],
            "HRD_ATP_Dist": meta['HRD_ATP_Dist'], "DFG_Mg_Dist": meta.get('DFG_Mg_Dist', 'N/A'), 
            "DFG_ATP_Dist": meta.get('DFG_ATP_Dist', 'N/A'), "PLoop_ATP_Dist": meta.get('PLoop_ATP_Dist', 'N/A'),
            "V104_RS2_Dist": meta['V104_RS2_Dist'], "I150_HRD_Dist": meta['I150_HRD_Dist'], "Shell_M118_M120_Dist": meta['Shell_M118_M120_Dist'],
            "Y156_N99_Dist": meta['Y156_N99_Dist'], "K105_E107_Dist": meta['K105_E107_Dist'], 
            "K105_E121_Dist": meta['K105_E121_Dist'], "K105_N99_Dist": meta['K105_N99_Dist'], "D220_HRD_Dist": meta['D220_HRD_Dist']
        })

        if data['landmarks']:
            lm = data['landmarks']
            res = data['residues']
            f_num = res[lm['f']].number
            d_num = f_num - 1
            hrd_d_num = res[lm['hrd'] + 2].number if (lm.get('hrd') is not None and lm['hrd'] + 2 < len(res)) else None
            spec = f"/{cid}"
            
            state_label = f"2dlabels create State_{cid} text 'Chain {cid} ({meta['Type']}): {meta['State']}' color white size 20 xpos 0.05 ypos {y_offset}"
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

            if all(lm.get(x) is not None for x in ['v104', 'm118', 'm120', 'rs2']):
                v104, m118, m120, rs2 = res[lm['v104']].number, res[lm['m118']].number, res[lm['m120']].number, res[lm['rs2']].number
                cxc_core.extend([f"color {spec}:{v104},{m118},{m120} teal", f"show {spec}:{v104},{m118},{m120} & sidechain",
                                 f"distance {spec}:{m118}@CA {spec}:{m120}@CA color teal radius 0.05", f"distance {spec}:{v104}@CA {spec}:{rs2}@CA color teal radius 0.05"])
                
            if lm.get('i150') is not None and lm.get('hrd') is not None:
                cxc_core.extend([f"color {spec}:{res[lm['i150']].number} tan", f"show {spec}:{res[lm['i150']].number} & sidechain",
                                 f"distance {spec}:{res[lm['i150']].number}@CA {spec}:{res[lm['hrd']].number}@CA color tan radius 0.05"])

            # --- POPULATE ALLOSTERIC VISUALIZER ---
            cxc_allo.extend([f"color {spec} #d3d3d3"])
            if all(lm.get(x) is not None for x in ['k105', 'n99', 'e107', 'e121']):
                k105, n99, e107, e121 = res[lm['k105']].number, res[lm['n99']].number, res[lm['e107']].number, res[lm['e121']].number
                cxc_allo.extend([
                    f"\n# aC-b4 Toggle Switch", f"color {spec}:{k105} hot pink", f"show {spec}:{k105} & sidechain",
                    f"color {spec}:{n99},{e107},{e121} pink", f"show {spec}:{n99},{e107},{e121} & sidechain",
                    f"distance {spec}:{k105}@CA {spec}:{n99}@CA color hot pink radius 0.05",
                    f"distance {spec}:{k105}@CA {spec}:{e107}@CA color hot pink radius 0.05",
                    f"distance {spec}:{k105}@CA {spec}:{e121}@CA color hot pink radius 0.05"
                ])
                
            if lm.get('y156') is not None and lm.get('n99') is not None:
                y156, n99 = res[lm['y156']].number, res[lm['n99']].number
                cxc_allo.extend([
                    f"\n# aE Anchor", f"color {spec}:{y156} lime green", f"show {spec}:{y156} & sidechain",
                    f"distance {spec}:{y156}@CA {spec}:{n99}@CA color lime green radius 0.05"
                ])
                
            if lm.get('d220') is not None and lm.get('hrd') is not None:
                d220, hrd = res[lm['d220']].number, res[lm['hrd']].number
                cxc_allo.extend([
                    f"\n# aF Scaffold", f"color {spec}:{d220} cyan", f"show {spec}:{d220} & sidechain",
                    f"distance {spec}:{d220}@CA {spec}:{hrd}@CA color cyan radius 0.05"
                ])

            # --- SHARED ELEMENTS (Catalytic & Ligand Coordination) ---
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

    run(session, "close all", log=False)
    return csv_rows

def main(session):
    base_dir = os.path.abspath(SEARCH_DIR)
    out_dir_core = os.path.join(base_dir, VIZ_OUT_DIR_CORE)
    out_dir_allo = os.path.join(base_dir, VIZ_OUT_DIR_ALLO)
    os.makedirs(out_dir_core, exist_ok=True)
    os.makedirs(out_dir_allo, exist_ok=True)
    
    chunk_list_file = os.environ.get("CHIMERAX_CHUNK")
    if not chunk_list_file or not os.path.exists(chunk_list_file): run(session, "quit"); return

    out_csv_name = f"{os.path.splitext(os.path.basename(chunk_list_file))[0]}_results.csv"
    landmarks_file = os.path.join(base_dir, LANDMARKS_JSON)
    if not os.path.exists(landmarks_file): run(session, "quit"); return
        
    with open(landmarks_file, 'r') as f: all_landmarks = json.load(f)
    with open(chunk_list_file, 'r') as f: files_to_process =[line.strip() for line in f if line.strip()]

    fasta_path = os.path.join(base_dir, "sequences.fasta")
    fasta_seqs = {}
    if os.path.exists(fasta_path):
        with open(fasta_path, 'r') as f:
            name = None
            seq_lines =[]
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if name: fasta_seqs[name] = "".join(seq_lines)
                    name = line[1:].strip().split()[0]
                    seq_lines =[]
                elif line: seq_lines.append(line)
            if name: fasta_seqs[name] = "".join(seq_lines)

    cols =["Simulation_ID", "Directory", "File", "Chain", "Type", "State", 
            "R_Spine", "C_Spine", "C_Helix", "Shell_State", "Spatial", "Dihedral", "ActLoop_NT", "ActLoop_CT", 
            "Phi_D", "Psi_D", "D1_Dist", "D2_Dist", "SB_Dist", 
            "HRD_ATP_Dist", "DFG_Mg_Dist", "DFG_ATP_Dist", "PLoop_ATP_Dist",
            "V104_RS2_Dist", "I150_HRD_Dist", "Shell_M118_M120_Dist",
            "Y156_N99_Dist", "K105_E107_Dist", "K105_E121_Dist", "K105_N99_Dist", "D220_HRD_Dist"]
           
    wt_dir = os.environ.get("CHIMERAX_WT_DIR")
    all_rows = list() 
    
    for full_cif_path in files_to_process:
        try:
            rows = process_model(session, full_cif_path, base_dir, out_dir_core, out_dir_allo, wt_dir, all_landmarks, fasta_seqs)
            all_rows.extend(rows)
        except Exception as e: print(f"Error processing {full_cif_path}: {e}")

    with open(os.path.join(base_dir, out_csv_name), 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=cols)
        writer.writeheader()
        writer.writerows(all_rows)

    run(session, "quit")

main(session)
