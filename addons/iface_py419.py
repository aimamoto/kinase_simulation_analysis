#!/usr/bin/env python3
"""Does SRC pY419 (PTR159) change the CSK-SRC docking interface?

Direct, assumption-free measurement on the AF3 ensembles: contacts, buried
surface area, polar bridges, and per-residue-pair contact frequency, compared
between SRC-unphosphorylated and SRC-pY419 conditions.
"""
import os, sys, json, collections
import numpy as np
import gemmi
from scipy.spatial import cKDTree

ROOT = os.environ.get("ALLOQUANT_SRC", "CSK-SRC_output")  # AlloQuant CSK-SRC output dir
CONDS = [
    'a-csk-wtcat-holo_b-src-wtcat-holo',        # unprimed, CSK+ATP  (primary contrast)
    'a-csk-wtcat-holo_b-src-wtcat-py159-holo',  # primed,   CSK+ATP
    'a-csk-wtcat-apo_b-src-wtcat-holo',         # unprimed, CSK apo  (replication)
    'a-csk-wtcat-apo_b-src-wtcat-py159-holo',   # primed,   CSK apo
]
RADII = {'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80, 'P': 1.80}
PROBE = 1.4
N_SPHERE = 92


def sphere_points(n=N_SPHERE):
    """Golden-spiral points on the unit sphere."""
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(theta) * np.sin(phi),
                     np.sin(theta) * np.sin(phi),
                     np.cos(phi)], axis=1)


SPHERE = sphere_points()


def sasa(coords, radii, score_idx=None):
    """Shrake-Rupley SASA (A^2). Returns per-atom SASA for score_idx atoms."""
    r_ext = radii + PROBE
    tree = cKDTree(coords)
    if score_idx is None:
        score_idx = np.arange(len(coords))
    out = np.zeros(len(score_idx))
    max_r = r_ext.max()
    for k, i in enumerate(score_idx):
        nb = tree.query_ball_point(coords[i], r_ext[i] + max_r)
        nb = [j for j in nb if j != i]
        pts = coords[i] + SPHERE * r_ext[i]
        if nb:
            d = np.linalg.norm(pts[:, None, :] - coords[nb][None, :, :], axis=2)
            buried = (d < r_ext[nb][None, :]).any(axis=1)
            acc = (~buried).sum()
        else:
            acc = len(pts)
        out[k] = 4 * np.pi * r_ext[i] ** 2 * acc / len(pts)
    return out


def load(path):
    """Return (chain_label -> dict of polymer heavy-atom arrays) for the two kinases."""
    st = gemmi.read_structure(path)
    st.setup_entities()
    st.remove_hydrogens()
    model = st[0]
    chains = {}
    for ch in model:
        poly = ch.get_polymer()
        if len(poly) < 50:          # skip ligand / ion chains
            continue
        xyz, rad, resid, resname = [], [], [], []
        for res in poly:
            for at in res:
                el = at.element.name
                if el == 'H':
                    continue
                xyz.append([at.pos.x, at.pos.y, at.pos.z])
                rad.append(RADII.get(el, 1.70))
                resid.append(res.seqid.num)
                resname.append(res.name)
        chains[ch.name] = dict(xyz=np.array(xyz), rad=np.array(rad),
                               resid=np.array(resid), resname=np.array(resname),
                               nres=len(poly))
    return chains


CHARGED_POS = {'ARG': ('NE', 'NH1', 'NH2'), 'LYS': ('NZ',), 'HIS': ('ND1', 'NE2')}
CHARGED_NEG = {'ASP': ('OD1', 'OD2'), 'GLU': ('OE1', 'OE2'), 'PTR': ('O1P', 'O2P', 'O3P')}


def polar_atoms(path):
    """Charged-group atoms per chain, for salt-bridge counting."""
    st = gemmi.read_structure(path)
    st.setup_entities()
    model = st[0]
    out = {}
    for ch in model:
        poly = ch.get_polymer()
        if len(poly) < 50:
            continue
        pos, neg = [], []
        for res in poly:
            for at in res:
                if res.name in CHARGED_POS and at.name in CHARGED_POS[res.name]:
                    pos.append([at.pos.x, at.pos.y, at.pos.z])
                if res.name in CHARGED_NEG and at.name in CHARGED_NEG[res.name]:
                    neg.append([at.pos.x, at.pos.y, at.pos.z])
        out[ch.name] = (np.array(pos).reshape(-1, 3), np.array(neg).reshape(-1, 3))
    return out


PHOS = {'P', 'O1P', 'O2P', 'O3P', 'OP1', 'OP2', 'OP3'}


def ptr_coords(path):
    """Whole pY419 residue and its phosphate group separately — they are ~8 A apart
    relative to CSK, so the distinction matters for what the text can claim."""
    st = gemmi.read_structure(path)
    for ch in st[0]:
        for res in ch:
            if res.name == 'PTR':
                allc = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in res])
                pho = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in res if a.name in PHOS])
                return ch.name, allc, pho
    return None, None, None


def analyse(path):
    ch = load(path)
    names = sorted(ch)
    assert len(names) == 2, f'expected 2 kinase chains, got {names} in {path}'
    csk, src = names[0], names[1]      # A = CSK, later letter = SRC
    A, B = ch[csk], ch[src]

    d = np.linalg.norm(A['xyz'][:, None, :] - B['xyz'][None, :, :], axis=2)
    n4 = int((d < 4.0).sum())
    n5 = int((d < 5.0).sum())
    mind = float(d.min())
    ires_a = sorted(set(A['resid'][(d < 4.5).any(axis=1)].tolist()))
    ires_b = sorted(set(B['resid'][(d < 4.5).any(axis=0)].tolist()))

    # buried surface area (only interface-vicinity atoms can change)
    near_a = np.where((d < 10.0).any(axis=1))[0]
    near_b = np.where((d < 10.0).any(axis=0))[0]
    allxyz = np.vstack([A['xyz'], B['xyz']])
    allrad = np.concatenate([A['rad'], B['rad']])
    off = len(A['xyz'])
    s_complex = (sasa(allxyz, allrad, near_a).sum()
                 + sasa(allxyz, allrad, near_b + off).sum())
    s_free = (sasa(A['xyz'], A['rad'], near_a).sum()
              + sasa(B['xyz'], B['rad'], near_b).sum())
    bsa = float(s_free - s_complex)

    # inter-chain salt bridges / polar contacts
    pa = polar_atoms(path)
    posA, negA = pa[csk]
    posB, negB = pa[src]
    sb = 0
    for P, N in ((posA, negB), (posB, negA)):
        if len(P) and len(N):
            sb += int((np.linalg.norm(P[:, None, :] - N[None, :, :], axis=2) < 4.0).sum())

    # is the pY mark itself anywhere near the interface?
    pch, pxyz, phos = ptr_coords(path)
    if pxyz is not None:
        other = B['xyz'] if pch == csk else A['xyz']
        ptr_min = float(np.linalg.norm(pxyz[:, None, :] - other[None, :, :], axis=2).min())
        phos_min = float(np.linalg.norm(phos[:, None, :] - other[None, :, :], axis=2).min())
    else:
        ptr_min = phos_min = np.nan

    pairs = [(int(a), int(b)) for a, b in zip(A['resid'][np.where(d < 4.5)[0]],
                                              B['resid'][np.where(d < 4.5)[1]])]
    return dict(n_contacts_4=n4, n_contacts_5=n5, min_dist=mind,
                n_iface_res_csk=len(ires_a), n_iface_res_src=len(ires_b),
                bsa=bsa, n_saltbridge=sb, ptr_min_dist=ptr_min, phos_min_dist=phos_min,
                com_dist=float(np.linalg.norm(A['xyz'].mean(0) - B['xyz'].mean(0))),
                pairs=sorted(set(pairs)))


def main():
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else None
    res = {}
    for cond in CONDS:
        cdir = os.path.join(ROOT, cond)
        seeds = sorted(os.listdir(cdir))
        if limit:
            seeds = seeds[:limit]
        rows = []
        for i, s in enumerate(seeds):
            p = os.path.join(cdir, s, 'model.cif')
            try:
                r = analyse(p)
            except Exception as e:                      # noqa: BLE001
                print(f'  !! {cond}/{s}: {e}', flush=True)
                continue
            r['seed'] = s
            rows.append(r)
            if (i + 1) % 20 == 0:
                print(f'  {cond}: {i+1}/{len(seeds)}', flush=True)
        res[cond] = rows
        print(f'{cond}: {len(rows)} models done', flush=True)
    out = os.environ.get("ALLOQUANT_IFACE_JSON", "analysis_iface_py419.json")
    with open(out, 'w') as fh:
        json.dump(res, fh)
    print('wrote', out)


if __name__ == '__main__':
    main()
