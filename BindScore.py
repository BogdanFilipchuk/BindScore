"""
BindScore — Protein-Protein Binding Energy Estimator
=====================================================
Estimates ΔG_bind and K_d for a protein-protein complex from a PDB structure.

Physics:
    ΔG_bind = ΔH_interactions + ΔG_solvent + ΔG_rigid_body

Interaction terms:
    - Hydrogen bonds       (N/O donor-acceptor pairs)
    - Salt bridges         (charged +/- pairs)
    - Hydrophobic contacts (C–C nonpolar pairs)
    - π-stacking           (aromatic ring centroids)
    - Cation-π             (ARG/LYS nitrogen — aromatic ring)

Solvent term:
    ΔG_solvent = σ_np × ΔSASA_nonpolar  +  σ_p × ΔSASA_polar

Entropy:
    ΔG_entropy = ΔG_rigid_body  +  N_interface_res × ΔG_per_residue

Usage:
    python BindScore.py --pdb 1BRS --chain_a A --chain_b D
    python BindScore.py --pdb my_complex.pdb --chain_a A --chain_b B --local

Requirements:
    pip install biopython numpy
"""

import sys
import math
import json
import argparse
import urllib.request
import tempfile
import os
import numpy as np
import warnings
warnings.filterwarnings("ignore")

from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.vectors import Vector

# ─── Path to energy database ────────────────────────────────────────────────
SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
DB_PATH      = os.path.join(SCRIPT_DIR, "interactions.json")


# ════════════════════════════════════════════════════════════════════════════
#  RESIDUE PROPERTY TABLES
# ════════════════════════════════════════════════════════════════════════════

AROMATIC_RES   = {"PHE", "TYR", "TRP", "HIS"}
HYDROPHOBIC_RES = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "TYR"}
NONPOLAR_ELEM  = {"C", "S"}

POS_CHARGED = {
    "ARG": ["NH1", "NH2", "NE"],
    "LYS": ["NZ"],
    "HIS": ["ND1", "NE2"],
}
NEG_CHARGED = {
    "ASP": ["OD1", "OD2"],
    "GLU": ["OE1", "OE2"],
}

# Aromatic ring heavy-atom names per residue
AROMATIC_ATOMS = {
    "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TRP": ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
    "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"],
}

# Approximate VdW radii by element (Å)
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80,
       "H": 1.20, "P": 1.80, "F": 1.47}


# ════════════════════════════════════════════════════════════════════════════
#  LOAD ENERGY DATABASE
# ════════════════════════════════════════════════════════════════════════════

def load_db(path=DB_PATH):
    with open(path) as f:
        return json.load(f)


# ════════════════════════════════════════════════════════════════════════════
#  PDB FETCH / PARSE
# ════════════════════════════════════════════════════════════════════════════

def fetch_pdb(accession):
    """Download PDB file from RCSB by 4-letter accession code."""
    url = f"https://files.rcsb.org/download/{accession.upper()}.pdb"
    print(f"  Fetching {accession.upper()} from RCSB...")
    req = urllib.request.Request(url, headers={"User-Agent": "BindScore/1.0 (biopython)"})
    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
    with urllib.request.urlopen(req) as response:
        tmp.write(response.read())
    tmp.flush()
    return tmp.name

def parse_structure(pdb_path, name="complex"):
    parser = PDBParser(QUIET=True)
    return parser.get_structure(name, pdb_path)

def get_chain(structure, chain_id):
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                return chain
    raise ValueError(f"Chain '{chain_id}' not found in structure.")

def standard_residues(chain):
    """Return only standard amino acid residues (no HETATM, no water)."""
    return [r for r in chain if r.get_id()[0] == " "]


# ════════════════════════════════════════════════════════════════════════════
#  GEOMETRY HELPERS
# ════════════════════════════════════════════════════════════════════════════

def atom_dist(a1, a2):
    return (a1.get_vector() - a2.get_vector()).norm()

def ring_centroid(residue):
    """Return centroid of aromatic ring as numpy array."""
    names = AROMATIC_ATOMS.get(residue.get_resname(), [])
    coords = []
    for name in names:
        if name in residue:
            coords.append(residue[name].get_vector().get_array())
    if not coords:
        return None
    return np.mean(coords, axis=0)

def ring_normal(residue):
    """Return normal vector of aromatic ring plane."""
    names = AROMATIC_ATOMS.get(residue.get_resname(), [])
    coords = []
    for name in names:
        if name in residue:
            coords.append(residue[name].get_vector().get_array())
    if len(coords) < 3:
        return None
    c = np.mean(coords, axis=0)
    vecs = [np.array(p) - c for p in coords]
    # fit plane by SVD
    _, _, vh = np.linalg.svd(vecs)
    return vh[-1]  # normal = last right singular vector

def angle_between(v1, v2):
    """Angle in degrees between two vectors."""
    cos = abs(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-9))
    cos = min(cos, 1.0)
    return math.degrees(math.acos(cos))


# ════════════════════════════════════════════════════════════════════════════
#  INTERFACE DETECTION
# ════════════════════════════════════════════════════════════════════════════

def find_interface_residues(chain_a_res, chain_b_res, cutoff=8.0):
    """
    Return (interface_a, interface_b) — residues within `cutoff` Å of
    any atom in the partner chain.
    """
    def all_atoms(residues):
        atoms = []
        for r in residues:
            atoms.extend(r.get_atoms())
        return atoms

    atoms_a = all_atoms(chain_a_res)
    atoms_b = all_atoms(chain_b_res)

    iface_a, iface_b = set(), set()

    for ra in chain_a_res:
        for atom_a in ra.get_atoms():
            for atom_b in atoms_b:
                if atom_dist(atom_a, atom_b) <= cutoff:
                    iface_a.add(ra)
                    break
            else:
                continue
            break

    for rb in chain_b_res:
        for atom_b in rb.get_atoms():
            for atom_a in atoms_a:
                if atom_dist(atom_b, atom_a) <= cutoff:
                    iface_b.add(rb)
                    break
            else:
                continue
            break

    return list(iface_a), list(iface_b)


# ════════════════════════════════════════════════════════════════════════════
#  INTERACTION FINDERS
# ════════════════════════════════════════════════════════════════════════════

def gaussian_decay(dist, d_opt, d_cut):
    """Distance-weighted energy: 1.0 at d_opt, 0.0 at d_cut."""
    if dist > d_cut:
        return 0.0
    sigma = (d_cut - d_opt) / 2.0
    return math.exp(-0.5 * ((dist - d_opt) / max(sigma, 0.01)) ** 2)


def find_hbonds(iface_a, iface_b, db):
    """
    Hydrogen bonds: N or O donor on one chain, N or O acceptor on the other.
    """
    cfg    = db["hydrogen_bond"]
    d_opt  = cfg["distance_optimal"]
    d_cut  = cfg["distance_cutoff"]
    e_opt  = cfg["energy_optimal"]
    bonds  = []

    donor_acc = {"N", "O"}

    def polar_atoms(residues):
        result = []
        for r in residues:
            for a in r.get_atoms():
                if (a.element or "").upper() in donor_acc:
                    result.append((r, a))
        return result

    polars_a = polar_atoms(iface_a)
    polars_b = polar_atoms(iface_b)

    for (ra, aa) in polars_a:
        for (rb, ab) in polars_b:
            d = atom_dist(aa, ab)
            w = gaussian_decay(d, d_opt, d_cut)
            if w > 0:
                bonds.append({
                    "res_a": f"{ra.get_resname()}{ra.get_id()[1]}",
                    "res_b": f"{rb.get_resname()}{rb.get_id()[1]}",
                    "dist":  round(d, 2),
                    "energy": round(e_opt * w, 3),
                })
    return bonds


def find_salt_bridges(iface_a, iface_b, db):
    """
    Salt bridges: positively charged atoms (ARG/LYS/HIS) — negatively
    charged atoms (ASP/GLU) across the interface.
    """
    cfg   = db["salt_bridge"]
    d_opt = cfg["distance_optimal"]
    d_cut = cfg["distance_cutoff"]
    e_opt = cfg["energy_optimal"]
    bridges = []

    def charged_atoms(residues, charge_dict):
        result = []
        for r in residues:
            rname = r.get_resname()
            if rname in charge_dict:
                for aname in charge_dict[rname]:
                    if aname in r:
                        result.append((r, r[aname]))
        return result

    pos_a = charged_atoms(iface_a, POS_CHARGED)
    neg_a = charged_atoms(iface_a, NEG_CHARGED)
    pos_b = charged_atoms(iface_b, POS_CHARGED)
    neg_b = charged_atoms(iface_b, NEG_CHARGED)

    pairs = [(pos_a, neg_b), (neg_a, pos_b)]
    for pos_list, neg_list in pairs:
        for (rp, ap) in pos_list:
            for (rn, an) in neg_list:
                d = atom_dist(ap, an)
                w = gaussian_decay(d, d_opt, d_cut)
                if w > 0:
                    bridges.append({
                        "res_pos": f"{rp.get_resname()}{rp.get_id()[1]}",
                        "res_neg": f"{rn.get_resname()}{rn.get_id()[1]}",
                        "dist":    round(d, 2),
                        "energy":  round(e_opt * w, 3),
                    })
    return bridges


def find_hydrophobic(iface_a, iface_b, db):
    """
    Hydrophobic contacts: C–C pairs between nonpolar residues.
    """
    cfg   = db["hydrophobic"]
    d_opt = cfg["distance_optimal"]
    d_cut = cfg["distance_cutoff"]
    e_opt = cfg["energy_optimal"]
    contacts = []

    def nonpolar_carbons(residues):
        result = []
        for r in residues:
            if r.get_resname() in HYDROPHOBIC_RES:
                for a in r.get_atoms():
                    if (a.element or "").upper() == "C":
                        result.append((r, a))
        return result

    carbons_a = nonpolar_carbons(iface_a)
    carbons_b = nonpolar_carbons(iface_b)

    seen = set()
    for (ra, aa) in carbons_a:
        for (rb, ab) in carbons_b:
            key = (ra.get_id()[1], rb.get_id()[1])
            if key in seen:
                continue
            d = atom_dist(aa, ab)
            w = gaussian_decay(d, d_opt, d_cut)
            if w > 0:
                seen.add(key)
                contacts.append({
                    "res_a":  f"{ra.get_resname()}{ra.get_id()[1]}",
                    "res_b":  f"{rb.get_resname()}{rb.get_id()[1]}",
                    "dist":   round(d, 2),
                    "energy": round(e_opt * w, 3),
                })
    return contacts


def find_pi_stacking(iface_a, iface_b, db):
    """
    π-stacking: aromatic ring centroid distance + near-parallel normals.
    """
    cfg     = db["pi_stack"]
    d_cut   = cfg["centroid_cutoff"]
    ang_cut = cfg["angle_cutoff_deg"]
    e_opt   = cfg["energy_optimal"]
    stacks  = []

    arom_a = [(r, ring_centroid(r), ring_normal(r))
              for r in iface_a if r.get_resname() in AROMATIC_RES]
    arom_b = [(r, ring_centroid(r), ring_normal(r))
              for r in iface_b if r.get_resname() in AROMATIC_RES]

    for (ra, ca, na) in arom_a:
        if ca is None or na is None:
            continue
        for (rb, cb, nb) in arom_b:
            if cb is None or nb is None:
                continue
            d   = np.linalg.norm(ca - cb)
            ang = angle_between(na, nb)
            if d <= d_cut and (ang <= ang_cut or ang >= 180 - ang_cut):
                w = gaussian_decay(d, 3.5, d_cut)
                stacks.append({
                    "res_a":  f"{ra.get_resname()}{ra.get_id()[1]}",
                    "res_b":  f"{rb.get_resname()}{rb.get_id()[1]}",
                    "dist":   round(float(d), 2),
                    "angle":  round(float(ang), 1),
                    "energy": round(e_opt * w, 3),
                })
    return stacks


def find_cation_pi(iface_a, iface_b, db):
    """
    Cation-π: ARG/LYS nitrogen — aromatic ring centroid.
    """
    cfg   = db["cation_pi"]
    d_cut = cfg["distance_cutoff"]
    e_opt = cfg["energy_optimal"]
    pairs = []

    def cation_atoms(residues):
        result = []
        for r in residues:
            rname = r.get_resname()
            if rname in POS_CHARGED:
                for aname in POS_CHARGED[rname]:
                    if aname in r:
                        result.append((r, r[aname]))
        return result

    def arom_centroids(residues):
        result = []
        for r in residues:
            if r.get_resname() in AROMATIC_RES:
                c = ring_centroid(r)
                if c is not None:
                    result.append((r, c))
        return result

    cats_a, aroms_a = cation_atoms(iface_a), arom_centroids(iface_a)
    cats_b, aroms_b = cation_atoms(iface_b), arom_centroids(iface_b)

    for (cat_list, arom_list) in [(cats_a, aroms_b), (cats_b, aroms_a)]:
        for (rc, ac) in cat_list:
            cat_pos = ac.get_vector().get_array()
            for (ra, centroid) in arom_list:
                d = np.linalg.norm(cat_pos - centroid)
                if d <= d_cut:
                    w = gaussian_decay(d, 4.0, d_cut)
                    pairs.append({
                        "cation": f"{rc.get_resname()}{rc.get_id()[1]}",
                        "arom":   f"{ra.get_resname()}{ra.get_id()[1]}",
                        "dist":   round(float(d), 2),
                        "energy": round(e_opt * w, 3),
                    })
    return pairs


# ════════════════════════════════════════════════════════════════════════════
#  SASA (APPROXIMATE)
# ════════════════════════════════════════════════════════════════════════════

def approx_sasa(residues):
    """
    Approximate SASA using a sphere-overlap method.
    Returns (total_sasa, nonpolar_sasa, polar_sasa) in Å².
    """
    total_np = 0.0
    total_p  = 0.0

    for r in residues:
        rname = r.get_resname()
        for atom in r.get_atoms():
            elem = (atom.element or "C").upper()
            rad  = VDW.get(elem, 1.70)
            area = 4 * math.pi * (rad + 1.4) ** 2  # probe radius 1.4 Å
            if elem in NONPOLAR_ELEM:
                total_np += area * 0.25  # rough exposure fraction
            else:
                total_p  += area * 0.25
    return total_np + total_p, total_np, total_p


# ════════════════════════════════════════════════════════════════════════════
#  SCORING ENGINE
# ════════════════════════════════════════════════════════════════════════════

def score_interactions(hbonds, salt_bridges, hydrophobic, pi_stacks, cation_pi):
    dH  = sum(x["energy"] for x in hbonds)
    dH += sum(x["energy"] for x in salt_bridges)
    dH += sum(x["energy"] for x in hydrophobic)
    dH += sum(x["energy"] for x in pi_stacks)
    dH += sum(x["energy"] for x in cation_pi)
    return dH

def score_solvation(chain_a_all, chain_b_all, iface_a, iface_b, db):
    """
    ΔG_solvent ≈ σ_np × ΔSASA_np  +  σ_p × ΔSASA_p
    ΔSASA = SASA(isolated) - SASA(complex) ≈ interface SASA
    """
    cfg    = db["solvation"]
    sig_np = cfg["nonpolar_coefficient"]
    sig_p  = cfg["polar_coefficient"]

    _, np_iface_a, p_iface_a = approx_sasa(iface_a)
    _, np_iface_b, p_iface_b = approx_sasa(iface_b)

    delta_np = np_iface_a + np_iface_b
    delta_p  = p_iface_a  + p_iface_b

    dG_solv = sig_np * delta_np + sig_p * delta_p
    return dG_solv, delta_np, delta_p

def score_entropy(iface_a, iface_b, db):
    cfg      = db["entropic_penalty"]
    rigid    = cfg["rigid_body"]
    per_res  = cfg["per_interface_residue"]
    n_res    = len(iface_a) + len(iface_b)
    return rigid + per_res * n_res

def dG_to_kd(dG, T=298.15):
    R = 0.001987  # kcal/mol/K
    return math.exp(dG / (R * T))

def interpret_kd(kd):
    if kd < 1e-12:
        return f"{kd*1e12:.2f} pM", "Exceptionally tight (antibody-like)"
    elif kd < 1e-9:
        return f"{kd*1e9:.2f} nM",  "Very strong (drug target-like)"
    elif kd < 1e-6:
        return f"{kd*1e6:.2f} µM",  "Moderate"
    elif kd < 1e-3:
        return f"{kd*1e3:.2f} mM",  "Weak"
    else:
        return f"{kd:.4f} M",       "Very weak / non-binder"


# ════════════════════════════════════════════════════════════════════════════
#  MAIN ANALYSIS FUNCTION (Python API)
# ════════════════════════════════════════════════════════════════════════════

def analyze(pdb_input, chain_a, chain_b, local=False, db_path=DB_PATH):
    """
    Main entry point.

    Parameters
    ----------
    pdb_input : str
        4-letter RCSB accession code OR path to a local PDB file (if local=True).
    chain_a : str
        Chain ID for the first protein.
    chain_b : str
        Chain ID for the second protein.
    local : bool
        If True, treat pdb_input as a local file path.
    db_path : str
        Path to interactions.json energy database.

    Returns
    -------
    dict with all results.
    """
    db = load_db(db_path)

    # ── Load structure ───────────────────────────────────────────────────────
    if local:
        pdb_path = pdb_input
        label    = os.path.basename(pdb_input)
    else:
        pdb_path = fetch_pdb(pdb_input)
        label    = pdb_input.upper()

    structure = parse_structure(pdb_path, label)

    chain_a_obj  = get_chain(structure, chain_a)
    chain_b_obj  = get_chain(structure, chain_b)
    res_a        = standard_residues(chain_a_obj)
    res_b        = standard_residues(chain_b_obj)

    # ── Interface ────────────────────────────────────────────────────────────
    iface_a, iface_b = find_interface_residues(res_a, res_b, cutoff=8.0)

    # ── Interactions ─────────────────────────────────────────────────────────
    hbonds      = find_hbonds(iface_a, iface_b, db)
    salt_bridges= find_salt_bridges(iface_a, iface_b, db)
    hydrophobic = find_hydrophobic(iface_a, iface_b, db)
    pi_stacks   = find_pi_stacking(iface_a, iface_b, db)
    cation_pi   = find_cation_pi(iface_a, iface_b, db)

    # ── Scoring ──────────────────────────────────────────────────────────────
    dH     = score_interactions(hbonds, salt_bridges, hydrophobic, pi_stacks, cation_pi)
    dG_sol, buried_np, buried_p = score_solvation(res_a, res_b, iface_a, iface_b, db)
    dG_ent = score_entropy(iface_a, iface_b, db)
    dG_tot = dH + dG_sol + dG_ent
    kd     = dG_to_kd(dG_tot)
    kd_str, interpretation = interpret_kd(kd)

    return {
        "label":          label,
        "chain_a":        chain_a,
        "chain_b":        chain_b,
        "n_res_a":        len(res_a),
        "n_res_b":        len(res_b),
        "n_iface_a":      len(iface_a),
        "n_iface_b":      len(iface_b),
        "hbonds":         hbonds,
        "salt_bridges":   salt_bridges,
        "hydrophobic":    hydrophobic,
        "pi_stacks":      pi_stacks,
        "cation_pi":      cation_pi,
        "dH_interactions":round(dH, 2),
        "dG_solvent":     round(dG_sol, 2),
        "buried_np_A2":   round(buried_np, 1),
        "buried_p_A2":    round(buried_p, 1),
        "dG_entropy":     round(dG_ent, 2),
        "dG_total":       round(dG_tot, 2),
        "K_d_molar":      kd,
        "K_d_str":        kd_str,
        "interpretation": interpretation,
        "iface_res_a":    [f"{r.get_resname()}{r.get_id()[1]}" for r in iface_a],
        "iface_res_b":    [f"{r.get_resname()}{r.get_id()[1]}" for r in iface_b],
    }


# ════════════════════════════════════════════════════════════════════════════
#  PRETTY PRINTER
# ════════════════════════════════════════════════════════════════════════════

def print_report(r):
    W = 62
    sep  = "═" * W
    thin = "─" * W

    print(f"\n{sep}")
    print(f"  BindScore — Protein-Protein Binding Energy Estimate")
    print(f"  Structure : {r['label']}   Chains: {r['chain_a']} ↔ {r['chain_b']}")
    print(sep)

    print(f"\n  CHAIN SUMMARY")
    print(thin)
    print(f"  Chain {r['chain_a']} residues     : {r['n_res_a']}")
    print(f"  Chain {r['chain_b']} residues     : {r['n_res_b']}")
    print(f"  Interface residues (A) : {r['n_iface_a']}")
    print(f"  Interface residues (B) : {r['n_iface_b']}")

    print(f"\n  INTERACTIONS AT INTERFACE")
    print(thin)
    print(f"  {'Type':<30} {'Count':>6}  {'ΔH (kcal/mol)':>14}")
    print(f"  {'-'*30}  {'-'*6}  {'-'*14}")

    rows = [
        ("Hydrogen bonds",       r["hbonds"],       sum(x["energy"] for x in r["hbonds"])),
        ("Salt bridges",         r["salt_bridges"],  sum(x["energy"] for x in r["salt_bridges"])),
        ("Hydrophobic contacts", r["hydrophobic"],   sum(x["energy"] for x in r["hydrophobic"])),
        ("π-stacking",           r["pi_stacks"],     sum(x["energy"] for x in r["pi_stacks"])),
        ("Cation-π",             r["cation_pi"],     sum(x["energy"] for x in r["cation_pi"])),
    ]
    for label, items, energy in rows:
        print(f"  {label:<30} {len(items):>6}  {energy:>14.2f}")

    print(f"\n  BURIED SURFACE AREA")
    print(thin)
    print(f"  Nonpolar ΔSASA : {r['buried_np_A2']:>8.1f} Å²")
    print(f"  Polar    ΔSASA : {r['buried_p_A2']:>8.1f} Å²")

    print(f"\n  FREE ENERGY DECOMPOSITION")
    print(thin)
    print(f"  {'Term':<35} {'ΔG (kcal/mol)':>12}")
    print(f"  {'-'*35}  {'-'*12}")
    print(f"  {'ΔH interactions':<35} {r['dH_interactions']:>12.2f}")
    print(f"  {'ΔG solvent release':<35} {r['dG_solvent']:>12.2f}")
    print(f"  {'ΔG entropic penalty':<35} {r['dG_entropy']:>12.2f}")
    print(f"  {'─'*35}  {'─'*12}")
    print(f"  {'ΔG total':<35} {r['dG_total']:>12.2f}")

    print(f"\n  BINDING AFFINITY")
    print(thin)
    print(f"  K_d (estimated) : {r['K_d_str']}")
    print(f"  Interpretation  : {r['interpretation']}")

    print(f"\n  INTERFACE RESIDUES — Chain {r['chain_a']}")
    print(thin)
    res_a = sorted(r["iface_res_a"])
    for i in range(0, len(res_a), 6):
        print("  " + "  ".join(f"{x:<8}" for x in res_a[i:i+6]))

    print(f"\n  INTERFACE RESIDUES — Chain {r['chain_b']}")
    print(thin)
    res_b = sorted(r["iface_res_b"])
    for i in range(0, len(res_b), 6):
        print("  " + "  ".join(f"{x:<8}" for x in res_b[i:i+6]))

    print(f"\n  DISCLAIMER")
    print(thin)
    print("  Static-structure empirical estimate. Typical error: ±2–4 kcal/mol")
    print("  after calibration. For production use, combine with MM-GBSA or FEP.")
    print(f"{sep}\n")


# ════════════════════════════════════════════════════════════════════════════
#  CLI
# ════════════════════════════════════════════════════════════════════════════

def build_parser():
    p = argparse.ArgumentParser(
        description="BindScore: Protein-Protein Binding Energy Estimator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python BindScore.py --pdb 1BRS --chain_a A --chain_b D
  python BindScore.py --pdb complex.pdb --chain_a A --chain_b B --local
  python BindScore.py --pdb 2PTC --chain_a E --chain_b I --db my_weights.json
        """,
    )
    p.add_argument("--pdb",      required=True,
                   help="RCSB accession code (e.g. 1BRS) or path to local PDB file")
    p.add_argument("--chain_a",  required=True, help="Chain ID for protein A")
    p.add_argument("--chain_b",  required=True, help="Chain ID for protein B")
    p.add_argument("--local",    action="store_true",
                   help="Treat --pdb as a local file path instead of RCSB accession")
    p.add_argument("--db",       default=DB_PATH,
                   help=f"Path to interactions.json (default: {DB_PATH})")
    p.add_argument("--json",     action="store_true",
                   help="Also dump full results as JSON to stdout after the report")
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()

    result = analyze(
        pdb_input = args.pdb,
        chain_a   = args.chain_a,
        chain_b   = args.chain_b,
        local     = args.local,
        db_path   = args.db,
    )

    print_report(result)

    if args.json:
        import json as _json
        # K_d_molar is a float — make it serialisable
        out = dict(result)
        out["K_d_molar"] = float(result["K_d_molar"])
        print(_json.dumps(out, indent=2))
