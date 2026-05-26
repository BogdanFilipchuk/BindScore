import math
from bindscore.parsing.pdb_utils_inter import Interaction

# ── Physical constants ────────────────────────────────────────────────────────
COULOMB_CONSTANT = 1389.4  
EPSILON_RELATIVE = 50.0   
MIN_ELECTROSTATIC_DISTANCE = 2.8  

# ── Condensed Database Lookups ────────────────────────────────────────────────
METAL_WELL_DEPTHS = {
    "ZN": 125.5, "CU": 125.5, "NI": 117.2, "CO": 104.6, "FE": 104.6,
    "MN": 83.7,  "CD": 75.3,  "HG": 75.3,  "PT": 146.4, "AU": 146.4,
    "MG": 62.8,  "CA": 50.2,  "NA": 20.9,  "K": 16.7
}

SALT_BRIDGE_CHARGES = {
    ("ARG", "NH1"): 1.0,  ("ARG", "NH2"): 1.0,  ("ARG", "NE"): 1.0,  ("LYS", "NZ"): 1.0,
    ("ASP", "OD1"): -1.0, ("ASP", "OD2"): -1.0, ("GLU", "OE1"): -1.0, ("GLU", "OE2"): -1.0
}

# Flattened lookup to eliminate nested residue-to-atom conditional forks
PARTIAL_CHARGES = {
    "N": -0.4157, "H": 0.2719, "CA": 0.0337, "HA": 0.0823, "C": 0.5973, "O": -0.5679,
    "ARG_CB": -0.0007, "ARG_HB2": 0.0327,  "ARG_HB3": 0.0327,  "ARG_CG": 0.0390,
    "ARG_HG2": 0.0285,  "ARG_HG3": 0.0285,  "ARG_CD": 0.0486,  "ARG_HD2": 0.0687,
    "ARG_HD3": 0.0687,  "ARG_NE": -0.5295,  "ARG_HE": 0.3456,  "ARG_CZ": 0.8076,
    "ARG_NH1": -0.8627, "ARG_HH11": 0.4478, "ARG_HH12": 0.4478, "ARG_NH2": -0.8627,
    "ARG_HH21": 0.4478, "ARG_HH22": 0.4478,
    "LYS_CB": -0.0094, "LYS_HB2": 0.0362,  "LYS_HB3": 0.0362,  "LYS_CG": 0.0187,
    "LYS_HG2": 0.0103,  "LYS_HG3": 0.0103,  "LYS_CD": -0.0479, "LYS_HD2": 0.0621,
    "LYS_HD3": 0.0621,  "LYS_CE": -0.0143,  "LYS_HE2": 0.1135, "LYS_HE3": 0.1135,
    "LYS_NZ": -0.3854,  "LYS_HZ1": 0.3400,  "LYS_HZ2": 0.3400,  "LYS_HZ3": 0.3400,
    "ASP_CB": -0.0303, "ASP_HB2": 0.0496,  "ASP_HB3": 0.0496,  "ASP_CG": 0.7994,
    "ASP_OD1": -0.8014, "ASP_OD2": -0.8014,
    "GLU_CB": 0.0560,  "GLU_HB2": -0.0173, "GLU_HB3": -0.0173, "GLU_CG": 0.0136,
    "GLU_HG2": -0.0425, "GLU_HG3": -0.0425, "GLU_CD": 0.8054,  "GLU_OE1": -0.8188,
    "GLU_OE2": -0.8188
}

WATER_HBOND_ENTHALPY = 4.25  
WATER_HBOND_CAPACITY = {
    "ARG_NE": 1, "ARG_NH1": 2, "ARG_NH2": 2, "LYS_NZ": 3, "ASP_OD1": 3, "ASP_OD2": 3,
    "GLU_OE1": 3, "GLU_OE2": 3, "HIS_ND1": 1, "HIS_NE2": 1, "ASN_OD1": 2, "ASN_ND2": 2,
    "GLN_OE1": 2, "GLN_NE2": 2, "SER_OG": 2,  "THR_OG1": 2, "TYR_OH": 2,  "TRP_NE1": 1,
    "CYS_SG": 1,  "MET_SD": 1, "N": 1, "O": 2, "OXT": 3
}

VDW_DEPTH = 0.50   
VDW_R0 = 3.80      
VDW_SIGMA = 0.70   

# ── Primitive potentials ──────────────────────────────────────────────────────
def _lj(r, epsilon, r_min):
    if r <= 0.0: return 0.0
    ratio = r_min / r
    return epsilon * (ratio**12 - 2.0 * ratio**6)

def _coulomb(q1, q2, r):
    if r <= 0.0: return 0.0
    return (COULOMB_CONSTANT * q1 * q2) / (EPSILON_RELATIVE * max(r, MIN_ELECTROSTATIC_DISTANCE))

def _gaussian_well(r, depth, r0, sigma):
    if r <= 0.0: return 0.0
    return -depth * math.exp(-((r - r0) ** 2) / (2.0 * sigma**2))

# ── Helper parsing functions ──────────────────────────────────────────────────
def _atom_symbol(atom):
    sym = atom.get("atom_symbol")
    if sym: return str(sym).upper().strip()
    name = str(atom.get("atom_name", "")).strip()
    return name[:1].upper() if name else ""

def get_partial_charge(atom):
    # Dynamic continuous search: Check sidechain-specific token first, fall back to universal backbone
    return PARTIAL_CHARGES.get(f"{atom['residue_name']}_{atom['atom_name']}", PARTIAL_CHARGES.get(atom["atom_name"], 0.0))

def _water_capacity(atom):
    return WATER_HBOND_CAPACITY.get(f"{atom['residue_name']}_{atom['atom_name']}", WATER_HBOND_CAPACITY.get(atom["atom_name"], 0))

def _residue_key(atom, side):
    return (side, atom.get("chain_id", atom.get("chain", "")), atom["residue_name"], atom["residue_seq"], atom.get("insertion_code", atom.get("icode", "")))

# ── Functional Dispatch Map ───────────────────────────────────────────────────
SCORING_MODULES = {
    "disulfide_bond": lambda a1, a2, r, **_: _lj(r, 251.0, 2.05),
    "metal_coordination": lambda a1, a2, r, **_: _lj(r, METAL_WELL_DEPTHS.get(_atom_symbol(a1) if _atom_symbol(a1) in METAL_WELL_DEPTHS else _atom_symbol(a2), 62.8), 2.1),
    "salt_bridge": lambda a1, a2, r, **_: _coulomb(SALT_BRIDGE_CHARGES.get((a1["residue_name"], a1["atom_name"]), 0.0), SALT_BRIDGE_CHARGES.get((a2["residue_name"], a2["atom_name"]), 0.0), r),
    "hydrogen_bond": lambda a1, a2, r, **_: _gaussian_well(r, 8.0, 3.0, 0.9),
    "halogen_bond": lambda a1, a2, r, **_: _lj(r, 8.4, 3.2),
    "hydrophobic_contact": lambda a1, a2, r, **_: _gaussian_well(r, VDW_DEPTH, VDW_R0, VDW_SIGMA) if r <= 5.0 else 0.0,
    
    "pi-pi_stacking": lambda a1, a2, r, **kw: _lj(
        r, 
        *[[10.5, 4.0], [11.3, 3.8], [11.6, 5.0], [10.5, 4.5]][
            0 if kw.get("angle") is None else 
            1 if (kw["angle"] % 180.0 <= 30.0 or kw["angle"] % 180.0 >= 150.0) else 
            2 if (60.0 <= kw["angle"] % 180.0 <= 120.0) else 3
        ]
    ),
    
    "dipole-dipole": lambda a1, a2, r, **_: (
        0.0 if (r <= 0.0 or r > 4.0 or get_partial_charge(a1) == 0.0 or get_partial_charge(a2) == 0.0) else
        (-0.50 if (get_partial_charge(a1) * get_partial_charge(a2) < 0.0) else 0.25) * (COULOMB_CONSTANT * abs(get_partial_charge(a1) * get_partial_charge(a2)) / (EPSILON_RELATIVE * max(r, MIN_ELECTROSTATIC_DISTANCE)))
    )
}

# ── Public API ────────────────────────────────────────────────────────────────
def interaction_energy(interaction, inter_obj=None):
    itype = interaction["type"] or "unclassified"
    r = interaction["distance"]

    if r <= 0.0 or r > 5.0: return 0.0

    if itype == "unclassified":
        if _atom_symbol(interaction["atom1"]) == "C" and _atom_symbol(interaction["atom2"]) == "C":
            return _gaussian_well(r, VDW_DEPTH, VDW_R0, VDW_SIGMA)
        return 0.0

    evaluator = SCORING_MODULES.get(itype)
    if not evaluator: return 0.0

    kwargs = {}
    if itype == "pi-pi_stacking" and inter_obj is not None:
        from pdb_utils_inter import get_angle_between_normals
        a1, a2 = interaction["atom1"], interaction["atom2"]
        r1 = inter_obj.rings1.get(a1["residue_seq"] if inter_obj.chain1 in inter_obj.protein.chains() else a1["atom_seq"])
        r2 = inter_obj.rings2.get(a2["residue_seq"] if inter_obj.chain2 in inter_obj.protein.chains() else a2["atom_seq"])
        if r1 and r2:
            kwargs["angle"] = get_angle_between_normals(r1["normal"], r2["normal"])

    return evaluator(interaction["atom1"], interaction["atom2"], r, **kwargs)


def binding_enthalpy(inter_obj):
    breakdown = {}
    salt_bridge_groups = {}

    for interaction in inter_obj.interactions:
        itype = interaction["type"] or "unclassified"
        if itype == "salt_bridge":
            key = (_residue_key(interaction["atom1"], "side1"), _residue_key(interaction["atom2"], "side2"))
            salt_bridge_groups.setdefault(key, []).append(interaction)
            continue

        energy = interaction_energy(interaction, inter_obj)
        interaction["energy"] = energy
        breakdown[itype] = breakdown.get(itype, 0.0) + energy

    for group in salt_bridge_groups.values():
        r_min = min(item["distance"] for item in group)
        q1 = 1.0 if group[0]["atom1"]["residue_name"] in {"ARG", "LYS"} else -1.0
        q2 = -1.0 if group[0]["atom2"]["residue_name"] in {"ASP", "GLU"} else 1.0
        
        group_energy = _coulomb(q1, q2, r_min)
        breakdown["salt_bridge"] = breakdown.get("salt_bridge", 0.0) + group_energy

        weights = [1.0 / max(item["distance"], 1e-6) for item in group]
        norm = sum(weights)
        for item, w in zip(group, weights):
            item["energy"] = group_energy * w / norm if norm else 0.0

    interface_count = {}
    atom_lookup = {}
    for interaction in inter_obj.interactions:
        if interaction["type"] in {"hydrogen_bond", "salt_bridge"}:
            for atom in (interaction["atom1"], interaction["atom2"]):
                aid = atom["atom_seq"]
                interface_count[aid] = interface_count.get(aid, 0) + 1
                atom_lookup[aid] = atom

    breakdown["desolvation"] = sum(min(count, _water_capacity(atom_lookup[aid])) * WATER_HBOND_ENTHALPY for aid, count in interface_count.items())
    breakdown["TOTAL"] = sum(breakdown.values())
    
    return breakdown