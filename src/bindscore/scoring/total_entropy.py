"""
binding_entropy.summary
=======================

Top-level entry point: combines the four entropy submodules into a single
estimate of the binding entropy contribution to Delta_G.

Modules used
------------
trans_rot   - Sackur-Tetrode translational + rigid-rotor rotational
              entropy of one chain relative to the other in the cage of
              the complex. Largest single contribution for most PPIs.
              References: Page & Jencks 1971; Finkelstein & Janin 1989.

hydrophobic - SASA-based estimate of the entropic part of the
              hydrophobic effect. The only typically FAVOURABLE component.
              References: Chothia 1974; Spolar, Livingstone & Record 1992.

sidechain   - Per-residue sidechain conformational entropy loss for
              interface residues, using the Pickett-Sternberg scale.
              References: Pickett & Sternberg 1993; Doig & Sternberg 1995;
              Halle & Nilsson 2002 (interface-specific validation).

backbone    - Coarse-grained backbone collective entropy from an
              Anisotropic Network Model normal-mode analysis via ProDy.
              References: Tirion 1996; Bahar et al. 1997; Atilgan et al.
              2001; Bakan, Meireles & Bahar 2011 (the ProDy package).

Explicitly omitted (see package __init__ docstring)
---------------------------------------------------
- Polarization entropy of solvent: requires solvent FEP, infeasible in 2 min.
- Bulk protein configurational entropy beyond the interface: assumed to
  cancel between bound and free states.

Sign convention
---------------
All values are reported as -T*Delta_S in kcal/mol. POSITIVE means the
component opposes binding (entropic penalty); NEGATIVE means it favours
binding. The total can be added directly to a Delta_H estimate to obtain
a Delta_G estimate.
"""

from __future__ import annotations
from dataclasses import dataclass

from . import trans_rot
from . import hydrophobic
from . import sidechain
from . import backbone



@dataclass
class Protein_Entropy:
    """Per-component binding entropy breakdown. All values in J/(mol·K)."""
    dS_trans_rot:   float
    dS_hydrophobic: float
    dS_sidechain:   float
    dS_backbone:    float
    dS_total:       float


# -----------------------------------------------------------------------------
# Main public function
# -----------------------------------------------------------------------------

def compute_total_entropy(
    complex_pdb: str,
    chain_a: str,
    chain_b: str,
    T: float = 300.0,
    return_breakdown: bool = False,
    HYDROPHOBIC_SETTING: str = "SASA",
) -> float | Protein_Entropy:
    """
    Compute the total binding entropy ΔS in J/(mol·K).

    Parameters
    ----------
    complex_pdb : str
        Path to a PDB file of the bound complex (both chains together).
    chain_a, chain_b : str
        Single-character chain IDs of the two binding partners.
    T : float
        Temperature in Kelvin (default 300).
    return_breakdown : bool
        If False (default), return the scalar total ΔS in J/(mol·K).
        If True, return an Protein_Entropy with each component exposed.
    HYDROPHOBIC_SETTING : str
        "SASA" (default) or "Sun" — selects the hydrophobic entropy method.

    Returns
    -------
    float (default) or Protein_Entropy — values in J/(mol·K)
    """
    # Each module is independent and can fail without bringing the others
    # down; we use try/except so a single module's failure doesn't kill
    # the whole estimate. Failed modules contribute 0 and emit a warning.

    import warnings
    pdb_path = str(complex_pdb)

    try:
        tr = trans_rot.compute(pdb_path, chain_a, chain_b, T=T)
        dS_tr = tr.total
    except Exception as exc:
        warnings.warn(f"trans_rot module failed: {exc}; contributing 0.")
        dS_tr = 0.0

    if HYDROPHOBIC_SETTING == "Sun":
        try:
            from . import protein_solvent_entropy
            import bindscore.pdb_file_treatment.protein_radius_estimation as radius_estimation
            coords_a   = radius_estimation.load_atoms(pdb_path, chain_id=chain_a)
            centroid_a = radius_estimation.find_centroid(coords_a)
            R_a        = radius_estimation.estimate_radius(coords_a, centroid_a)[0]
            coords_b   = radius_estimation.load_atoms(pdb_path, chain_id=chain_b)
            centroid_b = radius_estimation.find_centroid(coords_b)
            R_b        = radius_estimation.estimate_radius(coords_b, centroid_b)[0]
            dS_hp = (protein_solvent_entropy.dS_interfacial(R_a, T)
                   + protein_solvent_entropy.dS_bulk(R_a)
                   + protein_solvent_entropy.dS_interfacial(R_b, T)
                   + protein_solvent_entropy.dS_bulk(R_b))
        except Exception as exc:
            warnings.warn(f"protein solvent entropy module failed: {exc}; contributing 0.")
            dS_hp = 0.0
    else:  # "SASA"
        try:
            hp = hydrophobic.compute(pdb_path, chain_a, chain_b, T=T)
            dS_hp = hp.deltaS
        except Exception as exc:
            warnings.warn(f"hydrophobic module failed: {exc}; contributing 0.")
            dS_hp = 0.0

    try:
        sc = sidechain.compute(pdb_path, chain_a, chain_b, T=T)
        dS_sc = sc.deltaS
    except Exception as exc:
        warnings.warn(f"sidechain module failed: {exc}; contributing 0.")
        dS_sc = 0.0

    try:
        bb = backbone.compute(pdb_path, chain_a, chain_b, T=T)
        dS_bb = bb.deltaS
    except Exception as exc:
        warnings.warn(f"backbone (NMA) module failed: {exc}; contributing 0.")
        dS_bb = 0.0

    dS_total = dS_tr + dS_hp + dS_sc + dS_bb

    if return_breakdown:
        return Protein_Entropy(
            dS_trans_rot=dS_tr,
            dS_hydrophobic=dS_hp,
            dS_sidechain=dS_sc,
            dS_backbone=dS_bb,
            dS_total=dS_total,
        )
    return dS_total
