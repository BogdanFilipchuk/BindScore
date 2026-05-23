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
class EntropySummary:
    """Full per-component breakdown plus the total."""
    minusT_dS_trans_rot: float
    minusT_dS_hydrophobic: float
    minusT_dS_sidechain: float
    minusT_dS_backbone: float
    minusT_dS_total: float


# -----------------------------------------------------------------------------
# Main public function
# -----------------------------------------------------------------------------

def compute_total_entropy(
    complex_pdb: str,
    chain_a: str,
    chain_b: str,
    T: float = 300.0,
    return_breakdown: bool = False,
) -> float | EntropySummary:
    """
    Compute the total -T*Delta_S of binding (kcal/mol).

    Parameters
    ----------
    complex_pdb : str
        Path to a PDB file of the bound complex (both chains together).
    chain_a, chain_b : str
        Single-character chain IDs of the two binding partners.
    T : float
        Temperature in Kelvin (default 300).
    return_breakdown : bool
        If False (default), return the single scalar total. If True,
        return an EntropySummary dataclass with each component exposed.

    Returns
    -------
    float (default) or EntropySummary
    """
    # Each module is independent and can fail without bringing the others
    # down; we use try/except so a single module's failure doesn't kill
    # the whole estimate. Failed modules contribute 0 and emit a warning.

    import warnings

    try:
        tr = trans_rot.compute(complex_pdb, chain_a, chain_b, T=T)
        minusT_tr = tr.total
    except Exception as exc:
        warnings.warn(f"trans_rot module failed: {exc}; contributing 0.")
        minusT_tr = 0.0

    try:
        hp = hydrophobic.compute(complex_pdb, chain_a, chain_b)
        minusT_hp = hp.minusT_deltaS
    except Exception as exc:
        warnings.warn(f"hydrophobic module failed: {exc}; contributing 0.")
        minusT_hp = 0.0

    try:
        sc = sidechain.compute(complex_pdb, chain_a, chain_b)
        minusT_sc = sc.minusT_deltaS
    except Exception as exc:
        warnings.warn(f"sidechain module failed: {exc}; contributing 0.")
        minusT_sc = 0.0

    try:
        bb = backbone.compute(complex_pdb, chain_a, chain_b, T=T)
        minusT_bb = bb.minusT_deltaS
    except Exception as exc:
        warnings.warn(f"backbone (NMA) module failed: {exc}; contributing 0.")
        minusT_bb = 0.0

    total = minusT_tr + minusT_hp + minusT_sc + minusT_bb

    if return_breakdown:
        return EntropySummary(
            minusT_dS_trans_rot=minusT_tr,
            minusT_dS_hydrophobic=minusT_hp,
            minusT_dS_sidechain=minusT_sc,
            minusT_dS_backbone=minusT_bb,
            minusT_dS_total=total,
        )
    return total
