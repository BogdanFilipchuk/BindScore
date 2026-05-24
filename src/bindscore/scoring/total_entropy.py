"""
total_entropy.py
=======================

Top-level entry point: combines the four entropy submodules into a single estimate of the binding entropy contribution to Delta_G.

Modules used
------------
trans_rot, hydrophobic, sidechain and backbone submodules of the package.

Explicitly omitted (see package __init__ docstring and Jupyter)
---------------------------------------------------
- Polarization entropy of solvent: requires solvent FEP modelisation, infeasible in the goal time frame goal set by us in the design of the package.
- Bulk protein configurational entropy beyond the interface: assumed to cancel between bound and free states. (Would require using MD engines like OpenMM for simulations like RR-FEP. We tried, but for a 3 atom mo)

Sign convention
---------------
All values are reported as -T*Delta_S in kcal/mol. POSITIVE means the
component opposes binding (entropic penalty); NEGATIVE means it favours
binding. The total can be added directly to a Delta_H estimate to obtain
a Delta_G estimate.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional

from . import trans_rot
from . import protein_solvent_entropy
from . import sidechain
from . import backbone
from .trans_rot              import TransRotResult
from .protein_solvent_entropy import HydrophobicResult
from .sidechain              import SidechainResult
from .backbone               import BackboneResult


### MAIN DATACLASS
"""
Binding_Entropy_Summary will be returned by the module for final results. 
It contains all 4 contributions as proprietary objects with their own reuslts structure (see other modules),
as well as the 5 float variables taken from those objects giving the entropy dS results that were found.
"""
@dataclass
class Binding_Entropy_Summary:
    """Per-component binding entropy breakdown. All values in J/(mol·K)."""
    dS_trans_rot:   float
    dS_hydrophobic: float
    dS_sidechain:   float
    dS_backbone:    float
    dS_total:       float
    "List of submodule result objects — None if that module failed or was not run"
    trans_rot_detail:   Optional[TransRotResult]    = field(default=None, repr=False)
    hydrophobic_detail: Optional[HydrophobicResult] = field(default=None, repr=False)
    sidechain_detail:   Optional[SidechainResult]   = field(default=None, repr=False)
    backbone_detail:    Optional[BackboneResult]    = field(default=None, repr=False)

### MAIN FUNCTION
def compute_total_entropy(
    complex_pdb: str,
    chain_a: str,
    chain_b: str,
    T: float = 300.0,
    return_breakdown: bool = False,
) -> float | Binding_Entropy_Summary:
    """
    Compute the total binding entropy ΔS in J/(mol·K).

    Parameters
    ----------
    complex_pdb : str or pathlib path
        Path to a PDB file of the bound complex (both chains together).
    chain_a, chain_b : str
        Single-character chain IDs of the two binding partners.
    T : float
        Temperature in Kelvin (default 300).
    return_breakdown : bool
        If False (default), return the scalar total ΔS in J/(mol·K).
        If True, return a Binding_Entropy_Summary with each component exposed.

    Returns
    -------
    float (default) or Binding_Entropy_Summary — values in J/(mol·K)
    """
    # Each module is independent and can fail without bringing the others
    # down; we use try/except so a single module's failure doesn't kill
    # the whole estimate. Failed modules contribute 0 and emit a warning.

    import warnings
    pdb_path = str(complex_pdb)

    tr_result: Optional[TransRotResult] = None
    try:
        tr_result = trans_rot.compute(pdb_path, chain_a, chain_b, T=T)
        dS_tr = tr_result.total
    except Exception as exc:
        warnings.warn(f"trans_rot module failed: {exc}; contributing 0.")
        dS_tr = 0.0

    hp_result: Optional[HydrophobicResult] = None
    try:
        hp_result = protein_solvent_entropy.compute(pdb_path, chain_a, chain_b, T=T)
        dS_hp = hp_result.deltaS
    except Exception as exc:
        warnings.warn(f"hydrophobic (Sun) module failed: {exc}; contributing 0.")
        dS_hp = 0.0

    sc_result: Optional[SidechainResult] = None
    try:
        sc_result = sidechain.compute(pdb_path, chain_a, chain_b, T=T)
        dS_sc = sc_result.deltaS
    except Exception as exc:
        warnings.warn(f"sidechain module failed: {exc}; contributing 0.")
        dS_sc = 0.0

    bb_result: Optional[BackboneResult] = None
    try:
        bb_result = backbone.compute(pdb_path, chain_a, chain_b, T=T)
        dS_bb = bb_result.deltaS
    except Exception as exc:
        warnings.warn(f"backbone (NMA) module failed: {exc}; contributing 0.")
        dS_bb = 0.0

    dS_total = dS_tr + dS_hp + dS_sc + dS_bb

    if return_breakdown:
        return Binding_Entropy_Summary(
            dS_trans_rot=dS_tr,
            dS_hydrophobic=dS_hp,
            dS_sidechain=dS_sc,
            dS_backbone=dS_bb,
            dS_total=dS_total,
            trans_rot_detail=tr_result,
            hydrophobic_detail=hp_result,
            sidechain_detail=sc_result,
            backbone_detail=bb_result,
        )
    return dS_total
