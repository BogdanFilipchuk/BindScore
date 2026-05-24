"""
BROKEN MODULE.

binding_entropy.sidechain
=========================

Side-chain configurational entropy change on binding.

Method
------
Use the Pickett–Sternberg (1993) empirical scale: a single ΔS value per
residue type for full burial of its side chain.  For every residue that
loses solvent-accessible surface area (SASA) on binding, accumulate its
scale value.  The total ΔS is in J/(mol·K), negative because burial
removes rotameric freedom.

Sign convention
---------------
ΔS = S_bound − S_free, in J/(mol·K).
Burial only removes freedom, so deltaS ≤ 0.
Negative = entropy lost on binding = unfavourable contribution to ΔG_bind.

Expected magnitudes
-------------------
For typical protein–protein complexes with 20–30 interface residues,
ΔS_sidechain falls in the range  −200 to −800 J/(mol·K)
(i.e. −TΔS ≈ +14 to +57 kJ/mol at 300 K, or 3 to 14 kcal/mol).

References
----------
Pickett & Sternberg (1993), J. Mol. Biol. 231:825-839.
    Source of the per-residue −TΔS values (Table 1, at T_ref = 298 K).
Doig & Sternberg (1995), Protein Sci. 4:2247-2251.
    Confirmed the scale transfers from folding to binding contexts.
Halle & Nilsson (2002), Proteins 49:154-166.
    Validated rotamer-based estimates specifically for PPI interfaces.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Dict

from Bio.PDB import PDBParser

from .entropy_utils import identify_interface_residues


# -----------------------------------------------------------------------------
# Pickett–Sternberg per-residue scale
# -----------------------------------------------------------------------------
# Original table:  −T·ΔS in kcal/mol at T_ref = 298 K.
# We convert ONCE at import to ΔS in J/(mol·K) using:
#       ΔS [J/(mol·K)]  =  −(−T·ΔS) [kcal/mol] × 4184 / T_ref
# The leading minus sign flips the convention so ΔS is negative for an
# entropy LOSS. ALA / GLY / PRO have no side-chain rotameric freedom.

_PS_REF_T_K = 298.0
_KCAL_TO_J  = 4184.0

_MINUS_T_DELTAS_KCAL: Dict[str, float] = {
    "ALA": 0.00, "ARG": 2.13, "ASN": 0.81, "ASP": 0.61, "CYS": 0.55,
    "GLN": 1.29, "GLU": 1.06, "GLY": 0.00, "HIS": 0.99, "ILE": 0.84,
    "LEU": 0.78, "LYS": 1.56, "MET": 1.04, "PHE": 0.58, "PRO": 0.00,
    "SER": 0.55, "THR": 0.48, "TRP": 0.55, "TYR": 0.61, "VAL": 0.51,
}

_DELTAS_J_PER_MOL_K: Dict[str, float] = {
    res: -(v * _KCAL_TO_J) / _PS_REF_T_K
    for res, v in _MINUS_T_DELTAS_KCAL.items()
}


# -----------------------------------------------------------------------------
# Result container — simple: one number + one diagnostic count
# -----------------------------------------------------------------------------

@dataclass
class SidechainResult:
    """Result of the side-chain configurational entropy estimate."""
    deltaS: float = 0.0
    # ΔS_binding in J/(mol·K). Always ≤ 0.
    # Negative = entropy lost on binding (unfavourable for ΔG).
    n_interface_residues: int = 0
    # Total number of residues at the interface across both chains
    # (diagnostic, not used in any sum).


# -----------------------------------------------------------------------------
# Public entry point
# -----------------------------------------------------------------------------

def compute(
    complex_pdb: str,
    chain_a: str,
    chain_b: str,
    T: float = 300.0,
) -> SidechainResult:
    """
    Estimate the side-chain configurational entropy change on binding.

    Parameters
    ----------
    complex_pdb : str
        Path to PDB containing both chains in the bound pose.
    chain_a, chain_b : str
        Single-character chain IDs.
    T : float
        Temperature in K. Accepted for API consistency; not used in the
        formula because ΔS itself (a configurational entropy) is taken
        as T-independent over the relevant range.

    Returns
    -------
    SidechainResult
        ``deltaS`` is the total side-chain ΔS_binding in J/(mol·K).
        Always ≤ 0 (burial only removes rotameric freedom).
    """
    # Step 1: identify which residues lose SASA on binding.
    interface = identify_interface_residues(complex_pdb, chain_a, chain_b)

    # Step 2: build a {(chain_id, resnum) -> resname} map so we can look
    # up each interface residue's per-residue ΔS value.
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("s", complex_pdb)
    resname_map: Dict[tuple, str] = {}
    for model in structure:
        for chain in model:
            if chain.id not in (chain_a, chain_b):
                continue
            for res in chain:
                if res.id[0] != " ":   # skip HETATM / waters
                    continue
                resname_map[(chain.id, res.id[1])] = res.get_resname()
        break   # first model only

    # Step 3: sum per-residue contributions.
    total = 0.0
    n_residues = 0
    for chain_id, resnums in interface.items():
        for resnum in resnums:
            resname = resname_map.get((chain_id, resnum))
            if resname is None:
                continue   # interface residue we couldn't map back (rare)
            total += _DELTAS_J_PER_MOL_K.get(resname, 0.0)
            n_residues += 1

    return SidechainResult(deltaS=total, n_interface_residues=n_residues)
