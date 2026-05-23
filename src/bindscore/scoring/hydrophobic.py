"""
binding_entropy.hydrophobic
===========================

Hydrophobic contribution to the binding entropy from water reorganisation.

When a non-polar surface is buried at an interface, ordered water
molecules around it are released into bulk water. This is one of the
classic favourable contributions to binding (negative -T*Delta_S, i.e. it
*lowers* the free energy).

Rigorous calculation would require simulating water explicitly around the
free and bound surfaces (this is what Singh & Warshel do via Restraint-
Release FEP on the solvent). The standard fast approximation is to scale
the change in non-polar SASA by an empirical coefficient.

References
----------
Chothia (1974), Nature 248:338-339.
    Originally proposed that hydrophobic burial scales linearly with the
    change in non-polar surface area, with a coefficient of about
    25 cal/mol/A^2 for the *total* hydrophobic free energy.

Spolar, Livingstone & Record (1992), Biochemistry 31:3947-3955.
    Decomposed the hydrophobic effect into enthalpy and entropy. The
    entropic part of the hydrophobic effect at 300 K is approximately
    0.024 kcal/mol/A^2 of non-polar surface buried (this is the value
    we use as the default below).

Eisenberg & McLachlan (1986), Nature 319:199-203.
    Per-atom-type SASA coefficients (Eisenberg-McLachlan ASP scale). Our
    simpler implementation uses a single coefficient for all non-polar
    atoms, which is the standard further approximation.
"""

from __future__ import annotations
import os
from dataclasses import dataclass
from typing import Set

from .utils import split_chains_to_tempfiles


# Non-polar atom set: standard convention (C and S; H is implicit and
# excluded since most PDB files lack hydrogens, and freesasa handles this
# correctly by default). N and O are treated as polar.
_NONPOLAR_ELEMENTS: Set[str] = {"C", "S"}


# Default coefficient from Spolar et al. 1992 for the entropic component
# of the hydrophobic effect at 300 K. Sign convention: positive value
# means buried non-polar area FAVOURS binding (i.e. -T*Delta_S is
# *negative*).
_DEFAULT_HYDROPHOBIC_COEFF = 0.024  # kcal/mol/A^2


@dataclass
class HydrophobicResult:
    """Result of the SASA-based hydrophobic entropy estimate."""
    delta_nonpolar_sasa: float   # A^2 buried (positive = surface lost on binding)
    coefficient: float           # kcal/mol/A^2 used
    minusT_deltaS: float         # kcal/mol (typically negative = favourable)


# -----------------------------------------------------------------------------
# Per-atom-class SASA via freesasa
# -----------------------------------------------------------------------------

def _nonpolar_sasa(pdb_path: str) -> float:
    """
    Compute the total non-polar SASA of a PDB file (square Angstroms).

    Implementation note on the freesasa API:
      - `freesasa.calc(structure)` is a module-level function (NOT a
        Calc class). Returns a Result.
      - `Result.atomArea(i)` gives per-atom SASA in A^2.
      - `Structure.atomName(i)` gives the 4-character atom name; freesasa
        does NOT expose a separate element accessor, so we derive the
        element from the atom name following the standard PDB convention:
        the first non-space, non-digit character is the element symbol
        for the protein backbone and sidechain atoms (CA, CB, CG -> C;
        N, NZ -> N; O, OG -> O; SD -> S; etc.).
    """
    import freesasa

    structure = freesasa.Structure(pdb_path)
    result = freesasa.calc(structure)

    total_nonpolar = 0.0
    n_atoms = structure.nAtoms()
    for i in range(n_atoms):
        # Derive the element from the atom name. PDB atom names like " CA "
        # have the element in the first letter (skipping any leading digits).
        name = structure.atomName(i).strip()
        elem = ""
        for ch in name:
            if ch.isalpha():
                elem = ch.upper()
                break

        if elem in _NONPOLAR_ELEMENTS:
            total_nonpolar += result.atomArea(i)

    return total_nonpolar


# -----------------------------------------------------------------------------
# Public entry point
# -----------------------------------------------------------------------------

def compute(
    complex_pdb: str,
    chain_a: str,
    chain_b: str,
    coefficient: float = _DEFAULT_HYDROPHOBIC_COEFF,
) -> HydrophobicResult:
    """
    Estimate -T*Delta_S_hydrophobic from the change in non-polar SASA.

    Calculation
    -----------
    1. Compute non-polar SASA of the complex.
    2. Compute non-polar SASA of each chain in isolation (after splitting).
    3. Delta_SASA = (SASA_A_free + SASA_B_free) - SASA_complex
       Positive value = surface buried on binding.
    4. -T*Delta_S_hydrophobic = -coefficient * Delta_SASA
       Negative result (favourable) when surface is buried.

    The sign convention is the source of the most common bug in
    implementations of this formula: confirm against a known case
    (e.g. barnase-barstar) before trusting the absolute number.

    Parameters
    ----------
    coefficient : float
        kcal/mol/A^2. Default 0.024 from Spolar et al. 1992 represents
        the entropic part of the hydrophobic effect at 300 K.
    """
    # Free state SASA: each chain on its own
    path_a, path_b = split_chains_to_tempfiles(complex_pdb, chain_a, chain_b)
    try:
        sasa_a_free = _nonpolar_sasa(path_a)
        sasa_b_free = _nonpolar_sasa(path_b)
    finally:
        os.unlink(path_a)
        os.unlink(path_b)

    # Bound state SASA: the assembled complex
    sasa_complex = _nonpolar_sasa(complex_pdb)

    delta_sasa = (sasa_a_free + sasa_b_free) - sasa_complex
    # Sign: hydrophobic burial RELEASES ordered water -> favourable -> the
    # entropy CHANGE is POSITIVE -> -T*Delta_S is NEGATIVE.
    minusT_dS = -coefficient * delta_sasa

    return HydrophobicResult(
        delta_nonpolar_sasa=delta_sasa,
        coefficient=coefficient,
        minusT_deltaS=minusT_dS,
    )
