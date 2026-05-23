"""
solvation_entropy.py
────────────────────────────────────────────────────────────────────────────────
Solvent–solute entropy terms for a protein in water, following:
 
    Sun, Q. (2022). "The Hydrophobic Effects: Our Current Understanding."
    Molecules, 27(20), 7009.  https://doi.org/10.3390/molecules27207009
 
Water is partitioned into two populations at the protein surface:
  • Interfacial water  – weaker DA H-bonds  → higher entropy
  • Bulk water         – tetrahedral DDAA    → lower entropy (more ordered)
 
Public API
──────────
    solvent_interfacial_entropy(R, T) -> float   [J mol⁻¹ K⁻¹]
    solvent_bulk_entropy(R)           -> float   [J mol⁻¹ K⁻¹]
 
All other symbols are private (_-prefixed) implementation details.
"""
 
from __future__ import annotations
 
__all__ = ["solvent_interfacial_entropy", "solvent_bulk_entropy"]
 
# ── Water constants (Sun 2022, 293 K / 0.1 MPa) ──────────────────────────────
# Van't Hoff fit to Raman OH-stretch intensities (Fig. 5)
_DH_DDAA: float = 11.35e3      # J mol⁻¹   enthalpy of DDAA tetrahedral H-bond
_DS_DDAA: float = 29.66        # J mol⁻¹ K⁻¹  entropy of DDAA H-bond
 
# Effective radius of one water molecule (d = 3.8 Å → r = 1.9 Å, Sec. 3)
_R_H2O: float = 1.9            # Å
 
 
# ── Private helpers ───────────────────────────────────────────────────────────
 
def _interfacial_fraction(R_angstrom: float) -> float:
    """
    Fraction of water molecules in the interfacial layer for a spherical
    solute of radius *R_angstrom* (Å):
 
        f = 4 · r_H₂O / R          (Sun 2022, Sec. 3)
 
    Parameters
    ----------
    R_angstrom : sphere-equivalent protein radius in Ångströms
 
    Returns
    -------
    float : dimensionless interfacial fraction  (0 < f ≤ 1)
    """
    return 4.0 * _R_H2O / R_angstrom
 
 
# ── Public functions ──────────────────────────────────────────────────────────
 
def ds_interfacial(R_angstrom: float, T: float, IF: float) -> float:
    """
    Entropy contribution of the interfacial water layer (J mol⁻¹ K⁻¹).
 
    DA bonds at the protein surface are weaker than bulk DDAA bonds, so
    interfacial molecules carry higher configurational entropy.  The gain
    relative to bulk scales with the interfacial fraction *f*:
 
        ΔS_interfacial = (ΔH_DDAA / T) · f
 
    Parameters
    ----------
    R_angstrom : sphere-equivalent protein radius (Å)
    T          : temperature (K)
 
    Returns
    -------
    float : ΔS_interfacial in J mol⁻¹ K⁻¹  (positive — entropy gain)
    """
    f = IF
    return (_DH_DDAA / T) * f
 
 
def ds_bulk(R_angstrom: float, IF: float) -> float:
    """
    Entropy contribution of the bulk water surrounding the protein (J mol⁻¹ K⁻¹).
 
    Bulk water retains the ordered DDAA tetrahedral network; the ordering
    cost scales with the bulk fraction (1 − f):
 
        ΔS_bulk = −ΔS_DDAA · (1 − f)
 
    Parameters
    ----------
    R_angstrom : sphere-equivalent protein radius (Å)
 
    Returns
    -------
    float : ΔS_bulk in J mol⁻¹ K⁻¹  (negative — entropy cost)
    """
    f = IF
    return -_DS_DDAA * (1.0 - f)
