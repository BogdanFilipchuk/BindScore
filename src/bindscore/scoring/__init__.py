"""
binding_entropy
===============

Fast estimation of the entropic contribution to protein-protein binding.

This package decomposes the binding entropy into four physically distinct
components and computes each one with a fast approximation rather than the
rigorous (and expensive) Restraint-Release FEP scheme of Singh & Warshel
(Proteins, 2010). The goal is a sub-2-minute estimate suitable for
ranking and triage, NOT a publication-quality absolute value.

Components computed
-------------------
1. Translational + rotational entropy   (trans_rot module)
2. Hydrophobic solvation entropy        (hydrophobic module)
3. Sidechain conformational entropy     (sidechain module)
4. Backbone collective entropy via NMA  (backbone module)

Components explicitly NOT computed
----------------------------------
- Polarization entropy of the solvent (would require solvent FEP)
- Bulk protein configurational entropy beyond the interface
The output reports -T*Delta_S in kcal/mol so it can be added directly to
an enthalpy estimate to obtain Delta_G.

Public entry point
------------------
compute_total_entropy(complex_pdb, chain_a, chain_b)
    Returns a single number: -T*Delta_S_binding in kcal/mol.
"""

from .summary import compute_total_entropy

__all__ = ["compute_total_entropy"]
__version__ = "0.1.0"
