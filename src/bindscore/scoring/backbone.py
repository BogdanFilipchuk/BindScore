"""
binding_entropy.backbone
========================

Backbone (and collective) configurational entropy via Elastic Network
Model + Normal Mode Analysis.

The protein backbone fluctuates around its mean structure with a spectrum
of vibrational modes - from slow, large-amplitude global motions to fast,
small-amplitude local jiggling. When two chains associate, the modes of
the complex are generally stiffer (higher frequencies) than those of the
free chains, so the complex has LESS vibrational entropy. This loss is
what we estimate here.

Method
------
1. Build a coarse-grained Anisotropic Network Model (ANM) on the
   Cα atoms: every Cα-Cα pair within an 8-A cutoff is connected by a
   Hookean spring of unit force constant.
2. Diagonalise the resulting 3N-by-3N Hessian (ProDy does this).
3. Compute the classical harmonic-oscillator entropy of every non-zero
   mode and sum them: S = sum_k k_B * (1 - ln(hbar*omega_k / k_B T)).
4. Repeat for each free chain alone.
5. Delta_S = S_complex - S_A - S_B.

This is a coarse approximation (no atomic detail, no anharmonicity, no
solvent coupling, no explicit treatment of the inter-chain rigid-body
modes) but it captures the change in soft collective motion which is what
the textbook configurational entropy is meant to track.

References
----------
Tirion (1996), Phys. Rev. Lett. 77:1905-1908.
    The original elastic-network paper, showing that low-frequency modes
    of a protein can be reproduced with a single force constant.

Bahar, Atilgan & Erman (1997), Folding & Design 2:173-181.
    Gaussian Network Model: introduced the Cα-only coarse-graining.

Atilgan et al. (2001), Biophys. J. 80:505-515.
    Anisotropic Network Model: the 3D anisotropic version we use here.

Bakan, Meireles & Bahar (2011), Bioinformatics 27:1575-1577.
    The ProDy package, which implements the ANM in Python and is the
    primary computational dependency of this module.

Hsu et al. (2008), Proteins 71:455-466.
    Discusses the use of NMA-derived configurational entropies in binding
    free energy calculations. Also documents the limitations: NMA cannot
    capture anharmonic motions or solvent-driven fluctuations.

Note on the trans/rot double-counting
-------------------------------------
The complex has 6 more internal vibrational modes than the sum of the
two free chains (the 6 modes that were 'zero-frequency' rigid-body modes
of one chain relative to the other become finite-frequency inter-chain
modes in the bound state). These 6 modes carry the entropy of relative
trans+rot motion of the chains, which is also estimated by the trans_rot
module via the Sackur-Tetrode cage approximation. There is therefore some
degree of double-counting between this module and trans_rot. For a
2-minute triage estimate this is acceptable, but be aware: the totals
reported by this package may overstate the binding entropy penalty by
~5-15% from this overlap.
"""

from __future__ import annotations
import os
import math
import warnings
from dataclasses import dataclass

import numpy as np

from .utils import split_chains_to_tempfiles


# ----- Physical constants used in the harmonic-oscillator entropy -----
_R_KCAL = 1.987204e-3   # gas constant, kcal/mol/K
_KB_J   = 1.380649e-23  # Boltzmann constant, J/K
_HBAR_J = 1.054571817e-34  # reduced Planck, J*s


@dataclass
class BackboneResult:
    """Container for the NMA-based backbone entropy estimate."""
    S_complex: float    # vibrational entropy of the complex, kcal/mol/K
    S_chain_a: float    # vibrational entropy of chain A alone
    S_chain_b: float    # vibrational entropy of chain B alone
    delta_S: float      # S_complex - S_A - S_B, kcal/mol/K
    minusT_deltaS: float  # -T*Delta_S, kcal/mol


# -----------------------------------------------------------------------------
# Core ENM + NMA + entropy calculation
# -----------------------------------------------------------------------------

def _vibrational_entropy_from_eigvals(
    eigenvalues: np.ndarray,
    T: float,
    gamma: float,
) -> float:
    """
    Compute the classical harmonic-oscillator entropy for a set of normal
    mode eigenvalues, returning S in kcal/mol/K.

    Classical limit (kT >> hbar*omega), per mode:
        S_k = k_B * (1 - ln(hbar*omega_k / k_B T))

    The eigenvalues from a unit-mass, unit-spring ANM are in units of
    (force constant)/(mass). To get a physical frequency we need a
    conversion; the simplest assumption is to set the spring constant
    `gamma` to a representative value (default 1 kcal/mol/A^2) and
    masses to the average residue mass (110 Da). This gives entropies in
    the right order of magnitude. Because we are taking a DIFFERENCE
    between bound and free states with the same conversion factors, the
    absolute calibration partly cancels.
    """
    # Convert spring constant from kcal/mol/A^2 to SI (J/m^2 per molecule)
    gamma_SI = gamma * 4184.0 / 6.02214076e23 / 1e-20   # = J/(m^2)

    # Average residue mass (110 Da -> kg)
    mass_kg = 110.0 * 1.66053906660e-27

    # omega_k = sqrt(eigenvalue * gamma_SI / mass_kg)
    # Skip any non-positive eigenvalues (numerical noise near zero modes)
    valid = eigenvalues[eigenvalues > 1e-6]
    omegas = np.sqrt(valid * gamma_SI / mass_kg)   # rad/s

    # Classical harmonic-oscillator entropy per mode
    # S_k = k_B * (1 - ln(hbar*omega / kT))
    arg = _HBAR_J * omegas / (_KB_J * T)
    # Suppress log warnings for any residual near-zero modes
    with np.errstate(invalid="ignore", divide="ignore"):
        s_per_mode_J = _KB_J * (1.0 - np.log(arg))

    # Replace any non-finite from log of zero
    s_per_mode_J = s_per_mode_J[np.isfinite(s_per_mode_J)]

    # Sum and convert J/K per molecule -> kcal/mol/K per mole
    S_total_J_per_mol = s_per_mode_J.sum() * 6.02214076e23
    return S_total_J_per_mol / 4184.0


def _anm_entropy(pdb_path: str, T: float, cutoff: float, gamma: float) -> float:
    """
    Build an ANM on the Cα atoms of a PDB file and return its
    vibrational entropy in kcal/mol/K.

    Uses ProDy (Bakan et al. 2011, Bioinformatics 27:1575).
    """
    # We import ProDy lazily to avoid making it a hard dependency for
    # users who only want the trans/rot or sidechain estimates.
    import prody
    prody.confProDy(verbosity="none")

    atoms = prody.parsePDB(pdb_path).select("calpha")
    if atoms is None or len(atoms) < 4:
        # Too small to do meaningful NMA on
        warnings.warn(f"Too few Cα atoms in {pdb_path} for NMA; returning 0.")
        return 0.0

    anm = prody.ANM("nma")
    anm.buildHessian(atoms, cutoff=cutoff, gamma=gamma)

    # Diagonalise all 3N-6 internal modes. `zeros=False` makes ProDy drop
    # the 6 rigid-body modes automatically.
    n_modes = 3 * len(atoms) - 6
    anm.calcModes(n_modes=n_modes, zeros=False)

    eigvals = anm.getEigvals()
    return _vibrational_entropy_from_eigvals(eigvals, T, gamma)


# -----------------------------------------------------------------------------
# Public entry point
# -----------------------------------------------------------------------------

def compute(
    complex_pdb: str,
    chain_a: str,
    chain_b: str,
    T: float = 300.0,
    cutoff: float = 8.0,
    gamma: float = 1.0,
) -> BackboneResult:
    """
    Estimate -T*Delta_S for the backbone/collective configurational
    entropy change on binding.

    Parameters
    ----------
    complex_pdb : str
        Path to a PDB file containing both chains in their bound pose.
    chain_a, chain_b : str
        Single-character chain IDs.
    T : float
        Temperature in Kelvin (default 300).
    cutoff : float
        ANM contact cutoff in Angstroms (default 8.0, the standard ANM
        value from Atilgan et al. 2001).
    gamma : float
        Uniform spring force constant in kcal/mol/A^2. The numerical
        value mostly cancels in the bound-minus-free difference; we
        keep it as a tunable for transparency.
    """
    # Free-state entropies: each chain on its own
    path_a, path_b = split_chains_to_tempfiles(complex_pdb, chain_a, chain_b)
    try:
        S_a = _anm_entropy(path_a, T, cutoff, gamma)
        S_b = _anm_entropy(path_b, T, cutoff, gamma)
    finally:
        os.unlink(path_a)
        os.unlink(path_b)

    # Bound-state entropy: the complex as one assembly
    S_complex = _anm_entropy(complex_pdb, T, cutoff, gamma)

    delta_S = S_complex - S_a - S_b
    minusT_dS = -T * delta_S

    return BackboneResult(
        S_complex=S_complex,
        S_chain_a=S_a,
        S_chain_b=S_b,
        delta_S=delta_S,
        minusT_deltaS=minusT_dS,
    )
