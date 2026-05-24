"""
binding_entropy.backbone
========================

Backbone (collective) configurational entropy change on binding,
estimated from an Anisotropic Network Model + Normal Mode Analysis,
with Tidor–Karplus-style mode matching between the bound and free states.

Method
------
1. Build a Cα-level ANM for chain A alone, chain B alone, and the complex.
2. For each free-chain mode j, find the complex mode k that represents
   the same collective motion (highest overlap with V_a[:,j] in the
   chain-A subspace of the complex eigenvectors; analogously for B).
   Matching is done by HUNGARIAN ASSIGNMENT on the combined (A + B)
   set of free-chain modes against all complex modes, so each complex
   mode is matched to at most one free-chain mode.
3. For each matched pair (j, k), accumulate the log-ratio of eigenvalues:

        ΔS_vib  =  (R / 2) · Σ_paired  ln( λ_free / λ_complex )

   This is the closed-form classical-HO entropy difference between two
   sets of normal modes with matched frequencies (Tidor & Karplus 1994).
   Binding stiffens modes  ⇒  λ_complex > λ_free  ⇒  ln < 0  ⇒  ΔS < 0.

Why mode matching (and not a difference of absolute entropies)
--------------------------------------------------------------
S_complex from a sum of  k_B[1 − ln(ħω_k / k_BT)]  over ALL modes of the
complex is dominated by the softest 5–10 collective modes, and a bigger
elastic body always has softer low-frequency modes than a smaller one.
The "S_complex − S_A − S_B" difference therefore picks up a huge purely
geometric size bias (several × 10⁴ J/(mol·K)) that has nothing to do
with binding. Switching to Schlitter's QH formula does NOT fix this:
it only tames the stiff-mode tail, not the soft-mode growth with
system size.

Mode matching avoids the issue entirely. Each term in the sum compares
the SAME motion in two states, so size effects cancel by construction,
and the 6 inter-chain breathing modes of the complex are naturally
unmatched (no free-chain mode resembles them) and excluded — they
belong to the trans/rot cage term anyway.

Expected magnitudes
-------------------
For typical protein–protein complexes:
    deltaS ≈ −50 to −500 J/(mol·K)
giving −TΔS ≈ +0.4 to +4 kcal/mol at 300 K.
Compare Tidor & Karplus 1994 on insulin dimerization: ~ +10 kcal/mol.

References
----------
Tirion (1996), Phys. Rev. Lett. 77:1905-1908.
Bahar, Atilgan & Erman (1997), Folding & Design 2:173-181.
Atilgan et al. (2001), Biophys. J. 80:505-515.
Tidor & Karplus (1994), J. Mol. Biol. 238:405-414.  DOI:10.1006/jmbi.1994.1300
    Methodology paper: NMA binding entropy for protein–protein
    association, with the log-ratio formula used here.
Bakan, Meireles & Bahar (2011), Bioinformatics 27:1575-1577.  (ProDy)
"""

from __future__ import annotations
import os
import warnings
from dataclasses import dataclass
from typing import Tuple, Optional

import numpy as np
from scipy.optimize import linear_sum_assignment

from .entropy_utils import split_chains_to_tempfiles


# -----------------------------------------------------------------------------
# Physical constants
# -----------------------------------------------------------------------------
_NA   = 6.02214076e23   # Avogadro's number, 1/mol
_KB_J = 1.380649e-23    # Boltzmann constant, J/K
_R    = _KB_J * _NA     # universal gas constant, J/(mol·K) ≈ 8.314

# Minimum overlap (cosine) for a matched pair to be trusted. A free-chain
# mode whose best match has lower overlap than this is discarded: most
# likely its bound-state counterpart re-organised significantly and the
# log-ratio for that pair would be meaningless noise.
_OVERLAP_THRESHOLD = 0.30


# -----------------------------------------------------------------------------
# Result container — same minimalist style as SidechainResult
# -----------------------------------------------------------------------------

@dataclass
class BackboneResult:
    """Result of the NMA-based backbone configurational entropy estimate."""
    dS_backbone: float = 0.0
    # ΔS_binding in J/(mol·K). Typically negative (binding stiffens modes).
    n_modes_matched: int = 0
    # Number of (free-chain, complex) mode pairs that passed the overlap
    # quality filter and contributed to the sum. Diagnostic.


# -----------------------------------------------------------------------------
# ANM build helper
# -----------------------------------------------------------------------------

def _anm_build(
    pdb_path: str,
    cutoff: float,
    gamma: float,
) -> Tuple[Optional[object], np.ndarray, np.ndarray]:
    """
    Build ANM on Cα atoms of `pdb_path` and return:
        (atoms, eigenvalues, eigenvectors)

    Sorted by ascending eigenvalue.
    Eigenvectors have shape (3·N_atoms, n_modes); columns are unit vectors.
    The 6 rigid-body zero modes are excluded (ProDy `zeros=False`).
    """
    import prody
    prody.confProDy(verbosity="none")

    atoms = prody.parsePDB(pdb_path).select("calpha")
    if atoms is None or len(atoms) < 4:
        warnings.warn(f"Too few Cα atoms in {pdb_path}; skipping ANM.")
        return None, np.array([]), np.zeros((0, 0))

    anm = prody.ANM("nma")
    anm.buildHessian(atoms, cutoff=cutoff, gamma=gamma)
    n_modes = 3 * len(atoms) - 6
    anm.calcModes(n_modes=n_modes, zeros=False)

    eigvals = np.asarray(anm.getEigvals(), dtype=float)
    eigvecs = np.asarray(anm.getEigvecs(), dtype=float)  # (3N, n_modes)
    order = np.argsort(eigvals)
    return atoms, eigvals[order], eigvecs[:, order]


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
    Estimate ΔS_backbone in J/(mol·K) via mode-matched ANM.

    Parameters
    ----------
    complex_pdb : str
        Path to a PDB with both chains in the bound pose.
    chain_a, chain_b : str
        Single-character chain IDs.
    T : float
        Temperature in K. Cancels in the log-ratio formula, accepted
        only for API consistency with the other modules.
    cutoff : float
        ANM contact cutoff in Å (default 8.0, the Atilgan 2001 standard).
    gamma : float
        Uniform spring constant in kcal/mol/Å² (default 1.0). Numerically
        cancels in the eigenvalue ratio; kept for transparency.

    Returns
    -------
    BackboneResult
        ``deltaS``  in J/(mol·K). Typically −50 to −500 for a PPI.
    """

    # --- 1. ANM for each free chain ------------------------------------
    path_a, path_b = split_chains_to_tempfiles(complex_pdb, chain_a, chain_b)
    try:
        atoms_a, eigvals_a, eigvecs_a = _anm_build(path_a, cutoff, gamma)
        atoms_b, eigvals_b, eigvecs_b = _anm_build(path_b, cutoff, gamma)
    finally:
        os.unlink(path_a)
        os.unlink(path_b)

    # --- 2. ANM for the complex ----------------------------------------
    atoms_c, eigvals_c, eigvecs_c = _anm_build(complex_pdb, cutoff, gamma)

    if (atoms_c is None or eigvals_c.size == 0
            or atoms_a is None or eigvals_a.size == 0
            or atoms_b is None or eigvals_b.size == 0):
        return BackboneResult(dS_backbone=0.0, n_modes_matched=0)

    # --- 3. Identify chain-A and chain-B rows of the complex eigenvectors
    # getChids() returns a numpy array of chain IDs, one per Cα.
    chain_ids = atoms_c.getChids()
    idx_a = np.where(chain_ids == chain_a)[0]
    idx_b = np.where(chain_ids == chain_b)[0]

    # Sanity: the chain-by-chain ANMs must have the same number of Cα atoms
    # as the same chains do inside the complex. If a residue is parsed in
    # one but not the other (insertion code, alt-loc, etc.) the matching
    # breaks down. Bail out cleanly rather than producing garbage.
    if len(idx_a) != len(atoms_a) or len(idx_b) != len(atoms_b):
        warnings.warn(
            "Cα atom count mismatch between isolated chains and complex; "
            f"complex_A={len(idx_a)} vs free_A={len(atoms_a)}, "
            f"complex_B={len(idx_b)} vs free_B={len(atoms_b)}. "
            "Backbone ΔS skipped."
        )
        return BackboneResult(dS_backbone=0.0, n_modes_matched=0)

    # Each Cα contributes 3 rows (x, y, z) to the eigenvector matrix.
    rows_a = np.concatenate([[3*i, 3*i+1, 3*i+2] for i in idx_a])
    rows_b = np.concatenate([[3*i, 3*i+1, 3*i+2] for i in idx_b])

    # Chain blocks of the complex eigenvectors. Each column is the chain's
    # contribution to one complex mode. NOT renormalised: a complex mode
    # that is dominated by the OTHER chain will have small norm here and
    # will naturally lose the overlap competition for this chain's free
    # modes — which is exactly the behaviour we want.
    V_c_a = eigvecs_c[rows_a, :]   # (3·N_a, n_complex)
    V_c_b = eigvecs_c[rows_b, :]   # (3·N_b, n_complex)

    # --- 4. Combined overlap matrix (chain A modes stacked above chain B)
    # Row j < n_a:  free chain-A mode j, overlap measured against the
    #               chain-A subspace of complex modes.
    # Row j ≥ n_a:  free chain-B mode (j − n_a), overlap measured against
    #               the chain-B subspace of complex modes.
    n_a = eigvals_a.size
    n_b = eigvals_b.size
    n_c = eigvals_c.size

    overlap = np.zeros((n_a + n_b, n_c))
    overlap[:n_a, :] = np.abs(eigvecs_a.T @ V_c_a)
    overlap[n_a:, :] = np.abs(eigvecs_b.T @ V_c_b)

    # --- 5. Hungarian 1-to-1 assignment over the combined matrix --------
    # `linear_sum_assignment` minimises the cost, so feed it `-overlap`.
    # For a rectangular matrix it pairs every row with a distinct column
    # (free modes < complex modes, so all free modes get matched and the
    # 6 inter-chain breathing modes of the complex are left unpaired).
    row_idx, col_idx = linear_sum_assignment(-overlap)
    pair_overlap = overlap[row_idx, col_idx]

    # --- 6. Filter low-confidence matches --------------------------------
    keep = pair_overlap > _OVERLAP_THRESHOLD
    if keep.sum() == 0:
        return BackboneResult(dS_backbone=0.0, n_modes_matched=0)

    eigvals_free_concat = np.concatenate([eigvals_a, eigvals_b])
    lambda_free    = eigvals_free_concat[row_idx][keep]
    lambda_complex = eigvals_c[col_idx][keep]

    # Eigenvalues should all be > 0 (rigid-body modes already dropped).
    # Belt-and-braces guard against numerical zeros.
    safe = (lambda_free > 1e-12) & (lambda_complex > 1e-12)
    if safe.sum() == 0:
        return BackboneResult(dS_backbone=0.0, n_modes_matched=0)

    # --- 7. Tidor–Karplus log-ratio formula ------------------------------
    log_ratio = np.log(lambda_free[safe] / lambda_complex[safe])
    deltaS = float(0.5 * _R * log_ratio.sum())

    return BackboneResult(
        dS_backbone=deltaS,
        n_modes_matched=int(safe.sum()),
    )
