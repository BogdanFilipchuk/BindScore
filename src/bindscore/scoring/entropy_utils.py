"""
entropy_utils.py
Function with utilities for other entropy modules
========

Shared helpers used by the other entropy modules:
  - Parsing a complex PDB into per-chain temporary files
  - Computing phi/psi backbone dihedrals

We use Biopython for structural parsing because it is the de-facto standard
in the Python structural-biology ecosystem and its `Bio.PDB` module is
stable and well-documented.

References for the methods used here:
  - Hubbard & Thornton, NACCESS, Department of Biochemistry, University
    College London (1993). Defines the SASA-burial threshold for interface
    residue identification (Delta_SASA > 1 A^2).
  - Lee & Richards (1971) J. Mol. Biol. 55:379. Original SASA definition.
"""

from __future__ import annotations
import os
import tempfile
from typing import List, Tuple, Dict
import numpy as np

# Biopython 
from Bio.PDB import PDBParser, PDBIO, Select


# -----------------------------------------------------------------------------
# Chain splitting
# -----------------------------------------------------------------------------

class _ChainSelector(Select):
    """Biopython Select subclass that keeps only one chain ID."""
    def __init__(self, chain_id: str):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return 1 if chain.id == self.chain_id else 0


def split_chains_to_tempfiles(complex_pdb: str, chain_a: str, chain_b: str) -> Tuple[str, str]:
    """
    Write the two chains of a complex to separate PDB files in /tmp.

    Several downstream tools take file paths rather than in-memory objects like dicts or strings (such functions exist in pdb_utils_protein module),
    so we create each chain file on disk temporary.

    Returns
    -------
    (path_to_chain_a_pdb, path_to_chain_b_pdb)
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", complex_pdb)
    io = PDBIO()
    io.set_structure(structure)

    # Write chain A
    fd_a, path_a = tempfile.mkstemp(suffix=f"_chain{chain_a}.pdb")
    os.close(fd_a)
    io.save(path_a, _ChainSelector(chain_a))

    # Write chain B
    fd_b, path_b = tempfile.mkstemp(suffix=f"_chain{chain_b}.pdb")
    os.close(fd_b)
    io.save(path_b, _ChainSelector(chain_b))

    return path_a, path_b


# -----------------------------------------------------------------------------
# Backbone dihedral angles (phi/psi)
# -----------------------------------------------------------------------------

def _dihedral(p0, p1, p2, p3) -> float:
    """
    Compute a dihedral angle (degrees) defined by four points.

    Uses the standard cross-product formulation; this is the same math the
    Bio.PDB.vectors module uses internally, but written out so we don't have
    a hidden dependency on a less-stable submodule API.
    """
    b0 = p1 - p0
    b1 = p2 - p1
    b2 = p3 - p2

    # Normalise b1 so projections work cleanly
    b1 /= np.linalg.norm(b1)

    # Project b0 and b2 onto a plane perpendicular to b1
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))


def get_phi_psi_for_chain(chain) -> List[Tuple[int, str, float, float]]:
    """
    Compute (phi, psi) backbone dihedrals for every residue in a Biopython
    chain object.

    The first and last residue of a chain have undefined phi/psi (no
    preceding C or no following N atom respectively), so they are returned
    with NaN values and downstream code must handle this.

    Returns a list of tuples: (resseq, resname, phi_deg, psi_deg)
    """
    residues = [r for r in chain if r.id[0] == " "]  # standard residues only
    out: List[Tuple[int, str, float, float]] = []

    for i, res in enumerate(residues):
        try:
            n  = res["N"].get_vector().get_array()
            ca = res["CA"].get_vector().get_array()
            c  = res["C"].get_vector().get_array()
        except KeyError:
            # Residue missing backbone atoms (rare; sometimes happens at
            # disordered termini in low-resolution structures)
            continue

        # phi uses C_prev - N - CA - C
        phi = float("nan")
        if i > 0:
            try:
                c_prev = residues[i - 1]["C"].get_vector().get_array()
                phi = _dihedral(c_prev, n, ca, c)
            except KeyError:
                pass

        # psi uses N - CA - C - N_next
        psi = float("nan")
        if i < len(residues) - 1:
            try:
                n_next = residues[i + 1]["N"].get_vector().get_array()
                psi = _dihedral(n, ca, c, n_next)
            except KeyError:
                pass

        out.append((res.id[1], res.get_resname(), phi, psi))

    return out


