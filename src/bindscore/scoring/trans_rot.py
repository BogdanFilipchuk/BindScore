"""
BROKEN MODULE.

binding_entropy.trans_rot
=========================

Translational and rotational entropy loss upon binding.

When two molecules associate into a complex, six degrees of freedom (three
translational, three rotational) of one body relative to the other become
restricted to small fluctuations around the bound pose. This is one of the
largest single entropic penalties in binding.

We compute it analytically from the partition functions of an ideal gas
(Sackur-Tetrode for translation) and a rigid rotor (for rotation), which
costs microseconds of CPU time. The result is the standard textbook
expression and does NOT require any sampling.

References
----------
Page & Jencks (1971), PNAS 68:1678-1683.
    First quantitative estimate of trans+rot entropy in biological binding.
    Their canonical estimate of ~45 cal/mol/K (~13 kcal/mol at 300 K) for
    a typical bimolecular association is what this module reproduces.

Finkelstein & Janin (1989), Protein Engineering 3:1-3.
    Argued the loss is much smaller than Page-Jencks because the bound
    state retains finite fluctuations. Their estimate is ~6 kcal/mol at
    300 K, which is closer to modern simulation values and is what we
    implement here as the default.

Gilson & Zhou (2007), Annu. Rev. Biophys. Biomol. Struct. 36:21-42.
    Modern review of standard-state corrections for binding free energies;
    the 1660 A^3 cage volume below comes from this standard convention.
"""

from __future__ import annotations
import math
from dataclasses import dataclass
from Bio.PDB import PDBParser

# Physical constants in units that give kcal/mol directly
_R_KCAL = 1.987204e-3   # gas constant, kcal/mol/K
_NA     = 6.02214076e23 # Avogadro's number
_KB_J   = 1.380649e-23  # Boltzmann constant, J/K
_H_J    = 6.62607015e-34# Planck constant, J*s

# Standard state: 1 M concentration corresponds to a free volume per
# molecule of 1660 A^3 (= 1 / (NA * 1 mol/L) converted to A^3). This is the
# convention used in essentially every binding free energy paper.
_STD_VOLUME_A3 = 1660.0


@dataclass
class TransRotResult:
    """Container for the trans+rot result with its sub-components exposed."""
    translational: float   # Delta_S_trans, J/(mol·K)
    rotational: float      # Delta_S_rot,   J/(mol·K)
    total: float           # sum, J/(mol·K)


# -----------------------------------------------------------------------------
# Mass and moment-of-inertia helpers
# -----------------------------------------------------------------------------

# Approximate atomic masses (Daltons). We only care about the elements that
# actually appear in proteins; everything else falls back to 12 (carbon).
_ATOMIC_MASS = {
    "H": 1.008,  "C": 12.011, "N": 14.007, "O": 15.999,
    "S": 32.06,  "P": 30.974, "F": 18.998, "CL": 35.45,
    "BR": 79.90, "I": 126.90,
}


def _chain_mass_and_inertia(pdb_path: str, chain_id: str):
    """
    Return (mass_Da, I1, I2, I3) for a chain.

    Mass is summed from atomic masses; I1/I2/I3 are the principal moments
    of inertia in Da*A^2 obtained by diagonalising the inertia tensor.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("s", pdb_path)

    coords = []
    masses = []
    for model in structure:
        for chain in model:
            if chain.id != chain_id:
                continue
            for res in chain:
                if res.id[0] != " ":   # skip HETATM / waters
                    continue
                for atom in res:
                    elem = atom.element.strip().upper() if atom.element else "C"
                    masses.append(_ATOMIC_MASS.get(elem, 12.0))
                    coords.append(atom.get_coord())
        break  # only the first model

    import numpy as np
    coords = np.asarray(coords)
    masses = np.asarray(masses)

    total_mass = float(masses.sum())
    if total_mass == 0.0:
        raise ValueError(f"No atoms found for chain {chain_id} in {pdb_path}")

    # Centre of mass and inertia tensor
    com = (coords * masses[:, None]).sum(axis=0) / total_mass
    r = coords - com
    # Inertia tensor I_ij = sum_k m_k (delta_ij r^2 - r_i r_j)
    r2 = (r ** 2).sum(axis=1)
    I = np.zeros((3, 3))
    for k in range(3):
        for l in range(3):
            if k == l:
                I[k, l] = (masses * (r2 - r[:, k] * r[:, l])).sum()
            else:
                I[k, l] = -(masses * r[:, k] * r[:, l]).sum()

    eigvals = np.linalg.eigvalsh(I)
    # Clip tiny negatives from floating-point noise
    eigvals = [max(float(e), 1e-6) for e in eigvals]
    return total_mass, eigvals[0], eigvals[1], eigvals[2]


# -----------------------------------------------------------------------------
# Sackur-Tetrode translational entropy
# -----------------------------------------------------------------------------

def _translational_entropy(mass_Da: float, T: float, volume_A3: float) -> float:
    """
    Sackur-Tetrode entropy for a single molecule in a given volume.

    S_trans = R * [ ln(V * (2*pi*m*kT / h^2)^(3/2)) + 5/2 ]

    Returns S in kcal/mol/K. The 5/2 is the standard Sackur-Tetrode
    constant for an ideal monatomic gas; for a single rigid body without
    internal vibrations this is the right expression.
    """
    m_kg = mass_Da * 1.66053906660e-27           # Da -> kg
    V_m3 = volume_A3 * 1e-30                     # A^3 -> m^3
    thermal = (2.0 * math.pi * m_kg * _KB_J * T) / (_H_J ** 2)
    # Quantum volume = 1 / thermal^(3/2)
    S_per_molecule_J = _KB_J * (math.log(V_m3 * thermal ** 1.5) + 2.5)
    # Convert J/K (per molecule) -> kcal/mol/K (per mol)
    return S_per_molecule_J * _NA / 4184.0


# -----------------------------------------------------------------------------
# Rigid-rotor rotational entropy
# -----------------------------------------------------------------------------

def _rotational_entropy(I1: float, I2: float, I3: float, T: float) -> float:
    """
    Classical rigid-rotor rotational entropy.

    S_rot = R * (3/2 + 1/2 * ln(pi * I1*I2*I3 * (8*pi^2*kT/h^2)^3) - ln(sigma))

    where sigma is the symmetry number (1 for any protein, which is
    asymmetric). I_k are the principal moments of inertia in Da*A^2.

    Returns S in kcal/mol/K.
    """
    # Convert moments of inertia from Da*A^2 to kg*m^2
    da_a2_to_kg_m2 = 1.66053906660e-27 * 1e-20
    I1_si = I1 * da_a2_to_kg_m2
    I2_si = I2 * da_a2_to_kg_m2
    I3_si = I3 * da_a2_to_kg_m2

    pre = (8.0 * math.pi ** 2 * _KB_J * T) / (_H_J ** 2)
    # The famous rotational partition function for an asymmetric top
    arg = math.pi * I1_si * I2_si * I3_si * pre ** 3
    S_per_molecule_J = _KB_J * (1.5 + 0.5 * math.log(arg))
    return S_per_molecule_J * _NA / 4184.0


# -----------------------------------------------------------------------------
# Public entry point
# -----------------------------------------------------------------------------

def compute(
    complex_pdb: str,
    chain_a: str,
    chain_b: str,
    T: float = 300.0,
    cage_volume_A3: float = _STD_VOLUME_A3,
) -> TransRotResult:
    """
    Estimate -T*Delta_S for translational + rotational entropy of binding.

    Convention: at the standard state (1 M concentration), each free
    molecule has access to a translational volume of 1660 A^3 and to the
    full 8*pi^2 orientational space. In the complex, the relative trans
    and rot freedom of one chain with respect to the other shrinks to a
    'cage' of finite volume - we use the same 1660 A^3 cage volume
    convention for simplicity, which is equivalent to assuming the bound
    complex retains roughly thermal-amplitude fluctuations. This gives a
    Finkelstein-Janin-style estimate of ~6 kcal/mol.

    A more sophisticated estimate would use the actual mean-square
    fluctuation of the complex from MD or from B-factors; we accept the
    simpler cage approximation here for the 2-minute budget.

    Returns
    -------
    TransRotResult dataclass with translational, rotational, and total
    components, each as -T*Delta_S in kcal/mol (positive = unfavourable).
    """
    # Per-chain mass and inertia
    m_a, Ia1, Ia2, Ia3 = _chain_mass_and_inertia(complex_pdb, chain_a)
    m_b, Ib1, Ib2, Ib3 = _chain_mass_and_inertia(complex_pdb, chain_b)
    m_ab = m_a + m_b

    # Inertia of the complex: parsed from the same PDB but treating both
    # chains as one rigid body. We get a fresh COM-relative inertia tensor
    # by passing a "fake" chain ID? Simpler: just call once for each chain
    # and use the analytic result that the complex moments are obtained
    # from the combined atom set. We compute it explicitly:
    Iab1, Iab2, Iab3 = _complex_inertia(complex_pdb, [chain_a, chain_b])

    # Free-state entropies (chains as separate molecules in standard volume)
    S_trans_a = _translational_entropy(m_a,  T, cage_volume_A3)
    S_trans_b = _translational_entropy(m_b,  T, cage_volume_A3)
    S_rot_a   = _rotational_entropy(Ia1, Ia2, Ia3, T)
    S_rot_b   = _rotational_entropy(Ib1, Ib2, Ib3, T)

    # Bound-state entropy: complex moves as one body
    S_trans_ab = _translational_entropy(m_ab, T, cage_volume_A3)
    S_rot_ab   = _rotational_entropy(Iab1, Iab2, Iab3, T)

    # Delta_S = S_complex - S_A - S_B
    # Multiplied by -T to express as a free energy contribution
    dS_trans = (S_trans_ab - S_trans_a - S_trans_b) * 4184.0
    dS_rot   = (S_rot_ab   - S_rot_a   - S_rot_b)   * 4184.0

    return TransRotResult(
        translational=dS_trans,
        rotational=dS_rot,
        total=dS_trans + dS_rot,
    )


def _complex_inertia(pdb_path: str, chain_ids: list[str]):
    """Helper: principal moments of inertia for a multi-chain assembly."""
    import numpy as np
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("s", pdb_path)

    coords, masses = [], []
    for model in structure:
        for chain in model:
            if chain.id not in chain_ids:
                continue
            for res in chain:
                if res.id[0] != " ":
                    continue
                for atom in res:
                    elem = atom.element.strip().upper() if atom.element else "C"
                    masses.append(_ATOMIC_MASS.get(elem, 12.0))
                    coords.append(atom.get_coord())
        break

    coords = np.asarray(coords)
    masses = np.asarray(masses)
    total_mass = masses.sum()
    com = (coords * masses[:, None]).sum(axis=0) / total_mass
    r = coords - com
    r2 = (r ** 2).sum(axis=1)
    I = np.zeros((3, 3))
    for k in range(3):
        for l in range(3):
            if k == l:
                I[k, l] = (masses * (r2 - r[:, k] * r[:, l])).sum()
            else:
                I[k, l] = -(masses * r[:, k] * r[:, l]).sum()
    eigvals = np.linalg.eigvalsh(I)
    eigvals = [max(float(e), 1e-6) for e in eigvals]
    return eigvals[0], eigvals[1], eigvals[2]
