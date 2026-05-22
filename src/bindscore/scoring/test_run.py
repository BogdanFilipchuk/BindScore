import openmm as mm
import openmm.app as app
import openmm.unit as unit
import numpy as np
from conformational_entropy import compute_configurational_entropy

# --- Minimal system: 3 "carbon" atoms ---
system = mm.System()
for _ in range(3):
    system.addParticle(12.0)  # mass in amu

# Minimal topology
topology = app.Topology()
chain = topology.addChain()
residue = topology.addResidue("LIG", chain)
for i in range(3):
    topology.addAtom(f"C{i}", app.Element.getBySymbol("C"), residue)

# Positions: atoms spaced 0.15 nm apart
positions_bound = np.array([[0.0, 0.0, 0.0],
                             [0.15, 0.0, 0.0],
                             [0.30, 0.0, 0.0]]) * unit.nanometers

# Unbound: slightly spread out (more freedom)
positions_unbound = np.array([[0.0, 0.0, 0.0],
                               [0.20, 0.0, 0.0],
                               [0.40, 0.0, 0.0]]) * unit.nanometers

ligand_indices = [0, 1, 2]

result = compute_configurational_entropy(
    system, topology, positions_bound, positions_unbound, ligand_indices
)
print(f"-T*ΔS_conf = {result:.3f} kcal/mol")
