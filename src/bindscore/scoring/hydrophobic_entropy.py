### TEST. FULLY AI
# ── Core Python stdlib ──────────────────────────────────────────────
import os
import math
import warnings
from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict

# ── Numerical / Scientific ───────────────────────────────────────────
import numpy as np
import pandas as pd                     # for results tables
import scipy.optimize as opt            # for gamma calibration

# ── OpenMM ───────────────────────────────────────────────────────────
import openmm as mm
import openmm.app as app
import openmm.unit as unit

# ── FreeSASA ─────────────────────────────────────────────────────────
import freesasa

# ── MDTraj (trajectory handling + SASA cross-validation) ─────────────
import mdtraj as md

def compute_hydrophobic_entropy(system, topology, positions_water,
                                 positions_protein, ligand_indices):
    """
    Implements Eq. (12-13) and Fig. 4 from Singh & Warshel.
    
    Steps:
      1. Set ligand charges to zero → nonpolar ligand (ℓ')
      2. Annihilate ℓ' → "nothing" (ℓ") via softcore FEP
      3. Compare entropy in water vs. protein
    """

    def zero_ligand_charges(system, ligand_indices):
        """Create uncharged (nonpolar) ligand — the ℓ' state"""
        nb_force = [f for f in system.getForces()
                    if isinstance(f, mm.NonbondedForce)][0]
        for idx in ligand_indices:
            charge, sigma, epsilon = nb_force.getParticleParameters(idx)
            nb_force.setParticleParameters(idx, 0.0, sigma, epsilon)
        return system

    def annihilate_ligand(system, ligand_indices, n_windows=41, sim_time_ps=40):
        """
        Softcore FEP: shrink ℓ' → ℓ" (nothing)
        Uses λ-scaling of LJ epsilon for the nonpolar ligand atoms
        """
        # Add softcore custom nonbonded force
        softcore = mm.CustomNonbondedForce(
            "4*epsilon*lambda_sc*((sigma/r)^12 - (sigma/r)^6);"
            "sigma=0.5*(sigma1+sigma2); epsilon=sqrt(epsilon1*epsilon2)"
        )
        softcore.addGlobalParameter("lambda_sc", 1.0)
        softcore.addPerParticleParameter("sigma")
        softcore.addPerParticleParameter("epsilon")

        # Only apply softcore to ligand atoms
        softcore.addInteractionGroup(
            set(ligand_indices),
            set(range(system.getNumParticles())) - set(ligand_indices)
        )
        system.addForce(softcore)

        lambdas = np.linspace(1.0, 0.0, n_windows)
        dG_values = []

        integrator = mm.LangevinMiddleIntegrator(
            300*unit.kelvin, 1.0/unit.picosecond,
            0.001*unit.picoseconds
        )
        simulation = app.Simulation(topology, system, integrator)

        for lam_curr, lam_next in zip(lambdas[:-1], lambdas[1:]):
            simulation.context.setParameter("lambda_sc", lam_curr)
            simulation.step(int(sim_time_ps / 0.001))

            state_curr = simulation.context.getState(getEnergy=True)
            simulation.context.setParameter("lambda_sc", lam_next)
            state_next = simulation.context.getState(getEnergy=True)

            dE = (state_next.getPotentialEnergy() -
                  state_curr.getPotentialEnergy()
                 ).value_in_unit(unit.kilocalories_per_mole)
            beta = 1.0 / (0.001987 * 300)
            dG_values.append(-1/beta * np.log(np.exp(-beta * dE)))

        return sum(dG_values)  # ΔG' for ℓ' → ℓ"

    # --- In water ---
    sys_w = zero_ligand_charges(system, ligand_indices)
    # Entropy of ℓ' in water (release restraints on nonpolar ligand)
    dG_lprime_w  = run_rr_fep(sys_w, topology, positions_water, ligand_indices)
    # Entropy of ℓ" (nothing) in water
    dG_nothing_w = annihilate_ligand(sys_w, ligand_indices)
    minus_T_dS_hphob_w = -(dG_lprime_w - dG_nothing_w)  # Eq. 12

    # --- In protein ---
    sys_p = zero_ligand_charges(system, ligand_indices)
    dG_lprime_p  = run_rr_fep(sys_p, topology, positions_protein, ligand_indices)
    dG_nothing_p = annihilate_ligand(sys_p, ligand_indices)
    minus_T_dS_hphob_p = -(dG_lprime_p - dG_nothing_p)  # analogous in protein

    # Net hydrophobic entropy change (Eq. 13)
    minus_T_ddS_hphob = minus_T_dS_hphob_p - minus_T_dS_hphob_w

    return minus_T_ddS_hphob


# --- Final assembly (mirrors Table 6 in the paper) ---
def compute_total_binding_entropy(conf_entropy, polar_entropy, hphob_entropy,
                                   cage_correction):
    """
    -T*ΔS_bind = -T*ΔS_conf + (-T*ΔΔS_pol) + (-T*ΔΔS_hphob) + (-T*ΔS_cage)
    """
    total = conf_entropy + polar_entropy + hphob_entropy + cage_correction
    print(f"  Configurational:  {conf_entropy:.2f} kcal/mol")
    print(f"  Polar solvation:  {polar_entropy:.2f} kcal/mol")
    print(f"  Hydrophobic:      {hphob_entropy:.2f} kcal/mol")
    print(f"  Cage correction:  {cage_correction:.2f} kcal/mol")
    print(f"  TOTAL -T*ΔS_bind: {total:.2f} kcal/mol")
    return total