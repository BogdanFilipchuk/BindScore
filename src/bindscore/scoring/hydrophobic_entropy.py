"""
==================================
Corrected implementation of the hydrophobic entropy calculation from
Singh & Warshel (2010), Section II.2.3 and Figure 4.

KEY CORRECTION FROM PREVIOUS VERSION
------------------------------------
The hydrophobic entropy is NOT computed by alchemically annihilating the
ligand. It is computed by running TWO separate restraint-release (RR)
calculations in each environment:

    ℓ'  =  nonpolar ligand  (charges zeroed, LJ kept)
    ℓ"  =  "nothing"        (charges zeroed AND LJ zeroed)

For each environment (water and protein), we measure:

    -T·ΔS_hphob = -T·ΔS'_RR(ℓ")  -  -T·ΔS'_RR(ℓ')       (Eq. 12)

This is the entropy difference between having a phantom particle vs. a
ligand-shaped cavity in that solvent. The hydrophobic effect IS that
difference. Then (Eq. 13):

    -T·ΔΔS^{w→p}_hphob = -T·ΔS^p_hphob  -  -T·ΔS^w_hphob

No alchemical annihilation step is needed — the cycle is closed by the
paper's assumption that -T·ΔS'_1 ≈ 0 (states with strong restraints have
no entropy regardless of their LJ parameters).

PREREQUISITES
-------------
You need TWO completely separate systems:
  - system_water    : ligand in a water box
  - system_protein  : ligand in the protein active site
Each has its own Topology, System, and (eventually) Context.

The function `run_rr_fep` from conformational_entropy_commented.py is
reused. It expects restraints to be already attached to the system.
"""

import openmm as mm
import openmm.app as app
import openmm.unit as unit
import numpy as np
import copy

# Reuse from the configurational entropy file
from conformational_entropy_commented import (
    add_cartesian_restraints,
    run_rr_fep,
    N_REF_SETS,
    K_INITIAL,
    RT_KCAL,
)


# =============================================================================
# HELPER 1 — Switch the ligand to its nonpolar form (ℓ')
# =============================================================================

def make_nonpolar_ligand(system, ligand_indices, context=None):
    """
    Sets ligand atomic charges to zero, keeps LJ parameters intact.
    This produces the ℓ' state (nonpolar ligand) of Figure 4.

    Parameters
    ----------
    system         : OpenMM System (modified in place)
    ligand_indices : list of int, ligand atom indices
    context        : OpenMM Context, OPTIONAL. If a context already exists
                     for this system, we MUST call updateParametersInContext
                     to make the changes take effect. If you call this
                     function BEFORE creating the simulation/context, leave
                     it as None.

    Returns
    -------
    The same system (for convenience in chaining).
    """
    # Locate the NonbondedForce (there should be exactly one)
    nb_force = next(f for f in system.getForces()
                    if isinstance(f, mm.NonbondedForce))

    # Zero only the charges; sigma and epsilon are preserved so the
    # ligand still has its shape and excluded-volume effect (the cavity).
    for idx in ligand_indices:
        charge, sigma, epsilon = nb_force.getParticleParameters(idx)
        nb_force.setParticleParameters(idx, 0.0 * unit.elementary_charge,
                                            sigma, epsilon)

    # If a Context already exists, propagate changes to the GPU/CPU buffers.
    # Without this call, OpenMM silently ignores the parameter changes.
    if context is not None:
        nb_force.updateParametersInContext(context)

    return system


# =============================================================================
# HELPER 2 — Switch the ligand to its "nothing" form (ℓ")
# =============================================================================

def make_nothing_ligand(system, ligand_indices, context=None):
    """
    Sets ligand atomic charges AND LJ epsilon to zero. The atoms remain in
    the topology (so harmonic restraints can still act on them) but they
    no longer interact with the environment. This is the ℓ" state of
    Figure 4 — a "phantom" or "nothing" particle.

    NOTE on overlap: with epsilon = 0, the ligand atoms can occupy the same
    space as solvent atoms. This does NOT cause energy explosions because
    there is no interaction to compute. The atoms are inert tracking points
    that the restraint potential anchors to.

    Parameters and return value: same as make_nonpolar_ligand.
    """
    nb_force = next(f for f in system.getForces()
                    if isinstance(f, mm.NonbondedForce))

    for idx in ligand_indices:
        charge, sigma, epsilon = nb_force.getParticleParameters(idx)
        # Zero out charge AND epsilon. Keep sigma at original value
        # (it doesn't matter because epsilon = 0 makes the whole LJ term zero,
        # but some integrators dislike sigma = 0).
        nb_force.setParticleParameters(idx,
                                       0.0 * unit.elementary_charge,
                                       sigma,
                                       0.0 * unit.kilojoule_per_mole)

    if context is not None:
        nb_force.updateParametersInContext(context)

    return system


# =============================================================================
# HELPER 3 — Run the RR calculation with variational minimisation
# =============================================================================

def rr_with_minimisation(system, topology, positions, ligand_indices,
                          ref_coords_list=None,
                          n_ref_sets=N_REF_SETS):
    """
    Wrapper that runs the RR-FEP ladder N_REF_SETS times with different
    reference coordinates, then returns the minimum |ΔG_RR|. This is the
    same variational trick used in the configurational entropy calculation:
    the minimum across reference sets approximates the value where
    enthalpic contamination is zero, leaving pure -T·ΔS.

    Parameters
    ----------
    system          : OpenMM System (NO restraints attached yet — we add
                       them inside the loop). The ligand parameters should
                       already be set to ℓ' or ℓ" before calling this.
    topology        : OpenMM Topology
    positions       : initial positions (OpenMM Quantity)
    ligand_indices  : list of int
    ref_coords_list : OPTIONAL list of N_REF_SETS position arrays, one per
                       reference set. If None, the same `positions` is
                       reused for all sets (not ideal — see warning below).
    n_ref_sets      : int, number of reference sets

    Returns
    -------
    float : min |ΔG_RR| in kcal/mol, our best estimate of -T·ΔS_RR
    """
    # In production you should generate N_REF_SETS distinct snapshots by
    # running a short equilibration MD with a stiff restraint and saving
    # frames every few ps. Reusing the same positions gives the same R̄
    # every iteration and defeats the variational trick.
    if ref_coords_list is None:
        print("    WARNING: no distinct reference snapshots provided; "
              "reusing the same positions — variational minimisation "
              "will not work properly.")
        ref_coords_list = [positions] * n_ref_sets

    dG_list = []
    for i, ref_coords in enumerate(ref_coords_list):
        # Each call to add_cartesian_restraints bakes new R̄ anchors into
        # the system, so we need a fresh copy of the system each time.
        # Deepcopy is expensive but conceptually cleanest — for production,
        # rebuild the system from the force field each iteration instead.
        sys_copy = copy.deepcopy(system)
        add_cartesian_restraints(sys_copy, ref_coords, ligand_indices,
                                  K_kcal=K_INITIAL)

        dG = run_rr_fep(sys_copy, topology, ref_coords, ligand_indices)
        dG_list.append(dG)
        print(f"    Set {i+1}/{n_ref_sets}: ΔG_RR = {dG:.3f} kcal/mol")

    # Variational minimisation: take the run with smallest |ΔG_RR|
    return min(dG_list, key=abs)


# =============================================================================
# MAIN — Hydrophobic entropy of binding
# =============================================================================

def compute_hydrophobic_entropy(system_water, topology_water, positions_water,
                                  system_protein, topology_protein, positions_protein,
                                  ligand_indices_water,
                                  ligand_indices_protein,
                                  ref_coords_water_lprime=None,
                                  ref_coords_water_nothing=None,
                                  ref_coords_protein_lprime=None,
                                  ref_coords_protein_nothing=None):
    """
    Computes the hydrophobic contribution to the binding entropy:

        -T·ΔΔS^{w→p}_hphob = -T·ΔS^p_hphob  -  -T·ΔS^w_hphob   (Eq. 13)

    where each environment-specific term is (Eq. 12):

        -T·ΔS_hphob = -T·ΔS'_RR(ℓ")  -  -T·ΔS'_RR(ℓ')

    Note that ligand atom indices may differ between the water and protein
    systems because the topologies are different (different residues,
    different atom ordering). Pass both index lists explicitly.

    Parameters
    ----------
    system_water, topology_water, positions_water
        The ligand-in-water system.
    system_protein, topology_protein, positions_protein
        The ligand-in-protein system.
    ligand_indices_water, ligand_indices_protein
        Ligand atom indices in each system.
    ref_coords_*  (optional)
        Pre-generated reference coordinate snapshots for each of the four
        sub-calculations. Strongly recommended for proper variational
        minimisation.

    Returns
    -------
    float : -T·ΔΔS_hphob^{w→p} in kcal/mol
    """

    # ─────────────────────────────────────────────────────────────
    # PART A — In water
    # ─────────────────────────────────────────────────────────────
    print("\n=== HYDROPHOBIC ENTROPY: WATER ENVIRONMENT ===")

    # Step A1: nonpolar ligand ℓ' in water
    # Make a deep copy of the system so we don't mutate the caller's object
    sys_w_lprime = copy.deepcopy(system_water)
    make_nonpolar_ligand(sys_w_lprime, ligand_indices_water)

    print("  Running RR for ℓ' (nonpolar ligand) in water...")
    minus_T_dS_lprime_w = rr_with_minimisation(
        sys_w_lprime, topology_water, positions_water,
        ligand_indices_water,
        ref_coords_list=ref_coords_water_lprime,
    )

    # Step A2: "nothing" ℓ" in water
    sys_w_nothing = copy.deepcopy(system_water)
    make_nothing_ligand(sys_w_nothing, ligand_indices_water)

    print("  Running RR for ℓ\" (nothing) in water...")
    minus_T_dS_nothing_w = rr_with_minimisation(
        sys_w_nothing, topology_water, positions_water,
        ligand_indices_water,
        ref_coords_list=ref_coords_water_nothing,
    )

    # Equation 12: -T·ΔS^w_hphob = -T·ΔS'_ℓ",w  -  -T·ΔS'_ℓ',w
    minus_T_dS_hphob_w = minus_T_dS_nothing_w - minus_T_dS_lprime_w
    print(f"  -T·ΔS_hphob (water)   = {minus_T_dS_hphob_w:.3f} kcal/mol")

    # ─────────────────────────────────────────────────────────────
    # PART B — In protein
    # ─────────────────────────────────────────────────────────────
    print("\n=== HYDROPHOBIC ENTROPY: PROTEIN ENVIRONMENT ===")

    # Step B1: nonpolar ligand ℓ' in protein
    sys_p_lprime = copy.deepcopy(system_protein)
    make_nonpolar_ligand(sys_p_lprime, ligand_indices_protein)

    print("  Running RR for ℓ' (nonpolar ligand) in protein...")
    minus_T_dS_lprime_p = rr_with_minimisation(
        sys_p_lprime, topology_protein, positions_protein,
        ligand_indices_protein,
        ref_coords_list=ref_coords_protein_lprime,
    )

    # Step B2: "nothing" ℓ" in protein
    sys_p_nothing = copy.deepcopy(system_protein)
    make_nothing_ligand(sys_p_nothing, ligand_indices_protein)

    print("  Running RR for ℓ\" (nothing) in protein...")
    minus_T_dS_nothing_p = rr_with_minimisation(
        sys_p_nothing, topology_protein, positions_protein,
        ligand_indices_protein,
        ref_coords_list=ref_coords_protein_nothing,
    )

    # Equation 12 (analogue in protein):
    # -T·ΔS^p_hphob = -T·ΔS'_ℓ",p  -  -T·ΔS'_ℓ',p
    minus_T_dS_hphob_p = minus_T_dS_nothing_p - minus_T_dS_lprime_p
    print(f"  -T·ΔS_hphob (protein) = {minus_T_dS_hphob_p:.3f} kcal/mol")

    # ─────────────────────────────────────────────────────────────
    # PART C — Combine into binding contribution (Eq. 13)
    # ─────────────────────────────────────────────────────────────
    minus_T_ddS_hphob = minus_T_dS_hphob_p - minus_T_dS_hphob_w

    print("\n=== HYDROPHOBIC BINDING ENTROPY ===")
    print(f"  -T·ΔΔS_hphob (w→p) = {minus_T_ddS_hphob:.3f} kcal/mol")
    print(f"  (Compare to Table 6: T4-lysozyme expects -7.68 to +5.10)")

    return minus_T_ddS_hphob


# =============================================================================
# SANITY CHECK NOTES
# =============================================================================
# When you run this, expect (based on Table 6 of the paper):
#
#   T4 lysozyme/benzene:
#     -T·ΔS^w_hphob = -7.68  (water gives up entropy when ligand cavity forms)
#     -T·ΔS^p_hphob =  5.10  (protein structures water around empty cavity)
#
# The fact that the two terms have OPPOSITE signs is one of the paper's
# key findings — it's the "interesting compensation" they highlight in
# the Results section. A common mistake is to expect the hydrophobic
# contribution to be uniformly favourable; in fact it depends strongly
# on what the cavity looks like and how water orders around it.
#
# Two things that will go wrong if your implementation is buggy:
#   1. Both terms come out with the SAME sign → you've probably mixed up
#      ℓ' and ℓ" somewhere, or the LJ parameters aren't actually being
#      zeroed in the "nothing" state.
#   2. The "nothing" state gives a hugely larger |ΔG_RR| than ℓ' → the
#      phantom particle is being kicked around by leftover interactions.
#      Check that NonbondedForce.updateParametersInContext was called.
