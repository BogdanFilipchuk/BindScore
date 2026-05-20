"""
Conformational entropy module
=====================================
This module implements the harmonic Restraint Release Free Energy perturbation (RR-FEP) approach to compute the configurational
entropy contribution to protein-ligand binding. This module is following the Singh & Warshel article approach.

    Singh & Warshel (2010), Proteins 78(7):1724-1735
    "A Comprehensive Examination of the Contributions to Binding Entropy of Protein-Ligand Complexes"

SUMMARY 
------------------
We want to determine how much conformational freedom is lost by the ligand chain when it moves from an unbound state into a bound state at the protein active site.
We measure lost degrees of freedom indirectly (easier) by measuring the free energy cost of removing those dof's.
The dd_cartesian_harmonic_restraint_force() function adds those harmonic restraints, attached to a set of given reference coordinates ( TODO : EXPLAIN HOW ). 
The intensity of this force is given by K (analogous to a spring), and releasing K to zero gives the taken degrees of freedom back. 
The free energy of that release process, ΔG_RR, is what we aim to compute.

Therefore, the pipeline is:
  1. Attach restrains to all ligand atoms
  2. Slowly release the  restrains to K → 0 in both bound and unbound environments environments
  3. The free energy cost of this release ≈ -T*ΔS 
  4. Difference between bound and unbound gives the bindin entropy conformational penalty

Here are the references to key equations referenced from thepaper:
  Eq. (2): U_rest = (K/2) * Σ_i (R_i - R̄_i)²      [the spring potential]
  Eq. (3): FEP estimator via Zwanzig formula
  Eq. (4): -T*ΔS_conf_bind = min(ΔG_RR^bound) - min(ΔG_RR^unbound)
  Eq. (8): adds cage correction for 1M standard state
------------------

"""
# Note on the units:
# OpenMM works in kJ/mol and nm.
# The reference Reference paper (Singh and Warshel) uses kcal/mol and Å.
# Conversions: 1 kcal = 4.184 kJ  |  1 nm = 10 Å  |  so 1 kcal/mol/Å² = 418.4 kJ/mol/nm²
"TODO : Set unit conversions ?"

import openmm as mm
import openmm.app as app
import openmm.unit as unit
import numpy as np
import math


# =============================================================================
# TUNABLE PARAMETERS — collected here for easy adjustment
# =============================================================================

# --- harmonic_restraint schedule ---
K_INITIAL      = 10.0    # kcal/mol/Å²  Starting (stiff) harmonic constraint constant.
                          # Must be large enough that the ligand is nearly frozen.
                          # Reference paper uses 10.0. Increasing → safer but slower FEP.

K_FINAL        = 0.003   # kcal/mol/Å²  Final (loose) sharmonic constraint pring constant.
                          # Should be small enough that the ligand moves freely.
                          # Reference Reference paper uses 0.003. Decreasing → more complete release
                          # but harder to converge.

N_WINDOWS      = 41      # Number of FEP lambda windows between K_INITIAL and K_FINAL.
                          # Reference paper uses 41. More windows → better overlap between
                          # adjacent states → more accurate ΔG, but longer runtime.

SIM_TIME_PS    = 40      # Simulation time per window in picoseconds.
                          # Reference paper uses 40 ps. Increasing → better sampling per window.

TIMESTEP_FS    = 1.0     # MD timestep in femtoseconds. Reference paper uses 1 fs.
                          # Can try 2 fs with HBonds constraints but 1 fs is safer.

TEMPERATURE_K  = 300     # Simulation temperature in Kelvin. Reference paper uses 300 K.

# --- Sampling of reference coordinates ---
N_REF_SETS     = 8       # Number of different reference coordinate sets (R̄ sets).
                          # Reference paper uses 8. More sets → better chance of finding the
                          # true minimum (which removes enthalpic contamination).
                          # Each set = one full FEP ladder, so this multiplies runtime.

SPACING_PS     = 5       # How far apart (in ps) to space the reference snapshots.
                          # Should be large enough that snapshots are decorrelated.

# --- Cage harmonic_restraint for standard state correction ---
K_CAGE         = 0.22    # kcal/mol/Å²  Cage spring constant used during simulation.
                          # Limits how far the ligand's centre-of-mass can drift.
                          # Reference paper value: 0.22 (corresponds to ~4Å movement in protein)

K0             = 0.026   # kcal/mol/Å²  Cage spring constant at molar volume (1660 Å³).
                          # This is the analytical reference for 1M standard state.
                          # Reference paper value: 0.026. Do not change unless you change v0.

V0             = 1660.0  # Å³  Molar volume = 1 litre / Avogadro's number.
                          # Fixed by the definition of 1M standard state.

# --- Physical constants ---
kB_KCAL        = 0.001987  # Boltzmann constant in kcal/mol/K
RT_KCAL        = kB_KCAL * TEMPERATURE_K  # thermal energy in kcal/mol at 300K ≈ 0.596


# =============================================================================
# FUNCTION 1 — ADD HARMONIC CARTESIAN harmonic_restraintS
# =============================================================================

def add_cartesian_harmonic_restraint_force(system, positions, atom_indices, K_kcal) -> mm.CustomExternalForce:   
    """
    Attaches a harmonic restrain to each ligand atom, anchored at its current position.( Equation (2) of the paper). This is an essential step in realising the RR-FEP method.

    R̄_i are the reference (anchor) coordinates and K is simillar to a spring constant (Large K → ligand is frozen. K → 0 → ligand is free.)

    Parameters
    ----------
    system       : OpenMM System object (will be modified in-place)
    positions    : current atom positions (OpenMM Quantity with units)
    atom_indices : list of int, indices of ligand atoms to restrain
    K_kcal       : float, spring constant in kcal/mol/Å²

    Returns
    -------
    harmonic_restraint_force    : the CustomExternalForce object

    IMPORTANT: The reference coordinates R̄_i are SET HERE from `positions`.
    If you want different reference coordinates (different R̄ sets), call this
    function again with different positions — each call bakes in new anchors.
    """

    # Unit conversion: kcal/mol/Å² → kJ/mol/nm²
    # Factor 4.184: kcal → kJ
    # Factor 100:   Å² → nm²  (1 nm = 10 Å, so 1/Å² = 100/nm²)
    K_kJ_nm = K_kcal * 418.4

    # Define the energy expression as a string.
    # OpenMM evaluates this per particle - so x,y,z are current coords (in nm),
    # x0,y0,z0 are the per-particle reference point coords, will be set below.
    # K is a global parameter so we can change it during the simulation without rebuilding the force object.  ## TO DO : add an ability to change K (!!!!)
    harmonic_restraint_force = mm.CustomExternalForce("0.5 * K * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")   #That's some voodoo magic here

    # K is global: one value shared by all restrained atoms.
    # We add it as a global parameter so simulation.context.setParameter("K", ...) can update it on-the-fly during the FEP ladder.
    harmonic_restraint_force.addGlobalParameter("K", K_kJ_nm)

    # x0, y0, z0 are anchor point coordinate of each particle with respect to whom the harmonic potential will be calculate.
    # These are set once ig and do not change during a run. Honestly do not know what that
    harmonic_restraint_force.addPerParticleParameter("x0")
    harmonic_restraint_force.addPerParticleParameter("y0")
    harmonic_restraint_force.addPerParticleParameter("z0")

    # Register each ligand atom with its anchor coordinates
    for idx in atom_indices:
        # Extract position in nanometers (OpenMM's internal length unit)
        pos = positions[idx].value_in_unit(unit.nanometers)
        # pos is [x, y, z] in nm — these become the fixed anchor R̄_i
        harmonic_restraint_force.addParticle(idx, pos)

    # Add the harmonic_restraint to the system 
    system.addForce(harmonic_restraint_force)

    return harmonic_restraint_force   # return so caller can call setParameter("K", ...) later


# FUNCTION 2 — RUNNING THE FEP LADDER (releasing K from large to small) (God help me)


def run_rr_fep(system, topology, positions, atom_indices, K_initial=K_INITIAL, K_final=K_FINAL, n_windows=N_WINDOWS, sim_time_ps=SIM_TIME_PS):
    """
    Runs the harmonic_restraint_force Release FEP: gradually reduces K from K_initial to
    K_final in `n_windows` steps, computing the free energy change at each
    step via the Zwanzig (exponential averaging) formula (Equation 3).

    Conceptually:
      - At each window m, we equilibrate with spring constant K_m
      - Then we compute the energy difference if we instantly switched to K_{m+1}
      - The Zwanzig formula converts this to a free energy difference ΔΔG
      - Summing all ΔΔG gives the total ΔG_RR for this reference set

    The total ΔG_RR ≈ -T*ΔS ONLY when evaluated at the optimal reference
    coordinates R̄ (the ones that minimise |ΔG_RR|). See Eq. (4).

    Parameters
    ----------
    system       : OpenMM System (should already have harmonic_restraint_forces added via
                   add_cartesian_harmonic_restraints before calling this)
    topology     : OpenMM Topology
    positions    : initial positions (OpenMM Quantity)
    atom_indices : list of int, ligand atom indices (for reference only here)
    K_initial    : float, starting spring constant in kcal/mol/Å²
    K_final      : float, ending spring constant in kcal/mol/Å²
    n_windows    : int, number of FEP windows
    sim_time_ps  : float, simulation time per window in ps

    Returns
    -------
    dG_RR : float, total free energy of harmonic_restraint_force release in kcal/mol
            This equals -T*ΔS when computed at the optimal R̄.
    """

    # Build the K schedule: n_windows values from K_initial down to K_final.
    # Each consecutive pair (K_m, K_{m+1}) defines one FEP perturbation step.
    # NOTE: the Reference paper uses 4 separate "stages" with different spacing to ensure
    # good overlap. A simple linspace is a reasonable approximation; for
    # production runs consider log-spacing or the Reference paper's 4-stage schedule.
    K_values = np.linspace(K_initial, K_final, n_windows)

    # --- SET UP INTEGRATOR ---
    # Langevin dynamics at constant temperature.
    # friction coefficient 1/ps is a standard choice for aqueous systems.
    integrator = mm.LangevinMiddleIntegrator(
        TEMPERATURE_K * unit.kelvin,
        1.0 / unit.picosecond,           # friction coefficient
        TIMESTEP_FS * 0.001 * unit.picoseconds  # fs → ps conversion
    )

    # --- SET UP SIMULATION ---
    # Use CUDA if available for speed; fall back to CPU otherwise.
    # To force CPU: replace 'CUDA' with 'CPU'
    try:
        platform = mm.Platform.getPlatformByName('CUDA')
        simulation = app.Simulation(topology, system, integrator, platform)
    except Exception:
        platform = mm.Platform.getPlatformByName('CPU')
        simulation = app.Simulation(topology, system, integrator, platform)

    simulation.context.setPositions(positions)

    # Quick energy minimisation to remove any clashes before starting FEP
    simulation.minimizeEnergy()

    # Steps per window: sim_time_ps / timestep_ps
    steps_per_window = int(sim_time_ps / (TIMESTEP_FS * 0.001))

    # --- FEP LADDER ---
    dG_values = []  # will collect ΔΔG for each window

    for m, (K_curr, K_next) in enumerate(zip(K_values[:-1], K_values[1:])):

        # ── Step A: set current spring constant and equilibrate ──
        # Convert K from kcal/mol/Å² → kJ/mol/nm² for OpenMM
        simulation.context.setParameter("K", K_curr * 418.4)

        # Run MD at K_curr for sim_time_ps.
        # This equilibrates the system at the current harmonic_restraint_force strength.
        # The ligand explores configurations consistent with K_curr.
        simulation.step(steps_per_window)

        # ── Step B: compute energy at K_curr ──
        # This is U_m in Eq. (3): the potential energy with the current K
        state_curr = simulation.context.getState(getEnergy=True)
        E_curr = state_curr.getPotentialEnergy().value_in_unit(
                     unit.kilocalories_per_mole)

        # ── Step C: compute energy at K_next WITHOUT moving ──
        # We change only the spring constant (a global parameter),
        # keeping the coordinates frozen. This gives U_{m+1} at the same
        # configuration sampled under K_m. This is the "perturbation" in FEP.
        simulation.context.setParameter("K", K_next * 418.4)
        state_next = simulation.context.getState(getEnergy=True)
        E_next = state_next.getPotentialEnergy().value_in_unit(
                     unit.kilocalories_per_mole)

        # Restore K_curr so the next iteration's equilibration uses the right K
        simulation.context.setParameter("K", K_curr * 418.4)

        # ── Step D: Zwanzig formula (single-snapshot version) ──
        # Full Eq. (3) requires averaging exp(-β*ΔU) over many snapshots.
        # Here we use a single snapshot per window, which is the minimal
        # implementation. For production accuracy, collect ~100 snapshots
        # per window (run a sub-loop here) and average the exponentials.
        #
        # ΔΔG(m → m+1) = -kT * ln < exp(-ΔU/kT) >
        #
        # With one sample: ΔΔG ≈ -kT * ln( exp(-ΔU/kT) ) = ΔU
        # (reduces to simple energy difference — acceptable for small ΔK steps)
        delta_E = E_next - E_curr   # ΔU = U_{m+1} - U_m  (kcal/mol)
        beta    = 1.0 / RT_KCAL     # β = 1/kT  (mol/kcal)

        ddG = -RT_KCAL * np.log(np.exp(-beta * delta_E))
        # NOTE: with a single sample this simplifies to just delta_E.
        # Replace this line with an average over multiple snapshots for
        # production-quality calculations (see IMPROVEMENT NOTE below).

        dG_values.append(ddG)

    # Sum all window contributions → total ΔG_RR for this reference set
    dG_RR = float(np.sum(dG_values))

    return dG_RR


# =============================================================================
# FUNCTION 3 — MAIN CALCULATION: COMPUTE CONFIGURATIONAL BINDING ENTROPY
# =============================================================================

def compute_configurational_entropy(system, topology, positions_bound,positions_unbound, ligand_indices):
    """
    Implements Equations (4) and (8) from Singh & Warshel to compute:

        -T * ΔS_conf_bind = min(ΔG_RR^bound) - min(ΔG_RR^unbound) - T*ΔΔS_cage

    The three terms are:
      1. min(ΔG_RR^bound)   : cost of releasing springs in the protein
      2. min(ΔG_RR^unbound) : cost of releasing springs in water
                              (larger than bound because ligand is freer in water)
      3. T*ΔΔS_cage         : analytical correction to convert from the simulation
                              cage volume to the 1M standard state volume

    WHY WE TAKE THE MINIMUM (the variational trick):
      Every ΔG_RR value contains a mix of -T*ΔS (what we want) AND residual
      enthalpy (contamination from the specific reference coordinates chosen).
      At the optimal R̄, the enthalpy contribution goes to zero and
      ΔG_RR = -T*ΔS exactly. The minimum |ΔG_RR| across multiple R̄ sets
      is our best approximation to this optimum. This is the key insight of
      the Reference paper's variational minimisation (see Methods section II.2.1).

    Parameters
    ----------
    system            : OpenMM System (no harmonic_restraint_forces added yet)
    topology          : OpenMM Topology
    positions_bound   : ligand positions inside the protein (OpenMM Quantity)
    positions_unbound : ligand positions in water box (OpenMM Quantity)
    ligand_indices    : list of int, atom indices of the ligand

    Returns
    -------
    conf_entropy : float, -T*ΔS_conf_bind in kcal/mol
                  Positive value = entropy COST of binding (as expected,
                  since binding restricts the ligand's motion)
    """

    # ── STEP 1: Run N_REF_SETS FEP ladders for BOUND state ──
    # Each run uses a different reference coordinate set R̄ (snapshot from
    # equilibration). In practice you obtain these by running a short MD and
    # saving snapshots every SPACING_PS picoseconds (done externally, then
    # passed as a list of position arrays). Here we simulate this by running
    # the same starting positions N_REF_SETS times; for real use, pass in
    # N_REF_SETS distinct snapshots from a prior equilibration run.
    #
    # TUNABLE: N_REF_SETS (see top of file). More sets → better minimum estimate
    #           but N× longer runtime.

    dG_bound_list   = []
    dG_unbound_list = []

    for ref_set in range(N_REF_SETS):

        print(f"  Running reference set {ref_set + 1}/{N_REF_SETS}...")

        # --- BOUND STATE ---
        # Add springs anchored at the ligand's current position in the protein.
        # Each ref_set ideally uses a DIFFERENT snapshot from equilibration.
        # Passing the same positions every time gives the same R̄ (not ideal);
        # replace `positions_bound` with `snapshot_list_bound[ref_set]` in
        # production.
        harmonic_restraint_force_b = add_cartesian_harmonic_restraints(system, positions_bound, ligand_indices, K_kcal=K_INITIAL)

        # Run the FEP ladder: K goes from K_INITIAL → K_FINAL
        # Returns ΔG_RR for this reference set in kcal/mol
        dG_b = run_rr_fep(system, topology, positions_bound, ligand_indices)
        dG_bound_list.append(dG_b)

        # Remove the harmonic_restraint before adding a fresh one next iteration
        # (OpenMM doesn't let you modify forces after context creation,
        #  so we'd rebuild the system in a real implementation)
        # --- placeholder: in production, rebuild system each iteration ---

        # --- UNBOUND STATE (ligand in water) ---
        # Exactly the same procedure but with the ligand in a water box.
        # The ligand has MORE freedom here, so releasing the springs costs more,
        # giving a larger (more negative) ΔG_RR.
        harmonic_restraint_ub = add_cartesian_harmonic_restraints(
            system, positions_unbound, ligand_indices, K_kcal=K_INITIAL
        )
        dG_ub = run_rr_fep(system, topology, positions_unbound, ligand_indices)
        dG_unbound_list.append(dG_ub)

        print(f"    ΔG_RR bound   = {dG_b:.3f} kcal/mol")
        print(f"    ΔG_RR unbound = {dG_ub:.3f} kcal/mol")

    # ── STEP 2: Variational minimisation — take minimum |ΔG_RR| ──
    # We want the run where enthalpic contamination is smallest.
    # That is the run with the smallest absolute value of ΔG_RR,
    # NOT simply the most negative value. (Reference paper: "min" means min |ΔG|)
    min_dG_B  = min(dG_bound_list,   key=abs)   # Eq. (4): min(ΔG_RR^B)
    min_dG_UB = min(dG_unbound_list, key=abs)   # Eq. (4): min(ΔG_RR^UB)

    print(f"\n  Best bound   ΔG_RR = {min_dG_B:.3f} kcal/mol")
    print(f"  Best unbound ΔG_RR = {min_dG_UB:.3f} kcal/mol")

    # ── STEP 3: Cage correction — convert to 1M standard state ──
    # During simulation the ligand's centre-of-mass is confined to a small
    # "cage" volume v_cage (set by K_CAGE). We need to analytically add back
    # the entropy of expanding that cage to the molar volume v0 = 1660 Å³.
    #
    # Eq. (7): v_cage = (2π / (β * K_cage))^(3/2)
    # Eq. (6): -T*ΔS_cage = -kT * ln(v0 / v_cage)
    #
    # We do this separately for the cage used in the bound simulation (K_CAGE)
    # and the reference molar volume cage (K0), then take the difference.
    #
    # TUNABLE: K_CAGE and K0 (see top of file). Only change K_CAGE if you
    # used a different cage harmonic_restraint strength in your simulations.

    beta = 1.0 / RT_KCAL   # β = 1/kT in mol/kcal

    # Volume of the simulation cage (bound, using K_CAGE)
    # Units: (2π / (β * K)) is in Å³ when K is in kcal/mol/Å² and β in mol/kcal
    v_cage_B  = (2.0 * math.pi / (beta * K_CAGE)) ** 1.5   # Å³

    # Volume corresponding to molar volume cage (K0)
    v_cage_UB = (2.0 * math.pi / (beta * K0))    ** 1.5   # Å³

    # -T*ΔS_cage for bound state cage → molar volume
    minus_TdS_cage_B  = -RT_KCAL * math.log(V0 / v_cage_B)

    # -T*ΔS_cage for unbound reference cage → molar volume
    minus_TdS_cage_UB = -RT_KCAL * math.log(V0 / v_cage_UB)

    # Net cage correction: ΔΔS_cage (Eq. 8 correction term)
    T_delta_dS_cage = minus_TdS_cage_B - minus_TdS_cage_UB

    print(f"  Cage correction -T*ΔΔS_cage = {T_delta_dS_cage:.3f} kcal/mol")

    # ── STEP 4: Final answer — Equation (8) ──
    # -T*ΔS_conf_bind = min(ΔG_RR^B) - min(ΔG_RR^UB) - T*ΔΔS_cage
    #
    # Interpretation:
    #   min(ΔG_RR^B)   is small (protein already confines the ligand)
    #   min(ΔG_RR^UB)  is large (water gives lots of freedom to recover)
    #   Their difference is positive → binding costs entropy
    conf_entropy = min_dG_B - min_dG_UB - T_delta_dS_cage

    print(f"\n  -T*ΔS_conf_bind = {conf_entropy:.3f} kcal/mol")
    print(f"  (Positive = entropy cost of binding, as expected)")

    return conf_entropy


# =============================================================================
# IMPROVEMENT NOTE — Multi-snapshot Zwanzig averaging (production quality)
# =============================================================================
# In run_rr_fep(), the ΔΔG per window is currently computed from a SINGLE
# energy evaluation. For better accuracy, replace the single-snapshot block
# with a sub-loop that collects ~100 snapshots and averages:
#
#   delta_U_samples = []
#   for _ in range(100):                        # ← TUNABLE: samples per window
#       simulation.step(400)                    # short run between samples
#       E_curr = get_energy_at_K(K_curr)
#       E_next = get_energy_at_K(K_next)
#       delta_U_samples.append(E_next - E_curr)
#   delta_U = np.array(delta_U_samples)
#   ddG = -RT_KCAL * np.log(np.mean(np.exp(-delta_U / RT_KCAL)))
#
# This properly computes <exp(-βΔU)> and is much more accurate, especially
# for windows where K changes significantly.
#
# For even better accuracy, use the BAR (Bennett Acceptance Ratio) estimator
# instead of Zwanzig, collecting samples from both K_m and K_{m+1}.
