from openmm import CustomExternalForce
import numpy as np

def add_cartesian_restraints(system, positions, atom_indices, K_kcal):
    """
    Implements U_rest = (K/2) * sum_i (R_i - R_i_ref)^2
    Mirrors Eq. (2) in Singh & Warshel
    """
    K_kJ = K_kcal * 4.184  # convert kcal/mol/Å² → kJ/mol/nm²  (×4.184×100)
    K_kJ_nm = K_kcal * 418.4  

    restraint = CustomExternalForce(
        f"0.5 * K * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)"
    )
    restraint.addGlobalParameter("K", K_kJ_nm)
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")

    for idx in atom_indices:
        pos = positions[idx].value_in_unit(unit.nanometers)
        restraint.addParticle(idx, pos)

    system.addForce(restraint)
    return restraint


def run_rr_fep(system, topology, positions, atom_indices,
               K_initial=10.0, K_final=0.003,
               n_windows=41, sim_time_ps=40):
    """
    Restraint Release FEP: releases K from K_initial → K_final
    Mirrors the paper's RR-FEP protocol (41 windows, 40ps each)
    Returns -T*ΔS_RR in kcal/mol
    """
    # Lambda schedule: K interpolates from K_initial to K_final
    lambdas = np.linspace(0, 1, n_windows)
    K_values = K_initial * (1 - lambdas) + K_final * lambdas  

    integrator = mm.LangevinMiddleIntegrator(
        300*unit.kelvin,
        1.0/unit.picosecond,
        0.001*unit.picoseconds  # 1 fs timestep, matching paper
    )

    platform = mm.Platform.getPlatformByName('CUDA')
    simulation = app.Simulation(topology, system, integrator, platform)
    simulation.context.setPositions(positions)
    simulation.minimizeEnergy()

    dG_values = []
    steps = int(sim_time_ps / 0.001)  # 40ps at 1fs

    for i, (K_curr, K_next) in enumerate(zip(K_values[:-1], K_values[1:])):
        # Set current K
        simulation.context.setParameter("K", K_curr * 418.4)
        simulation.step(steps)

        # Collect energies for FEP (Eq. 3 in paper)
        state_curr = simulation.context.getState(getEnergy=True)
        E_curr = state_curr.getPotentialEnergy()

        simulation.context.setParameter("K", K_next * 418.4)
        state_next = simulation.context.getState(getEnergy=True)
        E_next = state_next.getPotentialEnergy()

        # ΔΔG'(m → m+1) = -β⁻¹ ln⟨exp(-β(U_{m+1} - U_m))⟩
        delta_E = (E_next - E_curr).value_in_unit(unit.kilocalories_per_mole)
        beta = 1.0 / (0.001987 * 300)  # 1/RT in kcal/mol
        dG_values.append(-1/beta * np.log(np.exp(-beta * delta_E)))

    return sum(dG_values)  # Total ΔG_RR ≈ -T*ΔS when at optimal restraint


def compute_configurational_entropy(system, topology, positions_bound,
                                     positions_unbound, ligand_indices):
    """
    Implements Eq. (4) and (8) from Singh & Warshel:
    -T*ΔS_conf = min(ΔG_RR^B) - min(ΔG_RR^UB) - T*ΔΔS_cage
    """
    # Run 8 sets with different restraint coordinates (paper protocol)
    n_sets = 8
    dG_bound_list = []
    dG_unbound_list = []

    for _ in range(n_sets):
        # Bound state
        sys_bound = add_cartesian_restraints(
            system, positions_bound, ligand_indices, K_kcal=10.0
        )
        dG_b = run_rr_fep(sys_bound, topology, positions_bound, ligand_indices)
        dG_bound_list.append(dG_b)

        # Unbound state (ligand in water)
        sys_unbound = add_cartesian_restraints(
            system, positions_unbound, ligand_indices, K_kcal=10.0
        )
        dG_ub = run_rr_fep(sys_unbound, topology, positions_unbound, ligand_indices)
        dG_unbound_list.append(dG_ub)

    # Take minimum |ΔG| values (variational minimization, Eq. 4)
    min_dG_B  = min(dG_bound_list,   key=abs)
    min_dG_UB = min(dG_unbound_list, key=abs)

    # Cage correction (Eq. 6-8): analytical term for 1M standard state
    import math
    K_cage = 0.22    # kcal/mol/Å²
    K0     = 0.026   # corresponds to molar volume (1660 Å³)
    kT     = 0.001987 * 300
    v_cage_B  = (2 * math.pi / (K_cage / kT)) ** 1.5
    v_cage_UB = (2 * math.pi / (K0    / kT)) ** 1.5
    v0        = 1660.0  # Å³
    minus_TdS_cage_B  = kT * math.log(v0 / v_cage_B)
    minus_TdS_cage_UB = kT * math.log(v0 / v_cage_UB)
    TdeltadS_cage = minus_TdS_cage_B - minus_TdS_cage_UB

    conf_entropy = min_dG_B - min_dG_UB - TdeltadS_cage
    return conf_entropy