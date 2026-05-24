import header

from pdb_utils_fetch import fetch_pdb_data
from pdb_utils_protein import Protein_Structure
from pdb_utils_inter import Interaction
from pdb_utils_enthalpy import binding_enthalpy, interaction_energy, get_partial_charge

# Reference ranges (kJ/mol)
reference_ranges = {
    'salt_bridge':         (-17.0,  -84.0),
    'hydrogen_bond':       ( -8.0,  -67.0),
    'hydrophobic_contact': ( -0.8,   -8.4),
    'pi-pi_stacking':      (-17.0,  -50.0),
    'halogen_bond':        ( -8.0,  -33.0),
    'dipole-dipole':       ( -2.0,  -21.0),
}

# Example proteins with known binding enthalpy values (kJ/mol) for validation
example_proteins = [
    {'pdb_id': '1DPU', 'chain_1': 'A', 'chain_2': 'B'},
    {'pdb_id': '1RST', 'chain_1': 'B', 'chain_2': 'P'},
    {'pdb_id': '2LQC', 'chain_1': 'A', 'chain_2': 'B'},
    {'pdb_id': '2MNU', 'chain_1': 'A', 'chain_2': 'B'},
    {'pdb_id': '2MWY', 'chain_1': 'A', 'chain_2': 'B'},
    {'pdb_id': '4F14', 'chain_1': 'A', 'chain_2': 'B'},
    {'pdb_id': '4Q6F', 'chain_1': 'A', 'chain_2': 'F'},
    {'pdb_id': '5E0M', 'chain_1': 'A', 'chain_2': 'C'},
    {'pdb_id': '5OVC', 'chain_1': 'A', 'chain_2': 'B'},
    {'pdb_id': '6EVO', 'chain_1': 'A', 'chain_2': 'C'},
    {'pdb_id': '6H8C', 'chain_1': 'A', 'chain_2': 'B'},
    {'pdb_id': '1BRS', 'chain_1': 'A', 'chain_2': 'D'},
    {'pdb_id': '3MXF', 'chain_1': 'A', 'chain_2': 'JQ1'},
    {'pdb_id': '3U5L', 'chain_1': 'A', 'chain_2': '08K'},
    {'pdb_id': '4LZR', 'chain_1': 'A', 'chain_2': 'LOC'},
    {'pdb_id': '4QB3', 'chain_1': 'A', 'chain_2': '30M'},
    {'pdb_id': '4XY9', 'chain_1': 'A', 'chain_2': '43U'},
    {'pdb_id': '5D0C', 'chain_1': 'A', 'chain_2': 'E0B'},
    {'pdb_id': '5D3S', 'chain_1': 'A', 'chain_2': '579'},
    {'pdb_id': '5DW2', 'chain_1': 'A', 'chain_2': '5GD'},
    {'pdb_id': '5FBX', 'chain_1': 'A', 'chain_2': '5W4'},
    {'pdb_id': '5IGK', 'chain_1': 'A', 'chain_2': 'BMF'},

]

exp_values = {'1DPU': -70.32, '1RST': -52.58, '2LQC': -28.92, 
              '2MNU': -19.25, '2MWY': -73.25, '4F14': -35.41, 
              '4Q6F': -41.06, '5E0M': -34.32, '5OVC': -28.46, 
              '6EVO': -36.42, '6H8C': -24.74, '1BRS': -77.85, 
              '3MXF': -35.24, '3U5L': -25.78, '4LZR': -37.67, 
              '4QB3': -27.71, '4XY9': -25.49, '5D0C': -42.70, 
              '5D3S': -40.90, '5DW2': -42.28, '5FBX': -65.18, 
              '5IGK': -46.42}

threshold = 5.0

comparison_table = []

for entry in example_proteins:
    pdb_id  = entry['pdb_id']
    chain_1 = entry['chain_1']
    chain_2 = entry['chain_2']

    pdb_data = fetch_pdb_data(pdb_id)
    
    # NMR Ensemble Shield: If the file contains stacked NMR models, truncate the text data right at the end of the first model block.
    if "MODEL" in pdb_data:
        pdb_data = pdb_data.split("ENDMDL")[0] + "ENDMDL\n"
        print("  → NMR ensemble detected. Extracted Model 1 coordinates exclusively.")

    protein = Protein_Structure(pdb_data)
    summary = protein.summary()

    print("\n--- Structure Summary ---")
    print(f"  Total atoms      : {summary['num_atoms_total']}")
    print(f"  Chains           : {summary['chain_ids']}")
    print(f"  Small molecules  : {summary['small_molecules']}")

    print(f"\nDetecting interactions between chain {chain_1} and chain {chain_2}...")
    inter = Interaction(protein, chain_1, chain_2, threshold=threshold)
    print(f"  Total contacts within {threshold} Å : {len(inter.interactions)}")

    type_counts: dict[str, int] = {}
    for i in inter.interactions:
        t = i['type'] or 'unclassified'
        type_counts[t] = type_counts.get(t, 0) + 1

    print("\n--- Interaction counts ---")
    for t, count in sorted(type_counts.items()):
        print(f"  {t:<22} {count}")

    print("\nCalculating binding enthalpy...")
    breakdown = binding_enthalpy(inter)

    print("\n--- Energy breakdown ---")
    for term, value in breakdown.items():
        if term != 'TOTAL':
            print(f"  {term:<22} {value:+.3f} kJ/mol")
    print(f"  {'─'*40}")
    print(f"  {'TOTAL':<22} {breakdown['TOTAL']:+.3f} kJ/mol")

    comparison_table.append({
    "pdb": pdb_id,
    "experimental": exp_values[pdb_id],
    "calculated": breakdown['TOTAL'],
    "error": breakdown['TOTAL'] - exp_values[pdb_id]
    })

    print("\n--- Top 5 interactions by |energy| ---")
    classified = [i for i in inter.interactions if i['type'] is not None]
    top5 = sorted(classified, key=lambda i: abs(i['energy']), reverse=True)[:5]
    for i in top5:
        a1, a2 = i['atom1'], i['atom2']
        print(
            f"  {a1['residue_name']}{a1['residue_seq']:<4} {a1['atom_name']:<5} --"
            f"  {a2['residue_name']}{a2['residue_seq']:<4} {a2['atom_name']:<5}"
            f"  [{i['type']:<22}]  r={i['distance']:.2f} Å  E={i['energy']:+.2f} kJ/mol"
        )

    print("\n--- Per-interaction averages vs reference ranges ---")
    by_type: dict[str, list[float]] = {}
    for i in inter.interactions:
        t = i['type']
        if t and t != 'unclassified':
            by_type.setdefault(t, []).append(i['energy'])

    for itype, energies in sorted(by_type.items()):
        avg      = sum(energies) / len(energies)
        lo, hi   = reference_ranges.get(itype, (None, None))
        if lo is not None and hi is not None:
            in_range = hi <= avg <= lo
            flag     = '✓' if in_range else '⚠ outside reference'
        else:
            flag = '(no reference)'
        ref_str = f"[{hi:.1f}, {lo:.1f}]" if lo is not None else "[ n/a ]"
        print(f"  {itype:<22}  n={len(energies):<3}  avg={avg:+.2f} kJ/mol  ref={ref_str}  {flag}")

print("\n\n Experimental vs Calculated Binding Enthalpy \n")

print(
    f"{'PDB':<8}"
    f"{'Experimental':>18}"
    f"{'Calculated':>18}"
    f"{'Error':>14}"
)

print("-" * 58)

for row in comparison_table:
    print(
        f"{row['pdb']:<8}"
        f"{row['experimental']:>18.2f}"
        f"{row['calculated']:>18.2f}"
        f"{row['error']:>14.2f}"
    )
