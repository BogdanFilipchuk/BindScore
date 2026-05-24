from pdb_utils_inter import *
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / 'scoring'))
from pdb_utils_enthalpy import *
from pdb_utils_protein import Protein_Structure

def get_dataset_interaction_list(pdb_id, chain_a, chain_b, threshold=5.0):
    """
    Iterates over the dataset and returns a flat list containing every single
    interaction dict across all proteins using the exact requested format.
    """
    all_extracted_interactions = []

    # 1. Fetch data and clean NMR ensembles if present
    pdb_data = fetch_pdb_data(pdb_id)
    if "MODEL" in pdb_data:
        pdb_data = pdb_data.split("ENDMDL")[0] + "ENDMDL\n"

    # 2. Build the structural objects
    protein = Protein_Structure(pdb_data)
    print('The entities in the protein structure are:', '\n Protein chains:', protein.chains(), '\n Small molecules:', protein.small_molecules())
    inter = Interaction(protein, chain_a, chain_b, threshold=threshold)

    # 3. Calculate enthalpies (populates the 'energy' keys inside inter.interactions)
    _ = binding_enthalpy(inter)

    # 4. Map the data directly to your exact dictionary schema
    for interaction in inter.interactions:
        itype = interaction['type'] or 'unclassified'
        energy = interaction.get('energy', 0.0)

        all_extracted_interactions.append({
            "atom1": interaction["atom1"],
            "atom2": interaction["atom2"],
            "distance": interaction["distance"],
            "interaction type": itype,
            "binding energy": energy
        })

    return all_extracted_interactions


# ── Execution ─────────────────────────────────────────────────────────────────
# Call the function directly with your configuration dataset

protein = '6PYH'
chain1 = 'A'
chain2 = 'B'
THRESHOLD = 5.0

dataset_values = get_dataset_interaction_list(protein, chain1, chain2, THRESHOLD)

# Example: Inspecting the first extracted dictionary entry from the master list
if dataset_values:
    print(f"Successfully extracted {len(dataset_values)} total interaction pairs.")
    print("Sample dictionary format:")
    print(dataset_values[0])