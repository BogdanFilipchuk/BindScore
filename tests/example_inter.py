

"""
Test script for protein interaction analysis
Tests specific entity pairs per protein structure
"""

from bindscore.pdb_file_treatment.pdb_utils_fetch import fetch_pdb_data
from bindscore.parsing.pdb_utils_protein import Protein_Structure
from bindscore.parsing.pdb_utils_inter import Interaction


def print_interactions(inter):
    interaction_types = {}
    for interaction in inter.interactions:
        itype = interaction['type']
        if itype:
            if itype not in interaction_types:
                interaction_types[itype] = []
            interaction_types[itype].append(interaction)

    if not interaction_types:
        print("  (no interactions found)")
        return

    for itype, interactions_list in sorted(interaction_types.items()):
        print(f"  {itype}: {len(interactions_list)} interactions")
        for i, interaction in enumerate(interactions_list):
            atom1 = interaction['atom1']
            atom2 = interaction['atom2']
            dist = interaction['distance']
            print(f"    [{i+1}] {atom1['residue_name']}{atom1['residue_seq']}:{atom1['atom_name']} <-> "
                  f"{atom2['residue_name']}{atom2['residue_seq']}:{atom2['atom_name']} "
                  f"({dist:.2f} A)")


def test_protein(pdb_id, description, pairs):
    print("\n" + "="*80)
    print(f"Testing: {pdb_id} - {description}")
    print("="*80)

    try:
        pdb_data = fetch_pdb_data(pdb_id)
        protein = Protein_Structure(pdb_data)

        summary = protein.summary()
        print(f"\nStructure Summary:")
        print(f"  Total atoms: {summary['num_atoms_total']}")
        print(f"  Chains: {summary['chain_ids']}")
        print(f"  Small molecules: {summary['small_molecules']}")
        print(f"  Metals: {len(protein.metals())} metal atoms")

        for entity1, entity2 in pairs:
            print(f"\n--- {entity1} - {entity2} ---")
            try:
                inter = Interaction(protein, entity1, entity2, threshold=5.0)
                print_interactions(inter)
            except Exception as e:
                print(f"  Error: {e}")

        print(f"\n[OK] {pdb_id} complete")
        return True

    except Exception as e:
        print(f"\n[FAIL] Error: {pdb_id}: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    test_proteins = [
        ("1BRS", "Barnase-Barstar - protein-protein interface", [
            ("A", "D"),       # barnase vs barstar
        ]),
        ("2PTC", "Trypsin-BPTI - serine protease inhibitor", [
            ("E", "I"),       # trypsin vs BPTI
        ]),
        ("4INS", "Insulin dimer - B-chain pi-pi stacking", [
            ("A", "B"),       # alpha + beta chain
        ]),
        ("1MLC", "Antibody-Lysozyme - CDR antigen interface", [
            ("H", "A"),       # heavy chain vs lysozyme
        ]),
        ("1C1Y", "Ras-Raf RBD - GTPase effector interface", [
            ("A", "B"),       # Ras vs Raf RBD
        ]),
    ]

    print("="*80)
    print("PROTEIN INTERACTION ANALYSIS TEST SUITE")
    print("="*80)
    print(f"Testing {len(test_proteins)} proteins (specific pairs only)")

    results = []
    for pdb_id, description, pairs in test_proteins:
        success = test_protein(pdb_id, description, pairs)
        results.append((pdb_id, success))

    print("\n" + "="*80)
    print("TEST SUMMARY")
    print("="*80)
    successful = sum(1 for _, success in results if success)
    print(f"Successful: {successful}/{len(results)}")
    for pdb_id, success in results:
        status = "[PASS]" if success else "[FAIL]"
        print(f"  {status} - {pdb_id}")


if __name__ == "__main__":
    main()