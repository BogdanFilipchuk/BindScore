from bindscore.scoring import *
from bindscore.pdb_file_treatment import pdb_utils_fetch as fetch
import bindscore.parsing
import pathlib as path
from bindscore.parsing.pdb_utils_protein import Protein_Structure
from bindscore.scoring.total_entropy import Binding_Entropy_Summary

"""
TEST SETUP
change bool variables to see each individual contribution
========
"""

DEBUG_OVERAL_ENTROPY:bool      = True
DEBUG_TRANS_ROT_ENTROPY:bool   = True
DEBUG_HYDROPHOBIC_ENTROPY:bool = True
DEBUG_SIDECHAIN_ENTROPY:bool   = True
DEBUG_BACKBONE_ENTROPY:bool    = True

T_KELVIN:float = 300.0

"Indicate the PDB file name if you want to test on another one"
# test_pdb_path:path.Path = path.Path(__file__).parent / "test_pdb_files" / "6PYH.pdb"
test_pdb_path:path.Path = path.Path(__file__).parent / "test_pdb_files" / "1BRS.pdb"
# test_pdb_path:path.Path = path.Path(__file__).parent / "test_pdb_files" / "1FQ9.pdb"

"Making a ProteinStructure object for your pdb"
myprotein = Protein_Structure(fetch.fetch_pdb_data(test_pdb_path))


def debug_entropy(result: Binding_Entropy_Summary, protein_id: str) -> None:
    if DEBUG_OVERAL_ENTROPY:
        print(f'[OVERALL]      {protein_id}: {result.dS_total:.4f} J/(mol·K)')
    if DEBUG_TRANS_ROT_ENTROPY:
        print(f'[TRANS_ROT]    {protein_id}: {result.dS_trans_rot:.4f} J/(mol·K)', end='')
        if result.trans_rot_detail:
            d = result.trans_rot_detail
            print(f'  (trans={d.translational:.4f}, rot={d.rotational:.4f})', end='')
        print()
    if DEBUG_HYDROPHOBIC_ENTROPY:
        print(f'[HYDROPHOBIC]  {protein_id}: {result.dS_hydrophobic:.4f} J/(mol·K)', end='')
        if result.hydrophobic_detail:
            d = result.hydrophobic_detail
            print(f'  (R_a={d.R_a:.1f}Å, R_b={d.R_b:.1f}Å)', end='')
        print()
    if DEBUG_SIDECHAIN_ENTROPY:
        print(f'[SIDECHAIN]    {protein_id}: {result.dS_sidechain:.4f} J/(mol·K)')
    if DEBUG_BACKBONE_ENTROPY:
        print(f'[BACKBONE]     {protein_id}: {result.dS_backbone:.4f} J/(mol·K)', end='')
        if result.backbone_detail:
            d = result.backbone_detail
            print(f'  (complex={d.S_complex:.4f}, A={d.S_chain_a:.4f}, B={d.S_chain_b:.4f})', end='')
        print()


"""
CHECKING THE OVERALL ENTROPY RESULTS
========
"""

total_entropy_result = compute_total_entropy(test_pdb_path, "A", "B", T=T_KELVIN, return_breakdown=True)
debug_entropy(total_entropy_result, myprotein.get_ID())