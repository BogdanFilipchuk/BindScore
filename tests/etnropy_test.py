from bindscore.scoring import *
from bindscore.pdb_file_treatment import pdb_utils_fetch as fetch
import bindscore.parsing
import pathlib as path
from bindscore.parsing.pdb_utils_protein import Protein_Structure

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

"Indicate the PDB file name if you want to test on another one"
test_pdb_path:path.Path = path.Path(__file__).parent / "test_pdb_files" / "6PYH.pdb"
# test_pdb_path:path.Path = path.Path(__file__).parent / "test_pdb_files" / "1FQ9.pdb"

"Making a ProteinStructure object for your pdb"
myprotein = Protein_Structure(fetch.fetch_pdb_data(test_pdb_path))


def debug_entropy(result: EntropySummary, protein_id: str, mode: str) -> None:
    if DEBUG_OVERAL_ENTROPY:
        print(f'[OVERALL]      {protein_id} ({mode}): {result.minusT_dS_total:.4f} kcal/mol')
    if DEBUG_TRANS_ROT_ENTROPY:
        print(f'[TRANS_ROT]    {protein_id} ({mode}): {result.minusT_dS_trans_rot:.4f} kcal/mol')
    if DEBUG_HYDROPHOBIC_ENTROPY:
        print(f'[HYDROPHOBIC]  {protein_id} ({mode}): {result.minusT_dS_hydrophobic:.4f} kcal/mol')
    if DEBUG_SIDECHAIN_ENTROPY:
        print(f'[SIDECHAIN]    {protein_id} ({mode}): {result.minusT_dS_sidechain:.4f} kcal/mol')
    if DEBUG_BACKBONE_ENTROPY:
        print(f'[BACKBONE]     {protein_id} ({mode}): {result.minusT_dS_backbone:.4f} kcal/mol')


"""
CHECKING THE OVERAL ENTROPY RESULTS
2 different modes of hydrophobic entropy calculations were implemented , default - SUN
========
"""

"Calculating the entropy result"
total_entropy_result = compute_total_entropy(test_pdb_path,"A","B", HYDROPHOBIC_SETTING="SASA", return_breakdown=True)
debug_entropy(total_entropy_result, myprotein.get_ID(), "SASA")

print()

"Calculating the entropy result"
total_entropy_result = compute_total_entropy(test_pdb_path,"A","B", HYDROPHOBIC_SETTING="Sun", return_breakdown=True)
debug_entropy(total_entropy_result, myprotein.get_ID(), "Sun")