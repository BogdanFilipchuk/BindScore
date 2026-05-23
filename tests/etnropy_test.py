from bindscore.scoring import *
from bindscore.pdb_file_treatment import pdb_utils_fetch as fetch
import bindscore.parsing
import pathlib as path
from bindscore.parsing.pdb_utils_protein import Protein_Structure

DEBUG_ENTROPY:bool = True


"Indicate the PDB file name if you want to test on another one"
test_pdb_path:path.Path = path.Path(__file__).parent / "test_pdb_files" / "6PYH.pdb"

"Making a ProteinStructure object for your pdb"
myprotein = Protein_Structure(fetch.fetch_pdb_data(test_pdb_path))

"Calculating the entropy result"
total_entropy_result = compute_total_entropy(test_pdb_path,"A","D")

print(f'The entropy value calculated for {myprotein.get_ID()} is {total_entropy_result} kcal/mol' )