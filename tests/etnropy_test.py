from bindscore.scoring import *
from bindscore.pdb_file_treatment import pdb_utils_fetch as fetch
import bindscore.parsing
import pathlib as path
from bindscore.parsing.pdb_utils_protein import Protein_Structure
test_pdb_path:path.Path = path.Path(__file__).parent / "test_pdb_files" / "6PYH.pdb"
# total_entropy_result = compute_total_entropy(test_pdb_path,"A","D")
# print("total_entropy_result)

myprotein = Protein_Structure(fetch.fetch_pdb_data(test_pdb_path))
print(myprotein.get_ID)