from bindscore.scoring import *
import bindscore.parsing
import pathlib as path
test_pdb_path:path.Path = path.Path(__file__).parent / "test_pdb_files" / "6PYH.pdb"
total_entropy_result = compute_total_entropy(test_pdb_path,"A","D")
print(total_entropy_result)