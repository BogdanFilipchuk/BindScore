import bindscore.scoring as sc
import bindscore.parsing
import pathlib as path
test_pdb_path:path.Path = path.Path(__file__).parent / "tests" / "6PYH.pdb"
total_entropy_result = sc.compute_total_entropy()