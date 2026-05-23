import freesasa as sasa
import pathlib 
import bindscore.web.utils
import bindscore.pdb_file_treatment.pdb_utils  as pdb_utils
import bindscore.pdbtime as pdbtime
"""
This module calculates the solvent available surface area contribution to entropy.
 It uses FreeSASA module, available at (add)
"""
 
current_directory=pathlib.Path('.')
pdb_file_path:pathlib.Path = (r"C:\Users\bogfi\VScode\6PYH.pdb")
MyProtein = pdbtime.Protein_Structure(pdb_file_path)

print(MyProtein.chains())

pdbname:str = "1A2B.pdb"     #FOR DEBUG, change PDB name if want to check something else. You are smart, figure it out.
# with open(pdbname,"r") as file:
#     pdbdata = file.read()
# # print(pdbdata)

# structure = sasa.Structure("1A2B.pdb")
# result = sasa.calc(structure)
# print(result.totalArea())
# classifier = sasa.Classifier()
# atom_class = classifier.classify("ALA", "CA")
# print(atom_class)

 
# def getsasa(pdb_path:pathlib.Path)->float:
#     ###We are running freesasa.calc on 3 things to calculate the delta SASA
#     sasa_two_chains_together  = sasa.calc(sasa.Structure(str(pdb_path)))       # both chains together
#     sasa_protein  = sasa.calc(sasa.Structure(str(pdb_utils.get_protein_from_pdb(pdb_path))))  # Protein alone
#     sasa_ligand   = sasa.calc(sasa.Structure(str(pdb_utils.get_ligand_from_pdb(pdb_path))))   # Ligand alone
#     buried_sasa = sasa_protein.totalArea() + sasa_ligand.totalArea() - sasa_complex.totalArea()
#     return buried_sasa

# print(getsasa(pathlib.Path(r"C:\Users\bogfi\VScode\BindScore Project\1A2B.pdb")))

# # get_ligand_from_pdb(pathlib.Path(r"C:\Users\bogfi\VScode\BindScore Project\1A2B.pdb"))
