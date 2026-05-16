import pathlib 
import bindscore.pdbtime
# def get_heteroatoms_from_pdb(pdb_path:pathlib.Path)->pathlib.Path:
#     """
#     Returns the path of a file containing only ligand part of the pdb.
#     """
#     lines = pdb_path.read_text().splitlines()
#     only_ligand_lines:list = [l for l in lines if l.startswith("HETATM")]
#     only_ligand_file_path = pdb_path.parent/f"{pdb_path.stem}_ligand_only.pdb"  #saves the atom only file. Old function. Idk if it will be usefull later.
#     with open(only_ligand_file_path,"w") as file:
#         file.write("\n".join(only_ligand_lines) + "\nEND\n")
#     return only_ligand_file_path

# def get_atoms_from_pdb(pdb_path:pathlib.Path)->pathlib.Path:
#     """
#     Returns the path of a file containing only protein part of the pdb.
#     """
#     with open(pdb_path,"r") as file:
#         lines = file.readlines()
#     only_protein_lines:list = [l for l in lines if l.startswith("ATOM")]
#     only_protein_file_path = pdb_path.parent/f"{pdb_path.stem}_protein_only.pdb"
#     with open(only_protein_file_path,"w") as file:
#         file.write("\n".join(only_protein_lines) + "\nEND\n")    

def fetch_pdb_data(pdb_id):
    '''
    Fetches the PDB file for the specified protein and saves it locally.
    Args:
        pdb_id (str): The 4-character PDB ID of the protein.
    Returns:
        str: The full path to the saved PDB file.
    '''

    pdb_id = pdb_id.upper()

    if len(pdb_id) != 4:
        raise ValueError("PDB ID must be exactly 4 characters long.")
    else:
        print(f"Fetching PDB file for {pdb_id}...")
        url = f'https://files.rcsb.org/download/{pdb_id}.pdb'

        with urllib.request.urlopen(url) as response:
            raw_pdb_data = response.read()
        
        pdb_data = raw_pdb_data.decode('utf-8')
    
        print(f"PDB file fetched successfully.")
        
        return pdb_data
    
# def fetch_pdb_file(pdb_id):
    
