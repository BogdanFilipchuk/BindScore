import os
from pathlib import Path
import requests

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

def fetch_pdb_data(pdb_id:str)->str:
    '''
    Fetches the PDB data for the specified protein from a local file path 
    of from the RCSB PDB database and returns the raw PDB data as a string.
    Args:
        pdb_id (str): The 4-character PDB ID or a local file path 
        to a .pdb file of the protein.
    Returns:
        str: The PDB data as a string.
    '''

    # Check if the input is a local file path
    if os.path.exists(pdb_id):
        print(f"Reading local PDB file: {pdb_id}")
        with open(pdb_id, 'r') as f:
            return f.read()

    # If not a local file, treat it as a PDB ID and fetch from RCSB PDB database
    pdb_id = pdb_id.upper()

    if len(pdb_id) != 4:
        raise ValueError("PDB ID must be exactly 4 characters long.")
    
    print(f"Fetching PDB file for {pdb_id}...")
    url = f'https://files.rcsb.org/download/{pdb_id}.pdb'

    raw_pdb_data = requests.get(url)
    raw_pdb_data.raise_for_status() # Error if the request is unsuccessful
    pdb_data = raw_pdb_data.text
        
    print(f"PDB data fetched successfully.")
        
    return pdb_data

# def fetch_pdb_file(pdb_id):
    
