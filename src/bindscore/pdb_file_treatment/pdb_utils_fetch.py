import os
import pathlib
from urllib import response
import gemmi
import requests
import pathlib

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
    If it detects a CIF format, it will convert it to PDB format automatically.
    Args:
        pdb_id (str): The 4-character PDB ID or a local file path 
        to a .pdb or .cif file of the protein.
    Returns:
        str: The PDB data as a string.
    '''

    # Check if the input is a local file path
    if os.path.exists(pdb_id):
        # Check if the file is in CIF or PDB format based on its extension
        ext = os.path.splitext(pdb_id)[1].lower()

        if ext in ('.cif', '.mmcif'):
            print(f"CIF file detected, converting to PDB...")
            doc = gemmi.cif.read_file(pdb_id)
            structure = gemmi.make_structure_from_block(doc.sole_block())
            return structure.make_pdb_string()
        elif ext in ('.pdb',):
            print(f"Reading local PDB file: {pdb_id}")
            with open(pdb_id, 'r') as f:
                return f.read()
        else:
            raise ValueError("Unsupported file format. Please provide a .pdb or .cif file.")

    # If not a local file, treat it as a PDB ID and fetch from RCSB PDB database
    pdb_id = pdb_id.upper()

    if len(pdb_id) != 4:
        raise ValueError("PDB or CIF ID must be exactly 4 characters long.")
    
    print(f"Fetching PDB file for {pdb_id}...")
    url = f'https://files.rcsb.org/download/{pdb_id}.pdb'

    raw_pdb_data = requests.get(url)

    # If it is a CIF file, convert it to PDB format
    if raw_pdb_data.status_code == 404:
        print(f"PDB not found, trying CIF...")
        raw_cif_data = requests.get(f'https://files.rcsb.org/download/{pdb_id}.cif')
        raw_cif_data.raise_for_status() # Error if the request is unsuccessful
        print("CIF file fetched, converting to PDB...")
        doc = gemmi.cif.read_string(raw_cif_data.text)
        structure = gemmi.make_structure_from_block(doc.sole_block())
        return structure.make_pdb_string()
    
    raw_pdb_data.raise_for_status() # Error if the request is unsuccessful AND the problem is not 404 (not found)
    pdb_data = raw_pdb_data.text
        
    print(f"PDB data fetched successfully.")
        
    return pdb_data


def fetch_pdb_file(pdb_id: str, save_dir: pathlib.Path = pathlib.Path('.')) -> pathlib.Path:
    '''
    Fetches the PDB data for the specified protein from a local file path
    or from the RCSB PDB database, saves it as a .pdb file, and returns the file path.
    If it detects a CIF format, it will convert it to PDB format automatically.
    Args:
        pdb_id (str): The 4-character PDB ID or a local file path
        to a .pdb or .cif file of the protein.
        save_dir (pathlib.Path): Directory to save the PDB file. Defaults to current directory.
    Returns:
        pathlib.Path: The path to the saved PDB file.
    '''
    pdb_path = pathlib.Path(pdb_id)

    if pdb_path.exists():
        ext = pdb_path.suffix.lower()

        if ext in ('.cif', '.mmcif'):
            print(f"CIF file detected, converting to PDB...")
            doc = gemmi.cif.read_file(str(pdb_path))
            structure = gemmi.make_structure_from_block(doc.sole_block())
            pdb_string = structure.make_pdb_string()
            base_name = pdb_path.stem
        elif ext == '.pdb':
            print(f"Local PDB file detected, returning path...")
            return pdb_path
        else:
            raise ValueError("Unsupported file format. Please provide a .pdb or .cif file.")

    else:
        pdb_id = pdb_id.strip().upper()
        if len(pdb_id) != 4:
            raise ValueError("PDB or CIF ID must be exactly 4 characters long.")

        print(f"Fetching PDB file for {pdb_id}...")
        response = requests.get(f'https://files.rcsb.org/download/{pdb_id}.pdb')

        if response.status_code == 404:
            print(f"PDB not found, trying CIF...")
            cif_response = requests.get(f'https://files.rcsb.org/download/{pdb_id}.cif')
            cif_response.raise_for_status()
            print("CIF file fetched, converting to PDB...")
            doc = gemmi.cif.read_string(cif_response.text)
            structure = gemmi.make_structure_from_block(doc.sole_block())
            pdb_string = structure.make_pdb_string()
        else:
            response.raise_for_status()
            pdb_string = response.text

        base_name = pdb_id

    save_path = save_dir / f'{base_name}.pdb'
    save_path.write_text(pdb_string)

    print(f"PDB file saved to {save_path}.")
    return save_path
    
