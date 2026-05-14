import urllib.request
from pathlib import Path
import numpy as np

# FETCHING PDB FILES FROM RCSB PDB DATABASE AND STORING THEIR RAW DATA AS STR, IT DOES NOT DOWNLOAD IT.

def fetch_not_save_pdb(pdb_id):
    '''
    Fetches the PDB file data for the specified protein.
    Args:
        pdb_id (str): The 4-character PDB ID of the protein.
    Returns:
        str: The raw PDB data as a string.
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
    