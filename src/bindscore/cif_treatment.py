import gemmi
import os

def cif_treatment(cif_id:str)->str:
    '''
    Converts CIF file to PDB format
    Args:
        cif_id (str): The 4-character CIF ID or a local file path to a .cif file.
    Returns:
        str: The PDB data as a string.
    '''

    # Check if the input is a local file path
    if os.path.exists(cif_id):
        print(f"Reading local CIF file: {cif_id}")
        cif = gemmi.read_cif_file(cif_id)
    
    # If not a local file, treat it as a CIF ID and fetch from RCSB PDB database
    cif_id = cif_id.upper()
    
    if len(cif_id) != 4:
        raise ValueError("CIF ID must be exactly 4 characters long.")
        
    print(f"Fetching CIF file for {cif_id}...")
    cif = gemmi.read_cif_file(f"https://files.rcsb.org/download/{cif_id}.cif")

    # Convert CIF to PDB format
    pdb_data = cif.to_pdb()
    return pdb_data