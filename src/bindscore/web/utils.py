import os
import requests
from requests.exceptions import RequestException

pdb_database_url:str = "https://files.rcsb.org"


def check_connection(url: str):
    try:
        response = requests.get(url)
        response.raise_for_status()
        #print(f"Connection with {url} established")
    except RequestException as err:
        raise ConnectionError(f"Cannot reach database: {err}")
#check_connection(database_url)


def fetch_pdb_file(database_url: str, pdb_id: str, save_dir: str = ".") -> str:
    """"
    Fetching a PDB file from the RCSB PDB database.
    Arguments:
        base_url: Base URL of the database 
        pdb_id:   The 4-character PDB ID of the target structure (Do we want to leave this like this ? To check later)
        save_dir: Directory to save the file (default: current directory)
    Returns:
        Path to the saved PDB file
    """
    check_connection(database_url)

    pdb_id = pdb_id.upper().strip()
    download_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    file_path = os.path.join(save_dir, f"{pdb_id}.pdb")
    
    if len(pdb_id) != 4:
        raise ValueError("PDB ID must be exactly 4 characters long.")
    

    print(f"Fetching {pdb_id} from {download_url} ...")
    response = requests.get(download_url)
    response.raise_for_status()
    with open(file_path, "wb") as f:
        f.write(response.content)
    print(f"Saved to {file_path}")

    return file_path

fetch_pdb_file(pdb_database_url,"6PYH" )
