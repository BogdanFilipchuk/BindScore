import os
import requests
from requests.exceptions import RequestException

database_url = "https://files.rcsb.org"


def check_connection(url: str):
    try:
        response = requests.get(url)
        response.raise_for_status()
    except RequestException as err:
        raise ConnectionError(f"Cannot reach database: {err}")


def fetch_pdb(base_url: str, pdb_id: str, save_dir: str = ".") -> str:
    """
    Fetch a PDB file from the RCSB PDB database.

    Args:
        base_url: Base URL of the database (used to verify connectivity)
        pdb_id:   The 4-character PDB ID of the target structure
        save_dir: Directory to save the file (default: current directory)

    Returns:
        Path to the saved PDB file
    """
    check_connection(base_url)
    pdb_id = pdb_id.upper().strip()
    download_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    file_path = os.path.join(save_dir, f"{pdb_id}.pdb")

    print(f"Fetching {pdb_id} from {download_url} ...")
    response = requests.get(download_url)
    response.raise_for_status()
    with open(file_path, "wb") as f:
        f.write(response.content)
    print(f"Saved to {file_path}")

    return file_path
