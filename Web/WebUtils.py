import requests 
import os
from requests.exceptions import RequestException

database_url = "https://files.rcsb.org"     #Should we just hardcode the url? It's not like we plan to add a choice of fetching url elsewhere. Remains to see. Will we add pdb file imports ? Who knows.


def check_connection(url:str):         #An utility function which raises an error if there are requests exceptions when accessing the database
    try:
        response = requests.get(url)
        response.raise_for_status()
    except RequestException as err:
      raise ConnectionError(f"Cannot reach database: {err}")

def fetch_pdb(url:str ,pdb_id: str, save_dir: str = ".") -> str:
    check_connection(url)
    """
    Fetch a PDB file from the RCSB PDB database.

    Args:
        pdb_id: The 4-character PDB ID of the target structure
        save_dir: Directory to save the file (default: current directory)

    Returns:
        Path to the saved PDB file
    """
    pdb_id = pdb_id.upper().strip()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    file_path = os.path.join(save_dir, f"{pdb_id}.pdb")

    print(f"Fetching {pdb_id} from {url} ...")
    urllib.request.urlretrieve(url, file_path)
    print(f"Saved to {file_path}")

    return file_path


# def parse_atoms(file_path: str, chains: list[str] | None = None) -> tuple[dict, list[str]]:
#     """
#     Parse ATOM records from a PDB file.

#     Args:
#         file_path: Path to the PDB file
#         chains: List of chain IDs to include (e.g. ['A', 'B']).
#                 If None, all chains are included.

#     Returns:
#         A tuple of:
#           - dict keyed by (chain, residue_number, res_name, atom_name) -> (x, y, z)
#           - sorted list of all chain IDs found in the file
#     """
#     coordinates = {}
#     all_chains = set()

#     with open(file_path, 'r') as f:
#         for line in f:
#             if not line.startswith("ATOM  "):
#                 continue

#             chain     = line[21].strip()
#             atom_name = line[12:16].strip()
#             res_name  = line[17:20].strip()
#             res_seq   = int(line[22:26].strip())
#             x         = float(line[30:38].strip())
#             y         = float(line[38:46].strip())
#             z         = float(line[46:54].strip())

#             all_chains.add(chain)  # always track, regardless of filter

#             if chains is not None and chain not in chains:
#                 continue

#             key = (chain, res_seq, res_name, atom_name)
#             coordinates[key] = (x, y, z)

#     return coordinates, sorted(all_chains)

# def get_chain_coords(file_path: str, chain_id: str) -> list[tuple]:
#     """
#     Return a list of (x, y, z) coordinates for all atoms in a given chain.
#     Mirrors the original coords() function, generalised to any chain.

#     Args:
#         file_path: Path to the PDB file
#         chain_id:  Single chain identifier (e.g. 'B')

#     Returns:
#         List of (x, y, z) tuples
#     """
#     coords_list = []
#     with open(file_path, 'r') as f:
#         for line in f:
#             if not line.startswith("ATOM  "):
#                 continue
#             if line[21].strip() != chain_id:
#                 continue
#             x = float(line[30:38].strip())
#             y = float(line[38:46].strip())
#             z = float(line[46:54].strip())
#             coords_list.append((x, y, z))
#     return coords_list


# def summarise(coordinates: dict) -> None:
#     """Print a brief summary of the parsed atoms."""
#     chains = sorted({k[0] for k in coordinates})
#     print(f"\nParsed {len(coordinates)} atoms across chain(s): {', '.join(chains)}")
#     for chain in chains:
#         n = sum(1 for k in coordinates if k[0] == chain)
#         print(f"  Chain {chain}: {n} atoms")


# # ---------------------------------------------------------------------------
# # Example usage
# # ---------------------------------------------------------------------------
# if __name__ == "__main__":
#     PDB_ID  = "1TIM"          # Change to any valid PDB ID
#     CHAINS  = ["A", "B"]      # Set to None to include all chains

#     # 1. Fetch the file
#     pdb_path = fetch_pdb(PDB_ID, save_dir=".")

#     # 2. Parse all atoms (optionally filtered by chain)
#     coordinates_with_ids = parse_atoms(pdb_path, chains=CHAINS)
#     summarise(coordinates_with_ids)

#     # 3. Example: grab just (x,y,z) list for chain B (original-style)
#     chain_b_coords = get_chain_coords(pdb_path, chain_id="B")
#     print(f"\nFirst 3 coords in chain B: {chain_b_coords[:3]}")

#     # 4. Example: look up a specific atom
#     sample_key = next(iter(coordinates_with_ids))
#     print(f"\nSample entry — key: {sample_key}, xyz: {coordinates_with_ids[sample_key]}")