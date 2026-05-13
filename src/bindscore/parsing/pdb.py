def parse_atoms(file_path: str, chains: list[str] | None = None) -> tuple[dict, list[str]]:
    """
    Parse ATOM records from a PDB file.

    Args:
        file_path: Path to the PDB file
        chains:    Chain IDs to include (e.g. ['A', 'B']). None = all chains.

    Returns:
        - dict keyed by (chain, residue_number, res_name, atom_name) -> (x, y, z)
        - sorted list of all chain IDs found in the file
    """
    coordinates = {}
    all_chains = set()

    with open(file_path, "r") as f:
        for line in f:
            if not line.startswith("ATOM  "):
                continue

            chain     = line[21].strip()
            atom_name = line[12:16].strip()
            res_name  = line[17:20].strip()
            res_seq   = int(line[22:26].strip())
            x         = float(line[30:38].strip())
            y         = float(line[38:46].strip())
            z         = float(line[46:54].strip())

            all_chains.add(chain)

            if chains is not None and chain not in chains:
                continue

            key = (chain, res_seq, res_name, atom_name)
            coordinates[key] = (x, y, z)

    return coordinates, sorted(all_chains)


def get_chain_coords(file_path: str, chain_id: str) -> list[tuple]:
    """
    Return a list of (x, y, z) coordinates for all atoms in a given chain.

    Args:
        file_path: Path to the PDB file
        chain_id:  Single chain identifier (e.g. 'B')

    Returns:
        List of (x, y, z) tuples
    """
    coords_list = []
    with open(file_path, "r") as f:
        for line in f:
            if not line.startswith("ATOM  "):
                continue
            if line[21].strip() != chain_id:
                continue
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            coords_list.append((x, y, z))
    return coords_list


def summarise(coordinates: dict) -> None:
    """Print a brief summary of parsed atoms."""
    chains = sorted({k[0] for k in coordinates})
    print(f"\nParsed {len(coordinates)} atoms across chain(s): {', '.join(chains)}")
    for chain in chains:
        n = sum(1 for k in coordinates if k[0] == chain)
        print(f"  Chain {chain}: {n} atoms")
