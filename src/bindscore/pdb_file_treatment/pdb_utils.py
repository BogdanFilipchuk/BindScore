import pathlib 

def get_ligand_from_pdb(pdb_path:pathlib.Path)->pathlib.Path:
    """
    Returns the path of a file containing only ligand part of the pdb
    """
    lines = pdb_path.read_text().splitlines()
    only_ligand_lines:list = [l for l in lines if l.startswith("HETATM")]
    only_ligand_file_path = pdb_path.parent/f"{pdb_path.stem}_ligand_only.pdb"  #saves the ligand only file. Check if we want to save in this format later
    with open(only_ligand_file_path,"w") as file:
        file.write("\n".join(only_ligand_lines))
    return only_ligand_file_path

def get_protein_from_pdb(pdb_path:pathlib.Path)->pathlib.Path:
    """
    Returns the path of a file containing only protein part of the pdb
    """
    with open(pdb_path,"r") as file:
        lines = file.readlines()
    only_protein_lines:list = [l for l in lines if l.startswit("ATOM")]
    only_protein_file_path = pdb_path.parent/f"{pdb_path.stem}_protein_only.pdb"
    with open(only_protein_file_path,"w") as file:
        file.write("\n".join(only_protein_lines))
    return only_protein_file_path
