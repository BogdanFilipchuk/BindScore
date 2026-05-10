import urllib.request
from pathlib import Path


def fetch_pdb(pdb_id):
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

        # Save the PDB data to a local file
        output_path = Path(f'{pdb_id}.pdb')
    
        with open(output_path, 'w') as f:
            f.write(pdb_data)
    
        print(f"PDB file saved to {output_path.resolve()}")
        
        return str(output_path.resolve())
    
class Protein_Structure:
    def __init__(self, pdb_file):
        '''
        Initializes the Protein_Structure object by parsing the PDB file.
        The chains (polypeptides) and small molecules are separated for easier access.
        Args:
            pdb_file (str): The path to the PDB file.
        Functions:
            pdb_file - The path to the PDB file.
            all_atoms() - Parses the PDB file and extracts relevant information.
            backbone() - Extracts the backbone atoms from the parsed PDB data.
            side_chain() - Extracts the side chain atoms from the parsed PDB data.
            hetero_atoms() - Extracts the hetero atoms from the parsed PDB data.
            chains() - Extracts the unique chains from the parsed PDB data.
            get_chain(chain_id) - Extracts the atoms for a specific chain from the parsed PDB data.
            small_molecules() - Extracts the unique small molecules from the parsed PDB data.
            get_small_molecule(small_molecule_id) - Extracts the atoms for a specific small molecule from the parsed PDB data.
            summary() - Provides a summary of the protein structure.
        '''

        self.pdb_file  = pdb_file
        self.all_atoms = self.parse_pdb()

    def parse_pdb(self):
        '''
        Parses the PDB file and extracts relevant information.
        Args:
            pdb_file (str): The path to the PDB file.
        Returns:
            list: A list of dictionaries containing the parsed information.
        '''

        parsed_data = []

        with open(self.pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    # Parse
                    info = {
                        'type' : line[0:6].strip(),
                        'atom_seq': line[6:11].strip(),
                        'atom_name': line[12:16].strip(),
                        'residue_name': line[17:20].strip(),
                        'chain': line[21].strip(),
                        'residue_seq': line[22:26].strip(),
                        'x': float(line[30:38].strip()),
                        'y': float(line[38:46].strip()),
                        'z': float(line[46:54].strip())
                    }
                    parsed_data.append(info)

        return parsed_data

    def backbone(self):
        '''
        Extracts the backbone atoms from the parsed PDB data.
        Returns:
            list: A list of dictionaries containing only the backbone atoms.
        '''

        #Backbone atoms are N, CA, C, O
        backbone_atoms = []
        for atom in self.all_atoms:
            if atom['type'] == 'ATOM' and atom['atom_name'] in {'N', 'CA', 'C', 'O'}:
                backbone_atoms.append(atom)

        return backbone_atoms

    def side_chain(self):
        '''
        Extracts the side chain atoms from the parsed PDB data.
        Returns:
            list: A list of dictionaries containing only the side chain atoms.
        '''

        # Side chain atoms are all atoms except the backbone atoms
        side_chain_atoms = []
        for atom in self.all_atoms:
            if atom['type'] == 'ATOM' and atom['atom_name'] not in {'N', 'CA', 'C', 'O'}:
                side_chain_atoms.append(atom)

        return side_chain_atoms
    
    def hetero_atoms(self):
        '''
        Extracts the hetero atoms from the parsed PDB data.
        Returns:
            list: A list of dictionaries containing only the hetero atoms.
        '''

        # Hetero atoms are all atoms that are not part of the standard amino acid residues
        hetero_atoms = []
        for atom in self.all_atoms:
            if atom['type'] == 'HETATM':
                hetero_atoms.append(atom)

        return hetero_atoms
    
    def chains(self):
        '''
        Extracts the unique chains from the parsed PDB data.
        Returns:
            set: A set of unique chains.
        '''

        chains = set()
        for atom in self.all_atoms:
            if atom['type'] == 'ATOM':
                chains.add(atom['chain'])

        return chains
    
    def get_chain(self, chain_id):
        '''
        Extracts the atoms for a specific chain from the parsed PDB data.
        Args:
            chain_id (str): The ID of the chain to extract atoms from.
        Returns:
            list: A list of dictionaries containing the atoms for the specified chain.
        '''

        residue_atoms = []
        for atom in self.all_atoms:
            if atom['type'] == 'ATOM' and atom['chain'] == chain_id:
                residue_atoms.append(atom)

        return residue_atoms
    
    def small_molecules(self):
        '''
        Extracts the unique small molecules from the parsed PDB data.
        Returns:
            set: A set of unique small molecules.
        '''

        small_molecules = set()
        for atom in self.all_atoms:
            if atom['type'] == 'HETATM':
                small_molecules.add(atom['residue_name'])

        return small_molecules
    
    def get_small_molecule(self, small_molecule_id):
        '''
        Extracts the atoms for a specific small molecule from the parsed PDB data.
        Args:
            small_molecule_id (str): The ID of the chain to extract atoms from.
        Returns:
            list: A list of dictionaries containing the atoms for the specified small molecule.
        '''

        #The water molecules are grouped together, not separate. So if the small_molecule_id is 'HOH', it will return all the water molecules.

        residue_atoms = []
        for atom in self.all_atoms:
            if atom['type'] == 'HETATM' and atom['residue_name'] == small_molecule_id:
                residue_atoms.append(atom)

        return residue_atoms

    def summary(self):
        '''
        Provides a summary of the protein structure.
        Returns:
            dict: A dictionary containing the summary.
        '''

        summary_dict = {
        'pdb_file'          : self.pdb_file,
        'num_atoms_total'   : len(self.all_atoms),
        'num_chains'        : len(self.chains()),
        'num_small_molecules' : len(self.small_molecules()),
        'chain_ids'         : sorted(self.chains()),
        'small_molecules'   : self.small_molecules()
        }

        # Can be extended

        return summary_dict
        

def example_usage():
    '''
    Example usage of the functions.
    '''

    fetch_pdb('1A2B')
    protein = Protein_Structure('1A2B.pdb')
    print(protein.summary())

example_usage()