import os
import pathlib
import sys
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / 'pdb_file_treatment'))
from pdb_utils_fetch import fetch_pdb_data

class Protein_Structure:
    def __init__(self, pdb_data:str):
        '''
        Initializes the Protein_Structure object by parsing the PDB file.
        The chains (polypeptides) and small molecules are separated for easier access.
        Args:
            pdb_data (str): The PDB data as a string.
        Functions:
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

        self.pdb_data  = pdb_data
        self.all_atoms = self.parse_pdb()

    def parse_pdb(self):
        '''
        Parses the PDB file and extracts relevant information.
        Args:
            pdb_data (str): The raw PDB data as a string.
        Returns:
            list: A list of dictionaries containing the parsed information.
        '''

        parsed_data = []

        for line in self.pdb_data.splitlines():
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # Every atom is now a dict as the one below, with all the relevant information extracted from the PDB file.
                # This will make it easier to access the information later on for the enthalpy calculations.
                info = {
                    'type' : line[0:6].strip(),
                    'atom_seq': line[6:11].strip(),
                    'atom_name': line[12:16].strip(),
                    'residue_name': line[17:20].strip(),
                    'chain': line[21].strip(),
                    'residue_seq': line[22:26].strip(),
                    'x': float(line[30:38].strip()),
                    'y': float(line[38:46].strip()),
                    'z': float(line[46:54].strip()),
                    'atom_symbol': line[76:78].strip()
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

        if chain_id not in self.chains():
            raise ValueError("Chain not found in the PDB data.")

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
            small_molecule_id (str): The ID of the small molecule to extract atoms from.
        Returns:
            list: A list of dictionaries containing the atoms for the specified small molecule.
        '''

        # The water molecules are grouped together, not separate. So if the small_molecule_id is 'HOH', 
        # it will return all the water molecules.

        if small_molecule_id not in self.small_molecules():
            raise ValueError("Small molecule not found in the PDB data.")

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

        print("Generating summary of the protein structure...")

        summary_dict = {
        'num_atoms_total'   : len(self.all_atoms),
        'num_chains'        : len(self.chains()),
        'num_small_molecules' : len(self.small_molecules()),
        'chain_ids'         : sorted(self.chains()),
        'small_molecules'   : self.small_molecules()
        }

        # Can be extended

        return summary_dict