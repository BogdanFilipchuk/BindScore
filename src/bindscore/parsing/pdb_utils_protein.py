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
            parse_pdb()                          - Parses the PDB data and extracts relevant information.
            backbone()                           - Extracts the backbone atoms (N, CA, C, O, OXT).
            side_chain()                         - Extracts the side chain atoms.
            ring_atoms()                         - Extracts the aromatic ring atoms (PHE, TYR, TRP, HIS).
            get_ring_centroid(ring_atoms)        - Calculates the centroid of a single ring.
            get_ring_normal(ring_atoms)          - Calculates the normal vector of a single ring.
            get_rings(chain)                     - Extracts all rings of a chain with their centroids and normals.
            hetero_atoms()                       - Extracts the hetero atoms.
            chains()                             - Extracts the unique chains.
            get_chain(chain_id)                  - Extracts the atoms for a specific chain.
            small_molecules()                    - Extracts the unique small molecules.
            get_small_molecule(small_molecule_id)- Extracts the atoms for a specific small molecule.
            get_entity(entity_id)                - Extracts the atoms for a chain or small molecule.
            summary()                            - Provides a summary of the protein structure.
            get_ID()                             - Returns the PDB ID.                           
        '''

        self.pdb_data  = pdb_data
        self.all_atoms = self.parse_pdb()
        self.conect = self.parse_conect()

    def parse_pdb(self):
        '''
        Parses the PDB data stored in self.pdb_data and extracts relevant information for each atom.
        Returns:
            list: A list of dictionaries, one per atom, containing type, atom_seq, atom_name,
                  residue_name, chain, residue_seq, x, y, z, and atom_symbol.
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
                    'chain': line[20:22].strip(),
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
        Extracts the backbone atoms (N, CA, C, O, OXT) from the parsed PDB data.
        Returns:
            list: A list of dictionaries containing only the backbone atoms.
        '''

        #Backbone atoms are N, CA, C, O, OXT
        backbone_atoms = []
        for atom in self.all_atoms:
            if atom['type'] == 'ATOM' and atom['atom_name'] in {'N', 'CA', 'C', 'O', 'OXT'}:
                backbone_atoms.append(atom)

        return backbone_atoms

    def side_chain(self):
        '''
        Extracts the side chain atoms (every ATOM except backbone N, CA, C, O, OXT) from the parsed PDB data.
        Returns:
            list: A list of dictionaries containing only the side chain atoms.
        '''

        # Side chain atoms are all atoms except the backbone atoms
        side_chain_atoms = []
        for atom in self.all_atoms:
            if atom['type'] == 'ATOM' and atom['atom_name'] not in {'N', 'CA', 'C', 'O', 'OXT'}:
                side_chain_atoms.append(atom)

        return side_chain_atoms
    
    def ring_atoms(self):
        '''
        Extracts the aromatic ring atoms from the parsed PDB data.
        Only atoms belonging to the aromatic ring of PHE, TYR, TRP, or HIS are included
        (side chain atoms outside the ring such as CB are excluded).
        Returns:
            list: A list of dictionaries containing only the aromatic ring atoms.
        '''

        aromatic_ring_atoms = {
        'PHE': {'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'},
        'TYR': {'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'},
        'TRP': {'CG', 'CD1', 'CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'},
        'HIS': {'CG', 'ND1', 'CD2', 'CE1', 'NE2'}
        }

        ring_atoms = []
        for atom in self.side_chain():
            if atom['residue_name'] in aromatic_ring_atoms and atom['atom_name'] in aromatic_ring_atoms[atom['residue_name']]:
                ring_atoms.append(atom)

        return ring_atoms
    
    def get_ring_centroid(self, ring_atoms):
        '''
        Calculates the geometric centroid of a ring by averaging the coordinates of its atoms.
        Args:
            ring_atoms (list): A list of dictionaries containing the information of the ring atoms.
        Returns:
            dict: A dictionary with keys 'x', 'y', 'z' giving the centroid coordinates.
        '''

        x_sum = sum(atom['x'] for atom in ring_atoms)
        y_sum = sum(atom['y'] for atom in ring_atoms)
        z_sum = sum(atom['z'] for atom in ring_atoms)
        num_atoms = len(ring_atoms)

        centroid = {
            'x': x_sum / num_atoms,
            'y': y_sum / num_atoms,
            'z': z_sum / num_atoms
        }

        return centroid
    
    def get_ring_normal(self, ring_atoms):
        '''
        Calculates a normal vector to the ring plane using the cross product of two
        vectors defined by the first three ring atoms. The vector is not normalized
        (magnitude depends on bond lengths) but its direction is what matters for
        downstream angle calculations.
        Args:
            ring_atoms (list): A list of dictionaries containing the information of the ring atoms.
        Returns:
            dict: A dictionary with keys 'x', 'y', 'z' giving the normal vector components.
        '''

        # Take 3 atoms to define 2 vectors in the ring plane
        atom1 = ring_atoms[0]
        atom2 = ring_atoms[1]
        atom3 = ring_atoms[2]

        # Two vectors in the ring plane
        v1 = {
            'x': atom2['x'] - atom1['x'],
            'y': atom2['y'] - atom1['y'],
            'z': atom2['z'] - atom1['z']
        }
        v2 = {
            'x': atom3['x'] - atom1['x'],
            'y': atom3['y'] - atom1['y'],
            'z': atom3['z'] - atom1['z']
        }

        # Cross product
        normal = {
            'x': v1['y'] * v2['z'] - v1['z'] * v2['y'],
            'y': v1['z'] * v2['x'] - v1['x'] * v2['z'],
            'z': v1['x'] * v2['y'] - v1['y'] * v2['x']
        }

        return normal

    def get_rings(self, chain):
        '''
        Extracts all aromatic rings of a chain, each with its constituent atoms, centroid and normal vector.
        Args:
            chain (str): The chain ID.
        Returns:
            list: A list of dictionaries, one per ring, with keys
                  'residue_name', 'residue_seq', 'chain', 'atoms', 'centroid', 'normal'.
        '''

        rings = []
        seen = set()

        for atom in self.ring_atoms():
            if atom['chain'] == chain:
                ring_id = (atom['residue_name'], atom['residue_seq'])
                if ring_id not in seen:
                    seen.add(ring_id)
                    ring_atoms = [
                        a for a in self.ring_atoms()
                        if a['residue_name'] == atom['residue_name'] and
                        a['residue_seq'] == atom['residue_seq'] and
                        a['chain'] == chain
                    ]
                    rings.append({
                        'residue_name': atom['residue_name'],
                        'residue_seq': atom['residue_seq'],
                        'chain': chain,
                        'atoms': ring_atoms,
                        'centroid': self.get_ring_centroid(ring_atoms),
                        'normal': self.get_ring_normal(ring_atoms)
                    })

        return rings
    
    def parse_conect(self):
        '''
        Parses the CONECT records from the PDB data.
        Returns:
            dict: A dictionary mapping atom_seq to a list of connected atom_seqs.
        '''

        conect = {}
        for line in self.pdb_data.splitlines():
            if line.startswith('CONECT'):
                atom_seq = line[6:11].strip()
                bonded = [
                    line[11:16].strip(),
                    line[16:21].strip(),
                    line[21:26].strip(),
                    line[26:31].strip()
                ]
                bonded = [b for b in bonded if b]  # remove empty strings
                if atom_seq not in conect:
                    conect[atom_seq] = []
                conect[atom_seq].extend(bonded)

        return conect

    def get_small_molecule_graph(self, small_molecule_id):
        '''
        Builds a connectivity graph for a specific small molecule using CONECT records.
        Args:
            small_molecule_id (str): The residue name of the small molecule.
        Returns:
            dict: A dictionary mapping atom_seq to a list of connected atom_seqs,
                filtered to only atoms belonging to the small molecule.
        '''

        sm_atoms = self.get_small_molecule(small_molecule_id)
        sm_atom_seqs = {atom['atom_seq'] for atom in sm_atoms}

        graph = {atom['atom_seq']: [] for atom in sm_atoms}

        for atom_seq in sm_atom_seqs:
            if atom_seq in self.conect:
                for bonded_seq in self.conect[atom_seq]:
                    if bonded_seq in sm_atom_seqs:  # only keep bonds within the molecule
                        graph[atom_seq].append(bonded_seq)

        return graph
    
    def get_small_molecule_rings(self, graph):
        '''
        Finds all 5 and 6 membered cycles in a connectivity graph.
        Args:
            graph (dict): A dictionary mapping atom_seq to a list of connected atom_seqs.
        Returns:
            list: A list of sets, each containing the atom_seqs of a ring.
        '''

        rings = []

        def dfs(start, current, path, visited):
            for neighbor in graph[current]:
                if neighbor == start and len(path) in {5, 6}:
                    ring = frozenset(path)
                    if ring not in rings:
                        rings.append(ring)
                elif neighbor not in visited and len(path) < 6:
                    visited.add(neighbor)
                    path.append(neighbor)
                    dfs(start, neighbor, path, visited)
                    path.pop()
                    visited.remove(neighbor)

        for start in graph:
            dfs(start, start, [start], {start})

        return rings
    
    def is_planar(self, ring_atoms, tolerance=0.1):
        '''
        Checks if a set of atoms is roughly coplanar by computing the ring normal
        and measuring each atom's deviation from the plane.
        Args:
            ring_atoms (list): A list of atom dicts belonging to the ring.
            tolerance (float): Maximum allowed deviation from the plane in Angstroms.
        Returns:
            bool: True if all atoms are within tolerance of the plane, False otherwise.
        '''

        # use first 3 atoms to define the plane
        a1, a2, a3 = ring_atoms[0], ring_atoms[1], ring_atoms[2]

        v1 = {'x': a2['x']-a1['x'], 'y': a2['y']-a1['y'], 'z': a2['z']-a1['z']}
        v2 = {'x': a3['x']-a1['x'], 'y': a3['y']-a1['y'], 'z': a3['z']-a1['z']}

        # normal vector of the plane
        normal = {
            'x': v1['y']*v2['z'] - v1['z']*v2['y'],
            'y': v1['z']*v2['x'] - v1['x']*v2['z'],
            'z': v1['x']*v2['y'] - v1['y']*v2['x']
        }

        magnitude = (normal['x']**2 + normal['y']**2 + normal['z']**2) ** 0.5

        # check deviation of each remaining atom from the plane
        for atom in ring_atoms[3:]:
            dx = atom['x'] - a1['x']
            dy = atom['y'] - a1['y']
            dz = atom['z'] - a1['z']
            deviation = abs(dx*normal['x'] + dy*normal['y'] + dz*normal['z']) / magnitude
            if deviation > tolerance:
                return False

        return True

    def get_small_molecule_aromatic_rings(self, small_molecule_id):
        '''
        Finds all aromatic rings in a small molecule by detecting planar 5 or 6
        membered cycles using CONECT records.
        Args:
            small_molecule_id (str): The residue name of the small molecule.
        Returns:
            list: A list of dictionaries in the same format as get_rings, each containing:
                residue_name, residue_seq, chain, atoms, centroid, normal.
        '''

        sm_atoms = self.get_small_molecule(small_molecule_id)
        atom_seq_to_atom = {atom['atom_seq']: atom for atom in sm_atoms}

        graph = self.get_small_molecule_graph(small_molecule_id)
        cycles = self.get_small_molecule_rings(graph)

        aromatic_rings = []
        for cycle in cycles:
            ring_atoms = [atom_seq_to_atom[seq] for seq in cycle]

            # filter to only C, N, O atoms
            if not all(atom['atom_symbol'] in {'C', 'N', 'O'} for atom in ring_atoms):
                continue

            if not self.is_planar(ring_atoms):
                continue

            aromatic_rings.append({
                'residue_name': ring_atoms[0]['residue_name'],
                'residue_seq' : ring_atoms[0]['residue_seq'],
                'chain'       : ring_atoms[0]['chain'],
                'atoms'       : ring_atoms,
                'centroid'    : self.get_ring_centroid(ring_atoms),
                'normal'      : self.get_ring_normal(ring_atoms)
            })

        return aromatic_rings

    def hetero_atoms(self):
        '''
        Extracts the hetero atoms (HETATM records) from the parsed PDB data.
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
        Extracts the unique chain IDs of the polypeptide chains from the parsed PDB data.
        Returns:
            set: A set of unique chain IDs.
        '''

        chains = set()
        for atom in self.all_atoms:
            if atom['type'] == 'ATOM':
                chains.add(atom['chain'])

        return chains
    
    def get_chain(self, chain_id):
        '''
        Extracts all ATOM records belonging to a specific chain.
        Args:
            chain_id (str): The ID of the chain to extract atoms from.
        Returns:
            list: A list of dictionaries containing the atoms for the specified chain.
        Raises:
            ValueError: If the chain is not present in the PDB data.
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
        Extracts the unique small molecule residue names from the parsed PDB data.
        Note: water molecules (HOH) are grouped under a single entry.
        Returns:
            set: A set of unique small molecule residue names.
        '''

        small_molecules = set()
        for atom in self.all_atoms:
            if atom['type'] == 'HETATM':
                small_molecules.add(atom['residue_name'])

        return small_molecules
    
    def get_small_molecule(self, small_molecule_id):
        '''
        Extracts all HETATM records belonging to a specific small molecule residue name.
        Note: if small_molecule_id is 'HOH', this returns the atoms of every water molecule together.
        Args:
            small_molecule_id (str): The residue name of the small molecule.
        Returns:
            list: A list of dictionaries containing the atoms for the specified small molecule.
        Raises:
            ValueError: If the small molecule is not present in the PDB data.
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
    
    def metals(self):
        metal_symbols = {'ZN', 'CA', 'MG', 'NA', 'K', 'FE', 'CU', 'MN', 'CO', 'NI', 'CD', 'HG', 'PT', 'AU'}
        metals = []
        for atom in self.all_atoms:
            if atom['type'] == 'HETATM' and atom['atom_symbol'].upper() in metal_symbols:
                metals.append(atom)
        return metals
    
    def get_entity(self, entity_id):
        '''
        Extracts the atoms for a specific entity, which can be either a chain or a small molecule.
        Args:
            entity_id (str): The ID of the chain (e.g. 'A') or small molecule (e.g. 'HOH', 'ATP').
        Returns:
            list: A list of dictionaries containing the atoms for the specified entity.
        Raises:
            ValueError: If the entity is not present in the PDB data.
        '''

        if entity_id in self.chains():
            return self.get_chain(entity_id)
        elif entity_id in self.small_molecules():
            return self.get_small_molecule(entity_id)
        else:
            raise ValueError(f"Entity '{entity_id}' not found in the PDB data.")

    def get_ID(self):
        '''
        Extracts the PDB ID from the PDB data.
        Returns:
            str: The PDB ID if found, otherwise 'Unknown'.
        '''

        for line in self.pdb_data.splitlines():
            if line.startswith('HEADER'):
                return line[62:66].strip()
        return 'Unknown'

    def summary(self):
        '''
        Provides a summary of the protein structure.
        Returns:
            dict: A dictionary containing the summary.
        '''

        print("Generating summary of the protein structure...")

        summary_dict = {
            'num_atoms_total'       : len(self.all_atoms),
            'num_chains'            : len(self.chains()),
            'num_small_molecules'   : len(self.small_molecules()),
            'chain_ids'             : sorted(self.chains()),
            'small_molecules'       : self.small_molecules(),
            'num_backbone_atoms'    : len(self.backbone()),
            'num_side_chain_atoms'  : len(self.side_chain()),
            'num_hetero_atoms'      : len(self.hetero_atoms()),
            'num_ring_atoms'        : len(self.ring_atoms()),
            'rings_per_chain'       : {chain: self.get_rings(chain) for chain in self.chains()}
        }

        return summary_dict
