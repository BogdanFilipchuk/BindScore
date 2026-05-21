import math
import os
import pathlib
from pdb_utils_protein import *


def get_angle_between_normals(normal1, normal2):
    '''
    Calculates the angle in degrees between two normal vectors using the dot product.
    The vectors do not need to be unit-normalized; magnitudes are accounted for.
    Args:
        normal1 (dict): A dictionary with keys 'x', 'y', 'z' for the first normal vector.
        normal2 (dict): A dictionary with keys 'x', 'y', 'z' for the second normal vector.
    Returns:
        float: The angle between the two normal vectors in degrees, in the range [0, 180].
    '''

    dot_product = normal1['x'] * normal2['x'] + normal1['y'] * normal2['y'] + normal1['z'] * normal2['z']

    magnitude1 = (normal1['x']**2 + normal1['y']**2 + normal1['z']**2) ** 0.5
    magnitude2 = (normal2['x']**2 + normal2['y']**2 + normal2['z']**2) ** 0.5

    cos_angle = dot_product / (magnitude1 * magnitude2)
    cos_angle = max(-1, min(1, cos_angle))  # clamp to avoid floating point errors in acos

    angle = math.acos(cos_angle) * (180 / math.pi)

    return angle

def calculate_distance(atom1, atom2):
    '''
    Calculates the Euclidean distance between two atoms from their 3D coordinates.
    Args:
        atom1 (dict): A dictionary with keys 'x', 'y', 'z' for the first atom.
        atom2 (dict): A dictionary with keys 'x', 'y', 'z' for the second atom.
    Returns:
        float: The distance in Angstroms.
    '''
    dx = atom1['x'] - atom2['x']
    dy = atom1['y'] - atom2['y']
    dz = atom1['z'] - atom2['z']
    return (dx**2 + dy**2 + dz**2) ** 0.5

class Interaction:
    def __init__(self, protein, chain1, chain2, threshold=5.0):
        '''
        Initializes the Interaction object between two entities of a protein structure.
        Each entity can be either a polypeptide chain or a small molecule.

        On initialization, atom pairs within `threshold` Angstroms are identified and
        categorized into interaction types (hydrogen bonds, pi-pi stacking). Aromatic
        rings of each chain entity are pre-cached for efficient pi-pi detection.
        Args:
            protein (Protein_Structure): The protein structure object.
            chain1 (str): The ID of the first entity (chain ID, e.g. 'A', or small
                        molecule residue name, e.g. 'ATP').
            chain2 (str): The ID of the second entity.
            threshold (float): The distance threshold in Angstroms for identifying
                            possible atom-atom interactions. Defaults to 5.0.
        Attributes:
            protein      - The Protein_Structure object.
            chain1       - The ID of the first entity.
            chain2       - The ID of the second entity.
            entity1      - The atoms of the first entity (list of dicts).
            entity2      - The atoms of the second entity (list of dicts).
            rings1       - Aromatic rings in entity1, keyed by residue_seq (empty for small molecules).
            rings2       - Aromatic rings in entity2, keyed by residue_seq (empty for small molecules).
            interactions - List of detected atom-pair interactions with distance and type.
        Functions:
            find_interactions(threshold)         - Identifies atom pairs within threshold and categorizes them.
            categorize_interaction(atom1, atom2) - Categorizes a single atom pair into an interaction type.
        '''

        self.protein = protein
        self.chain1 = chain1
        self.chain2 = chain2
        self.entity1 = protein.get_entity(chain1)
        self.entity2 = protein.get_entity(chain2)

        # cache rings for fast lookup during categorize_interaction
        if chain1 in self.protein.chains():
            self.rings1 = {r['residue_seq']: r for r in self.protein.get_rings(self.chain1)}
        else:
            self.rings1 = {}

        if chain2 in self.protein.chains():
            self.rings2 = {r['residue_seq']: r for r in self.protein.get_rings(self.chain2)}
        else:
            self.rings2 = {}

        self.interactions = self.find_interactions(threshold)

    def find_interactions(self, threshold):
        '''
        Iterates over all atom pairs (one from each entity) and identifies those within
        the distance threshold. Each candidate pair is then passed to categorize_interaction
        to determine the interaction type.
        Args:
            threshold (float): The distance threshold in Angstroms.
        Returns:
            list: A list of dictionaries, one per atom pair within threshold, with keys
                'atom1', 'atom2', 'distance' (in Angstroms), and 'type' (str or None).
        '''

        interactions = []
        for atom1 in self.entity1:
            for atom2 in self.entity2:
                distance = calculate_distance(atom1, atom2)
                if distance <= threshold:
                    interaction_type = self.categorize_interaction(atom1, atom2)
                    interactions.append({
                        'atom1': atom1,
                        'atom2': atom2,
                        'distance': distance,
                        'type': interaction_type
                    })
        return interactions
    
    def categorize_interaction(self, atom1, atom2):
        '''
        Categorizes the interaction between two atoms based on the entity types,
        atom types, and (for pi-pi stacking) ring geometry.

        Detection rules:
            - Hydrogen bond (protein-protein): both atoms are N/O/F; one must be a known
            donor (or backbone N) and the other a known acceptor (or backbone O).
            - Hydrogen bond (involving a small molecule): the small molecule N/O/F is
            treated as both potential donor and acceptor; the protein side must be
            capable of either role.
            - Pi-pi stacking (protein-protein): triggered when both atoms are the CG
            of an aromatic residue (PHE, TYR, TRP, HIS). Ring centroid distance must
            be <= 5.5 A and the angle between ring normals must fall in the parallel
            (<=30 deg or >=150 deg) or T-shaped (60-120 deg) range.
        Args:
            atom1 (dict): A dictionary containing the information of the first atom.
            atom2 (dict): A dictionary containing the information of the second atom.
        Returns:
            str or None: 'hydrogen_bond', 'pi-pi stacking', or None if no interaction is detected.
        '''

        h_bond_donor = {('ARG', 'NE'), ('ARG', 'NH1'), ('ARG', 'NH2'), ('ASN', 'ND2'), 
                        ('GLN', 'NE2'), ('HIS', 'ND1'), ('HIS', 'NE2'), ('HYP', 'OD'), 
                        ('LYS', 'NZ'), ('SER', 'OG'), ('THR', 'OG1'), ('TRP', 'NE1'), 
                        ('TYR', 'OH')}
        h_bond_acceptor = {('ASN', 'OD1'), ('ASP', 'OD1'), ('ASP', 'OD2'), ('GLU', 'OE1'), 
                           ('GLU', 'OE2'), ('GLN', 'OE1'), ('HIS', 'ND1'), ('HIS', 'NE2'), 
                           ('SER', 'OG'), ('THR', 'OG1'), ('TYR', 'OH')}
        listNO = {'N', 'O', 'F'}
        aromatic_residues = {'PHE', 'TYR', 'TRP', 'HIS'}

        if self.chain1 in self.protein.chains():
            if self.chain2 in self.protein.chains():
                # protein-protein interaction

                if (atom1['atom_symbol'] in listNO) and (atom2['atom_symbol'] in listNO):
                    if ((atom1['residue_name'], atom1['atom_name']) in h_bond_donor or atom1['atom_name'] =='N') and ((atom2['residue_name'], atom2['atom_name']) in h_bond_acceptor or atom2['atom_name'] =='O'):
                        return 'hydrogen_bond'
                    elif ((atom2['residue_name'], atom2['atom_name']) in h_bond_donor or atom2['atom_name'] =='N') and ((atom1['residue_name'], atom1['atom_name']) in h_bond_acceptor or atom1['atom_name'] =='O'):
                        return 'hydrogen_bond'
                elif (atom1['atom_name'] == 'CG' and atom2['atom_name'] == 'CG' and atom1['residue_name'] in aromatic_residues and atom2['residue_name'] in aromatic_residues):

                    ring1 = self.rings1.get(atom1['residue_seq'])
                    ring2 = self.rings2.get(atom2['residue_seq'])

                    if ring1 is not None and ring2 is not None:
                        # centroid-centroid distance
                        dx = ring1['centroid']['x'] - ring2['centroid']['x']
                        dy = ring1['centroid']['y'] - ring2['centroid']['y']
                        dz = ring1['centroid']['z'] - ring2['centroid']['z']
                        distance = (dx**2 + dy**2 + dz**2) ** 0.5

                        if distance <= 5.5:
                            angle = get_angle_between_normals(ring1['normal'], ring2['normal'])
                            if angle <= 30 or angle >= 150 or (60 <= angle <= 120):
                                return 'pi-pi stacking'
                
            elif self.chain2 in self.protein.small_molecules():
                # protein-small molecule interaction

                if (atom1['atom_symbol'] in listNO) and (atom2['atom_symbol'] in listNO):
                    if ((atom1['residue_name'], atom1['atom_name']) in h_bond_donor or atom1['atom_name'] =='N') or ((atom1['residue_name'], atom1['atom_name']) in h_bond_acceptor or atom1['atom_name'] =='O'):
                        return 'hydrogen_bond'
                    
        elif self.chain1 in self.protein.small_molecules():
            if self.chain2 in self.protein.chains():
                # small molecule-protein interaction

                if (atom1['atom_symbol'] in listNO) and (atom2['atom_symbol'] in listNO):
                    if ((atom2['residue_name'], atom2['atom_name']) in h_bond_donor or atom2['atom_name'] =='N') or ((atom2['residue_name'], atom2['atom_name']) in h_bond_acceptor or atom2['atom_name'] =='O'):
                        return 'hydrogen_bond'
                    
            elif self.chain2 in self.protein.small_molecules():
                # small molecule-small molecule interaction

                if (atom1['atom_symbol'] in listNO) and (atom2['atom_symbol'] in listNO):
                    return 'hydrogen_bond'

if __name__ == '__main__':
    import sys

    pdb_id    = sys.argv[1] if len(sys.argv) > 1 else '1brs'
    threshold = 5.0

    pdb_data = fetch_pdb_data(pdb_id)
    protein  = Protein_Structure(pdb_data)

    print(f"\nStructure       : {pdb_id}")
    print(f"Total atoms     : {len(protein.all_atoms)}")
    print(f"Chains          : {sorted(protein.chains())}")
    print(f"Small molecules : {sorted(protein.small_molecules())}")
    print(f"Ring atoms      : {len(protein.ring_atoms())}")

    chains = sorted(protein.chains())

    print(f"\nChain-chain interactions (threshold = {threshold} A)")
    print(f"{'pair':>8} | {'pairs<thr':>10} | {'h-bonds':>8} | {'pi-pi':>6}")
    print("-" * 50)

    summary = []
    for i, c1 in enumerate(chains):
        for c2 in chains[i+1:]:
            inter = Interaction(protein, c1, c2, threshold=threshold)
            hb    = sum(1 for x in inter.interactions if x['type'] == 'hydrogen_bond')
            pp    = sum(1 for x in inter.interactions if x['type'] == 'pi-pi stacking')
            total = len(inter.interactions)
            print(f"{c1}-{c2:>4} | {total:>10} | {hb:>8} | {pp:>6}")
            summary.append((c1, c2, hb, pp, inter))

    # detail the most interactive pair
    if summary:
        top = max(summary, key=lambda s: s[2] + s[3])
        c1, c2, hb, pp, inter = top
        if hb + pp > 0:
            print(f"\nTop pair: {c1}-{c2}  (categorized interactions only)")
            for x in inter.interactions:
                if x['type'] is not None:
                    a1, a2 = x['atom1'], x['atom2']
                    print(f"  {a1['residue_name']} {a1['residue_seq']:>4} {a1['atom_name']:>4} "
                          f"<-> "
                          f"{a2['residue_name']} {a2['residue_seq']:>4} {a2['atom_name']:>4}  "
                          f"d = {x['distance']:.2f} A  type = {x['type']}")