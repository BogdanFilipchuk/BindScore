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
    cos_angle = max(-1, min(1, cos_angle))

    angle = math.acos(cos_angle) * (180 / math.pi)

    return angle

def calculate_distance(atom1, atom2):
    '''
    Calculates the distance between two atoms from their 3D coordinates.
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

def ring_pair_key(ring1, ring2):
    '''
    Builds an order-independent unique key for a pair of rings, using the lowest
    atom_seq of each ring as its identifier. Suitable for deduplication when
    multiple atom-pair iterations refer to the same ring stacking event.
    Args:
        ring1, ring2 (dict): Ring dictionaries with an 'atoms' list.
    Returns:
        tuple: A sorted tuple of two atom_seq strings.
    '''
    id1 = min(a['atom_seq'] for a in ring1['atoms'])
    id2 = min(a['atom_seq'] for a in ring2['atoms'])
    return tuple(sorted([id1, id2]))

class Interaction:
    def __init__(self, protein, chain1, chain2, threshold=5.0):
        '''
        Initializes the Interaction object between two entities of a protein
        structure. Each entity can be either a polypeptide chain or a small
        molecule (including metals and water).

        Aromatic rings for each entity are pre-cached as dictionaries keyed by
        residue_seq (protein chains) or as lists (small molecules) to avoid
        redundant computation during the atom-pair loop. All interactions within
        the distance threshold are detected and stored on initialization.
        Args:
            protein (Protein_Structure): The parsed protein structure object.
            chain1 (str): ID of the first entity — chain letter (e.g. 'A') or
                small molecule residue name (e.g. 'ATP', 'ZN').
            chain2 (str): ID of the second entity, same format as chain1.
            threshold (float): Distance threshold in Angstroms for detecting
                atom-atom contacts. Defaults to 5.0.
        Attributes:
            protein      - The Protein_Structure object.
            chain1       - ID of the first entity.
            chain2       - ID of the second entity.
            entity1      - Atoms of the first entity (list of dicts).
            entity2      - Atoms of the second entity (list of dicts).
            rings1       - Aromatic rings of entity1: dict keyed by residue_seq
                        for protein chains, list for small molecules, empty
                        dict if no rings are present.
            rings2       - Same as rings1 for entity2.
            interactions - List of all detected atom-pair interactions within
                        threshold, each with keys atom1, atom2, distance, type.
        '''

        self.protein = protein
        self.chain1 = chain1
        self.chain2 = chain2
        self.entity1 = protein.get_entity(chain1)
        self.entity2 = protein.get_entity(chain2)

        # cache rings for fast lookup during categorize_interaction
        if chain1 in self.protein.chains():
            self.rings1 = {r['residue_seq']: r for r in self.protein.get_rings(self.chain1)}
        elif chain1 in self.protein.small_molecules():
            self.rings1 = list(self.protein.get_small_molecule_aromatic_rings(self.chain1))
        else:
            self.rings1 = {}
            
        if chain2 in self.protein.chains():
            self.rings2 = {r['residue_seq']: r for r in self.protein.get_rings(self.chain2)}
        elif chain2 in self.protein.small_molecules():
            self.rings2 = list(self.protein.get_small_molecule_aromatic_rings(self.chain2))
        else:
            self.rings2 = {}

        self.interactions = self.find_interactions(threshold)

    def find_interactions(self, threshold):
        '''
        Iterates over all atom pairs (one from each entity) and identifies those
        within the distance threshold. Each candidate pair is passed to
        categorize_interaction to determine the interaction type. All pairs within
        threshold are recorded, including those with type None (unclassified contacts).

        Deduplication sets for hydrophobic contacts, salt bridges, and pi-pi
        stacking are initialized here and passed through to categorize_interaction
        so that multi-atom-pair interactions are reported only once per residue
        pair or ring pair.
        Args:
            threshold (float): Distance threshold in Angstroms. Only atom pairs
                within this distance are considered.
        Returns:
            list: A list of dictionaries, one per atom pair within threshold, with
                keys 'atom1', 'atom2', 'distance' (float, Angstroms), and 'type'
                (str or None).
        '''

        interactions = []
        seen_hydrophobic = set()
        seen_salt_bridge = set()
        seen_pi_pi = set()

        for atom1 in self.entity1:
            for atom2 in self.entity2:
                distance = calculate_distance(atom1, atom2)
                if distance <= threshold:
                    interaction_type = self.categorize_interaction(atom1, atom2, seen_hydrophobic, seen_salt_bridge, seen_pi_pi)
                    interactions.append({
                        'atom1': atom1,
                        'atom2': atom2,
                        'distance': distance,
                        'type': interaction_type
                    })
        return interactions
    
    def categorize_interaction(self, atom1, atom2, seen_hydrophobic=None, seen_salt_bridge=None, seen_pi_pi=None):
        '''
        Categorizes the interaction between two atoms based on the entity types
        (chain or small molecule), atom identity, geometry, and chemical rules.
        Checks are applied in priority order so that each atom pair is assigned
        the most specific interaction type that applies.

        Detection rules (in order of priority):
            1. Disulfide bond: both atoms are CYS SG, separated by 1.8-2.5 A.
            2. Metal coordination: one atom belongs to a metal residue (ZN, CA, MG,
            NA, K, FE, CU, MN, CO, NI, CD, HG, PT, AU); the other is a known
            coordinating atom (HIS ND1/NE2, CYS SG, ASP OD1/OD2, GLU OE1/OE2,
            SER OG, THR OG1, MET SD), a backbone O/OXT, a water oxygen, or any S.
            3. Salt bridge (protein-protein only): one atom is from a positively
            charged side chain (ARG NH1/NH2/NE or LYS NZ), the other from a
            negatively charged side chain (ASP OD1/OD2 or GLU OE1/OE2).
            Reported once per residue pair via seen_salt_bridge.
            4. Hydrogen bond (protein-protein): both atoms are N/O/F; one must be a
            known donor (or backbone N, excluding PRO/HYP) and the other a known
            acceptor (or backbone O/OXT). Pairing is checked in both directions.
            5. Hydrogen bond (protein-SM or SM-protein): both atoms are N/O/F; the
            protein atom must be capable of either donor or acceptor role. The SM
            atom is treated as both.
            6. Hydrogen bond (SM-SM): any N/O/F pair qualifies.
            7. Pi-pi stacking (all entity combinations): for protein chains, triggered
            by CG of PHE/TYR/TRP/HIS; for small molecules, by membership in a
            detected aromatic ring. Ring centroid distance must be <= 5.5 A and
            the normal-to-normal angle must be in the parallel (<=30 or >=150 deg)
            or T-shaped (60-120 deg) range. Reported once per ring pair via
            seen_pi_pi.
            8. Hydrophobic contact (protein-protein only): both residues in
            {LEU, ILE, VAL, PHE, TYR, TRP, MET, ALA, PRO}, both atoms are
            non-backbone carbons (excluding C and CA). Reported once per residue
            pair via seen_hydrophobic.
            9. Halogen bond (all entity combinations): distance < 4.0 A; one atom is
            a halogen (F, CL, BR, I) and the other is polar (known donor/acceptor,
            backbone N/O/OXT, or S).
            10. Dipole-dipole (all entity combinations): distance < 4.0 A; both atoms
            are polar (known donor/acceptor, backbone N/O/OXT, or S) and neither
            is a halogen. Catches polar contacts that do not satisfy stricter
            h-bond donor/acceptor pairing (e.g. carbonyl-carbonyl, donor-donor,
            sulfur-mediated contacts).

        Args:
            atom1 (dict): Parsed atom dict with keys type, atom_seq, atom_name,
                residue_name, chain, residue_seq, x, y, z, atom_symbol.
            atom2 (dict): Same format as atom1.
            seen_hydrophobic (set or None): Residue-seq pairs already reported as
                hydrophobic_contact; further pairs for the same residues return None.
                Mutated in place on first detection.
            seen_salt_bridge (set or None): Residue-seq pairs already reported as
                salt_bridge; further pairs for the same residues return None.
                Mutated in place on first detection.
            seen_pi_pi (set or None): Ring-pair keys already reported as
                pi-pi_stacking; further pairs for the same rings return None.
                Mutated in place on first detection.
        Returns:
            str or None: One of 'disulfide_bond', 'metal_coordination', 'salt_bridge',
                'hydrogen_bond', 'pi-pi_stacking', 'hydrophobic_contact',
                'halogen_bond', 'dipole-dipole', or None if no interaction is detected.
        '''

        if atom1['residue_seq'] == atom2['residue_seq'] and atom1['chain'] == atom2['chain']:
            return None

        # Hydrogen Bonds
        h_bond_donor = {('ARG', 'NE'), ('ARG', 'NH1'), ('ARG', 'NH2'), ('ASN', 'ND2'), 
                        ('GLN', 'NE2'), ('HIS', 'ND1'), ('HIS', 'NE2'), ('HYP', 'OD'), 
                        ('LYS', 'NZ'), ('SER', 'OG'), ('THR', 'OG1'), ('TRP', 'NE1'), 
                        ('TYR', 'OH')}
        h_bond_acceptor = {('ASN', 'OD1'), ('ASP', 'OD1'), ('ASP', 'OD2'), ('GLU', 'OE1'), 
                           ('GLU', 'OE2'), ('GLN', 'OE1'), ('HIS', 'ND1'), ('HIS', 'NE2'), 
                           ('SER', 'OG'), ('THR', 'OG1'), ('TYR', 'OH')}
        listNO = {'N', 'O', 'F'}
        
        # Hydrophobic Contacts
        hydrophobic_residues = {'LEU', 'ILE', 'VAL', 'PHE', 'TYR', 'TRP', 'MET', 'ALA', 'PRO'}
        
        # Pi-Pi Stacking
        aromatic_residues = {'PHE', 'TYR', 'TRP', 'HIS'}

        # Salt Bridges
        pos_charged = {('ARG', 'NH1'), ('ARG', 'NH2'), ('ARG', 'NE'), ('LYS', 'NZ')}
        neg_charged = {('ASP', 'OD1'), ('ASP', 'OD2'), ('GLU', 'OE1'), ('GLU', 'OE2')}

        # Halogens
        halogens = {'F', 'CL', 'BR', 'I'}

        # Metals
        coordinating_atoms = {('HIS', 'ND1'), ('HIS', 'NE2'), ('CYS', 'SG'), ('ASP', 'OD1'), 
                              ('ASP', 'OD2'), ('GLU', 'OE1'), ('GLU', 'OE2'), ('SER', 'OG'), 
                              ('THR', 'OG1'), ('MET', 'SD')}
        metal_symbols = {'ZN', 'CA', 'MG', 'NA', 'K', 'FE', 'CU', 'MN', 'CO', 'NI', 'CD', 'HG', 'PT', 'AU'}


        def rings_are_stacking(ring1, ring2):
            """Return True if two rings satisfy centroid distance and angle criteria."""
            
            dx = ring1['centroid']['x'] - ring2['centroid']['x']
            dy = ring1['centroid']['y'] - ring2['centroid']['y']
            dz = ring1['centroid']['z'] - ring2['centroid']['z']
            dist = (dx**2 + dy**2 + dz**2) ** 0.5
            if dist > 5.5:
                return False
            angle = get_angle_between_normals(ring1['normal'], ring2['normal'])
            return angle <= 30 or angle >= 150 or (60 <= angle <= 120)

        def find_ring_for_atom(atom, rings):
            rings_iter = rings.values() if isinstance(rings, dict) else rings
            for ring in rings_iter:
                for ratom in ring['atoms']:
                    if ratom['atom_seq'] == atom['atom_seq']:
                        return ring
            return None

        res_at1 = (atom1['residue_name'], atom1['atom_name'])
        res_at2 = (atom2['residue_name'], atom2['atom_name'])

        # Disulfide Bonds
        if (atom1['atom_name'] == 'SG' and atom2['atom_name'] == 'SG' and atom1['residue_name'] == 'CYS' and atom2['residue_name'] == 'CYS'):
            distance = calculate_distance(atom1, atom2)
            if 1.8 <= distance <= 2.5:
                return 'disulfide_bond'
            
        # Metal Coordination
        if atom1['atom_symbol'].upper() in metal_symbols:
            if res_at2 in coordinating_atoms or atom2['atom_name'] in {'O', 'OXT'} or atom2['residue_name'] == 'HOH' or atom2['atom_symbol'] == 'S':
                return 'metal_coordination'

        if atom2['atom_symbol'].upper() in metal_symbols:
            if res_at1 in coordinating_atoms or atom1['atom_name'] in {'O', 'OXT'} or atom1['residue_name'] == 'HOH' or atom1['atom_symbol'] == 'S':
                return 'metal_coordination'

        if self.chain1 in self.protein.chains():
            if self.chain2 in self.protein.chains():
                # protein-protein interaction

                # Salt Bridges
                if (res_at1 in pos_charged and res_at2 in neg_charged) or (res_at1 in neg_charged and res_at2 in pos_charged):

                        if seen_salt_bridge is not None:
                            res_key = (atom1['residue_seq'], atom2['residue_seq'])
                            if res_key in seen_salt_bridge:
                                return None
                            seen_salt_bridge.add(res_key)
                                
                        return 'salt_bridge'

                # Hydrogen Bonds
                elif (atom1['atom_symbol'] in listNO) and (atom2['atom_symbol'] in listNO):
                    if (res_at1 in h_bond_donor or (atom1['atom_name'] =='N' and atom1['residue_name'] not in {'PRO', 'HYP'})) and (res_at2 in h_bond_acceptor or atom2['atom_name'] =='O' or atom2['atom_name'] =='OXT'):
                        return 'hydrogen_bond'
                    elif (res_at2 in h_bond_donor or (atom2['atom_name'] =='N' and atom2['residue_name'] not in {'PRO', 'HYP'})) and (res_at1 in h_bond_acceptor or atom1['atom_name'] =='O' or atom1['atom_name'] =='OXT'):
                        return 'hydrogen_bond'
                
                # Pi-Pi Stacking                    
                elif (atom1['atom_name'] == 'CG' and atom2['atom_name'] == 'CG' and atom1['residue_name'] in aromatic_residues and atom2['residue_name'] in aromatic_residues):

                    ring1 = self.rings1.get(atom1['residue_seq'])
                    ring2 = self.rings2.get(atom2['residue_seq'])

                    if ring1 and ring2 and rings_are_stacking(ring1, ring2):
                        
                        if seen_pi_pi is not None:
                            key = ring_pair_key(ring1, ring2)
                            if key in seen_pi_pi:
                                return None
                            seen_pi_pi.add(key)

                        return 'pi-pi_stacking'
                
                # Hydrophobic Contact
                elif atom1['residue_name'] in hydrophobic_residues and atom2['residue_name'] in hydrophobic_residues:
                    if atom1['atom_symbol'] == 'C' and atom2['atom_symbol'] == 'C':
                        if atom1['atom_name'] not in {'C', 'CA'} and atom2['atom_name'] not in {'C', 'CA'}:
                            
                            if seen_hydrophobic is not None:
                                res_key = (atom1['residue_seq'], atom2['residue_seq'])
                                if res_key in seen_hydrophobic:
                                    return None
                                seen_hydrophobic.add(res_key)
                                
                            return 'hydrophobic_contact'
                
            elif self.chain2 in self.protein.small_molecules():
                # protein-small molecule interaction

                # Hydrogen Bonds
                if (atom1['atom_symbol'] in listNO) and (atom2['atom_symbol'] in listNO):
                    if (res_at1 in h_bond_donor or (atom1['atom_name'] =='N' and atom1['residue_name'] not in {'PRO', 'HYP'})) or (res_at1 in h_bond_acceptor or atom1['atom_name'] =='O' or atom1['atom_name'] =='OXT'):
                        return 'hydrogen_bond'
                    
                # Pi-Pi Stacking
                elif atom1['atom_name'] == 'CG' and atom1['residue_name'] in aromatic_residues:
                    
                    ring1 = self.rings1.get(atom1['residue_seq'])     
                    ring2 = find_ring_for_atom(atom2, self.rings2)     
                    
                    if ring1 and ring2 and rings_are_stacking(ring1, ring2):

                        if seen_pi_pi is not None:
                            key = (min(a['atom_seq'] for a in ring1['atoms']), atom2['residue_seq'])
                            if key in seen_pi_pi:
                                return None
                            seen_pi_pi.add(key)

                        return 'pi-pi_stacking'
            
        elif self.chain1 in self.protein.small_molecules():
            if self.chain2 in self.protein.chains():
                # small molecule-protein interaction

                # Hydrogen Bonds
                if (atom1['atom_symbol'] in listNO) and (atom2['atom_symbol'] in listNO):
                    if (res_at2 in h_bond_donor or (atom2['atom_name'] =='N' and atom2['residue_name'] not in {'PRO', 'HYP'})) or (res_at2 in h_bond_acceptor or atom2['atom_name'] =='O' or atom2['atom_name'] =='OXT'):
                        return 'hydrogen_bond'
                
                # Pi-Pi Stacking
                elif atom2['atom_name'] == 'CG' and atom2['residue_name'] in aromatic_residues:
                    
                    ring1 = find_ring_for_atom(atom1, self.rings1) 
                    ring2 = self.rings2.get(atom2['residue_seq']) 
                    
                    if ring1 and ring2 and rings_are_stacking(ring1, ring2):
                        
                        if seen_pi_pi is not None:
                            key = (min(a['atom_seq'] for a in ring2['atoms']), atom1['residue_seq'])
                            if key in seen_pi_pi:
                                return None
                            seen_pi_pi.add(key)
                        
                        return 'pi-pi_stacking'
                    
            elif self.chain2 in self.protein.small_molecules():
                # small molecule-small molecule interaction

                # Hydrogen Bonds
                if (atom1['atom_symbol'] in listNO) and (atom2['atom_symbol'] in listNO):
                    return 'hydrogen_bond'
                
                # Pi-Pi Stacking
                ring1 = find_ring_for_atom(atom1, self.rings1)
                ring2 = find_ring_for_atom(atom2, self.rings2)
                
                if ring1 and ring2 and rings_are_stacking(ring1, ring2):
                    
                    if seen_pi_pi is not None:
                            key = ring_pair_key(ring1, ring2)
                            if key in seen_pi_pi:
                                return None
                            seen_pi_pi.add(key)
                    
                    return 'pi-pi_stacking'
        
        # Halogen Bonds
        if calculate_distance(atom1, atom2) < 4.0 and ((atom1['atom_symbol'] in halogens and (res_at2 in h_bond_donor or res_at2 in h_bond_acceptor or atom2['atom_name'] in {'N', 'O', 'OXT'} or atom2['atom_symbol'] == 'S')) or (atom2['atom_symbol'] in halogens and (res_at1 in h_bond_donor or res_at1 in h_bond_acceptor or atom1['atom_name'] in {'N', 'O', 'OXT'} or atom1['atom_symbol'] == 'S'))):
                return 'halogen_bond'

        # Dipole-Dipole
        if calculate_distance(atom1, atom2) < 4.0 and (res_at1 in h_bond_donor or res_at1 in h_bond_acceptor or atom1['atom_name'] in {'N', 'O', 'OXT'} or atom1['atom_symbol'] == 'S') and (res_at2 in h_bond_donor or res_at2 in h_bond_acceptor or atom2['atom_name'] in {'N', 'O', 'OXT'} or atom2['atom_symbol'] == 'S'):
                return 'dipole-dipole'
