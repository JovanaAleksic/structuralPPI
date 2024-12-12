import numpy as np
from Bio.PDB import *
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings
from scipy.spatial.distance import pdist, squareform
import pandas as pd

class PDockQ:
    def __init__(self):
        # Model coefficients from the paper
        self.weights = {
            'bias': -6.262,
            'log_interface_score': 0.533,
            'normalized_interface_energy': -0.492,
            'log_interface_size': 0.463,
            'log_num_interface_residues': 0.367,
            'ratio_polar_interface': -0.184,
            'normalized_surface_comp': 0.206
        }

    def get_accessible_residues(self, chain, neighbor_radius=10.0):
        """
        Identify surface-accessible residues
        """
        atoms = [atom for residue in chain for atom in residue 
                if atom.get_name() in ['CB', 'CA']]
        ns = NeighborSearch(atoms)
        
        accessible_residues = set()
        for atom in atoms:
            neighbors = ns.search(atom.get_coord(), neighbor_radius)
            if len(neighbors) < 15:  # Fewer neighbors indicates surface exposure
                accessible_residues.add(atom.get_parent())
        
        return accessible_residues

    def get_interface_residues(self, chain1, chain2, cutoff=10.0):
        """
        Identify interface residues between two chains
        Returns interface residues and their contacts
        """
        interface_residues1 = set()
        interface_residues2 = set()
        contacts = []

        # Get CB atoms (CA for GLY)
        atoms1 = []
        atoms2 = []
        for residue in chain1:
            if 'CB' in residue:
                atoms1.append(residue['CB'])
            elif 'CA' in residue:
                atoms1.append(residue['CA'])
                
        for residue in chain2:
            if 'CB' in residue:
                atoms2.append(residue['CB'])
            elif 'CA' in residue:
                atoms2.append(residue['CA'])

        ns = NeighborSearch(atoms2)
        for atom1 in atoms1:
            res1 = atom1.get_parent()
            neighbors = ns.search(atom1.get_coord(), cutoff)
            
            if neighbors:
                interface_residues1.add(res1)
                for neighbor in neighbors:
                    res2 = neighbor.get_parent()
                    interface_residues2.add(res2)
                    contacts.append((res1, res2, atom1 - neighbor))

        return interface_residues1, interface_residues2, contacts

    def calculate_interface_energy(self, contacts):
        """Calculate interface energy using a knowledge-based potential"""
        energy = 0.0
        for res1, res2, distance in contacts:
            # Distance-dependent energy term
            distance_term = -1.0 / (1.0 + np.exp(distance - 8.0))
            
            # Simple residue-type based term
            res1_name = res1.get_resname()
            res2_name = res2.get_resname()
            
            # Group residues into hydrophobic, polar, and charged
            hydrophobic = set(['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP'])
            polar = set(['SER', 'THR', 'ASN', 'GLN', 'TYR'])
            charged = set(['LYS', 'ARG', 'HIS', 'ASP', 'GLU'])
            
            # Add residue-type specific energy term
            if res1_name in hydrophobic and res2_name in hydrophobic:
                energy += -0.5 * distance_term
            elif res1_name in charged and res2_name in charged:
                energy += -0.3 * distance_term
            elif res1_name in polar and res2_name in polar:
                energy += -0.2 * distance_term
            else:
                energy += -0.1 * distance_term

        return energy

    def calculate_surface_composition(self, residues):
        """Calculate surface amino acid composition"""
        aa_groups = {
            'polar': ['SER', 'THR', 'ASN', 'GLN', 'TYR'],
            'charged': ['LYS', 'ARG', 'HIS', 'ASP', 'GLU'],
            'hydrophobic': ['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP']
        }
        
        composition = {group: 0 for group in aa_groups}
        total = 0
        
        for res in residues:
            resname = res.get_resname()
            for group, aas in aa_groups.items():
                if resname in aas:
                    composition[group] += 1
                    total += 1
        
        if total > 0:
            for group in composition:
                composition[group] /= total
                
        return composition

    def calculate_pdockq(self, cif_file):
        """
        Calculate pDockQ score for a protein complex
        
        Parameters:
        cif_file (str): Path to the CIF file
        
        Returns:
        float: pDockQ score between 0 and 1, dict: feature values
        """
        # Suppress Bio.PDB warnings
        warnings.filterwarnings('ignore', category=PDBConstructionWarning)
        
        # Parse structure
        parser = MMCIFParser()
        structure = parser.get_structure('complex', cif_file)
        
        # Get chains
        chains = [chain for chain in structure[0]]
        if len(chains) < 2:
            raise ValueError("Structure must contain at least 2 chains for interface analysis")

        # Calculate features across all interfaces
        all_interface_sizes = []
        all_interface_energies = []
        all_interface_residues = set()
        total_contacts = 0
        surface_compositions = []

        # Analyze all chain pairs
        for i in range(len(chains)):
            for j in range(i+1, len(chains)):
                interface_res1, interface_res2, contacts = self.get_interface_residues(chains[i], chains[j])
                
                if contacts:
                    # Interface size
                    interface_size = len(interface_res1) + len(interface_res2)
                    all_interface_sizes.append(interface_size)
                    
                    # Interface energy
                    energy = self.calculate_interface_energy(contacts)
                    all_interface_energies.append(energy)
                    
                    # Keep track of all interface residues
                    all_interface_residues.update(interface_res1)
                    all_interface_residues.update(interface_res2)
                    
                    # Count contacts
                    total_contacts += len(contacts)
                    
                    # Surface composition
                    surface_comp = self.calculate_surface_composition(interface_res1.union(interface_res2))
                    surface_compositions.append(surface_comp)

        if not all_interface_sizes:
            return 0.0, {}

        # Calculate features
        features = {}
        features['log_interface_size'] = np.log(np.mean(all_interface_sizes))
        features['log_interface_score'] = np.log(abs(np.mean(all_interface_energies)))
        features['normalized_interface_energy'] = np.mean(all_interface_energies) / len(all_interface_residues)
        features['log_num_interface_residues'] = np.log(len(all_interface_residues))
        
        # Calculate ratio of polar interface residues
        if surface_compositions:
            avg_composition = {k: np.mean([comp[k] for comp in surface_compositions]) 
                             for k in surface_compositions[0]}
            features['ratio_polar_interface'] = avg_composition['polar']
            features['normalized_surface_comp'] = (avg_composition['charged'] + 
                                                 avg_composition['polar']) / (avg_composition['hydrophobic'] + 0.1)

        # Calculate pDockQ score
        score = self.weights['bias']
        for feature, value in features.items():
            score += self.weights[feature] * value
        
        # Convert to probability
        pdockq = 1.0 / (1.0 + np.exp(-score))
        
        return min(max(pdockq, 0.0), 1.0), features

# Example usage
if __name__ == "__main__":

    pdockq_calculator = PDockQ()
    score, features = pdockq_calculator.calculate_pdockq(sys.argv[1])
    print(f"\npDockQ score: {score:.3f}")
    print("\nFeature values:")
    for feature, value in features.items():
        print(f"{feature}: {value:.3f}")
