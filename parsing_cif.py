from Bio.PDB import *
import numpy as np
from Bio.PDB.SASA import ShrakeRupley

def identify_interface_regions(chain1, chain2, distance_cutoff=8.0):
    """
    Identify continuous interface regions between two chains.
    Also returns detailed residue-level interaction data.
    """
    # Find interacting residues with detailed information
    interacting_res1 = {}  # residue_number: (resname, [interacting_partners])
    interacting_res2 = {}
    
    atoms2 = list(chain2.get_atoms())
    ns = NeighborSearch(atoms2)
    
    # Find all interacting residues and their partners
    for atom1 in chain1.get_atoms():
        res1 = atom1.get_parent()
        res1_num = res1.id[1]
        neighbors = ns.search(atom1.get_coord(), distance_cutoff)
        
        if neighbors:
            if res1_num not in interacting_res1:
                interacting_res1[res1_num] = (res1.get_resname(), set())
            
            for atom2 in neighbors:
                res2 = atom2.get_parent()
                res2_num = res2.id[1]
                interacting_res1[res1_num][1].add(res2_num)
                
            if res2_num not in interacting_res2:
                    interacting_res2[res2_num] = (res2.get_resname(), set())
                interacting_res2[res2_num][1].add(res1_num)
    
    # Convert sets to sorted lists for consistent output
    for res_num in interacting_res1:
        interacting_res1[res_num] = (
            interacting_res1[res_num][0],
            sorted(interacting_res1[res_num][1])
        )
    for res_num in interacting_res2:
        interacting_res2[res_num] = (
            interacting_res2[res_num][0],
            sorted(interacting_res2[res_num][1])
        )
    
    # Find continuous regions (keeping existing functionality)
    def find_continuous_regions(residue_numbers):
        if not residue_numbers:
            return []
        residue_numbers = sorted(residue_numbers)
        regions = []
        start = residue_numbers[0]
        prev = start
        for num in residue_numbers[1:]:
            if num - prev > 1:
                regions.append((start, prev))
                start = num
            prev = num
        regions.append((start, prev))
        return regions
    
    regions1 = find_continuous_regions(interacting_res1.keys())
    regions2 = find_continuous_regions(interacting_res2.keys())
    
    def get_residue_names(chain, region):
        residues = {res.id[1]: res.get_resname() for res in chain.get_residues()}
        start_res = residues.get(region[0], "UNK")
        end_res = residues.get(region[1], "UNK")
        return start_res, end_res
    
    regions_info1 = [(r[0], r[1], *get_residue_names(chain1, r)) for r in regions1]
    regions_info2 = [(r[0], r[1], *get_residue_names(chain2, r)) for r in regions2]
    
    return {
        'chain1_regions': regions_info1,
        'chain2_regions': regions_info2,
        'chain1_interacting_residues': len(interacting_res1),
        'chain2_interacting_residues': len(interacting_res2),
        'detailed_interactions': {
            'chain1': interacting_res1,
            'chain2': interacting_res2
        }
    }

def export_interacting_residues(results, output_file):
    """
    Export interacting residues to a file.
    
    Parameters:
    results: Results dictionary from calculate_interface_area
    output_file: Path to output file
    """
    with open(output_file, 'w') as f:
        for pair, data in results.items():
            chain1_id = data['chain1_id']
            chain2_id = data['chain2_id']
            interactions = data['interface_regions']['detailed_interactions']
            
            f.write(f"Interface between chains {chain1_id} and {chain2_id}\n")
            f.write(f"Interface area: {data['interface_area']:.2f} Å²\n\n")
            
            # Write chain 1 interacting residues
            f.write(f"Chain {chain1_id} interacting residues:\n")
            f.write("ResNum\tResType\tInteracting_with\n")
            for res_num in sorted(interactions['chain1'].keys()):
                resname, partners = interactions['chain1'][res_num]
                partners_str = ','.join(map(str, partners))
                f.write(f"{res_num}\t{resname}\t{partners_str}\n")
            
            f.write(f"\nChain {chain2_id} interacting residues:\n")
            f.write("ResNum\tResType\tInteracting_with\n")
            for res_num in sorted(interactions['chain2'].keys()):
                resname, partners = interactions['chain2'][res_num]
                partners_str = ','.join(map(str, partners))
                f.write(f"{res_num}\t{resname}\t{partners_str}\n")
            f.write("\n" + "-"*50 + "\n")

# Keep the existing calculate_interface_area function but update it to include the new detailed interaction data
def calculate_interface_area(cif_file, chain1_id=None, chain2_id=None):
    """
    Calculate interface area and identify interaction regions between chains.
    """
    # [Previous implementation remains the same until the results dictionary]
    parser = MMCIFParser()
    structure = parser.get_structure('complex', cif_file)
    model = structure[0]
    
    chains = list(model.get_chains())
    if len(chains) < 2:
        raise ValueError("Structure must contain at least 2 chains for interface analysis")
    
    sr = ShrakeRupley()
    
    def get_chain_sasa(chain):
        new_structure = Structure.Structure('temp')
        new_model = Model.Model(0)
        new_structure.add(new_model)
        new_model.add(chain.copy())
        sr.compute(new_structure, level="R")
        return sum(residue.sasa for residue in new_model.get_residues())
    
    results = {}
    
    if chain1_id and chain2_id:
        chain_pairs = [(chain1_id, chain2_id)]
    else:
        chain_pairs = [(c1.id, c2.id) 
                      for i, c1 in enumerate(chains) 
                      for c2 in chains[i+1:]]
    
    for c1_id, c2_id in chain_pairs:
        try:
            chain1 = model[c1_id]
            chain2 = model[c2_id]
            
            interface_regions = identify_interface_regions(chain1, chain2)
            
            sasa1 = get_chain_sasa(chain1)
            sasa2 = get_chain_sasa(chain2)
            
            complex_structure = Structure.Structure('complex')
            complex_model = Model.Model(0)
            complex_structure.add(complex_model)
            complex_model.add(chain1.copy())
            complex_model.add(chain2.copy())
            sr.compute(complex_structure, level="R")
            
            complex_sasa = sum(residue.sasa for residue in complex_model.get_residues())
            interface_area = (sasa1 + sasa2 - complex_sasa) / 2
            
            results[f"{c1_id}_{c2_id}"] = {
                'interface_area': interface_area,
                'chain1_sasa': sasa1,
                'chain2_sasa': sasa2,
                'complex_sasa': complex_sasa,
                'chain1_id': c1_id,
                'chain2_id': c2_id,
                'interface_regions': interface_regions
            }
            
        except Exception as e:
            print(f"Error analyzing chains {c1_id}-{c2_id}: {str(e)}")
            continue
    
    return results


cif_file = "PredictionChai_MDM2-FLT1/pred.model_idx_2.rank_0.cif"
output_file = "intResChai_FLT1_MDM2.txt"
# Analyze the interface
results = calculate_interface_area(cif_file)

# Export detailed residue information
export_interacting_residues(results, output_file)