import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def parse_interface_file(filename):
    """Parse the interface file and extract interaction data."""
    interactions = []
    current_chain = None
    reading_residues = False
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            
            if not line:
                continue
                
            if 'Chain A interacting residues:' in line:
                current_chain = 'A'
                reading_residues = False
                continue
            elif 'Chain B interacting residues:' in line:
                current_chain = 'B'
                reading_residues = False
                continue
            elif line.startswith('ResNum'):
                reading_residues = True
                continue
            
            if reading_residues and current_chain and not line.startswith('Interface'):
                try:
                    parts = line.split('\t')
                    if len(parts) >= 3:
                        res_num = int(parts[0])
                        interacting_with = [int(x) for x in parts[2].split(',')]
                        
                        if current_chain == 'A':
                            for partner in interacting_with:
                                interactions.append((res_num, partner))
                except ValueError:
                    continue
    
    return interactions

def create_interaction_heatmap(filename, output_file='interface_heatmapChai.png'):
    """Create a heatmap visualization of protein-protein interactions."""
    # Set the style
    plt.style.use('seaborn-white')
    sns.set_style("ticks")
    
    # Parse the data
    interactions = parse_interface_file(filename)
    
    # Create full range of residue numbers
    chek2_nums = list(range(1, 492))  # 1-543
    mdm2_nums = list(range(1, 1338))   # 1-491
    
    # Create matrix for full protein sequences
    matrix = np.zeros((len(chek2_nums), len(mdm2_nums)))
    
    # Fill in interactions
    for a_res, b_res in interactions:
        try:
            a_idx = a_res - 1
            b_idx = b_res - 1
            if 0 <= a_idx < len(chek2_nums) and 0 <= b_idx < len(mdm2_nums):
                matrix[a_idx, b_idx] = 1
        except IndexError:
            continue
    
    # Create figure with specific size and DPI
    plt.figure(figsize=(20, 15), dpi=100)
    
    # Create the heatmap with enhanced styling and flipped y-axis
    ax = sns.heatmap(matrix,
                    xticklabels=50,
                    yticklabels=50,
                    cmap=['white', '#FF1F1F'],
                    cbar=False,
                    square=True,
                    linewidths=0,
                    rasterized=True)
    
    # Flip the y-axis
    ax.invert_yaxis()
    
    # Customize the plot
    plt.xlabel('MDM2 Residue Number', fontsize=12, fontweight='bold')
    plt.ylabel('FLT1 Residue Number', fontsize=12, fontweight='bold')
    plt.title('FLT1-MDM2 Interface Interactions', fontsize=14, fontweight='bold', pad=20)
    
    # Add grid lines at major ticks
    ax.grid(False)
    
    # Customize tick labels
    plt.xticks(np.arange(0, len(mdm2_nums), 50), 
               np.arange(1, len(mdm2_nums)+1, 50),
               rotation=0)
    plt.yticks(np.arange(0, len(chek2_nums), 50), 
               np.arange(1, len(chek2_nums)+1, 50),
               rotation=0)
    
    # Add spines (axis lines)
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)
    
    # Make tick labels larger
    ax.tick_params(axis='both', which='major', labelsize=10)
    
    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    
    # Save the plot with high quality
    plt.savefig(output_file, 
                dpi=300, 
                bbox_inches='tight',
                facecolor='white',
                edgecolor='none')
    plt.show()

# Example usage:
if __name__ == "__main__":
    input_file = 'intResChai_FLT1_MDM2.txt'
    try:
        create_interaction_heatmap(input_file)
        print("Heatmap created successfully!")
    except Exception as e:
        print(f"Error: {str(e)}")