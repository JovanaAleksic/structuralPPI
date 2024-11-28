import pandas as pd
import os



def translate_dna(dna_sequence):
    # Define the codon to amino acid dictionary
    codon2aa = {
        'TCA':'S','TCC':'S','TCG':'S','TCT':'S','TTC':'F','TTT':'F','TTA':'L',
        'TTG':'L','TAC':'Y','TAT':'Y','TAA':'*','TAG':'*','TGC':'C','TGT':'C',
        'TGA':'*','TGG':'W','CTA':'L','CTC':'L','CTG':'L','CTT':'L','CCA':'P',
        'CCC':'P','CCG':'P','CCT':'P','CAC':'H','CAT':'H','CAA':'Q','CAG':'Q',
        'CGA':'R','CGC':'R','CGG':'R','CGT':'R','ATA':'I','ATC':'I','ATT':'I',
        'ATG':'M','ACA':'T','ACC':'T','ACG':'T','ACT':'T','AAC':'N','AAT':'N',
        'AAA':'K','AAG':'K','AGC':'S','AGT':'S','AGA':'R','AGG':'R','GTA':'V',
        'GTC':'V','GTG':'V','GTT':'V','GCA':'A','GCC':'A','GCG':'A','GCT':'A',
        'GAC':'D','GAT':'D','GAA':'E','GAG':'E','GGA':'G','GGC':'G','GGG':'G',
        'GGT':'G'
    }
    
    # Convert the DNA sequence to uppercase
    dna_sequence = dna_sequence.upper()
    
    # Initialize an empty protein sequence
    protein_sequence = ""
    
    # Iterate through the DNA sequence in steps of 3 (codon by codon)
    for i in range(0, len(dna_sequence), 3):
        # Extract the codon
        codon = dna_sequence[i:i+3]
        
        # If the codon is incomplete (less than 3 nucleotides), break the loop
        if len(codon) < 3:
            break
        
        # Translate the codon to an amino acid and add it to the protein sequence
        amino_acid = codon2aa.get(codon, 'X')  # 'X' for unknown codons
        protein_sequence += amino_acid
        
        # # Stop translation if we encounter a stop codon
        # if amino_acid == '*':
        #     break
    
    return protein_sequence[:-1]




def create_separate_fasta_files(input_xlsx, output_dir):
    """
    Create separate FASTA files, each containing MDM2 and one protein from the list
    
    Parameters:
    input_xlsx (str): Path to input Excel file
    output_dir (str): Directory to store output FASTA files
    """
    try:
        
        # Read Excel file
        df = pd.read_excel(input_xlsx, sheet_name="Sheet1")
        
        # Get MDM2 sequence
        mdm2_row = df[df['SYMBOL'].str.contains('MDM2', case=False, na=False)]
        if mdm2_row.empty:
            raise ValueError("MDM2 sequence not found in the Excel file")
            
        mdm2_seq = translate_dna(mdm2_row['Nucleotides'].iloc[0])
        
        # Process each protein
        for idx, row in df.iterrows():
            protein_name = row['SYMBOL']
            if 'MDM2' in protein_name.upper():
                continue  # Skip MDM2 as we're combining it with others
            
            # Translate nucleotide sequence
            protein_seq = translate_dna(row['Nucleotides'])
            if protein_seq:
                # Create output file name
                output_file = os.path.join(output_dir, f"{protein_name}_MDM2.fasta")
                
                # Write to FASTA file
                with open(output_file, 'w') as fasta_out:
                    # Write both sequences
                    fasta_out.write(f">{protein_name}\n{protein_seq}\n")
                    fasta_out.write(f">MDM2\n{mdm2_seq}\n")
                
                print(f"Created FASTA file for {protein_name}: {output_file}")
                
    except Exception as e:
        print(f"Error processing file: {e}")

# Example usage
if __name__ == "__main__":
    input_file = "cancerPathwayGenesTab.xlsx"  # Your input Excel file
    output_directory = "fasta_files"  # Directory to store output files
    create_separate_fasta_files(input_file, output_directory)