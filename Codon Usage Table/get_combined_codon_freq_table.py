import os
import matplotlib.pyplot as plt
from Bio import SeqIO
from get_codon_freq_table import calculate_codon_frequencies, plot_codon_frequencies

def render_combined_codon_table(fasta_filepath):
    """
    Parses a given FASTA file, aggregates the overall codon frequencies from every 
    sequence record, and graphically renders a combined frequency distribution chart.
    
    Args:
        fasta_filepath (str): The absolute or relative path to the target FASTA file.
    """
    if not os.path.exists(fasta_filepath):
        print(f"Error: Target FASTA file '{fasta_filepath}' does not exist.")
        return
        
    # Initialize the global dictionary for codon frequencies
    global_freq = {}
    
    records_parsed = 0
    # Parse the FASTA file using Biopython's standard IO
    try:
        for record in SeqIO.parse(fasta_filepath, "fasta"):
            records_parsed += 1
            # Calculate codon frequencies for the current sequence
            local_freq = calculate_codon_frequencies(record.seq)
            
            # Aggregate the sequence-level dictionary manually
            for codon, count in local_freq.items():
                if codon in global_freq:
                    global_freq[codon] += count
                else:
                    global_freq[codon] = count
                    
        print(f"Successfully processed {records_parsed} sequence records from '{fasta_filepath}'.")
        print(f"Found {len(global_freq)} unique aggregate codons globally.\nGenerating Matplotlib Histogram...")
        
        # Plotting the combined codon frequency table
        plot_codon_frequencies(global_freq)
        
    except Exception as e:
        print(f"An error occurred while parsing the FASTA file: {e}")

# =============================================================================
# Stress Testing & Execution Demonstration
# =============================================================================
#
# if __name__ == "__main__":
#     # Using the `exFasta3Rec.fasta` test file
#     test_fasta = os.path.join(os.path.dirname(os.path.dirname(__file__)), "datasets", "exFasta3Rec.fasta")
#
#     print(f"Generating Combined Codon Frequency table for: {test_fasta}")
#     render_combined_codon_table(test_fasta)
