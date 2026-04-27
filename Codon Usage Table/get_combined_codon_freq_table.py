import os
import matplotlib.pyplot as plt
from Bio import SeqIO
from get_codon_freq_table import calculate_codon_frequencies, plot_codon_frequencies

def render_combined_codon_table(fasta_filepath):
    """
    Parses a given FASTA file, aggregates the overall codon frequencies from every 
    sequence record natively, and graphically renders a combined frequency distribution chart.
    
    Args:
        fasta_filepath (str): The absolute or relative path to the target FASTA file.
    """
    if not os.path.exists(fasta_filepath):
        print(f"Error: Target FASTA file '{fasta_filepath}' does not exist.")
        return
        
    # Initialize the global dictionary for tally accumulation
    global_freq = {}
    
    records_parsed = 0
    # Parse the FASTA file robustly utilizing Biopython's standard IO interface
    try:
        for record in SeqIO.parse(fasta_filepath, "fasta"):
            records_parsed += 1
            # Retrieve isolated tallies structurally from the existing robust logic
            local_freq = calculate_codon_frequencies(record.seq)
            
            # Aggregate the sequence-level dictionary manually natively
            for codon, count in local_freq.items():
                if codon in global_freq:
                    global_freq[codon] += count
                else:
                    global_freq[codon] = count
                    
        print(f"Successfully processed {records_parsed} sequence records from '{fasta_filepath}'.")
        print(f"Found {len(global_freq)} unique aggregate codons globally.\nGenerating Matplotlib Histogram...")
        
        # Reuse the established, beautifully constructed plotting interface natively
        plot_codon_frequencies(global_freq)
        
    except Exception as e:
        print(f"An error occurred while parsing the FASTA file: {e}")

# =============================================================================
# Stress Testing & Execution Demonstration
# =============================================================================
#
# if __name__ == "__main__":
#     # We utilize the requested `exFasta3Rec.fasta` test file path mathematically relative here:
#     test_fasta = os.path.join(os.path.dirname(os.path.dirname(__file__)), "datasets", "exFasta3Rec.fasta")
#
#     print(f"Initiating Combined Codon Frequency analysis on targeting file: {test_fasta}")
#     render_combined_codon_table(test_fasta)
