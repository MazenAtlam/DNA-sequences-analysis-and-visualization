import matplotlib.pyplot as plt
from Bio.Seq import Seq

def calculate_codon_frequencies(sequence):
    """
    Calculates the frequency of each codon in a given DNA sequence.
    
    The sequence is linearly processed in non-overlapping increments of 3 bases.
    Any trailing sequences shorter than 3 bases at the end are ignored.
    
    Args:
        sequence (str or Bio.Seq.Seq): The DNA sequence to be analyzed.
        
    Returns:
        dict: A Python dictionary safely mapping each complete 3-letter 
              codon string to its absolute integer frequency tally.
    """
    # Enforce uppercase string mapping
    seq_str = str(sequence).upper()
    
    # Traverse sequence in non-overlapping increments of 3
    codons = [seq_str[i:i+3] for i in range(0, len(seq_str) - len(seq_str) % 3, 3)]
    
    # Calculate distribution frequencies
    codon_freq = {}
    for codon in codons:
        if codon in codon_freq:
            codon_freq[codon] += 1
        else:
            codon_freq[codon] = 1
    
    return codon_freq

def plot_codon_frequencies(codon_freq):
    """
    Renders a bar chart of codon frequencies using matplotlib.
    
    Args:
        codon_freq (dict): A dictionary of codon frequencies 
                           originating from calculate_codon_frequencies().
    """
    if not codon_freq:
        print("Empty or invalid sequence detected. No codons to plot.")
        return
        
    # Sort the codon keys alphabetically for better visualization
    sorted_codons = sorted(codon_freq.keys())
    frequencies = [codon_freq[codon] for codon in sorted_codons]
    
    # Adjust plot dimensions based on the number of unique elements
    fig_width = min(max(len(sorted_codons) * 0.4, 8), 16)
    plt.figure(figsize=(fig_width, 6))
    
    # Initialize the bar chart geometry
    plt.bar(sorted_codons, frequencies, color='teal', edgecolor='black', alpha=0.8)
    
    # Annotate the plot
    plt.title('Codon Usage Frequency Distribution', fontsize=16)
    plt.xlabel('Codons', fontsize=14)
    plt.ylabel('Absolute Frequency Count', fontsize=14)
    
    # Rotate X-axis coordinates to prevent clustering and overlap
    plt.xticks(rotation=90)
    
    # Set Y-ticks to integers
    max_freq = max(frequencies) if frequencies else 0
    plt.yticks(range(0, max_freq + 2))
    
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.show()

# =============================================================================
# Stress Testing & Execution Demonstration
# =============================================================================
#
# if __name__ == "__main__":
#     # Dataset to force frequency overlaps
#     test_sequence = Seq("ATGCGTAAATGATAGTATAAAAATGCCCTGA")
#     
#     print(f"Codon Processing Test on sequence of length: {len(test_sequence)}")
#     print("-" * 50)
#     
#     # Applying codon frequency counting
#     freq_extracted = calculate_codon_frequencies(test_sequence)
#     
#     print(f"Total Unique Codons Extracted: {len(freq_extracted)}")
#     for structural_codon, struct_count in freq_extracted.items():
#         print(f"  {structural_codon}: {struct_count}")
#         
#     print("\nRendering corresponding Matplotlib Bar Chart...")
#     plot_codon_frequencies(freq_extracted)
