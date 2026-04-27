import matplotlib.pyplot as plt
from Bio.Seq import Seq

def calculate_codon_frequencies(sequence):
    """
    Calculates the frequency of each codon in a given DNA sequence.
    
    The sequence is linearly processed in non-overlapping increments of 3 bases.
    Any trailing sequences strictly shorter than 3 bases at the end are automatically ignored.
    
    Args:
        sequence (str or Bio.Seq.Seq): The DNA sequence to be analyzed.
        
    Returns:
        dict: A Python dictionary safely mapping each complete 3-letter 
              codon string to its absolute integer frequency tally.
    """
    # Enforce uppercase string mapping to establish standardized codon keys computationally
    seq_str = str(sequence).upper()
    
    # Traverse sequence robustly in non-overlapping structural increments of 3
    # The condition effectively omits the trailing fragments via slicing mechanics
    codons = [seq_str[i:i+3] for i in range(0, len(seq_str) - len(seq_str) % 3, 3)]
    
    # Calculate exact distribution frequencies manually from scratch
    codon_freq = {}
    for codon in codons:
        if codon in codon_freq:
            codon_freq[codon] += 1
        else:
            codon_freq[codon] = 1
    
    return codon_freq

def plot_codon_frequencies(codon_freq):
    """
    Renders a bar chart of clustered codon frequencies graphically using matplotlib.
    
    Args:
        codon_freq (dict): A dictionary of structural codon frequencies 
                           originating from calculate_codon_frequencies().
    """
    if not codon_freq:
        print("Empty or invalid sequence detected. No codons to plot.")
        return
        
    # Organically sort the extracted codon keys alphabetically for a uniform and 
    # stable visual baseline representation
    sorted_codons = sorted(codon_freq.keys())
    frequencies = [codon_freq[codon] for codon in sorted_codons]
    
    # Calibrate plot dimensions adaptively based on the density of unique elements
    fig_width = min(max(len(sorted_codons) * 0.4, 8), 16)
    plt.figure(figsize=(fig_width, 6))
    
    # Initialize the bar chart geometry directly mapping frequencies accurately
    plt.bar(sorted_codons, frequencies, color='teal', edgecolor='black', alpha=0.8)
    
    # Annotate and clarify aesthetic constraints
    plt.title('Codon Usage Frequency Distribution', fontsize=16)
    plt.xlabel('Codons', fontsize=14)
    plt.ylabel('Absolute Frequency Count', fontsize=14)
    
    # Rotate X-axis coordinates strictly to alleviate clustering and native overlap visual bugs
    plt.xticks(rotation=90)
    
    # Enforce integer constraints natively on the Y-ticks since counts are wholly integer
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
#     # Comprehensive structural dataset to force frequency overlaps
#     test_sequence = Seq("ATGCGTAAATGATAGTATAAAAATGCCCTGA")
#     
#     print(f"Executing Codon Processing Test on sequence of length: {len(test_sequence)}")
#     print("-" * 50)
#     
#     # Subvert logic sequentially to tally codon dictionary structures natively
#     freq_tally = calculate_codon_frequencies(test_sequence)
#     
#     print(f"Total Unique Codons Extracted: {len(freq_tally)}")
#     for structural_codon, struct_count in freq_tally.items():
#         print(f"  {structural_codon}: {struct_count}")
#         
#     print("\nRendering corresponding Matplotlib UI Bar Chart...")
#     plot_codon_frequencies(freq_tally)
