import numpy as np
import matplotlib.pyplot as plt
from Bio.Seq import Seq

def generate_dot_plot(seq1, seq2):
    """
    Generates and displays a dot plot to visualize structural alignments 
    between two DNA sequences.
    
    The dot plot compares the sequences nucleotide by nucleotide. Matches
    are represented as '1' (dots) in a 2D matrix, while mismatches are '0',
    creating diagonal lines where extended homologies exist.
    
    Args:
        seq1 (str or Bio.Seq.Seq): The first DNA sequence corresponding to the Y-axis.
        seq2 (str or Bio.Seq.Seq): The second DNA sequence corresponding to the X-axis.
    """
    
    seq1_str = str(seq1).upper()
    seq2_str = str(seq2).upper()
    
    len1 = len(seq1_str)
    len2 = len(seq2_str)
    
    # Initialize a 2D matrix with zeroes representing mismatches
    dot_matrix = np.zeros((len1, len2))
    
    # Populate the matrix with 1s where nucleotides natively match
    for i in range(len1):
        for j in range(len2):
            if seq1_str[i] == seq2_str[j]:
                dot_matrix[i, j] = 1

    # Render the dot plot using Matplotlib
    plt.figure(figsize=(8, 8))
    
    # Find all coordinates where the matrix is 1 (matches)
    y_coords, x_coords = np.where(dot_matrix == 1)
    
    # Plot as simple black small open circles ('o') instead of solid rectangles
    plt.scatter(x_coords, y_coords, marker='o', facecolors='none', edgecolors='black', s=30)
    
    # Adjust axes limits to encompass the matrix dimensions fully
    # Setting Y-axis descending so that y=0 is positioned at the top like standard imshow
    plt.xlim(-0.5, len2 - 0.5)
    plt.ylim(len1 - 0.5, -0.5)
    
    # Contextual plot aesthetics and labeling
    plt.title('Sequence Alignment Dot Plot', fontsize=16)
    plt.xlabel('Sequence 2', fontsize=14)
    plt.ylabel('Sequence 1', fontsize=14)
    
    # Display the nucleotides directly on the axes if the sequences are sufficiently short
    if len1 <= 50 and len2 <= 50:
        plt.yticks(range(len1), list(seq1_str))
        plt.xticks(range(len2), list(seq2_str))

        # Style the tick labels to be bigger, have a background of a brown square,
        # and use a monospace font so that bounding boxes are uniformly sized.
        bbox_props = dict(boxstyle="square,pad=0.3", fc="brown", ec="brown")
        
        for label in plt.gca().get_yticklabels():
            label.set_fontsize(14)
            label.set_family("monospace")
            label.set_weight("bold")
            label.set_color("white")
            label.set_bbox(bbox_props)
            
        for label in plt.gca().get_xticklabels():
            label.set_fontsize(14)
            label.set_family("monospace")
            label.set_weight("bold")
            label.set_color("white")
            label.set_bbox(bbox_props)
            
    plt.tight_layout()
    plt.show()

# =============================================================================
# Stress Testing & Execution Demonstration
# =============================================================================
#
# if __name__ == "__main__":
#     # Test Case 1: Exact self-alignment. Provides a solid diagonal line.
#     test_seq_a = Seq("ATGCGTAAATGATAGTATAAAA")
#     test_seq_b = Seq("ATGCGTAAATGATAGTATAAAA")
#
#     print("Executing Test Case 1: Exact Self-Alignment.")
#     generate_dot_plot(test_seq_a, test_seq_b)
#
#     # Test Case 2: Partial overlap and repeats mapping.
#     test_seq_c = Seq("ATGATAGTAT")
#
#     print("Executing Test Case 2: Partial overlap and shifting alignments.")
#     generate_dot_plot(test_seq_a, test_seq_c)
