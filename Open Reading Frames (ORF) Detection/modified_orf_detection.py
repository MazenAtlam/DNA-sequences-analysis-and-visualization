from orf_detection import get_orfs
from Bio.Seq import Seq

def get_orfs_with_min_length(sequence, min_length=100):
    """
    Scans a given DNA sequence for complete ORFs and filters them 
    based on a minimum length threshold.
    
    This function wraps the base `get_orfs` logic and removes any 
    ORFs whose nucleotide length is less than the `min_length`.

    Args:
        sequence (str or Bio.Seq.Seq): The DNA sequence to be analyzed.
        min_length (int): The minimum acceptable length for an ORF to be included.
                          Default is 100 nucleotides.
        
    Returns:
        list[dict]: A list of dictionaries containing the filtered ORF details.
    """
    
    # Fetch all valid ORFs across the 3 forward frames
    all_orfs = get_orfs(sequence)
    
    # Filter out any ORF that its length is less than min_length
    filtered_orfs = [orf for orf in all_orfs if orf["length"] >= min_length]
    
    return filtered_orfs

# =============================================================================
# Stress Testing & Execution Demonstration
# =============================================================================
#
# if __name__ == "__main__":
#     # Same sequence from orf_detection.py
#     # Expected ORFs:
#     # 1. ATGCGTAAATGA (length: 12)
#     # 2. ATGCCCTGA    (length: 9)
#     # 3. ATGATAGTATAA (length: 12)
#     test_sequence = Seq("ATGCGTAAATGATAGTATAAAATGCCCTGA")
#
#     print(f"Testing DNA Sequence: {test_sequence}")
#     print("=" * 40)
#
#     # Test with min_length = 10 to filter out the 3rd ORF exclusively
#     threshold = 10
#     filtered_results = get_orfs_with_min_length(test_sequence, min_length=threshold)
#
#     print(f"Filtered ORFs found (min_length={threshold} nt): {len(filtered_results)}\n")
#
#     for count, orf in enumerate(filtered_results, start=1):
#         print(f"ORF #{count}:")
#         for key, val in orf.items():
#             print(f"  {key}: {val}")
#         print("-" * 25)
