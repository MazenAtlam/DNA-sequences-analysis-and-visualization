from Bio.Seq import Seq

def get_orfs(sequence):
    """
    Scans a given DNA sequence for Open Reading Frames (ORFs) through 
    all three forward reading frames (offsets 0, 1, and 2).
    
    An ORF is defined as a sequence of codons that begins with the start codon 'ATG' 
    and ends with one of the stop codons ('TAA', 'TAG', 'TGA') in the same reading frame.

    Args:
        sequence (str or Bio.Seq.Seq): The DNA sequence to be analyzed.
        
    Returns:
        list[dict]: A list of dictionaries containing the extracted ORF details:
            - 'sequence' (str): The nucleotide sequence of the ORF.
            - 'start' (int): The 0-based start index (inclusive) of the ORF in the parent sequence.
            - 'end' (int): The 0-based end index (exclusive) of the ORF in the parent sequence.
            - 'length' (int): The total length of the ORF in nucleotides.
            - 'frame' (int): The forward reading frame (0, 1, or 2).
    """
    
    # Convert sequence to an upper-case string
    seq_str = str(sequence).upper()
    
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}
    orfs = []
    
    # Scan all 3 forward reading frames
    for frame in range(3):
        i = frame
        while i <= len(seq_str) - 3:
            codon = seq_str[i:i+3]
            
            # Found a start codon
            if codon == start_codon:
                stop_found = False
                j = i + 3
                
                # Scan the following codons in increments of 3 for the first stop codon
                while j <= len(seq_str) - 3:
                    stop_candidate = seq_str[j:j+3]
                    
                    if stop_candidate in stop_codons:
                        # Found a complete, valid ORF
                        orf_seq = seq_str[i:j+3]
                        orfs.append({
                            "sequence": orf_seq,
                            "start": i,
                            "end": j + 3,
                            "length": len(orf_seq),
                            "frame": frame
                        })
                        
                        # Update the search index `i` to `j + 3`
                        # to search the remainder of that specific reading frame
                        i = j + 3
                        stop_found = True
                        break
                    
                    j += 3
                
                # If the rest of the sequence has no stop codon,
                # break to advance to the next frame.
                if not stop_found:
                    break
            else:
                # Advance if the current codon was not a start codon
                i += 3
                
    return orfs

# =============================================================================
# Stress Testing & Execution Demonstration
# =============================================================================
#
# if __name__ == "__main__":
#     # Small test DNA sequence containing multiple ORFs of different frames
#     # Frame 0: ATG...TAA 
#     # Frame 1: .ATG...TAG
#     # Frame 2: ..ATG...TGA
#     test_sequence = Seq("ATGCGTAAATGATAGTATAAAATGCCCTGA")
#
#     print(f"Testing DNA Sequence: {test_sequence}")
#     print("=" * 40)
#
#     found_orfs = get_orfs(test_sequence)
#     print(f"Total ORFs found: {len(found_orfs)}\n")
#
#     for count, orf in enumerate(found_orfs, start=1):
#         print(f"ORF #{count}:")
#         for key, val in orf.items():
#             print(f"  {key}: {val}")
#         print("-" * 25)
