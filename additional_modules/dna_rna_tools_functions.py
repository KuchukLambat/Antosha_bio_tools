def is_nucleic_acid(seq: str) -> bool:
    """
    Check if the sequence matches the specified alphabet.
    
    Args:
        data (string): nucleotide sequence
    
    Return:
        boolean result
    """
    bases = set(seq.upper())
    
    valid_chars_dna = {'A', 'G', 'C', 'T'}
    valid_chars_rna = {'A', 'G', 'C', 'U'}
    
    return bases <= valid_chars_dna or bases <= valid_chars_rna

def transcribe(seq: str) -> str:
    """
    This function transcribes the nucleotide sequence.

    Args:
        data (string): nucleotide sequence
    Return:
        string: transcribed sequence.
    """
    if ("U" not in seq) and ("u" not in seq):
        return seq.replace("T", "U").replace("t", "u")

    else:
        return f"Your sequence has already been transcribed: {seq}"


def reverse(seq: str) -> str:
    """
    This function reverses the nucleotide sequence.
    
    Args:
        data (string): nucleotide sequence.
    Return:
        string: reversed sequence.
    """
    return seq[::-1]


def complement(seq: str) -> str:
    """
    This function creates a complementary nucleotide sequence.
    
    Args:
        data (string): nucleotide sequence.
    Return:
        string: complementary sequence.
    """
    if ("T" in seq) or ("t" in seq):
        complement_rules = {
            "A": "T",
            "a": "t",
            "T": "A",
            "t": "a",
            "G": "C",
            "g": "c",
            "C": "G",
            "c": "g",
        }
    elif ("U" in seq) or ("u" in seq):
        complement_rules = {
            "A": "U",
            "a": "u",
            "U": "A",
            "u": "a",
            "G": "C",
            "g": "c",
            "C": "G",
            "c": "g",
        }
    return "".join([complement_rules[base] for base in seq])


def reverse_complement(seq: str) -> str:
    """
    This function creates a reverse complementary nucleotide sequence. 
    
    Args:
        data (string): nucleotide sequence.
    Return:
        string: reverse complementary sequence.
    """
    complement_seq = complement(seq)
    complement_reverse_seq = reverse(complement__seq)

    return complement_reverse_seq