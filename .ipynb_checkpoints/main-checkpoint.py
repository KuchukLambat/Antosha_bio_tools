from additional_modules.dna_rna_tools_functions import is_nucleic_acid, transcribe, reverse, complement, reverse_complement
from additional_modules.filtering_functions import get_bounds, filter_gc, filter_length, filter_quality
import os

def run_dna_rna_tools(*seqs: str) -> list[str]:
    """

    The function accepts an arbitrary number of arguments with DNA or RNA sequences,
    as well as the name of the procedure to be performed (this is always the last argument, procedure).
    After that, it performs the specified action on all the transmitted sequences and returns the result.
    
    Args:
        data (strings): string or strings of data, where the last argument is a procedure.
        
    Return:
        list or string: processed sequences
    """

    *seqs, procedure = seqs

    if procedure not in (
        "is_nucleic_acid",
        "transcribe",
        "reverse",
        "complement",
        "reverse_complement",
    ):
        return None

    outputs = list()

    for seq in seqs:

        output = is_nucleic_acid(seq)

        if output:
            dict_procedures = {
                "is_nucleic_acid": is_nucleic_acid,
                "transcribe": transcribe,
                "reverse": reverse,
                "complement": complement,
                "reverse_complement": reverse_complement
            }

            output = dict_procedures[procedure](seq)
            
        else:
            output = f"Your sequence has extraneous characters: {seq}"

        outputs.append(output)

    if len(outputs) == 1:
        return outputs[0]
    return outputs


def filter_fastq(input_fastq: str, gc_bounds: tuple = (0, 100), 
                 length_bounds: tuple = (0, 2**32), quality_threshold: int = 0,
                output_fastq: str = 'filtered/output_fastq.fastq'):
    """
    The function filters readings by a set of conditions: GC composition, length and quality.
    
    Args:
        input_fastq: a string with path to the file; gc_bounds, length_bounds, quality_threshold: 
        integer or tuple of integers.
        
    Return:
        output_fastq: file with filtered readings
    """
    if not os.path.isfile(input_fastq):
        return None

    if not os.path.isdir('filtered'):
        os.mkdir('filtered')
    
    with (open(input_fastq, 'r') as input_fastq, open(output_fastq, 'w') as output_fastq):
        
        for line in input_fastq:
            key = line
            sequence = input_fastq.readline().strip()
            plus_string = input_fastq.readline()
            sequence_quality = input_fastq.readline().strip()
            print(sequence, sequence_quality)

            if (filter_gc(sequence, gc_bounds) and filter_length(sequence, length_bounds) and 
            filter_quality(sequence_quality, quality_threshold)):
                
                output_fastq_record = f"{key}{sequence}\n{plus_string}{sequence_quality}\n"
                output_fastq.write(output_fastq_record)


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


def get_bounds(bounds: int | tuple) -> bool:
    """
    The function outputs a range of values.
    
    Args:
        data: integer or tuple of integers.
    
    Return:
        two integers.
    """
    if type(bounds) == int:
        min_bound = 0
        max_bound = bounds
    
    else:
        min_bound = min(bounds)
        max_bound = max(bounds)

    return min_bound, max_bound


def filter_gc(sequence: str, gc_bounds: int | tuple=(0, 100)) -> bool:
    """
    The function filters readings by GC composition.
    
    Args:
        sequence: string, gc_bounds: integer or tuple of integers.
    
    Return:
        boolean result
    """
    gc = (sequence.count('G') + sequence.count('C')) / len(sequence) * 100

    min_gc_bound, max_gc_bound = get_bounds(gc_bounds)
    
    return min_gc_bound <= gc <= max_gc_bound


def filter_length(sequence: str, length_bounds: int | tuple=(0, 2**32)) -> bool:
    """
    The function filters readings by length.
    
    Args:
        sequence: string, length_bounds: integer or tuple of integers.
    
    Return:
        boolean result
    """
    sequence_length = len(sequence)
    
    min_length_bound, max_length_bound = get_bounds(length_bounds)
    
    return  min_length_bound <= sequence_length <= max_length_bound


def filter_quality(sequence_quality: str, quality_threshold: int=0) -> bool:
    """
    The function filters readings by quality
    
    Args:
        sequence_quality: string, quality_threshold: integer or tuple of integers.
    
    Return:
        boolean result
    """
    
    sequence_quality_score = [(ord(chr_quality) - 33) for chr_quality in sequence_quality]
    sequence_quality_score = sum(sequence_quality_score) / len(sequence_quality_score)

    return quality_threshold <= sequence_quality_score