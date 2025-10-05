from additional_modules import dna_rna_tools_functions as drtf
from additional_modules import filtering_functions as ff

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

        output = drtf.is_nucleic_acid(seq)

        if output:
            dict_procedures = {
                "is_nucleic_acid": drtf.is_nucleic_acid(seq),
                "transcribe": drtf.transcribe(seq),
                "reverse": drtf.reverse(seq),
                "complement": drtf.complement(seq),
                "reverse_complement": drtf.reverse_complement(seq)
            }

            output = dict_procedures[procedure]
            
        else:
            output = f"Your sequence has extraneous characters: {seq}"

        outputs.append(output)

    if len(outputs) == 1:
        return outputs[0]
    return outputs


def filter_fastq(seqs: dict, gc_bounds: tuple = (0, 100), length_bounds: tuple = (0, 2**32), quality_threshold: int = 0) -> dict:
    """
    The function filters readings by a set of conditions: GC composition, length and quality.
    
    Args:
        seqs: dictionary fastq readings; gc_bounds, length_bounds, quality_threshold: 
        integer or tuple of integers.
        
    Return:
        dictionary of filtred readings
    """
    filtered_fastq = dict()
    
    for key, value in seqs.items():
        sequence, sequence_quality = value[0], value[1]
        
        if (ff.filter_gc(sequence, gc_bounds) and ff.filter_length(sequence, length_bounds) and 
            ff.filter_quality(sequence_quality, quality_threshold)):
            
            filtered_fastq[key] = value
    
    return filtered_fastq