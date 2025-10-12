from additional_modules import dna_rna_tools_functions as is_nucleic_acid, transcribe, reverse, complement, reverse_complement
from additional_modules import filtering_functions as get_bounds, filter_gc, filter_length, filter_quality
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
            sequence = input_fastq.readline()
            plus_string = input_fastq.readline()
            sequence_quality = input_fastq.readline()

            if (filter_gc(sequence, gc_bounds) and filter_length(sequence, length_bounds) and 
            filter_quality(sequence_quality, quality_threshold)):
                
                output_fastq_record = f"{key}{sequence}{plus_string}{sequence_quality}"
                output_fastq.write(output_fastq_record)