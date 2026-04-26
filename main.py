from abc import ABC, abstractmethod
import os
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


class BiologicalSequence(ABC):
    """
    An abstract class for working with biological sequences.
    """
    
    @abstractmethod
    def __len__(self):
        pass

    @abstractmethod
    def __getitem__(self, key):
        pass

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def __repr__(self):
        pass

    @abstractmethod
    def check_alphabet(self):
        pass


class NucleicAcidSequence(BiologicalSequence):
    """
    A class for working with nucleotide sequences.
    
    Attributes:
        seq (str): nucleotide sequence.
    """
    
    def __init__(self, seq):
        self.seq = seq

    def __len__(self):
        return len(self.seq)

    def __getitem__(self, key):
        return self.seq[key]

    def __str__(self):
        return self.seq

    def __repr__(self):
        return self.seq

    def check_alphabet(self):
        """
        Check if the sequence matches the specified alphabet.
    
        Args:
            str: nucleotide sequence
    
        Return:
            boolean result
        """
        return set(self.seq.upper()) <= {'A', 'G', 'C', 'T', 'U'}

    @property
    @abstractmethod
    def _complement_nucleotides(self):
        """
        Return the complementarity rule of the inheritor object.
        """
        pass
        
    def complement(self):
        """
        Create a complementary nucleotide sequence.
    
        Args:
            str : nucleotide sequence.
            
        Return:
            object : complementary sequence.
        """
        complement_seq = self.seq.translate(self._complement_nucleotides)
        return self.__class__(complement_seq)

    def reverse(self):
        """ 
        Reverse the nucleotide sequence.
        
        Args:
            str : nucleotide sequence.
            
        Return:
            object : reversed sequence.
        """
        return self.__class__(self.seq[::-1])

    def reverse_complement(self):
        """ 
        Create a reverse complementary nucleotide sequence.

        Args:
            str : nucleotide sequence.
            
        Return:
            object : reverse complementary sequence.
        """
        return self.complement().reverse()
        

class DNASequence(NucleicAcidSequence):
    """
    A class for working with DNA sequences.
    """
    
    @property
    def _complement_nucleotides(self):
        """ 
        Define the complementarity rule.
        
        Return:
            tuple
        """
        return str.maketrans('ATGCatgc', 'TACGtacg')

    def transcribe(self):
        """ 
        Transcribes DNA into RNA.

        Args:
            str : DNA sequence

        Return:
            RNASequence : RNA sequence
        """
        return RNASequence(self.seq.replace("T", "U").replace("t", "u"))


class RNASequence(NucleicAcidSequence):
    """
    A class for working with RNA sequences.
    """
    
    @property
    def _complement_nucleotides(self):
        """ 
        Define the complementarity rule.

        Return:
            tuple
        """
        return str.maketrans('AUGCaugc', 'UACGuacg')
    

class AminoAcidSequence(BiologicalSequence):
    """
    A class for working with amino acid sequences.
    """

    def __init__(self, seq):
        self.seq = seq

    def __len__(self):
        return len(self.seq)

    def __getitem__(self, key):
        return self.seq[key]

    def __str__(self):
        return self.seq

    def __repr__(self):
        return self.seq

    def check_alphabet(self):
        """
        Check if the sequence matches the specified alphabet.
    
        Args:
            str: amino acid sequence.
    
        Return:
            boolean result
        """
        return set(self.seq.upper()) <= set('ACDEFGHIKLMNPQRSTVWY')
    
        
    def get_approximate_peptide_mass(self):
        """ 
        Calculates the approximate mass of the peptide.

        Args:
            str : Amino acid sequence.

        Return:
            float : Approximate mass of the peptide.
        """
        return len(self.seq) * 110


def fit(value, bounds: int | tuple) -> bool:
    """
    The function outputs True if value is in a range.
    
    Args:
        value: number
        data: number or tuple of numbers.
    
    Return:
        bool
    """
    if type(bounds) == int:
        min_bound = 0
        max_bound = bounds
    
    else:
        min_bound = min(bounds)
        max_bound = max(bounds)

    return min_bound <= value <= max_bound


def filter_fastq(input_fastq_file: str, gc_bounds: tuple = (0, 100), 
                 length_bounds: tuple = (0, 2**32), quality_threshold: int = 0,
                output_fastq_file: str = 'filtered/output_fastq.fastq'):
    """
    The function filters readings by a set of conditions: GC composition, length and quality.
    
    Args:
        input_fastq_file: a string with path to the file; 
        gc_bounds, length_bounds, quality_threshold: number or tuple of numbers.
        
    Return:
        output_fastq_file: file with filtered readings
    """
    if not os.path.isfile(input_fastq_file):
        return None

    if not os.path.isdir('filtered'):
        os.mkdir('filtered')
        
    filtered_reads = list()

    for rec in SeqIO.parse(input_fastq_file, "fastq"):
        rec_gc_fraction = 100 * gc_fraction(rec.seq)
        rec_length = len(rec)
        rec_quality_score = sum(rec.letter_annotations["phred_quality"]) / rec_length
        
        if (fit(rec_gc_fraction, gc_bounds) and fit(rec_length, length_bounds) and 
            quality_threshold <= rec_quality_score):
            filtered_reads.append(rec)

    SeqIO.write(filtered_reads, output_fastq_file, "fastq")