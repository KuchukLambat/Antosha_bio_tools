# Antosha_bio_tools :)

**Antosha_bio_tools** this is a set of functions that allows performing certain operations with nucleotide sequences.

This set is divided into two parts: a function that allows you to change the nucleotide sequences (run_dna_rna_tools) 

and a function for filtering sequencing reads (filter_fastq).


## Installation

To work with functions, they must be cloned:

~~~
git clone https://github.com/KuchukLambat/Antosha_bio_tools
~~~


## Working with the tools

**`run_dna_rna_tools`** to work with this function, it is necessary to pass string variables, 
where the last argument indicates the name of the operation to be performed.

Implemented functions:

-`is_nucleic_acid` - returns a Boolean result to check if the sequence matches the specified alphabet.

-`transcribe` - returns the transcribed sequence.

-`reverse` - return the reversed sequence.

-`complement` - return the complementary sequence.

-`reverse_complement` - return the reverse complementary sequence.<br>

______
**`def filter_fastq`** this function filters readings by a set of conditions: 

-`filter_gc`  - function filters readings by GC composition.

-`filter_length` - function filters readings by length.

-`filter_quality` - function filters readings by quality.

-adding function `get_bounds` to set a range of values.

## Module bio_files_processor.py
It contains three functions:
`filter_fastq` - The function filters readings by a set of conditions: GC composition, length and quality.
`convert_multiline_fasta_to_oneline` - This function converts multiline fasta to oneline
parse_blast_output - This function parse blast output and writes the name of the proteins in alphabetical order to the file

Athour: Shaposhnikov Anton
