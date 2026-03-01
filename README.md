# Biological Sequence and FastQ Filter 
Antosha_bio_tools :)
--- 

**This toolkit for managing biological sequences (DNA, RNA, Proteins) and filtering sequencing data using object-oriented principles and the Biopython library.**

**Install:**
```
git clone https://github.com/KuchukLambat/Antosha_bio_tools
```

## 🧬 Class Architecture
The project implements a class hierarchy based on Abstract Base Classes (abc.ABC), ensuring a strict interface and high polymorphism.

1. `BiologicalSequence` **(Abstract Base Class)**  
The foundational interface for all sequence types. It enforces the implementation of:

+ `len()`: get sequence length.
+ `[]`: support for indexing and slicing.
+ `str()` / `repr()`: string representation.
+ `check_alphabet()`: check against a specific alphabet.

2. `NucleicAcidSequence`  
An abstract subclass for nucleotide sequences. It implements shared logic:

+ `complement()`: generates a complementary sequence using the polymorphic _complement_nucleotides property.
+ `reverse()`: reverses the sequence direction.
+ `reverse_complement()`: a combined method for the reverse-complementary strand.

3. `DNA/RNA Sequence`  
A class designed for DNA and RNA sequences.
+ DNASequence: Defines DNA complementarity rules. Includes the `transcribe()` method, which returns an RNASequence object.
+ RNASequence: defines RNA complementarity rules.

5. `AminoAcidSequence`  
A class designed for protein sequences.  
Includes `get_approximate_peptide_mass()`: calculates the estimated mass of the peptide based on the average amino acid weight (110 Da).

## 🛠 Filtering Tools
`filter_fastq`

A function for filtering reads from FASTQ files. It uses the Biopython (SeqIO and gc_fraction) library for parsing.

**Filtering Parameters:**
+ gc_bounds: GC-content range in percentages.
+ length_bounds: acceptable read length range.
+ quality_threshold: minimum average Phred quality score.

## Examples

```
my_dna = DNASequence('ATGCatgc')
print(f'Check alphabet: {my_dna.check_alphabet()}')
print(f'Sequence:\t\t\t {my_dna}')
print(f'Complement sequence:\t\t {my_dna.complement()}')
print(f'Reverse sequence:\t\t {my_dna.reverse()}')
print(f'Reverse complement sequence:\t {my_dna.reverse_complement()}')
print(f'Transcribe sequence:\t\t {my_dna.transcribe()}')
```
    #Check alphabet: True
    #Sequence:			 ATGCatgc
    #Complement sequence:		 TACGtacg
    #Reverse sequence:		 cgtaCGTA
    #Reverse complement sequence:	 gcatGCAT
    #Transcribe sequence:		 AUGCaugc


```
my_rna = RNASequence('AUGCaugc')
print(f'Check alphabet: {my_rna.check_alphabet()}')
print(f'Sequence:\t\t\t {my_rna}')
print(f'Complement sequence:\t\t {my_rna.complement()}')
print(f'Reverse sequence:\t\t {my_rna.reverse()}')
print(f'Reverse complement sequence:\t {my_rna.reverse_complement()}')
```
     #Check alphabet: True
     #Sequence:			 AUGCaugc
     #Complement sequence:		 UACGuacg
     #Reverse sequence:		 cguaCGUA
     #Reverse complement sequence:	 gcauGCAU


```
my_protein = AminoAcidSequence('MALWMRLLPLLALLALWGPDPAAA')
print(f'Check alphabet: {my_protein.check_alphabet()}')
print(f'Sequence: {my_protein}')
print(f'Length of sequence: {len(my_protein)}')
print(f'Complement sequence: {my_protein.get_approximate_peptide_mass()}')
```
     #Check alphabet: True
     #Sequence: MALWMRLLPLLALLALWGPDPAAA
     #Length of sequence: 24
     #Complement sequence: 2640


```
filter_fastq('example_fastq.fastq', gc_bounds=(55,60), length_bounds=(30,80), quality_threshold=30)
```
        Output fastq file:
        
        @SRX079804:1:SRR292678:1:1101:190845:190845 1:N:0:1 BH:changed:1
        CCTCAGCGTGGATTGCCGCTCATGCAGGAGCAGATAATCCCTTCGCCATCCCATTAAGCGCCGTTGTCGGTATTCC
        +
        FF@FFCFEECEBEC@@BBBBDFBBFFDFFEFFEB8FFFFFFFFEFCEB/>BBA@AFFFEEEEECE;ACD@DBBEEE
        @SRX079804:1:SRR292678:1:1101:475517:475517 1:N:0:1 BH:failed
        TGTAGCGGGAGGGTGGAAGCAGTGGGCCCTACCACCTACACAACCTGTTTGCTCAAGAT
        +
        GEFEEGGGG<BEEBDE@E;EBB;@C?@?C<EDE=EGBBBBF8<FFEBE>BDBDEBCBCA
        @SRX079804:1:SRR292678:1:1101:601307:601307 1:N:0:1 BH:changed:1
        TTGGCGTGCTGATGATTATCGGTATCTTCAAAGGCGCGCAGCCTGCGGGCTG
        +
        GGGEGFGGEGE:EE>GFFGGGGDCGEBFFF>G=EBFFEC?DFGAD?DDECBE
        @SRX079804:1:SRR292678:1:1101:667761:667761 1:N:0:1 BH:failed
        CAGCCTTTTGAGGTCGTCTATCGCAGCGTGTCCGCGACGTTTTGTTGCG
        +
        GGFGGG=GGGG@GFGGG@GGEGFGGGEGGGFGGEG@EEEDE8EE=E=DE


Athour: Shaposhnikov Anton