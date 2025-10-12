def convert_multiline_fasta_to_oneline(input_file: str, output_file: str):
    """
    This function converts multiline fasta to oneline

    Args:
        input_file: file with sequences in multiline fasta format
    Return:
        output_file: file with sequences in oneline fasta format 
    """
    with (open(input_file, 'r') as input_file,  
        open(output_file, 'w') as output_file):
        
        description = input_file.readline()
        sequence = str()
            
        for line in input_file:
            if line.startswith('>'):
                fasta_record = f'\n{description}{sequence}'
                output_file.write(fasta_record)
                    
                description = line
                sequence = str()
                    
            else:
                sequence += line.strip()
 
 
def lower_chr(string: str):
    """
    This function converts string to lowercase for key of sort().

    Args:
        string of protein name
    Return:
        lowercase string
    """
    return string.lower()


def parse_blast_output(input_file: str, output_file: str):
    """
    This function parse blast output and 
    writes the name of the proteins in 
    alphabetical order to the file

    Args:
        input_file: blast output
    Return:
        output_file: sorted protein names
    """
    protein_records = list()
    
    with open(input_file, 'r') as input_file, open(output_file, 'w') as output_file:
        
        for line in input_file:
            if line.startswith('Description'):
                colum_names_line = line
                protein_name_end = colum_names_line.find('Name')
                accession_number = colum_names_line.find('Accession')
    
                target_line = input_file.readline()
                protein_record = f'{target_line[:protein_name_end].strip()} ({target_line[accession_number:].strip()})'
                protein_records.append(protein_record)
    
        protein_records.sort(key=lower_chr)
        output_file.write('\n'.join(protein_records))


def search_target_genes(genes: str | tuple, line: str) -> bool:
    """
    This function searches target gene in selected line from gbk file
    Args:
        genes: string or tuple with target genes; line form gbk file
    Return:
        boolean result 
    """
    gene =  line.strip()
    if gene[gene.find('"')+1:-1] in genes:
        return True
    return False

    
def parsin_gene_record(input_file: str, line: str) -> str:
    """
    This function selects genes record from gbk file

    Args:
        input_file: gbk file; line: string from input_file
    Return:
        string with gene name and amino acid sequence
    """
    if '/gene=' in line:
        gene_name = line.strip()
        gene_name = gene_name[gene_name.find('"')+1:-1]

    while not '/translation=' in line:
        line = input_file.readline()

    amino_seq = line[line.find('"')+1:].strip()
    line = input_file.readline()
            
    while line.startswith(' '*21):
        amino_seq += line.strip()
        line = input_file.readline()
        
    return f'>{gene_name}\n{amino_seq[:-1]}'


def select_genes_from_gbk_to_fasta(input_file: str, genes: str | tuple, n_before: int, n_after: int, output_fasta: str):
    """
    This function selects genes from gbk file and writes them to fasta format

    Args:
        input_file: gbk file; genes: string or tuple of target genes; 
        n_before and n_after: number of selected genes; output_fasta: name of output file
    Return:
        output_fasta: file in fasta format
    """
    before_after_target_genes = list()
    n_after_counter = 0
    
    with open(input_file, 'r') as input_file:
        
        for line in input_file:

            if '/gene=' not in line:
                continue

            if search_target_genes(genes, line):
                n_after_counter = n_after

            
            gene_record = parsin_gene_record(input_file, line)
            print(before_after_target_genes)
            before_after_target_genes.append(gene_record)
    
            if len(before_after_target_genes) > (n_before + n_after):
                del before_after_target_genes[0]

            
            if n_after_counter == 0:
                continue
                
            elif n_after_counter > 1:
                n_after_counter -= 1
                
            elif n_after_counter == 1:
                output_file = open(output_fasta, 'w')
                record = '\n'.join(before_after_target_genes)
                output_file.write(record)
                output_file.close()