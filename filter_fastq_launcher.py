import argparse
from main import filter_fastq
from loguru import logger

parser = argparse.ArgumentParser(description='Filtering sequencing reads')
    
parser.add_argument('--input_fastq_file', '-i', help='Path to the input FASTQ file')
    
parser.add_argument('--gc_bounds', '-g', nargs='+', type=int, default=(0, 100), 
                    help='GC composition: one number (upper threshold) or two (range)\n Example: (40,60)')
parser.add_argument('--length_bounds', '-l', nargs='+', type=int, default=(0, 2**32), 
                    help='Length: one number (upper threshold) or two (range)\n Example: (40,60)')
parser.add_argument('--quality_threshold', '-q', type=int, default=0, 
                    help='Minimum quality threshold (Phred score)')
    
parser.add_argument('--output_fastq_file', '-o', default='filtered/output_fastq.fastq', 
                    help='Path to the output FASTQ file')
       
args = parser.parse_args()

logger.remove()
logger.add('filter_fastq.log', rotation='5 MB', level='INFO')
logger.info(f"The beginning of file filtering: {args.input_fastq_file}")

filter_fastq(args.input_fastq_file, 
            args.gc_bounds, 
            args.length_bounds, 
            args.quality_threshold,
            args.output_fastq_file)

logger.info(f'End filtering: {args.input_fastq_file}')