import os
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from main import fit, filter_fastq 


def test_fit_int_bound():
    """
    1. Checking with integer (0 <= 5 <= 10)
    """
    assert fit(5, 10) is True
        
def test_fit_tuple_bound():
    """
    2. Checking with a tuple (range)
    """
    assert fit(15, (10, 20)) is True
        
def test_fit_out_of_range():
    """
    3. Checking the value out of range
    """
    assert fit(25, (10, 20)) is False


@pytest.fixture
def create_sample_fastq(tmp_path):
    """
    Creates a FASTQ file for tests
    """
    file_path = str(tmp_path / 'test.fastq')
    records = [
        SeqRecord(Seq('GCGC'), id='high_gc', description='', 
                  letter_annotations={'phred_quality': [40, 40, 40, 40]}),
        SeqRecord(Seq('ATATATAT'), id='long_low_gc', description='', 
                  letter_annotations={'phred_quality': [10, 10, 10, 10, 10, 10, 10, 10]})
        ]
    SeqIO.write(records, file_path, 'fastq')
    return str(file_path)


def test_filter_gc_bounds(create_sample_fastq, tmp_path):
    """
    4. Checking filtration by GC-composition
    """
    output = str(tmp_path / 'gc_out.fastq')
    filter_fastq(create_sample_fastq, gc_bounds=(80, 100), output_fastq_file=output)
        
    results = list(SeqIO.parse(output, 'fastq'))
    assert len(results) == 1
    assert results[0].id == 'high_gc'

def test_filter_length(create_sample_fastq, tmp_path):
    """
    5. Length filtering check
    """
    output = str(tmp_path / 'len_out.fastq')
    filter_fastq(create_sample_fastq, length_bounds=(5, 10), output_fastq_file=output)

    results = list(SeqIO.parse(output, 'fastq'))
    assert len(results) == 1
    assert results[0].id == 'long_low_gc'

def test_filter_quality(create_sample_fastq, tmp_path):
    """
    6. Quality control of the filter (threshold)
    """
    output = str(tmp_path / 'qual_out.fastq')
    filter_fastq(create_sample_fastq, quality_threshold=30, output_fastq_file=output)
        
    results = list(SeqIO.parse(output, 'fastq'))
    assert len(results) == 1
    assert results[0].id == 'high_gc'

def test_file_creation(create_sample_fastq, tmp_path):
    """
    7. Read/write test: check that the file is actually being created
    """
    output = str(tmp_path / 'exists.fastq')
    filter_fastq(create_sample_fastq, output_fastq_file=output)
    assert os.path.exists(output)

def test_error_invalid_path():
    """
    8. Test for a non-existent file
    """
    result = filter_fastq('non_existent.fastq')
    assert result is None