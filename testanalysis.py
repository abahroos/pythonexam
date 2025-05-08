import pytest
from kmer_analysis import read_sequences, generate_kmer_counts, generate_next_char_counts, write_output
import os

# Sample sequences used for testing
sample_sequences = [
    "ATGTCTGTCTGAA",
    "TCTGAA"
]

def test_read_sequences(tmp_path):
    """
    Test reading sequences from a file, skipping headers like > or @.
    """
    # Create a temporary FASTA-style file with headers
    file = tmp_path / "test.fa"
    file.write_text(">header1\nATGTCTGTCTGAA\n>header2\nTCTGAA\n")

    # Use the read_sequences function
    sequences = read_sequences(str(file))

    # Check that the function returns the correct sequences
    assert sequences == ["ATGTCTGTCTGAA", "TCTGAA"]

def test_generate_kmer_counts_typical():
    """
    Test k-mer counting on typical input sequences.
    """
    k = 2
    result = generate_kmer_counts(sample_sequences, k)

    # Expected k-mer counts calculated manually
    expected = {
        'AT': 1, 'TG': 3, 'GT': 2, 'TC': 3, 'CT': 3, 'GA': 3, 'AA': 2
    }

    # Assert each expected k-mer count matches the result
    for kmer in expected:
        assert result[kmer] == expected[kmer]

def test_generate_kmer_counts_edge_empty():
    """
    Test k-mer counting with an empty sequence list (edge case).
    """
    k = 3
    result = generate_kmer_counts([], k)

    # An empty list should return an empty dictionary
    assert result == {}

def test_generate_next_char_counts_typical():
    """
    Test next character frequency counting on typical input sequences.
    """
    k = 2
    result = generate_next_char_counts(sample_sequences, k)

    # Check that specific k-mers lead to expected next characters
    assert result['AT']['G'] == 1
    assert result['TG']['T'] == 2
    assert result['CT']['G'] == 3

def test_generate_next_char_counts_edge_empty():
    """
    Test next character counting with an empty sequence list (edge case).
    """
    k = 3
    result = generate_next_char_counts([], k)

    # An empty list should return an empty dictionary
    assert result == {}

def test_write_output(tmp_path):
    """
    Test writing output file for k-mer and next character counts.
    """
    # Prepare dummy k-mer counts and next character counts
    kmer_counts = {'AT': 1, 'TG': 2}
    next_char_counts = {'AT': {'G': 1}, 'TG': {'T': 2}}

    # Define a temporary output file path
    output_file = tmp_path / "output.txt"

    # Write the output
    write_output(kmer_counts, next_char_counts, str(output_file))

    # Read back the file
    with open(output_file, 'r') as f:
        lines = f.read().strip().split('\n')

    # Verify the expected lines are present in the output
    assert "AT: 1 total, next {G: 1}" in lines
    assert "TG: 2 total, next {T: 2}" in lines
