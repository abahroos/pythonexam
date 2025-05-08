"""
Unit tests for kmer_analysis.py using pytest.

Each function in the k-mer analysis pipeline is tested independently to verify:
- Correct handling of DNA sequence input
- Accurate frequency calculation of k-mers and their next characters
- Handling of edge cases (e.g., empty input, short sequences)
- Correct file writing format and output structure

These tests are designed to reflect realistic input patterns while maintaining reproducibility and clarity.
"""

import pytest
from kmer_analysis import (
    read_sequences,
    generate_kmer_counts,
    generate_next_char_counts,
    write_output
)
import os

# Sample input used across multiple tests
sample_sequences = [
    "ATGTCTGTCTGAA",
    "TCTGAA"
]

def test_read_sequences(tmp_path):
    """
    Test whether read_sequences correctly parses DNA sequence data
    from a FASTA-style file, skipping any header lines starting with '>'.
    """
    file = tmp_path / "test.fa"
    file.write_text(">header1\nATGTCTGTCTGAA\n>header2\nTCTGAA\n")

    sequences = read_sequences(str(file))

    assert sequences == ["ATGTCTGTCTGAA", "TCTGAA"], "Sequences should match content after header lines"

def test_generate_kmer_counts_typical():
    """
    Test generate_kmer_counts on standard input sequences with k=2.
    Manually verify output matches expected counts of each overlapping k-mer.
    """
    k = 2
    result = generate_kmer_counts(sample_sequences, k)

    expected = {
        'AT': 1,
        'TG': 3,
        'GT': 2,
        'TC': 3,
        'CT': 3,
        'GA': 3,
        'AA': 2
    }

    assert result == expected, f"Expected {expected}, but got {result}"

def test_generate_kmer_counts_edge_empty():
    """
    Test generate_kmer_counts with an empty list of sequences.
    Should return an empty dictionary instead of erroring or miscounting.
    """
    k = 3
    result = generate_kmer_counts([], k)
    assert result == {}, "Empty input should yield empty k-mer counts"

def test_generate_next_char_counts_typical():
    """
    Test generate_next_char_counts for standard inputs with k=2.
    Confirms the frequency of characters that follow each k-mer is counted correctly.
    """
    k = 2
    result = generate_next_char_counts(sample_sequences, k)

    assert result['AT']['G'] == 1, "'AT' should be followed by 'G' once"
    assert result['TG']['T'] == 2, "'TG' should be followed by 'T' twice"
    assert result['CT']['G'] == 3, "'CT' should be followed by 'G' three times"

def test_generate_next_char_counts_edge_empty():
    """
    Test generate_next_char_counts with empty input.
    Ensures that no exceptions are raised and an empty dictionary is returned.
    """
    k = 3
    result = generate_next_char_counts([], k)
    assert result == {}, "Empty input should return empty next-character dictionary"

def test_write_output(tmp_path):
    """
    Verify that write_output correctly formats and writes k-mer and
    next character counts to a file.
    """
    kmer_counts = {'AT': 1, 'TG': 2}
    next_char_counts = {'AT': {'G': 1}, 'TG': {'T': 2}}

    output_file = tmp_path / "output.txt"
    write_output(kmer_counts, next_char_counts, str(output_file))

    with open(output_file, 'r') as f:
        lines = f.read().strip().split('\n')

    assert "AT: 1 total, next {G: 1}" in lines, "Output should contain formatted count for AT"
    assert "TG: 2 total, next {T: 2}" in lines, "Output should contain formatted count for TG"
