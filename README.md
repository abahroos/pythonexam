# Asha Bahroos

# K-mer Analysis Project

## Overview

This project explores k-mer analysis—a foundational technique in bioinformatics and genome assembly. A **k-mer** is a substring of length `k` found within a larger DNA sequence, which is composed of the four nucleotides: A, C, G, and T.

The goal of this script is to:
- Count how often each unique k-mer appears across all sequences in a file
- Track the frequency of each **next character** that immediately follows a given k-mer

This type of k-mer context tracking is a key starting point for reconstructing full genomes from fragmented sequencing data, though actual assembly is beyond the scope of this project.

> Think of this as building a map of all short patterns in the DNA and how they typically continue—essentially prepping the data for assembly down the road.

## Repository Structure

```
439539-kmer-project/
├── kmer_analysis.py         # Main script: reads input, extracts k-mers, writes output
├── test_kmer_analysis.py    # Test suite using Pytest to validate each function
├── README.md                # You're reading it!
├── DOCUMENTATION.md         # Detailed design decisions, edge case handling, and k-mer logic
```

## How to Run the Script

### Requirements
- Pytest (for running tests)

Install pytest if you don’t already have it to be able to run this:
```bash
pip install pytest
```

### Command-Line Usage

To run the analysis:
```bash
python kmer_analysis.py <input_file> <k> [output_file]
```

**Arguments:**
- `<input_file>`: Path to the input FASTA/FASTQ file containing DNA reads  
- `<k>`: The k-mer length (an integer, e.g., 3, 5, etc.)  
- `[output_file]` (optional): File to write results to. If not provided, defaults to `kmer_output.txt`

**Example:**
```bash
python kmer_analysis.py ../shared/439539/assignment_bash_data/BalReg/BalReg_SRR952792_Trim_Left.fastq 3 output.txt
```

## Output Format

The output file will list each k-mer, how many times it appeared, and a breakdown of which characters typically follow it. Example:

```
AT: 5 total, next {G: 3, T: 2}
TG: 7 total, next {G: 5, A: 2}
GT: 4 total, next {C: 4}
```

This tells us, for instance, that the k-mer "AT" occurred 5 times, and was followed by "G" in 3 of those cases and "T" in 2.

## How to Run the Tests

The script is modularized and all major functions are covered by unit tests using **Pytest**, including edge cases like:
- Very short sequences
- Repetitive patterns
- Mixed-case input (converted to uppercase internally)

To run all tests:
```bash
pytest
```
Make sure you’re in the root directory of the project when you do this.

## Notes

- Assumes DNA input (A, C, G, T only). Lowercase is converted to uppercase.
- Ignores FASTA/FASTQ headers and non-sequence lines:
  - Lines starting with `>`, `@`, or `+` are skipped
- Can be used with standard `.fa`, `.fasta`, or `.fastq` files
- Designed to handle large input files efficiently with minimal memory footprint

## Why This Matters

Understanding how short patterns repeat and link together in DNA is what makes large-scale sequencing projects possible. By tracking not just the k-mers but what follows them, this project simulates a critical component of genome assembly pipelines—without diving into graph theory or actual assembly logic.
