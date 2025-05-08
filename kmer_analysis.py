import sys
from collections import defaultdict

def read_sequences(filename):
    """
    Reads DNA sequences from a file, ignoring headers commonly found
    in FASTA and FASTQ formats.

    Parameters:
    filename (str): Path to the sequence file.

    Returns:
    List[str]: A list of DNA sequence strings (all uppercase).
    """
    sequences = []

    try:
        # Open the input file for reading
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue  # Skip empty lines

                # Skip header or metadata lines found in common formats
                if line.startswith(('>', '@', '+')):
                    continue

                # Convert to uppercase to standardize input (ACGT)
                sequences.append(line.upper())

    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        sys.exit(1)  # Exit with error code if file does not exist

    return sequences

def generate_kmer_counts(sequences, k):
    """
    Counts how often each k-mer (substring of length k) occurs across
    all sequences using a sliding window approach.

    Parameters:
    sequences (List[str]): List of DNA sequences.
    k (int): Desired length of each k-mer.

    Returns:
    Dict[str, int]: A dictionary mapping each k-mer to its count.
    """
    kmer_counts = defaultdict(int)

    for seq in sequences:
        # Slide a window of size k across each sequence
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            kmer_counts[kmer] += 1  # Increment count for this k-mer

    return dict(kmer_counts)

def generate_next_char_counts(sequences, k):
    """
    Builds a nested dictionary to track which characters follow each
    k-mer and how often.

    Parameters:
    sequences (List[str]): List of DNA sequences.
    k (int): Length of each k-mer.

    Returns:
    Dict[str, Dict[str, int]]: A dictionary where each k-mer maps to
    another dictionary of character frequencies that follow it.
    """
    next_char_counts = defaultdict(lambda: defaultdict(int))

    for seq in sequences:
        # We stop at len(seq) - k to ensure there's a next character
        for i in range(len(seq) - k):
            kmer = seq[i:i+k]              # Current k-mer
            next_char = seq[i+k]           # Character immediately after k-mer
            next_char_counts[kmer][next_char] += 1  # Count transition

    # Convert nested defaultdicts to normal dicts for clean output
    return {kmer: dict(next_chars) for kmer, next_chars in next_char_counts.items()}

def write_output(kmer_counts, next_char_counts, output_filename):
    """
    Writes a formatted output file listing each k-mer, its total count,
    and the distribution of characters that follow it.

    Parameters:
    kmer_counts (Dict[str, int]): Frequency of each k-mer.
    next_char_counts (Dict[str, Dict[str, int]]): Characters that follow each k-mer and their counts.
    output_filename (str): Path to the output file.
    """
    with open(output_filename, 'w') as f:
        for kmer in sorted(kmer_counts.keys()):
            total = kmer_counts[kmer]

            # Get the next-character mapping if it exists
            next_chars = next_char_counts.get(kmer, {})

            # Format the next-char dictionary into a readable string
            next_chars_str = ', '.join(f'{char}: {count}' for char, count in sorted(next_chars.items()))

            # Final line format: AT: 5 total, next {G: 3, T: 2}
            f.write(f"{kmer}: {total} total, next {{{next_chars_str}}}\n")

def main():
    """
    Main function that parses command-line arguments, processes the
    input file, performs k-mer analysis, and writes output.
    """
    # Basic usage check
    if len(sys.argv) < 3:
        print("Usage: python kmer_analysis.py <input_file> <k> [output_file]")
        sys.exit(1)

    input_file = sys.argv[1]

    try:
        k = int(sys.argv[2])
    except ValueError:
        print("Error: k must be an integer.")
        sys.exit(1)

    # Optional output file name, defaulting to 'kmer_output.txt'
    output_file = sys.argv[3] if len(sys.argv) > 3 else "kmer_output.txt"

    # Step 1: Read and clean the sequences
    sequences = read_sequences(input_file)

    # Step 2: Count all valid k-mers across all sequences
    kmer_counts = generate_kmer_counts(sequences, k)

    # Step 3: Track what characters follow each k-mer
    next_char_counts = generate_next_char_counts(sequences, k)

    # Step 4: Write the results to a file in readable format
    write_output(kmer_counts, next_char_counts, output_file)

if __name__ == "__main__":
    # If this file is run directly, call main()
    main()
