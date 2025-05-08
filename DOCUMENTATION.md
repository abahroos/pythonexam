# DOCUMENTATION.md

## 1. Data Structures Used

**Main structures used:**

- `Dict[str, int]` – for storing **total counts** of each k-mer  
  - Each unique k-mer (like `"AT"` or `"TG"`) maps to the number of times it appears across all sequences.  
  - Example: `{ 'AT': 5, 'TG': 7 }`  
  - This lets you efficiently track how frequently each substring of length `k` occurs.

- `Dict[str, Dict[str, int]]` – for tracking **what character comes after** each k-mer  
  - Each k-mer maps to another dictionary that stores how often specific characters follow it.  
  - Example: `{ 'AT': {'G': 3, 'T': 2}, 'TG': {'G': 5, 'A': 2} }`  
  - This structure reflects biological reality: the same k-mer can be followed by different nucleotides in different reads.

**Why dictionaries?**
- Dictionary lookups and updates are fast (average O(1) time), which is important when processing large datasets.
- They also make the logic clean and readable—especially with nested mappings.
- We used `collections.defaultdict` to simplify counting. This avoids writing repeated `if key in dict` checks.

**Example using `defaultdict`:**
```python
from collections import defaultdict

# Automatically initializes any new key to 0
kmer_counts = defaultdict(int)
kmer_counts["AT"] += 1  # No need to check if "AT" exists first
```

Using this pattern in both k-mer counts and next-character counts reduced potential bugs and made the code cleaner.

---

## 2. How Edge Cases Are Handled

Several edge cases were considered in both code and testing:

- **Empty sequences or files:**  
  - If a sequence is empty, the functions return empty dictionaries.
  - The code handles this naturally without crashing, and the tests confirm it.

- **Sequences shorter than `k`:**  
  - If a sequence is shorter than the given `k`, the `range(len(seq) - k + 1)` logic just skips over it.  
  - This prevents any invalid k-mers from being counted and avoids index errors.

- **Sequences shorter than `k+1`:**  
  - These sequences might produce a valid k-mer but don’t have a following character.
  - We still count the k-mer, but we don’t record a "next" character. This is intentional—there’s no character to follow.

- **Last k-mer in a sequence:**  
  - If the k-mer is at the very end of the string (so nothing comes after it), it still contributes to total k-mer counts.  
  - But no next-character frequency is recorded. This keeps the output biologically accurate.

- **Headers in FASTA or FASTQ files:**  
  - Lines that start with `>`, `@`, or `+` are treated as metadata and skipped.  
  - This makes the function compatible with common bioinformatics file formats.

---

## 3. Avoiding Overcounting and Losing Data

To make sure we neither overcount k-mers nor miscount what comes after them:

- **Sliding Window Approach:**  
  - The code uses `seq[i:i+k]` to extract each k-mer in order. This is the standard method for getting all overlapping substrings.
  - For next-character analysis, we only count the character at `seq[i + k]` if it exists. This avoids going out of bounds or adding non-existent characters.

- **No duplication between functions:**  
  - `generate_kmer_counts()` and `generate_next_char_counts()` are completely separate.  
  - This makes sure each function has a clear responsibility and avoids accidental double-counting.

- **Handling the End of Sequences:**  
  - We intentionally avoid recording next-character data if there isn’t a character after the k-mer.
  - This helps maintain the integrity of the results and ensures no "ghost" characters are added.

- **Tests Back Everything Up:**  
  - Edge cases like empty input, short sequences, and output formatting are all covered in unit tests.
  - This helps confirm that our approach doesn’t silently fail or overcount under special conditions.

---

## Conclusion

This k-mer analysis tool was built with clarity, performance, and biological accuracy in mind. It uses efficient data structures (`dict` and `defaultdict`), handles edge cases cleanly, and separates counting logic to avoid confusion or overlap. The functions are modular, tested with pytest, and compatible with standard DNA formats like FASTA and FASTQ. With these foundations, the script is robust enough to scale up to larger genomic datasets while remaining easy to understand and maintain.
