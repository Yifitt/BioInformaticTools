# bio_seq Python Toolkit 🧬

A comprehensive Python class for DNA, RNA, and protein sequence analysis with common bioinformatics functionalities.

---

## Features

- Sequence validation (DNA, RNA, Protein)
- Random sequence generation (DNA/RNA)
- Transcription (DNA → RNA)
- Complement and Reverse Complement of DNA/RNA
- GC content calculation (whole sequence and subsequences)
- Translation of DNA/RNA sequences into amino acids
- Codon usage frequency for given amino acids
- Generation of six reading frames
- Protein extraction from reading frames
- Motif and k-mer searching with position and frequency counts
- Finding most frequent k-mers in a sequence
- Protein molecular mass calculation
- Reading and writing sequences from/to text and FASTA files

---

## Installation

Clone the repository or copy the `bio_seq` class and utility functions into your project.

---

## Usage Example

```python
from bio_structure import bio_seq, read_FASTA, writeTextFile

# Create a DNA sequence object
my_dna = bio_seq("ATGCGTACGTTAG", "DNA", "Example Sequence")

# Get sequence info
print(my_dna.get_seq_info())

# Transcribe DNA to RNA
print(my_dna.transcription())

# Translate DNA to protein sequence (amino acids)
print(my_dna.translate_seq())

# Calculate GC content
print(my_dna.gc_content())

# Get reverse complement of the DNA sequence
print(my_dna.reverseComplement())

# Generate random DNA sequence of length 15
my_dna.generate_rnd_seq(length=15, seq_type="DNA")
print(my_dna.seq)

# Reading a FASTA file
fasta_data = read_FASTA("example.fasta")
for label, seq in fasta_data.items():
    print(label, seq)

# Write a sequence to a text file
writeTextFile("output.txt", my_dna.seq)
****
```
## Functions Overview

| Function                      | Description                                           |
|------------------------------|-----------------------------------------------------|
| `__validateSeq()`             | Validates sequence type and content                  |
| `get_seq_info()`              | Returns sequence label, type, length and sequence    |
| `generate_rnd_seq(length, type)` | Creates a random DNA or RNA sequence               |
| `countNucFrequency()`         | Counts nucleotide occurrences                         |
| `transcription()`             | Converts DNA to RNA                                   |
| `Complement()`                | Returns complement sequence                           |
| `reverseComplement()`         | Returns reverse complement sequence                   |
| `gc_content()`                | Calculates GC content percentage                       |
| `gc_content_subsec(k)`        | GC content for subsequences of length k               |
| `translate_seq(init_pos)`     | Translates DNA/RNA to amino acid sequence             |
| `codon_usage(aminoacid)`      | Codon frequency for an amino acid                      |
| `gen_reading_frames()`        | Generates 6 reading frames                             |
| `proteins_from_rf(aa_seq)`    | Extracts proteins from amino acid sequences            |
| `all_proteins_from_orfs()`    | Finds all proteins in all open reading frames         |
| `find_motif_positions(seq, kmer)` | Finds motif start positions                       |
| `count_kmer(seq, kmer)`       | Counts occurrences of k-mer (including overlaps)      |
| `find_most_frequent_kmers(seq, k_len)` | Finds most frequent k-mers                     |
| `calculate_protein_mass(seq)` | Calculates molecular mass of a protein sequence        |

---

## File Handling Utilities

| Function                          | Description                                   |
|----------------------------------|-----------------------------------------------|
| `readTextFile(filePath)`          | Read sequence data from a text file            |
| `writeTextFile(filePath, seq, mode='w')` | Write sequence to a file (default overwrite)  |
| `read_FASTA(filePath)`             | Parse FASTA files and return `{label: sequence}` dictionary |
