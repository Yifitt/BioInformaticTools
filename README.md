# bio_seq Python Toolkit 🧬

A comprehensive Python class for DNA, RNA, and protein sequence analysis, featuring common bioinformatics functionalities.

---

## Features

- **Sequence validation** for DNA, RNA, and Protein sequences  
- **Random sequence generation** (DNA/RNA)  
- **Transcription** (DNA → RNA) and **back-transcription** (RNA → DNA)  
- **Complement** and **reverse complement** computation for DNA/RNA  
- **GC content calculation** (whole sequence and in subsequences)  
- **Translation** of DNA/RNA sequences into amino acid sequences  
- **Codon usage frequency** for given amino acids  
- Generation of **six reading frames** (including reverse complement frames)  
- Extraction of **possible proteins** from reading frames  
- **Motif and k-mer searching** with positional information and counts (including overlapping kmers)  
- Finding **most frequent k-mers** of a specified length  
- **Protein molecular mass** calculation  
- **Melting temperature** calculation for DNA/RNA sequences  
- Nucleotide frequency **visualization** via matplotlib bar plots  
- Utility methods for **reading and writing sequences** from/to text and FASTA files

---

## Installation

Clone the repository or copy the `bio_seq` class along with supporting modules (`bio_structure.py` and utility functions) into your project.

---

## Usage Example

```python
from bio_structure import *
from bio_seq import *
from utilities import *

# Create a DNA sequence object
my_dna = bio_seq("ATGCGTACGTTAG", "DNA", "Example Sequence")

# Print sequence info
print(my_dna.get_seq_info())

# Transcribe DNA to RNA
print(my_dna.transcription())

# Translate DNA to amino acids (protein)
print(my_dna.translate_seq())

# Calculate GC content percentage
print(my_dna.gc_content)

# Get reverse complement sequence
print(my_dna.reverseComplement())

# Generate a random DNA sequence of length 15
my_dna.generate_rnd_seq(length=15, seq_type="DNA")
print(my_dna.seq)

# Find motif positions (k-mer) in the sequence
motif_positions = my_dna.find_motif_positions("CGT")
print(motif_positions)

# Find most frequent 3-mers (trimers)
print(my_dna.find_most_frequent_kmers(3))

# Plot nucleotide frequency as horizontal bar chart
my_dna.plot_nuc_frequency(bar_direction="Horizontal")

# Read sequences from a FASTA file
fasta_data = read_FASTA("example.fasta")
for label, seq in fasta_data.items():
    print(f"{label}: {seq}")

# Write sequence to a text file
writeTextFile("output.txt", my_dna.seq)
```
## Functions Overview

| Function                          | Description                                             |
|----------------------------------|---------------------------------------------------------|
| `__validateSeq()`                 | Validates sequence type and content                      |
| `get_seq_info()`                  | Returns sequence label, type, length, and sequence       |
| `generate_rnd_seq(length, type)` | Creates a random DNA or RNA sequence                      |
| `nucleotide_frequency`            | Counts nucleotide occurrences, returns a dictionary      |
| `transcription()`                 | Converts DNA to RNA                                       |
| `backtranscription()`             | Converts RNA to DNA                                       |
| `Complement()`                   | Returns complement sequence                               |
| `reverseComplement()`             | Returns reverse complement sequence                       |
| `gc_content`                     | Calculates GC content percentage                          |
| `gc_content_subsec(k)`            | GC content for subsequences of length k                   |
| `translate_seq(init_pos)`         | Translates DNA/RNA to amino acid sequence                 |
| `codon_usage(aminoacid)`          | Codon frequency for a given amino acid                    |
| `gen_reading_frames()`            | Generates six reading frames (3 forward + 3 reverse)      |
| `proteins_from_rf(aa_seq)`        | Extracts possible proteins from amino acid sequences      |
| `all_proteins_from_orfs()`        | Finds all proteins in all open reading frames             |
| `find_motif_positions(kmer)`      | Finds motif start positions in sequence                    |
| `count_kmer(kmer)`                | Counts occurrences of k-mer (including overlaps)          |
| `find_most_frequent_kmers(k_len)`| Finds most frequent k-mers of given length                |
| `calculate_protein_mass()`        | Calculates molecular mass of a protein sequence            |
| `melting_temp()`                  | Calculates melting temperature of DNA/RNA sequence         |
| `plot_nuc_frequency(bar_direction)` | Plots nucleotide frequency as bar chart                 |

## File Handling Utilities

| Function                             | Description                                      |
|------------------------------------|------------------------------------------------|
| `readTextFile(filePath)`            | Reads sequence data from a plain text file       |
| `writeTextFile(filePath, seq, mode='w')` | Writes a sequence to a text file (default overwrite mode) |
| `read_FASTA(filePath)`              | Parses a FASTA file and returns a dictionary `{label: sequence}` |

