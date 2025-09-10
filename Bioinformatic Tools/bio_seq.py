import functools
from itertools import chain
from typing import List, Callable

from Biyoinformatik.Bio_Seq.exceptions import *
from bio_structure import *


import random
import collections
import re
import matplotlib.pyplot as plt
from Bio import Align
import warnings
import time



def measure_time(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        print(f"'{func.__name__}' metodu {end_time - start_time:.4f} saniye sürdü.")
        return result

    return wrapper


class BioSeq:
    """DNA sequence class.Default value: ATCG, DNA, No label"""
    def __init__(self, seq="ATCG", seq_type="DNA", label='No Label'):
        """Sequence initialization, validation."""
        self.seq = seq.upper()
        self.label = label
        self.seq_type = seq_type.upper()
        self.is_valid = self.__validate_seq()
        assert self.is_valid, f"Provided data does not seem to be a correct {self.seq_type} sequence."

    def __validate_seq(self):
        """Check the sequence to make sure it is a valid sequence."""
        if len(self.seq) != 0:
            if self.seq_type == "RNA" or self.seq_type == "DNA":
                return set(NUCLEOTIDE_BASE[self.seq_type]).issuperset(self.seq)
            elif self.seq_type == "PROTEIN":
                return AminoAcids.issuperset(self.seq)
            else:
                return False
        else:
            return False

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        return f"[Label]: {self.label}, [Biotype]: {self.seq_type}, [Length]: {len(self.seq)}"

    def __bool__(self):
        return bool(self.seq)

    def __iter__(self):
        return iter(self.seq)
    
    def get_seq_biotype(self):
        """Returns sequence type"""
        return self.seq_type

    def get_seq_info(self):
        """Returns 4 strings. Full sequence information"""
        return f"[Label]: {self.label}\n[Sequence]: {self.seq}\n[Biotype]: {self.seq_type}\n[Length]: {len(self.seq)}"

    def generate_rnd_seq(self, length=10, seq_type="DNA"):
        """Generate a random DNA sequence, provided the length"""
        seq = ''.join([random.choice(NUCLEOTIDE_BASE[seq_type])
                       for _ in range(length)])
        self.__init__(seq.upper(), seq_type, "Randomly generated sequence")


    @property
    @functools.lru_cache()
    def nucleotide_frequency(self):
        """Count nucleotides in a given sequence. Return a dictionary"""
        return dict(collections.Counter(self.seq))

    @measure_time
    def transcription(self):
        """DNA -> RNA Transcription.Replacing Thymine with Uracil"""
        mapping  = str.maketrans('Tt','Uu')
        return self.seq.translate(mapping) if self.seq_type == "DNA" else None

    @measure_time
    def backtranscription(self):
        """RNA -> DNA Replacing Uracil with Thymine"""
        mapping = str.maketrans('Uu', 'Tt')
        return self.seq.translate(mapping) if self.seq_type == "RNA" else None

    @measure_time
    def reverse_complement(self):
        """Swapping Adenine->Thymine and Guanine->Cytosine.Reversing newly generated string."""
        mapping = None
        if self.seq_type == "DNA" or self.seq_type == "RNA":
            if self.seq_type == "DNA":
                mapping = str.maketrans('ATCGatcg','TAGCtagc')
            elif self.seq_type == "RNA":
                mapping = str.maketrans('AUCGaucg','UAGCuacg')

            return self.seq.translate(mapping)[::-1]
        else:
            raise InvalidSequenceError(f"Invalid Sequence Type {self.seq_type}!\nOnly DNA and RNA is allowed")

    @measure_time
    def complement(self):
        """Swapping Adenine->Thymine and Guanine->Cytosine. in  DNA sequence
           Swapping Adenine->Uracil and Guanine->Cytosine. in  RNA sequence

           returns: Complemented version of DNA OR RNA input sequence"""
        mapping = None
        if self.seq_type == "DNA" or self.seq_type == "RNA":
            if self.seq_type == "DNA":
                mapping = str.maketrans('ATCGatcg','TAGCtacg')
            elif self.seq_type == "RNA":
                mapping = str.maketrans('AUCGaucg','UAGCuacg')
            return self.seq.translate(mapping)
        else:
            raise InvalidSequenceError(f"Invalid Sequence Type {self.seq_type}!\nOnly DNA and RNA is allowed")

    @property
    def gc_content(self):
        """ GC Content in a DNA/RNA sequence"""
        return round((self.seq.count('C') + self.seq.count('G'))/ len(self.seq)*100)
    
    def gc_content_subseq(self, k=20):
        """GC Content in a DNA/RNA sub-sequence length k. k=20 by default"""
        res = []
        for i in range(0, len(self.seq) - k + 1, k):
            subseq = self.seq[i:i + k]
            res.append(
                round((subseq.count('C') + subseq.count('G')) / len(subseq) * 100))
        return res

    @measure_time
    def translate_seq(self,init_pos = 0):
        """Translates a DNA sequence into an aminoacid sequence."""
        if self.seq_type == "DNA":
            return[DNA_Codons[self.seq[pos:pos+3]]for pos in range(init_pos,len(self.seq)-2,3)]
        elif self.seq_type == "RNA":
            return[RNA_Codons[self.seq[pos:pos+3]]for pos in range(init_pos,len(self.seq)-2,3)]
        else:
           raise InvalidSequenceError(f"Provided data does not seem to be appropriate sequence type ({self.seq_type})")

    @measure_time
    def codon_usage(self, aminoacid):
        """Provides the frequency of each codon encoding a given aminoacid in a DNA sequence"""
        tmplist = []
        if self.seq_type == "DNA":
            for i in range(0, len(self.seq) - 2, 3):
                if DNA_Codons[self.seq[i:i + 3]] == aminoacid:
                    tmplist.append(self.seq[i:i + 3])

        elif self.seq_type == "RNA":
            for i in range(0, len(self.seq) - 2, 3):
                if RNA_Codons[self.seq[i:i + 3]] == aminoacid:
                    tmplist.append(self.seq[i:i + 3])

        freqDict = dict(collections.Counter(tmplist))
        totalW = sum(freqDict.values())
        for seq in freqDict:
            freqDict[seq] = round(freqDict[seq] / totalW, 2)
        return freqDict

    @measure_time
    def gen_reading_frames(self):
        """Generate the six reading frames of a DNA sequence, including reverse complement"""
        frames = [self.translate_seq(0), self.translate_seq(1), self.translate_seq(2)]
        tmp_seq = BioSeq(self.reverse_complement(), self.seq_type)
        frames.append(tmp_seq.translate_seq(0))
        frames.append(tmp_seq.translate_seq(1))
        frames.append(tmp_seq.translate_seq(2))
        del tmp_seq
        return frames

    @staticmethod
    def proteins_from_rf(seq):
        """Compute all possible proteins in an aminoacid seq and return a list of possible proteins"""
        current_prot = []
        proteins = []
        # aa resembles amino-acid here
        for aa in seq:
            if aa == "_":
                # STOP accumulating amino acids if _ - STOP was found
                if current_prot:
                    for p in current_prot:
                        proteins.append(p)
                    current_prot = []
            else:
                # START accumulating amino acids if M - START was found
                if aa == "M":
                    current_prot.append("")
                for i in range(len(current_prot)):
                    current_prot[i] += aa
        return proteins

    @measure_time
    def all_proteins_from_orfs(self, start=0, end=0, ordered=False,largest= False):
        """Compute all possible proteins for all open reading frames.
        ordered' and 'Largest' are set to False by default.
        Largest' returns the largest protein in that sequence when it's set to 'True'"""
        if end > start:
            tmp_seq = BioSeq(
                self.seq[start: end], self.seq_type)
            rfs = tmp_seq.gen_reading_frames()
        else:
            rfs = self.gen_reading_frames()
        res = []
        for rf in rfs:
            prots = self.proteins_from_rf(rf)
            for p in prots:
                res.append(p)
        if ordered:
            if largest:
                return sorted(res, key=len, reverse=True)[0]
            return sorted(res, key=len, reverse=True)
        return res

    @measure_time
    def find_motif_positions(self,kmer):
        """Finds the starting positions of given kmer in a sequence.
        Parameters:
            kmer (str): The specific k-mer to search for in the sequence.
        """
        return[match.start() for match in re.finditer(f'(?=({kmer}))', self.seq)]

    @staticmethod
    def __find_kmers(seq: str, k: int) -> List[str]:
        """ Find k-mers in string """

        n = len(seq) - k + 1
        return [] if n < 1 else [seq[i:i + k] for i in range(n)]

    def count_kmer(self,kmer):
        """
        Counts the number of times a specific k-mer appears in a given sequence,including overlapping k-mers.
            
        Parameters:
        kmer (str): The specific k-mer to search for in the sequence.

        Returns:
        int: The number of times the k-mer appears in the sequence.
        
        """
        return len(re.findall(f'(?=({kmer}))',self.seq)) #RegEx

    def find_most_frequent_kmers(self, k_len):
        """
        Finds the most frequent k-mers of a given length in a DNA string.

        Parameters:
            k_len (int): The length of the k-mers to search for.

        Returns:
            list: A list of the most frequent k-mers in the DNA string.
        """
        # 1. A dictionary to store k-mer frequencies.
        kmer_frequencies = {}

        # 2. A loop to iterate through the DNA string and extract k-mers of a given length,
        # while also incrementing the frequency of each k-mer
        for i in range(len(self.seq) - k_len + 1):
            kmer = self.seq[i:i+k_len]
            if kmer in kmer_frequencies:
                kmer_frequencies[kmer] += 1
            else:
                kmer_frequencies[kmer] = 1
        
        # 3. A variable to store the maximum value in our dictionary
        highest_frequency = max(kmer_frequencies.values())

        # 4. List Comprehension to scan the k-mer dictionary and return a list of k-mers 
        # with the highest frequency
        return [
            kmer for kmer, frequency in kmer_frequencies.items()
            if frequency == highest_frequency
        ]
    def calculate_molecular_mass(self):
        """Calculates the mass of a sequence and returns it.
         Returns:
            mass: The rounded mass of the molecule .
            """
        #We assign a mass variable in order to store the mass of the amino-acids
        mass = 0
        if self.seq_type == "PROTEIN":
            for i in self.seq:
                mass += amino_acid_weights[i]
        elif self.seq_type == "DNA":
            for i in self.seq:
                mass += nucleotide_weights[i]
        elif self.seq_type == "RNA":
            for i in self.seq:
                mass += rna_nucleotide_weights[i]

        mass -= (len(self.seq) - 1) * 18
        return round(mass,3)

    def melting_temp(self):
        """Calculates the melting temperature of a given DNA sequence and returns it."""
        tm = 0
        if self.seq_type == "DNA" or self.seq_type == "RNA":
            a = self.seq.count('A')
            g = self.seq.count('G')
            c = self.seq.count('C')
            if self.seq_type == "DNA":
                t = self.seq.count('T')
                tm = 2 * (a + t) + 4 * (g + c)
            else:
                u = self.seq.count('U')
                tm = 2 * (a + u) + 4 * (g + c)

        return tm

    def plot_nuc_frequency(self,bar_direction = "Horizontal"):
        """Plots nucleotide frequency of a given DNA sequence and plots it using matplotlib bar graph
        Parameters:
            bar_direction (str): The direction of bar graph. Defaults to 'Horizontal'.
            """
        frequency_dict = self.nucleotide_frequency
        seq_len = len(self.seq)
        frequency_list = [(frequency_dict.get(nuc, 0) / seq_len)*100 for nuc in NUCLEOTIDE_BASE["DNA"]]

        fig, ax = plt.subplots()
        if bar_direction.lower() == "horizontal":
            ax.barh(NUCLEOTIDE_BASE["DNA"], frequency_list)
        elif bar_direction.lower() == "vertical":
            ax.bar(NUCLEOTIDE_BASE["DNA"], frequency_list)
        ax.set_xlabel('Nucleotides')
        ax.set_ylabel('Frequency')
        ax.set_title('Bar Plot')

    @staticmethod
    def pairwise_alignment(seq1, seq2,alignerMode = "local"):
        """
         Perform pairwise sequence alignment between two sequences using Biopython's PairwiseAligner.

    Parameters:
        seq1 (str): First sequence to be aligned.
        seq2 (str): Second sequence to be aligned.
        alignerMode (str, optional): Alignment mode to use;
            "local" for local alignment (Smith-Waterman),
            "global" for global alignment (Needleman-Wunsch).
            Default is "local".

    Returns:
        alignments (Alignments): An iterable of alignment objects representing the optimal alignments.
        score (float): The alignment score between seq1 and seq2.
        """
        aligner = Align.PairwiseAligner()
        aligner.mode = alignerMode
        alignments = aligner.align(seq1, seq2)
        score = aligner.score(seq1, seq2)

        return alignments, score

    @staticmethod
    def calculate_hamming_distance(seq1,seq2):
        if len(seq1) != len(seq2):
            warnings.warn("Strings have different lengths, computing up to shortest length", UserWarning)
        return len([(n1, n2) for n1, n2 in zip(seq1, seq2) if n1 != n2])

    @staticmethod
    def find_longest_subsequence(sequences:List[str],)-> str:

        shortest = min(map(len, sequences))

        common = functools.partial(BioSeq.common_kmers, sequences)
        start = BioSeq.binary_search(common, 1, shortest)

        if start >= 0:
            # Hill climb to find max
            candidates = []
            for k in range(start, shortest + 1):
                if kmers := common(k):
                    candidates.extend(kmers)
                else:
                    break

            return str(max(candidates, key=len))
        else:
            raise NoCommonKmersError("There are not any common kmers in sequences.")

    @staticmethod
    def binary_search(f: Callable, low: int, high: int) -> int:
        """ Binary search """

        hi, lo = f(high), f(low)
        mid = (high + low) // 2

        if hi and lo:
            return high

        if lo and not hi:
            return BioSeq.binary_search(f, low, mid)

        if hi and not lo:
            return BioSeq.binary_search(f, mid, high)

        return -1

    @staticmethod
    def common_kmers(seqs: List[str], k: int) -> List[str]:
        """ Find k-mers common to all elements """

        kmers = [set(BioSeq.__find_kmers(seq, k)) for seq in seqs]
        counts = collections.Counter(chain.from_iterable(kmers))
        n = len(seqs)
        return [kmer for kmer, freq in counts.items() if freq == n]

    def plot_amino_acid_frequency(self, bar_direction="horizontal"):
        """Plots amino acid frequency of a given Sequence"""
        if self.seq_type != "PROTEIN":
            aa_seq = self.translate_seq()
            aa_seq_frequency = dict(collections.Counter(aa_seq))
            seq_len = len(aa_seq)
        else:
            aa_seq_frequency = dict(collections.Counter(self.seq))
            seq_len = len(self.seq)

        frequency_list = [(aa_seq_frequency.get(aa, 0) / seq_len) * 100 for aa in AMINOACID_BASE]

        fig, ax = plt.subplots()
        if bar_direction.lower() == "horizontal":
            ax.barh(sorted(AMINOACID_BASE), frequency_list)
        elif bar_direction.lower() == "vertical":
            ax.bar(sorted(AMINOACID_BASE), frequency_list)
        ax.set_xlabel('AminoAcids')
        ax.set_ylabel('Frequency')
        ax.set_title('Bar Plot')




    













    
