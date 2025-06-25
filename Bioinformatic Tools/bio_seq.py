from bio_structure import *
import random
import collections
import re
import matplotlib.pyplot as plt

class bio_seq:
    """DNA sequence class.Default value: ATCG, DNA, No label"""
    def __init__(self, seq="ATCG", seq_type="DNA", label='No Label'):
        """Sequence initialization, validation."""
        self.seq = seq.upper()
        self.label = label
        self.seq_type = seq_type.upper()
        self.is_valid = self.__validateSeq()
        assert self.is_valid, f"Provided data does not seem to be a correct {self.seq_type} sequence."

    #DNA Toolkit Functions

    def __validateSeq(self): # basina 2 tane '_' ekleyince private olur
        """Check the sequence to make sure it is a valid DNA sequence."""
        if self.seq_type == "RNA" or self.seq_type == "DNA":
            return set(NUCLEOTIDE_BASE[self.seq_type]).issuperset(self.seq) # bool dondurur
        elif self.seq_type == "PROTEIN":
            return AminoAcids.issuperset(self.seq)
        else:
            return False
    
    def get_seq_biotype(self):
        """Returns sequence type"""
        return self.seq_type

    def get_seq_info(self):
        """Returns 4 strings. Full sequence information"""
        return f"[Label]: {self.label}\n[Sequence]: {self.seq}\n[Biotype]: {self.seq_type}\n[Length]: {len(self.seq)}"

    def generate_rnd_seq(self, length=10, seq_type="DNA"):
        """Generate a random DNA sequence, provided the length"""
        seq = ''.join([random.choice(NUCLEOTIDE_BASE[seq_type])
                       for x in range(length)])
        self.__init__(seq, seq_type, "Randomly generated sequence")

    @property
    def nucleotide_frequency(self):
        """Count nucleotides in a given sequence. Return a dictionary"""
        return dict(collections.Counter(self.seq))
    
    def transcription(self):
        """DNA -> RNA Transcription.Replacing Thymine with Uracil"""
        return self.seq.replace("T", "U") if self.seq_type == "DNA" else None

    def backtranscription(self):
        """RNA -> DNA Replacing Uracil with Thymine"""
        return self.seq.replace("U", "T") if self.seq_type == "RNA" else None

    def reverseComplement(self):
        """Swapping Adenine->Thymine and Guanine->Cytosine.Reversing newly generated string."""
        if self.seq_type == "DNA":
            mapping = str.maketrans('ATCG','TAGC') #bir ceviri tablosu yaratir.Her harfi eslestirir
        elif self.seq_type == "RNA":
            mapping = str.maketrans('AUCG','UAGC')
        return self.seq.translate(mapping)[::-1] #her harfi eslesenle degisir
    
    def Complement(self):
        """Swapping Adenine->Thymine and Guanine->Cytosine."""
        if self.seq_type == "DNA":
            mapping = str.maketrans('ATCG','TAGC')
        else:
            mapping = str.maketrans('AUCG','UAGC')
        return self.seq.translate(mapping)

    @property
    def gc_content(self):
        """ GC Content in a DNA/RNA sequence"""
        return round((self.seq.count('C') + self.seq.count('G'))/ len(self.seq)*100)
    
    def gc_content_subsec(self, k=20):
        """GC Content in a DNA/RNA sub-sequence length k. k=20 by default"""
        res = []
        for i in range(0, len(self.seq) - k + 1, k):
            subseq = self.seq[i:i + k]
            res.append(
                round((subseq.count('C') + subseq.count('G')) / len(subseq) * 100))
        return res
    def translate_seq(self,init_pos = 0):
        """Translates a DNA sequence into an aminoacid sequence."""
        if self.seq_type == "DNA":
            return[DNA_Codons[self.seq[pos:pos+3]]for pos in range(init_pos,len(self.seq)-2,3)]
        elif self.seq_type == "RNA":
            return[RNA_Codons[self.seq[pos:pos+3]]for pos in range(init_pos,len(self.seq)-2,3)]
        else:
           return f"Provided data does not seem to be appropriate sequence."

    def codon_usage(self, aminoacid):
        """Provides the frequency of each codon encoding a given aminoacid in a DNA sequence"""
        tmpList = []
        if self.seq_type == "DNA":
            for i in range(0, len(self.seq) - 2, 3):
                if DNA_Codons[self.seq[i:i + 3]] == aminoacid:
                    tmpList.append(self.seq[i:i + 3])

        elif self.seq_type == "RNA":
            for i in range(0, len(self.seq) - 2, 3):
                if RNA_Codons[self.seq[i:i + 3]] == aminoacid:
                    tmpList.append(self.seq[i:i + 3])

        freqDict = dict(collections.Counter(tmpList))
        totalWight = sum(freqDict.values())
        for seq in freqDict:
            freqDict[seq] = round(freqDict[seq] / totalWight, 2)
        return freqDict
    
    def gen_reading_frames(self):
        """Generate the six reading frames of a DNA sequence, including reverse complement"""
        frames = []
        frames.append(self.translate_seq(0))
        frames.append(self.translate_seq(1))
        frames.append(self.translate_seq(2))
        tmp_seq = bio_seq(self.reverseComplement(), self.seq_type)
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


    def all_proteins_from_orfs(self, start=0, end=0, ordered=False,largest= False):
        """Compute all possible proteins for all open reading frames."""
        """'ordered' and 'Largest' are set to False by default."""
        """'Largest' returns the largest protein in that sequence when it's set to 'True'"""
        if end > start:
            tmp_seq = bio_seq(
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

    def find_motif_positions(self,kmer):
        """Finds the starting positions of given kmer in a sequence.
        Parameters:
            kmer (str): The specific k-mer to search for in the sequence.
        """
        return[match.start() for match in re.finditer(f'(?=({kmer}))', self.seq)]


    def count_kmer(self,kmer):
        """
        Counts the number of times a specific k-mer appears in a given sequence,including overlapping k-mers.
            
        Parameters:
        kmer (str): The specific k-mer to search for in the sequence.

        Returns:
        int: The number of times the k-mer appears in the sequence.
        
        """
        return len(re.findall(f'(?=({kmer}))',self.seq)) #Kisa yol regular expressions kullaniliyor

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

        # 2. A loop to iterate through the DNA string and exract k-mers of a given length,
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
    def calculate_protein_mass(self):
        """Calculates the mass of a protein and returns it.
         Returns:
            mass: The rounded mass of the protein.
            """
        #We assign a mass variable in order to store the mass of the amino-acids
        mass = 0
        if self.seq_type == "PROTEIN":
            for i in self.seq:
                mass += amino_acid_weights[i]
        else:
            return f"This method is for protein sequences only"
        return round(mass,3)

    def melting_temp(self):
        """Calculates the melting temperature of a given DNA sequence and returns it."""
        tm = 0
        if self.seq_type == "DNA" or self.seq_type == "RNA":
            A = self.seq.count('A')
            G = self.seq.count('G')
            C = self.seq.count('C')
            if self.seq_type == "DNA":
                T = self.seq.count('T')
                tm = 2 * (A + T) + 4 * (G + C)
            else:
                U = self.seq.count('U')
                tm = 2 * (A + U) + 4 * (G + C)

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
        plt.show()
