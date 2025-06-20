from Python.myclasses.bio_structure import *
import random
import collections
import re

class bio_seq:
    """DNA sequence class.Default value: ATCG, DNA, No label"""
    def __init__(self, seq="ATCG", seq_type="DNA", label='No Label'):
        """Sequence initialization, validation."""
        self.seq = seq.upper()
        self.label = label
        self.seq_type = seq_type
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

    def countNucFrequency(self):
        """Count nucleotides in a given sequence. Return a dictionary"""
        return dict(collections.Counter(self.seq))
    
    def transcription(self):
        """DNA -> RNA Transcription.Replacing Thymine with Uracil"""
        if self.seq_type == "DNA":
            return self.seq.replace("T", "U")
        return "Not a DNA sequence"
    
    def reverseComplement(self):
        """Swapping Adenine->Thymine and Guanine->Cytosine.Reversing newly generated string."""
        if self.seq_type == "DNA":
            mapping = str.maketrans('ATCG','TAGC') #bir ceviri tablosu yaratir.Her harfi eslestirir
        else:
            mapping = str.maketrans('AUCG','UAGC')
        return self.seq.translate(mapping)[::-1] #her harfi eslesenle degisir
    
    def Complement(self):
        """Swapping Adenine->Thymine and Guanine->Cytosine."""
        if self.seq_type == "DNA":
            mapping = str.maketrans('ATCG','TAGC')
        else:
            mapping = str.maketrans('AUCG','UAGC')
        return self.seq.translate(mapping) 
    
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
    
    def proteins_from_rf(self, aa_seq,):
        """Compute all possible proteins in an aminoacid seq and return a list of possible proteins"""
        current_prot = []
        proteins = []
        for aa in aa_seq:
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

    def all_proteins_from_orfs(self, startReadPos=0, endReadPos=0, ordered=False,Largest= False):
        """Compute all possible proteins for all open reading frames."""
        """'ordered' and 'Largest' are set to False by default."""
        """'Largest' returns the largest protein in that sequence when it's set to 'True'"""
               
        if endReadPos > startReadPos:
            tmp_seq = bio_seq(
                self.seq[startReadPos: endReadPos], self.seq_type)
            rfs = tmp_seq.gen_reading_frames()
        else:
            rfs = self.gen_reading_frames()

        res = []
        for rf in rfs:
            prots = self.proteins_from_rf(rf)
            for p in prots:
                res.append(p)

        if ordered:
            if Largest:
                return sorted(res, key=len, reverse=True)[0]
            return sorted(res, key=len, reverse=True)
        return res

    def find_motif_positions(self,sequence,kmer):
        """With this function we can find the position of mutations or other kmer start positions"""

        positions = ''
        
        for i in range(len(sequence)):
            if sequence[i] == kmer[0]:
                if sequence[i:i+len(kmer)] == kmer:
                    positions += str(i+1)+' '
        return positions
        
    def count_kmer(self,sequence,kmer):
        """
        Counts the number of times a specific k-mer appears in a given sequence,including overlapping k-mers.
            
        Parameters:
        sequence (str): The DNA sequence to search in.
        kmer (str): The specific k-mer to search for in the sequence.

        Returns:
        int: The number of times the k-mer appears in the sequence.
        
        """
        #return len(re.findall(f'(?=({kmer}))', sequence)) #Kisa yol regular expressions kullaniliyor
        
        kmer_count = 0
        for position in range(len(sequence) - (len(kmer)-1)):
            print(sequence[position:position+len(kmer)] + "=" + kmer)
            if sequence[position:position+len(kmer)] == kmer:
               kmer_count+=1
        return kmer_count



        
    def find_most_frequent_kmers(self, sequence, k_len):
        """
        Finds the most frequent k-mers of a given length in a DNA string.

        Parameters:
            sequence (str): The DNA string to search.
            k_len (int): The length of the k-mers to search for.

        Returns:
            list: A list of the most frequent k-mers in the DNA string.
        """
        # 1. A dictionary to store k-mer frequencies.
        kmer_frequencies = {}

        # 2. A loop to iterate through the DNA string and exract k-mers of a given length,
        # while also incrementing the frequency of each k-mer
        for i in range(len(sequence) - k_len + 1):
            kmer = sequence[i:i+k_len]
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
    

    def calculate_protein_mass(self,sequence):
        """Calculates the mass of a protein and prints it."""
        #We assign a mass variable in order to store the mass of the amino-acids
        mass = 0

        for i in sequence:
            mass += amino_acid_weights[i]
            
        return round(mass,3)


    