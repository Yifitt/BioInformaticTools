NUCLEOTIDE_BASE = {
    "DNA":["A","T","G","C","a","t","g","c"],
    "RNA":["A","U","G","C","a","u","g","c"],
                   }

AMINOACID_BASE = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y",
                  "_"}

DNA_Codons = {
    # 'M' = START, '_' = STOP
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A", # A = Alanine
    "TGT": "C", "TGC": "C", # C = Cysteine
    "GAT": "D", "GAC": "D", # D = Aspartic Acid
    "GAA": "E", "GAG": "E", # E = Glutamic Acid
    "TTT": "F", "TTC": "F", # F = Phenylalanine
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G", # G = Glycine
    "CAT": "H", "CAC": "H", # H = Histidine
    "ATA": "I", "ATT": "I", "ATC": "I", # I = Isoleucine
    "AAA": "K", "AAG": "K", # K = Lysine
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", #L = Leucine
    "ATG": "M", # M = Methionine(Start Codon)
    "AAT": "N", "AAC": "N", # N = Asparagine
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P", # P = Proline
    "CAA": "Q", "CAG": "Q", # Q = Glutamine
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R", # R = Arginine
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S", # S = Serine
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", # T = Threonine
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", # V = Valine
    "TGG": "W", # W = Tryptophan
    "TAT": "Y", "TAC": "Y", # Y = Tyrosine
    "TAA": "_", "TAG": "_", "TGA": "_"

}
RNA_Codons = {
    # 'M' = START, '_' = STOP
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",#Alanine
    "UGU": "C", "UGC": "C", #Cysteine
    "GAU": "D", "GAC": "D",#Aspartic Acid
    "GAA": "E", "GAG": "E", #Glutamic Acid
    "UUU": "F", "UUC": "F", #Phenylalanine
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G", #Glycine
    "CAU": "H", "CAC": "H", #Histidine
    "AUA": "I", "AUU": "I", "AUC": "I", # I = Isoleucine
    "AAA": "K", "AAG": "K", # K = Lysine
    "UUA": "L", "UUG": "L", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L", #L = Leucine
    "AUG": "M", # M = Methionine(Start Codon)
    "AAU": "N", "AAC": "N", # N = Asparagine
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P", # P = Proline
    "CAA": "Q", "CAG": "Q", # Q = Glutamine
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R", # R = Arginine
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "AGU": "S", "AGC": "S", # S = Serine
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T", # T = Threonine
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V", # V = Valine
    "UGG": "W", # W = Tryptophan
    "UAU": "Y", "UAC": "Y", # Y = Tyrosine
    "UAA": "_", "UAG": "_", "UGA": "_"

}

amino_acid_weights = {
    "A": 71.03711,
    "C": 103.00919,
    "D": 115.02694,
    "E": 129.04259,
    "F": 147.06841,
    "G": 57.02146,
    "H": 137.05891,
    "I": 113.08406,
    "K": 128.09496,
    "L": 113.08406,
    "M": 131.04049,
    "N": 114.04293,
    "P": 97.05276,
    "Q": 128.05858,
    "R": 156.10111,
    "S": 87.03203,
    "T": 101.04768,
    "V": 99.06841,
    "W": 186.07931,
    "Y": 163.06333
}

nucleotide_weights = {
    "A": 313.21,
    "T": 304.2,
    "C": 289.18,
    "G": 329.21,
}


rna_nucleotide_weights = {
    "A": 329.21,
    "U": 306.17,
    "C": 305.18,
    "G": 345.21,
}

AminoAcids = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
"M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"}
    
