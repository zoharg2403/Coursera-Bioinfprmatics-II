### GeneticCode() ###
# genetic code represented as an array containing 64 elements
def GeneticCode():
    return {'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAU': 'N', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
            'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGU': 'S', 'AUA': 'I', 'AUC': 'I', 'AUG': 'M', 'AUU': 'I',
            'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAU': 'H', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R', 'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
            'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAU': 'D', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G', 'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
            'UAA': '*', 'UAC': 'Y', 'UAG': '*', 'UAU': 'Y', 'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
            'UGA': '*', 'UGC': 'C', 'UGG': 'W', 'UGU': 'C', 'UUA': 'L', 'UUC': 'F', 'UUG': 'L', 'UUU': 'F'}


# ---------------------------------------------------------------------------------------------
### complementDNAstrand(s) ###
# Input: DNA sequence string s
# Outputn: complement DNA string

def ReversedComplementDNA(s):
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    letters = list(s)
    letters = [basecomplement[base] for base in letters]
    return ''.join(letters[::-1])


# ---------------------------------------------------------------------------------------------
### Translation(RNAseq, genetic_code) ###
# Protein Translation Problem: Translate an RNA string into an amino acid string
#       Input: An RNA sequence string and the dictionary genetic_code
#              (Default Argument: genetic_code = GeneticCode())
#       Output: The translation of Pattern into an amino acid string Peptide
# Sample Input: 'AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'
# Sample Output: 'MAMAPRTEINSTRING'

def Translation(RNAseq, genetic_code=GeneticCode()):
    AAseq = ''
    for i in range(0, len(RNAseq) - 2, 3):
        AAseq += genetic_code[RNAseq[i:i + 3]]
    return AAseq


# ---------------------------------------------------------------------------------------------
### PeptideEncoding(DNAseq, Peptide) ###
# Peptide Encoding Problem: Find substrings of a genome encoding a given amino acid sequence
#       Input: A DNA string DNAseq, an amino acid string Peptide, and the array genetic_code
#              (Default Argument: genetic_code = GeneticCode())
#       Output: All substrings of DNAseq encoding Peptide (if any such substrings exist)
# Sample Input:
#               ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA
#               MA
# Sample Output: ['ATGGCC', 'GGCCAT', 'ATGGCC']

def PeptideEncoding(DNAseq, Peptide):
        substrings = list()
        for i in range(len(DNAseq) - len(Peptide) * 3 + 1):
                sub_seq = DNAseq[i:i+len(Peptide) * 3]
                if Translation(sub_seq.replace('T', 'U')) == Peptide or Translation(ReversedComplementDNA(sub_seq).replace('T', 'U')) == Peptide:
                        substrings.append(sub_seq)
        return substrings

# ---------------------------------------------------------------------------------------------
### TheoreticalSpectrum(Peptide, AAmass) ###
# Generating Theoretical Spectrum Problem: Generate the theoretical spectrum of a cyclic peptide
#       Input: An amino acid string Peptide (represents a cyclic peptide)
#       Output: The theoretical spectrum Peptide = the collection of all of the masses of its subpeptides
#       (including 0 and the mass of the entire peptide) ordered from smallest to largest (list)
# Sample Input: 'LEQN'
# Sample Output: [0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370, 371, 484]

def AminoAcidMass():
    # keys are the symbols of Alphabet (amino acids), values are the masses of each amino acid
    return {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114,
            'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}

def CyclicSpectrum(peptide, AAmass = AminoAcidMass()):
    n = len(peptide)
    PrefixMass = [0]
    for i in range(n):
        PrefixMass.append(PrefixMass[i] + AAmass[peptide[i]])
    peptideMass = PrefixMass[n]
    cSpectrum = [0]
    for i in range(n):
        for j in range(i+1, n+1):
            cSpectrum.append(PrefixMass[j]-PrefixMass[i])
            if i > 0 and j < n:
                cSpectrum.append(peptideMass-(PrefixMass[j]-PrefixMass[i]))
    return sorted(cSpectrum)

# ---------------------------------------------------------------------------------------------
### BFCyclopeptideSequencing(Spectrum) ###
# Counting Peptides with Given Mass Problem: Compute the number of peptides of given mass
#    Input: An integer n
#    Output: The number of linear peptides having integer mass n
# Sample Input: 1024
# Sample Output: 14712706211

def AminoAcidMass_Table():
    '''
    G	A	S	P	V	T	C	I/L	N	D	K/Q	E	M	H	F	R	Y	W
    57	71	87	97	99	101	103	113	114	115	128	129	131	137	147	156	163	186
    '''
    mass = '57 71 87 97 99 101 103 113 114 115 128 129 131 137 147 156 163 186'
    return [int(m) for m in mass.split()]

def CountingPeptides(n):
    massTalbe = AminoAcidMass_Table()
    m = len(massTalbe)
    table = [0] * (n+1)
    table[0] = 1
    for i in range(n+1):
        currSum = 0
        for j in range(m):
            if i - massTalbe[j] >= 0:
                currSum += table[i-massTalbe[j]]
        table[i] += currSum
    return table[n]




