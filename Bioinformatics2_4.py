import Bioinformatics2_3 as BI2_3
### CyclopeptideScor(peptide, spectrum) ###
# Cyclopeptide Scoring Problem: Compute the score of a cyclic peptide against a spectrum.
#   Input: An amino acid string peptide and a collection of integers spectrum.
#   Output: The score of Peptide against Spectrum, Score(Peptide, Spectrum).
# Sample Input:
#    NQEL
#    [0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484]
# Sample Output: 11

def CyclicSpectrumDict(peptide, massDict = BI2_3.AminoAcidMass()):
    n = len(peptide)
    PrefixMass = [0]
    for i in range(n):
        PrefixMass.append(PrefixMass[i] + massDict[peptide[i]])
    peptideMass = PrefixMass[n]
    cSpectrumDict = {0: 1}
    for i in range(n):
        for j in range(i + 1, n + 1):
            s = PrefixMass[j] - PrefixMass[i]
            cSpectrumDict[s] = cSpectrumDict.get(s, 0) + 1
            if i > 0 and j < n:
                s = peptideMass - (PrefixMass[j] - PrefixMass[i])
                cSpectrumDict[s] = cSpectrumDict.get(s, 0) + 1
    return (cSpectrumDict)

def CyclopeptideScor(peptide, spectrum):
    theoSpectrumDict = CyclicSpectrumDict(peptide)
    score = 0
    spectrumDict = dict()
    for s in spectrum:
        spectrumDict[s] = spectrumDict.get(s, 0) + 1
    for s, v in theoSpectrumDict.items():
        v0 = spectrumDict.get(s, 0)
        if v0 >= v:
            score += v
        else:
            score += v0
    return score

# ---------------------------------------------------------------------------------------------
