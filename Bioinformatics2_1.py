### Composition ###

# Code Challenge: Solve the String Composition Problem.
# Input: An integer k and a string Text.
# Output: Composition(k, Text) (the k-mers can be provided in any order).
# Sample Input:
#     k = 5
#     Text = 'CAATCCAAC'
# Sample Output:
#     kmers =
#         'CAATC
#         AATCC
#         ATCCA
#         TCCAA
#         CCAAC'

def Composition(k, Text):
    kmersList = list()
    for i in range(len(Text)-k+1):
        kmersList.append(Text[i:i+k])
    return kmersList

# ---------------------------------------------------------------------------------------------
### PathToGenome(Path) ###
# Code Challenge: Solve the String Spelled by a Genome Path Problem.
# Sample Input:
#   ACCGA
#   CCGAA
#   CGAAG
#   GAAGC
#   AAGCT
# Sample Output:
#   ACCGAAGCT

def PathToGenome(Path):
    Genome = Path[0]
    for i in range(1,len(Path)):
        Genome += Path[i][-1]
    return Genome

# ---------------------------------------------------------------------------------------------
### PathToGenome(Path) ###
# Code Challenge: Solve the Overlap Graph Problem (restated below).
#      Input: A collection Patterns of k-mers.
#      Output: The overlap graph Overlap(Patterns), in the form of an adjacency list. (You may return the nodes and
#      their edges in any order.)
#
# Sample Input:
#     ATGCG
#     GCATG
#     CATGC
#     AGGCA
#     GGCAT
#     GGCAC
# Sample Output:
#     ATGCG -> []
#     CATGC -> [ATGCG]
#     GCATG -> [CATGC]
#     GGCAT -> [GCATG]
#     AGGCA -> [GGCAC, GGCAT]
#     GGCAC -> []

def OverlapGraph(kmers):
    kmerDict = dict() # create empty dictionary of kmers
    for kmer in kmers: kmerDict.setdefault(kmer, [])
    for i in range(len(kmers)):
        suffix = kmers[i][1:]
        for j in range( len(kmers) ) :
            prefix = kmers[j][:-1]
            if prefix == suffix : kmerDict[ kmers[i] ].append( kmers[j] )
    return kmerDict

# ---------------------------------------------------------------------------------------------
### DeBruijnGraph(k, Text) ###
# Code Challenge: Solve the De Bruijn Graph from a String Problem.
#      Input: An integer k and a string Text.
#     Output: DeBruijn(k, Text), in the form of an adjacency list.

# Sample Input:
#   4
#   AAGATTCTCTAAGA
# Sample Output:
#   AAG -> AGA,AGA
#   AGA -> GAT
#   ATT -> TTC
#   CTA -> TAA
#   CTC -> TCT
#   GAT -> ATT
#   TAA -> AAG
#   TCT -> CTA,CTC
#   TTC -> TCT

def DeBruijnGraph(k, Text):

    PathGraph = dict() # Adjacency list - empty dict
    EdgesList = list() # Edge = k_mer
    for i in range(len(Text) - k + 1): EdgesList.append(Text[i:i+k])
    for kmer in EdgesList:
        if kmer[0:k-1] in PathGraph:
            PathGraph[kmer[0:k - 1]].append(kmer[1:k])
        else:
            PathGraph[kmer[0:k-1]] = [] # prefix = key in PathGraph
            PathGraph[kmer[0:k-1]].append(kmer [1:k])
    return PathGraph

# ---------------------------------------------------------------------------------------------
### DeBruijnGraph_Anotherway(kmers) ###
# DeBruijn Graph from k-mers Problem: Construct the de Bruijn graph from a set of k-mers.
#      Input: A collection of k-mers Patterns.
#      Output: The adjacency list of the de Bruijn graph DeBruijn(Patterns).

# Sample Input: ['GAGG', 'CAGG', 'GGGG', 'GGGA', 'CAGG', 'AGGG', 'GGAG']
# Sample Output:
#     AGG -> GGG
#     CAG -> AGG,AGG
#     GAG -> AGG
#     GGA -> GAG
#     GGG -> GGA,GGG

def DeBruijnGraph_Anotherway(kmers):
    Graph = dict() # empty dict
    for kmer in kmers:
        if kmer[:-1] in Graph:  # if prefix is key in Graph
            Graph[kmer[-1]].append(kmer[1:])  # add suffix
        else: # if prefix is not a key in Graph
            Graph[kmer[:-1]] = {kmer[1:]}  ## add prefix (key) and suffix
    return Graph

