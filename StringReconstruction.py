
### StringReconstruction() ###
# Code Challenge: Solve the String Reconstruction Problem.
    # Input: An integer k followed by a list of k-mers Patterns.
    # Output: A string Text with k-mer composition equal to Patterns. (If multiple answers exist, you may return any one.)
# pseudocode:
#     StringReconstruction(Patterns)
#         dB <- DeBruijn(Patterns)
#         path <- EulerianPath(dB)
#         Text <- PathToGenome(path)
#         return Text
#
# Sample Input:
#     4
#    kmers = ['CTTA', 'ACCA', 'TACC', 'GGCT', 'GCTT', 'TTAC']
# Sample Output:
#     GGCTTACCA

import random

### Part A - data

# with open('dataset_203_7.txt', 'r') as file:
#     k = int(file.readline().strip())
#     kmers = [x.strip() for x in file.readlines()]


kmers = ['AAAT', 'AATG', 'ACCC', 'ACGC', 'ATAC', 'ATCA', 'ATGC', 'CAAA', 'CACC','CATA', 'CATC', 'CCAG',
         'CCCA', 'CGCT', 'CTCA', 'GCAT', 'GCTC', 'TACG', 'TCAC', 'TCAT', 'TGCA']



edges = []
for kmer1 in kmers:
    for kmer2 in [x for x in kmers if x != kmer1]:
        if kmer1[:-1] == kmer2[1:]:
            edges.append([kmer1, kmer2])

def get_adjacent_edges(pos, edges):
    choices = []
    for edge in edges:
        if pos == edge[0]:
            choices.append(edge[1])
    return choices


def check_edges(edges):
    start = [x[0] for x in edges]
    end = [x[1] for x in edges]
    vertices = set(start).union(set(end))
    in_dict = {x: 0 for x in vertices}
    out_dict = {x: 0 for x in vertices}
    for edge in edges:
        out_dict[edge[0]] += 1
        in_dict[edge[1]] += 1
    for vertex in vertices:
        if out_dict[vertex] > in_dict[vertex]:
            x = vertex
    if x:
        return x
    else:
        return random.sample(vertices, 1)


stack = []
circuit = []

pos = check_edges(edges)

while len(edges) > 0:
    choices = get_adjacent_edges(pos, edges)
    while choices == []:
        circuit.append(pos)
        pos = stack.pop()
        choices = get_adjacent_edges(pos, edges)
    stack.append(pos)
    pos = random.choice(choices)
    edges.remove([stack[-1], pos])
stack.append(pos)

while len(stack) > 0:
    circuit.append(stack.pop())

circuit = list(map(str, circuit))

s = ''
s = s + circuit[0]
for x in circuit[1:]:
    s = s + x[-1]

output = open('output.txt', 'w')
output.write(s)