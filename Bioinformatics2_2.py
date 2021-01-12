import random
import copy
from functools import reduce
from itertools import product

### HamiltonianCycle(Graph) ###
# Input: The adjacency list of a directed graph.
# Output: An Hamiltonian cycle in this graph.
#
# Sample Input:
#     0 -> 3
#     1 -> 0
#     2 -> 1,6
#     3 -> 2
#     4 -> 2
#     5 -> 4
#     6 -> 5,8
#     7 -> 9
#     8 -> 7
#     9 -> 6
# Sample Output:
# [8, 7, 9, 6, 5, 4, 2, 1, 0, 3]

def HamiltonianCycle(Graph):
    nodes = list(Graph.keys())
    Cycle = list()

    all_nodes_flag = 0
    while not all_nodes_flag:  # while there are unexplored edges in Graph
        # initialization - starting node
        node_visit_mark = [0]*len(nodes)  # mark which nodes were visited -> 1 = visited
        Cycle.append(random.choice(nodes))  # random start node
        node_visit_mark[nodes.index(Cycle[0])] = 1
        # other nodes
        next_node = random.choice(Graph[Cycle[0]])
        while not node_visit_mark[nodes.index(next_node)]:  # while next was not visited
            Cycle.append(next_node)
            node_visit_mark[nodes.index(Cycle[-1])] = 1
            next_node = random.choice(Graph[Cycle[-1]])
        if all(node_visit_mark):
            all_nodes_flag = 1
        else:
            Cycle = list()
    return Cycle

# ---------------------------------------------------------------------------------------------
### EulerianCycle(Graph) ###
# Code Challenge: Solve the Eulerian Cycle Problem.
#      Input: The adjacency list of an Eulerian directed graph.
#      Output: An Eulerian cycle in this graph.
#
# Sample Input:
#     0 -> 3
#     1 -> 0
#     2 -> 1,6
#     3 -> 2
#     4 -> 2
#     5 -> 4
#     6 -> 5,8
#     7 -> 9
#     8 -> 7
#     9 -> 6
#    ->     Graph = {0:[3], 1:[0], 2:[1. 6], 3:[2], 4:[2], 5:[4], 6:[5, 8], 7:[9], 8:[7], 9:[6]}
# Sample Output:
#     6->8->7->9->6->5->4->2->1->0->3->2->6

def EulerianCycle (Graph):
    # create edges list
    edges = []
    for n in list(Graph.keys()):
        for e in Graph[n]:
            edges.append((n, e))

    Cycle = []
    all_edges_flag = 0
    while not all_edges_flag:  # while there are unexplored edges in Graph
        # initialization - starting edge
        edge_visit_mark = [0]*len(edges)  # mark which edges were visited -> 1 = visited
        start_edge = random.choice(edges)  # random start edge
        Cycle.append(start_edge[0])
        Cycle.append(start_edge[1])
        edge_visit_mark[edges.index(start_edge)] = 1
        # other edges
        next_edge = (Cycle[-1], random.choice(Graph[Cycle[-1]]))
        while not edge_visit_mark[edges.index(next_edge)]:  # while next edge was not visited
            Cycle.append(next_edge[1])
            edge_visit_mark[edges.index(next_edge)] = 1
            next_edge = (Cycle[-1], random.choice(Graph[Cycle[-1]]))
        if all(edge_visit_mark):  # and Graph[Cycle[-1]] == Cycle[0]:
            all_edges_flag = 1
        else:
            Cycle = list()
    return Cycle

## another version
def EulerianCycle_v2(Graph):
    edge_dict = copy.deepcopy(Graph)

    # initializing - choosing random start node
    cur_node = random.choice(list(edge_dict.keys()))
    Cycle = [cur_node]

    flag = 1
    while flag:
        Cycle.append(edge_dict[cur_node][0])  # add next node to Cycle

        # deleting the edge visited from edge_dict, if the node has no more edges -> deleting the node
        if len(edge_dict[cur_node]) == 1:
            del edge_dict[cur_node]
        else:
            edge_dict[cur_node] = edge_dict[cur_node][1:]

        # checking the cycle:
        if Cycle[-1] in edge_dict:  # if next node has more edges (its in edge_dict), update cur_node and continue
            cur_node = Cycle[-1]
        else:  # if next node has no more edges
            if len(edge_dict) > 0:  # if there are unvisited edges, restart new Cycle
                edge_dict = copy.deepcopy(Graph)
                cur_node = random.choice(list(edge_dict.keys()))
                Cycle = [cur_node]
                print('B')
            else: # if there are no unvisited edges
                flag = 0
    return Cycle

## fast version
def EulerianCycle_fast(edge_dict):
    '''Generates an Eulerian cycle from the given edges.'''
    current_node = list(edge_dict.keys())[0]
    path = [current_node]
    # Get the initial cycle.
    while True:
        path.append(edge_dict[current_node][0])

        if len(edge_dict[current_node]) == 1:
            del edge_dict[current_node]
        else:
            edge_dict[current_node] = edge_dict[current_node][1:]

        if path[-1] in edge_dict:
            current_node = path[-1]
        else:
            break

    # Continually expand the initial cycle until we're out of edge_dict.
    while len(edge_dict) > 0:
        for i in range(len(path)):
            if path[i] in edge_dict:
                current_node = path[i]
                cycle = [current_node]
                while True:
                    cycle.append(edge_dict[current_node][0])

                    if len(edge_dict[current_node]) == 1:
                        del edge_dict[current_node]
                    else:
                        edge_dict[current_node] = edge_dict[current_node][1:]

                    if cycle[-1] in edge_dict:
                        current_node = cycle[-1]
                    else:
                        break

                path = path[:i] + cycle + path[i+1:]
                break
    return path

# ---------------------------------------------------------------------------------------------
### EulerianPath(Graph) ###
#Code Challenge: Solve the Eulerian Path Problem.
#      Input: The adjacency list of a directed graph that has an Eulerian path.
#      Output: An Eulerian path in this graph.

# Sample Input:
#     0 -> 2
#     1 -> 3
#     2 -> 1
#     3 -> 0,4
#     6 -> 3,7
#     7 -> 8
#     8 -> 9
#     9 -> 6
#    ->     Graph = {0:[2], 1:[3], 2:[1], 3:[0,4], 6:[3, 7], 7:[8], 8:[9], 9:[6]}
# Sample Output:
#     6->7->8->9->6->3->0->2->1->3->4
#    ->     [6, 7, 8, 9, 6, 3, 0, 2, 1, 3, 4]

def EulerianPath(edge_dict):
    # Determine the unbalanced edges.
    out_values = reduce(lambda a, b: a + b, edge_dict.values())
    for node in set(out_values + list(edge_dict.keys())):
        out_value = out_values.count(node)
        if node in edge_dict:
            in_value = len(edge_dict[node])
        else:
            in_value = 0

        if in_value < out_value:
            unbalanced_from = node
        elif out_value < in_value:
            unbalanced_to = node

    # Add an edge connecting the unbalanced edges.
    if unbalanced_from in edge_dict:
        edge_dict[unbalanced_from].append(unbalanced_to)
    else:
        edge_dict[unbalanced_from] = [unbalanced_to]

    # Get the Eulerian Cycle from the edges, including the unbalanced edge.
    cycle = EulerianCycle_fast(edge_dict)

    # Find the location of the unbalanced edge in the eulerian cycle.
    divide_point = list(filter(lambda i: cycle[i:i + 2] == [unbalanced_from, unbalanced_to], range(len(cycle) - 1)))[0]

    # Remove the unbalanced edge, and shift appropriately, overlapping the head and tail.
    return cycle[divide_point + 1:] + cycle[1:divide_point + 1]

# ---------------------------------------------------------------------------------------------
### universal_circular_string(k) ###
# Code Challenge: Solve the k-Universal Circular String Problem.
#      Input: An integer k.
#      Output: A k-universal circular string.
# Sample Input: 4
# Sample Output: 0000110010111101

def universal_circular_string (k):
    # Create edges dict for binary k mers (2^k k mers)
    universal_dict = {}
    for kmer in [''.join(item) for item in product('01', repeat=k)]:
        if kmer[:-1] in universal_dict:
            universal_dict[kmer[:-1]].append(kmer[1:])
        else:
            universal_dict[kmer[:-1]] = [kmer[1:]]

    # Get the cycle:
    path = EulerianCycle_fast(universal_dict)
    # remove the repeated last entry for the associated path (last and first entry are overlapping):
    path_string = ''.join([item[0] for item in path[:-1]])
    return path_string

# ---------------------------------------------------------------------------------------------
### PairedComposition(k, d, Text) ###
#      Input: k and d integers and string Text
#      Output:(k,d)-mer composition of Text (list of lists) -> the collection of all (k,d)- mers in Text
#      (including repeated (k,d)-mers)
# Given a string Text, a (k,d)-mer is a pair of k-mers in Text separated by distance d.
# For example, [AAT,TGG] is a (3,4)-mer in TAATGCCATGGGATGTT.

# Sample Input: 3, 1, 'TAATGCCATGGGATGTT'
# Sample Output:
# [[AAT,CCA], [ATG,CAT], [ATG,GAT], [CAT,GGA], [CCA,GGG], [GCC,TGG], ->
#                                          -> [GGA,GTT, [GGG,TGT], [TAA,GCC], [TGC,ATG], [TGG,ATG]]

def PairedComposition(k, d, Text):
    composition = []
    for i in range(len(Text)-2*k-d+1):
        composition.append([Text[i:i+k], Text[i+k+d:i+d+2*k]])
    sorted_composition = sorted(composition, key=lambda x: x[0])
    return sorted_composition

# ---------------------------------------------------------------------------------------------
### StringReconstructionFromPairs(k, d) ###
# Challenge: Solve the String Reconstruction from Read-Pairs Problem.
#      Input: Integers k and d followed by a collection of paired k-mers PairedReads.
#      Output: A string Text with (k, d)-mer composition equal to PairedReads.

# Sample Input: k = 4 ; d = 2 ;
# [[GAGA,TTGA], [TCGT,GATG], [CGTG,ATGT], [TGGT,TGAG], [GTGA,TGTT], ->
#                                        -> [GTGG,GTGA], [TGAG,GTTG], [GGTC,GAGA], [GTCG,AGAT]]
# Sample Output:
#     GTGGTCGTGAGATGTTGA

def StringReconstructionFromPairs(k, d, paired_kmers):

    ## Construct a dictionary of edges from the paired reads:
    # Graph: edge = paired kmers ((kmer1, kmer2)); node = prefix (start node) / suffix (end node)
    # the graph is represented by dictionery: key = prefix (start node of an edge); value =  suffix (ending node of an edge)

    paired_kmers_dict = {}
    for pair in paired_kmers:
        if (pair[0][:-1], pair[1][:-1]) in paired_kmers_dict:  # if paired prefix is in dict
            paired_kmers_dict[(pair[0][:-1], pair[1][:-1])].append((pair[0][1:], pair[1][1:]))  # add paired suffix
        else:  # if paired prefix is not in dict
            paired_kmers_dict[(pair[0][:-1], pair[1][:-1])] = [(pair[0][1:], pair[1][1:])]  # add paired prefix (= key) and paired sufffix (= value)

    ## Get an eulerian path from the paired edges:
    paired_path = EulerianPath(paired_kmers_dict)

    ## Recombine the paths, accounting for their overlaps.
    strings = [paired_path[0][i] + ''.join(map(lambda x: x[i][-1], paired_path[1:])) for i in range(2)]
    return strings[0][:k + d] + strings[1]

# ---------------------------------------------------------------------------------------------
### ContigGeneration(kmers) ###
# Contig Generation Problem: Generate the contigs from a collection of reads (with imperfect coverage).
#      Input: A collection of k-mers kmers.
#      Output: All contigs in DeBruijn(Patterns).

# Sample Input: ['ATG', 'ATG', 'TGT', 'TGG', 'CAT', 'GGA', 'GAT', 'AGA']
# Sample Output:['AGA', 'ATG', 'ATG', 'CAT', 'GAT', 'TGGA', 'TGT']

def flatten (l):
    return flatten(l[0]) + (flatten(l[1:]) if len(l) > 1 else []) if type(l) is list else [l]

def ContigGeneration (kmers):
    # Construct a dictionary of edges.
    edges = {}
    for kmer in kmers:
        if kmer[:-1] in edges:
            edges[kmer[:-1]].append(kmer[1:])
        else:
            edges[kmer[:-1]] = [kmer[1:]]

    # Determine the balanced and unbalanced edges.
    balanced, unbalanced = [], []
    out_values = reduce(lambda a, b: a + b, edges.values())
    for node in set(out_values + list(edges.keys())):
        out_value = out_values.count(node)
        if node in edges:
            in_value = len(edges[node])
        else:
            in_value = 0

        if in_value == out_value == 1:
            balanced.append(node)
        else:
            unbalanced.append(node)

    # Generate the contigs.
    get_contigs = lambda s, c: flatten(
        [c + e[-1] if e not in balanced else get_contigs(e, c + e[-1]) for e in edges[s]])
    contigs = sorted(flatten([get_contigs(start, start) for start in set(unbalanced) & set(edges.keys())]))
    return contigs


