import copy
from collections import Counter

from debruijn import get_graph_from_reads



def construct_naive_debruijn_graph(reads,k):
    vertices, edges = get_graph_from_reads(reads, k)
    for edge in edges:
        edges[edge] = list(Counter(edges[edge]).keys())
    return (vertices, edges)

def output_contigs(g, branch_kmer, already_pull_out):
    """ Perform searching for Eulerian path in the graph to output genome assembly"""
    V = g[0]
    E = g[1]
    # Pick starting node (the vertex with zero in degree)
    starts = []
    for k in list(V.keys()):
        if V[k].indegree == 0:
            starts.append(k)

    print('Number of kmers have no income edges: ', len(starts))
    contig = []
    for i in range(len(starts)):
        start = starts[i]
        current = start
        vec = []
        output = []
        contig_copy = []
        DFS(current, E, vec, output, contig_copy, branch_kmer, already_pull_out)
        contig.extend(contig_copy)

    return contig

def DFS(current, E, vec, output, contig_copy, branch_kmer, already_pull_out):
    if current in vec:
        return
    vec.append(current)
    if current in already_pull_out:
        if len(vec) == 1:
            vec.pop()
            return
        vec.pop()
        if vec not in output:
            result = vec[0]
            for i in range(1, len(vec)):
                result += vec[i][-1]
            output.append(copy.deepcopy(vec))
            contig_copy.append(result)
        return
    if current in branch_kmer or len(E[current]) == 0:
        if vec not in output:
            result = vec[0]
            for i in range(1, len(vec)):
                result += vec[i][-1]
            # print(result,len(result))
            output.append(copy.deepcopy(vec))
            contig_copy.append(result)
        vec.pop()
        return
    for i in range(len(E[current])):
        DFS(E[current][i], E, vec, output, contig_copy, branch_kmer, already_pull_out)
    vec.pop()

