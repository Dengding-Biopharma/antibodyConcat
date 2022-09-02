import copy
from collections import Counter

from tqdm import trange

from debruijn import get_graph_from_reads



def construct_naive_debruijn_graph(reads,k,pruning):
    vertices, edges = get_graph_from_reads(reads, k)

    if pruning:
        edges = pruningEdges(edges,2)
    else:
        for edge in edges:
            edges[edge] = list(Counter(edges[edge]).keys())

    # for edge in edges:
    #     print(edge,edges[edge])
    # quit()
    branch_kmer = []
    count = 0
    for edge in list(edges):
        if len(edges[edge]) > 1:
            count += 1
            branch_kmer.append(edge)
    print('branch number: ', count)
    return (vertices, edges)

def output_contigs(g):
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
    for i in trange(len(starts)):
        start = starts[i]
        current = start
        vec = []
        output = []
        contig_copy = []
        DFS(current, E, vec, output, contig_copy)
        contig.extend(contig_copy)

    return contig

def DFS(current, E, vec, output, contig_copy):
    if current in vec:
        return
    vec.append(current)
    print(111,vec)
    if len(E[current]) == 0:
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
        DFS(E[current][i], E, vec, output, contig_copy)
    vec.pop()

def pruningEdges(edges, threshold):
    for edge in edges:
        previous = edges[edge]
        counter = Counter(previous)
        if len(counter) == 0:
            continue
        if len(counter) == 1:
            edges[edge] = list(counter)
        else:
            maxCountKmer = [counter.most_common(1)[0][0]]
            maxCount = counter.most_common(1)[0][1]
            for i in range(1, len(counter)):
                nextCount = counter.most_common()[i][1]
                if nextCount >= maxCount / threshold:
                    maxCountKmer += [counter.most_common()[i][0]]
            edges[edge] = maxCountKmer
    return edges