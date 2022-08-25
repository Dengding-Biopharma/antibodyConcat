import copy
from collections import Counter
import sys
from pprint import pprint
from tqdm import trange


class Node:
    """ Class Node to represent a vertex in the de bruijn graph """

    def __init__(self, lab):
        self.label = lab
        self.indegree = 0
        self.outdegree = 0


class Edge:
    def __init__(self, lab):
        self.label = lab


def read_reads(fname):
    """ Read short reads in FASTA format. It is assumed that one line in the input file correspond to one read. """
    f = open(fname, 'r')
    lines = f.readlines()
    f.close()
    reads = []

    for line in lines:
        if line[0] != '>':
            reads = reads + [line.rstrip()]
    return reads


def get_kmers(sequences, k):
    """
    Returns dictionary with keys representing all possible kmers in a sequence
    and values counting their occurrence in the sequence.
    """
    short_sequences = []
    for sequence in sequences:
        if len(sequence) < k:
            short_sequences.append(sequence)
    for sequence in short_sequences:
        sequences.remove(sequence)
    for short_sequence in short_sequences:
        concat = False
        for i in range(len(sequences)):
            sequence = sequences[i]
            if sequence[len(sequence) - 3:] == short_sequence[:3]:
                sequences[i] = sequence + short_sequence[3:]
                concat = True
            if short_sequence[len(short_sequence) - 3:] == sequence[:3]:
                sequences[i] = short_sequence + sequence[3:]
                concat = True
        if concat:
            short_sequences.remove(short_sequence)
    # print('111111111111inputs: ',sequences)
    # dict to store kmers
    kmers = {}
    for sequence in sequences:
        # count how many times each occurred in this sequence (treated as cyclic)
        for i in range(len(sequence)):
            kmer = sequence[i:i + k]

            # for cyclic sequence get kmers that wrap from end to beginning
            if len(kmer) != k:
                continue
            # count occurrence of this kmer in sequence
            if kmer in kmers:
                kmers[kmer] += 1
            else:
                kmers[kmer] = 1

    return list(kmers.keys())


def get_graph_from_kmers(kmers, k):
    edges = dict()
    vertices = dict()
    for index in range(len(kmers)):
        kmer = kmers[index]
        vertices[kmer] = Node(kmer)
        edges[kmer] = []
        for edge in edges:
            if kmer[1:] == edge[:k - 1]:  # kmer -> edge
                edges[kmer] += [edge]
                vertices[kmer].outdegree += 1
                vertices[edge].indegree += 1
            if kmer[:k - 1] == edge[1:]:  # edge -> kmer
                edges[edge] += [kmer]
                vertices[edge].outdegree += 1
                vertices[kmer].indegree += 1

    return vertices, edges


def get_graph_from_reads(reads, k):
    # a = k-1
    # short_sequences = []
    # for sequence in reads:
    #     if len(sequence) <= k:
    #         short_sequences.append(sequence)
    # for sequence in short_sequences:
    #     reads.remove(sequence)
    #
    # for short_sequence in short_sequences:
    #     concat = False
    #     for i in range(len(reads)):
    #         sequence = reads[i]
    #         if sequence[len(sequence)-a:] == short_sequence[:a]:
    #             reads[i] = sequence + short_sequence[a:]
    #             concat = True
    #         elif short_sequence[len(short_sequence)-a:] == sequence[:a]:
    #             reads[i] = short_sequence + sequence[a:]
    #             concat = True
    #     if concat:
    #         short_sequences.remove(short_sequence)

    # print('inputs before: ',reads)
    edges = dict()
    vertices = dict()
    for index in range(len(reads)):
        read = reads[index]
        i = 0
        while i + k < len(read):
            v1 = read[i:i + k]
            v2 = read[i + 1:i + k + 1]
            if v1 in edges.keys():
                if v2 not in edges[v1]:
                    vertices[v1].outdegree += 1
                edges[v1] += [v2]
            else:
                vertices[v1] = Node(v1)
                vertices[v1].outdegree += 1
                edges[v1] = [v2]
            if v2 in edges.keys():
                if v2 not in edges[v1]:
                    vertices[v2].indegree += 1
            else:
                vertices[v2] = Node(v2)
                vertices[v2].indegree += 1
                edges[v2] = []
            i += 1
    # print('input after',edges)
    # print('wasted reads: ',wasted)
    return vertices, edges


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


def pruningErrorContigFromBranch(current, edges, vertices, vec, output, depth, already_pull_out):
    if depth == 0:
        return
    vec.append(current)
    if vertices[current].outdegree == 0:
        if vec not in output:
            output.append(copy.deepcopy(vec))
        vec.pop()
        return
    elif vertices[current].indegree > 1:
        vec.pop()
        if vec not in output:
            output.append(copy.deepcopy(vec))
        return
    for i in range(len(edges[current])):
        if edges[current][i] not in already_pull_out:
            pruningErrorContigFromBranch(edges[current][i], edges, vertices, vec, output, depth - 1, already_pull_out)
    vec.pop()


def pruningErrorContigFromHead(current, edges, vertices, vec, output, depth, already_pull_out, edge_count_table):
    if depth == 0 or current in already_pull_out:
        return
    vec.append(current)
    if vertices[current].indegree > 1:
        name = vec[-2] + vec[-1][-1]
        vec.pop()
        if vec not in output:
            if edge_count_table[name] < 2:
                output.append(copy.deepcopy(vec))
        return
    for i in range(len(edges[current])):
        pruningErrorContigFromHead(edges[current][i], edges, vertices, vec, output, depth - 1, already_pull_out,
                                   edge_count_table)
    vec.pop()


def construct_graph(reads, k, threshold=3, final=False):
    """ Construct de bruijn graph from sets of short reads with k length word"""
    pull_out_read = []
    vertices, edges = get_graph_from_reads(reads, k)
    # vertices,edges = get_graph_from_kmers(get_kmers(reads,k),k)
    # print(1,edges)

    edge_count_table = dict()
    for edge in edges:
        counter = Counter(edges[edge])
        for i in list(counter.keys()):
            edge_name = edge + i[-1]
            count = counter[i]
            if edge_name not in edge_count_table.keys():
                edge_count_table[edge_name] = count
            else:
                edge_count_table[edge_name] = count + edge_count_table[edge_name]

    print('number of {}mer: '.format(k), len(vertices))

    edges = pruningEdges(edges, threshold)
    # for edge in edges:
    #     print(edge,edges[edge])

    branch_kmer = []
    count = 0
    for edge in list(edges):
        if len(edges[edge]) > 1:
            count += 1
            branch_kmer.append(edge)
    print('branch number: ', count)

    already_pull_out = []
    # print(branch_kmer)

    for kmer in branch_kmer:
        vec = []
        output = []
        if kmer in already_pull_out:
            continue
        pruningErrorContigFromBranch(kmer, edges, vertices, vec, output, 5, already_pull_out)
        # print(kmer,edges[kmer])
        # print(output)
        for o in output:
            for item in o:
                if item not in already_pull_out and item not in branch_kmer:
                    # print(item)
                    already_pull_out.append(item)
                    edges.pop(item)
        # print('----------------')
    # print(edges)


    starts = []
    for k in list(vertices.keys()):
        if vertices[k].indegree == 0:
            starts.append(k)

    for start in starts:
        vec = []
        output = []
        pruningErrorContigFromHead(start, edges, vertices, vec, output, 5, already_pull_out, edge_count_table)
        for o in output:
            for item in o:
                if item not in already_pull_out and item not in branch_kmer:
                    already_pull_out.append(item)
                    edges.pop(item)

    for read in reads:
        for kmer in branch_kmer:
            if kmer in read:
                pull_out_read.append(read)
                break

    # pprint(pull_out_read)
    if final:
        pull_out_read = []
        branch_kmer = []

    return (vertices, edges), pull_out_read, branch_kmer, already_pull_out, edge_count_table


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


def printPath(vec):
    # Print elements in the vector
    for ele in vec:
        print(ele, end=" ")
    print()


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


def print_graph(g):
    """ Print the information in the graph to be (somewhat) presentable """
    V = g[0]
    E = g[1]
    for k in V.keys():
        print("name: ", V[k].label, ". indegree: ", V[k].indegree, ". outdegree: ", V[k].outdegree)
        print("Edges: ")
        for e in E[k]:
            print(e)
        print()
