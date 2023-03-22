from collections import Counter
 
sample_genome = (["CTTA",
                  "ACCA",
                  "TACC",
                  "GGCT",
                  "GCTT",
                  "TTAC"])
sample_k = 4

 
 
def read2B(filename):
    """Read rosalind.info 2B input data"""
    with open(filename, "r") as file:
        k= file.readline().split()
        
        kmers = [line.strip() for line in file.readlines()]
        return k, kmers
    
def construct_debruijn(kmers):
    pn_kmers = kmers
    pn_kmers.sort()
    
    graph = {}
    for i in range(0, len(pn_kmers)):
        
        #turn every prefix into a key in the dictionary and add suffix
        prefix = pn_kmers[i][:-1]
        suffix = pn_kmers[i][1:]
        if prefix in graph:
            graph[prefix].append(suffix)
        else:
            graph[prefix] = [suffix]
            
    return graph

def eulerian_cycle(graph):
    #form a cycle Cycle by randomly walking in Graph and without visiting the same edge twice
    cycle = [next(iter(graph.keys()))]
    
    while len(graph) > 0:
        node = cycle[-1]
        
        edges = graph.get(node, [])
        
        if len(edges) > 0:
            cycle.append(edges[0])
            edges.pop(0)
            if len(edges) == 0:
                graph.pop(node)
        else:
            for i, newstart in enumerate(cycle):
                if newstart in graph:
                    cycle = cycle[i:] + cycle[1:i+1]
                    break
    
    return cycle

def eulerian_path(graph):
    balanced = Counter()
    
    #calculate in-edges
    for node, edges in graph.items():
        balanced[node] += len(edges)
        for edge in edges:
            balanced[edge] += -1
    
    for node in balanced:
        if balanced[node] == 1:
            to_node = node
        elif balanced[node] == -1:
            from_node = node
    
    #add edge
    if from_node in graph:
        graph[from_node].append(to_node)
    else:
        graph[from_node] = [to_node]  
    
    cycle = eulerian_cycle(graph)
    
    #break cycle at added edge, w and v become newEnd and newStart
    for i in range(0, len(cycle)-1):
        if cycle[i] == from_node and cycle[i+1] == to_node:
            #form a path by traversing cycle from newStart to newEnd
                path = cycle[i+1:] + cycle[1:i+1]
    return path

def reconstruct_kmer_comp(kmers, k):
    graph = construct_debruijn(kmers)
    path = eulerian_path(graph)
    string = path[0]
    
    for i in range(1, len(path)):
        string += path[i][-1]
    
    return string

### TEST
sample_output = reconstruct_kmer_comp(sample_genome, sample_k)
print(sample_output)

### OUTPUT
# k, kmers = read2B("rosalind_ba3h.txt")
# output = reconstruct_kmer_comp(kmers, k)
# print(output)