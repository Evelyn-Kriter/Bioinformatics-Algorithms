from collections import Counter

sample_kmers = ["ATG",
                "ATG",
                "TGT",
                "TGG",
                "CAT",
                "GGA",
                "GAT",
                "AGA"]
                #"CTT",
               # "TTC",
                #"TCT"]

def read2B(filename):
    """Read rosalind.info 2B input data"""
    with open(filename, "r") as file:
        kmers = [line.strip() for line in file.readlines()]
        return kmers
    
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

def generate_contigs(graph):
    num_incoming = Counter()
    num_outgoing = Counter()
    
    #calculate in-edges
    for node, edges in graph.items():
        num_outgoing[node] += len(edges)
        contigs = []
        for edge in edges:
            num_incoming[edge] += 1
            
        for v in graph:#range(0, len(graph[node])): #length of KEYS in graph
            if num_incoming[v] != 1 or num_outgoing[v] != 1: #v is not a 1-in-1-out node
                if num_outgoing[v] > 0:
                    for w in graph[v]: #range(0, len(num_outgoing)): #for each outgoing edge (v, w) from v
                        #nonbranchingpath = the path consisting of the single edge (v, w)
                        paths = [v, w]
                        while num_incoming[w] == 1 and num_outgoing[w] == 1: #w is a 1-in-1-out node
                            #extend NonBranchingPath by the outgoing edge (w, u) from w
                            w = graph[w][0]
                            paths.append(w)  
                        contigs.append(paths)
                        
        for v in graph:#range(0, len(graph[node])): #length of KEYS in graph
            for w in graph[v]: #range(0, len(num_outgoing)): #for each outgoing edge (v, w) from v
                if num_incoming[v] == 1 or num_outgoing[v] == 1:   
                    cycle = [v]
                while (num_incoming[v] == 1 and num_outgoing[v] == 1 and cycle[0] != cycle[-1]): #w is a 1-in-1-out node
                    #extend NonBranchingPath by the outgoing edge (w, u) from w
                    num_outgoing[v] += -1
                    v = graph[v].pop(0)
                    cycle.append(v)
                    #print("I got here!")
                    #print(cycle)
                    if cycle[0] == cycle[-1]:
                        contigs.append(cycle)
                        break

                
    return contigs

def clean_answer(contigs):
    #trimmed_contigs = [len(contigs)]
    cleaned_contigs = []
    
    for i in range(0, len(contigs)):
        contig_list = []
        contig_string = contigs[i][0]
        for j in range(1, len(contigs[i])):
            #contig_string = ""
            #stripped_contigs = "".join(contigs[i])
            contig_string += contigs[i][j][-1]
            
        cleaned_contigs.append(contig_string)
        
    return cleaned_contigs
    

### TEST
graph = (construct_debruijn(sample_kmers))

### OUTPUT
#kmers = read2B("rosalind_ba3k.txt")
#graph = (construct_debruijn(kmers))
contigs = generate_contigs(graph)

output = clean_answer(contigs)
print(" ".join(output))


