
sample_genome = (["GAGG",
                  "CAGG",
                  "GGGG",
                  "GGGA",
                  "CAGG",
                  "AGGG",
                  "GGAG"])


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

### TEST OUTPUT
# sample_graph = construct_debruijn(sample_genome)
# for key,value in sample_graph.items():
#     print(key, "->", ",".join(value))

### OUTPUT
list_of_kmers = read2B("rosalind_ba3e.txt")
graph = construct_debruijn(list_of_kmers)
with open("output_ba3e.txt", "w") as file:
    for key,value in graph.items():
       print(key, "->", ",".join(value), file = file)