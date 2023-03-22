
sample_graph = {0:[3],
                1:[0],
                2:[1,6],
                3:[2],
                4:[2],
                5:[4],
                6:[5,8],
                7:[9],
                8:[7],
                9:[6]}
 
def read2B(filename):
    """Read rosalind.info 2B input data"""
    with open(filename, "r") as file:
        graph = {}
        
        for line in file:
            key, value = line.strip().split("->")
            value = list(map(str.strip, value.split(",")))
            graph[key.strip()] = value
            
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

### TEST OUTPUT
#sample_output = eulerian_cycle(sample_graph)
#print(sample_output) 

### OUTPUT
graph = read2B("rosalind_ba3f.txt")
output = eulerian_cycle(graph)
with open("output_ba3f.txt", "w") as file:
   print("->".join(output), file=file)