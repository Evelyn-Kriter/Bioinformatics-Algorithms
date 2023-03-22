from collections import Counter


sample_graph = {0:[2],
                1:[3],
                2:[1],
                3:[0, 4],
                6:[3, 7],
                7:[8],
                8:[9],
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


### TEST
# sample_path = eulerian_path(sample_graph)
# print("->".join(map(str, sample_path)))

### OUTPUT
graph = read2B("rosalind_ba3g.txt")
path = eulerian_path(graph)
print("->".join(map(str, path)))
with open("output_ba3g.txt", "w") as file:
   print("->".join(map(str, path)), file=file)
   
