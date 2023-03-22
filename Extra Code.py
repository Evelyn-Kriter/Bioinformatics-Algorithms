LinearSpaceAlignment(v, w, top, bottom, left, right)
  if left = right
    output path formed by bottom − top vertical edges
  if top = bottom
    output path formed by right − left horizontal edges
  middle ← ⌊ (left + right)/2⌋
  midEdge ← MiddleEdge(v, w, top, bottom, left, right)
  midNode ← vertical coordinate of the initial node of midEdge
  LinearSpaceAlignment(v, w, top, midNode, left, middle)
  output midEdge
  if midEdge = "→" or midEdge = "↘"
    middle ← middle + 1
  if midEdge = "↓" or midEdge ="↘"
    midNode ← midNode + 1
  LinearSpaceAlignment(v, w, midNode, bottom, middle, right)
  
def SouthOrEast(n, m):
    if n =0 and m = 0:
        return 0
    x = -inf
    y = -inf
    if n >0:
        x = SouthOrEast(n-1, m) + weight of "down" into (n,m)
    if m >0:
        y = SouthOrEast(n, m-1) + weight of the "east" edge into (n,m)
    return max(x,y)

def OutputLCS (backtrack, v, i, j):
    if i =0 or j=0:
        return v
    if backtrack(i, j) = "down":
        OutputLCS(backtrack, v, i-1, j)
    else if backtrack(i, j) = "east":
        OutputLCS(backtrac, v, i, j-1)
    else:
        OutputLCS(backtrack, v, i-1, j-1)
        return v
    
def LongestPath(Graph, source, sink):
    for each node a in Graph:
        s[a] = -inf
        source = 0
    #topologically order Graph
        for each node a: #(from source to sink in the topological order)
            s[a] = max{s[b]+weight of edge from b to a}
        return s[sink]
    
def construct_debruijn(kmers):
    pn_kmers = kmers
    pn_kmers.sort()
    k = len(pn_kmers[0]) # may want to use later to not hardcode this
    
    graph = []*len(pn_kmers)
    for i in range(0, len(pn_kmers)):
    
        #find possible next nodes in graph
        next_nodes = []
        kmer = pn_kmers[i][:-1] #hardcoded, needs changing
        overlap = pn_kmers[i][1:-1] #hardcoded, needs changing
        print(pn_kmers)
        print(overlap)
        print("~")

        for j in range(0, len(pn_kmers)):
            print(pn_kmers[j][:-2])
            if (pn_kmers[j][:-2] == overlap):
                print('found one!')
                next_nodes.append(pn_kmers[j][:-1])
            
        pn_nodes = list(set(next_nodes)) #remove duplicates
        
        if len(pn_nodes) < 1:
            graph.append(kmer + " -> ")
        elif len(pn_nodes) == 1:    
            graph.append(kmer + " -> " + pn_nodes[0])
        else:
            graph.append(kmer + " -> " + pn_nodes[0])
            for m in range(1, len(pn_nodes)):
                graph[i] += ", " + pn_nodes[m]
        
    return graph


    