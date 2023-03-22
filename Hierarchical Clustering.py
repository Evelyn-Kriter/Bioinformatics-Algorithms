import numpy as np

sample_D = 7
sample_n = np.array([[.00, 0.74, 0.85, 0.54, 0.83, 0.92, 0.89],
                    [0.74, 0.00, 1.59, 1.35, 1.20, 1.48, 1.55],
                    [0.85, 1.59, 0.00, 0.63, 1.13, 0.69, 0.73],
                    [0.54, 1.35, 0.63, 0.00, 0.66, 0.43, 0.88],
                    [0.83, 1.20, 1.13, 0.66, 0.00, 0.72, 0.55],
                    [0.92, 1.48, 0.69, 0.43, 0.72, 0.00, 0.80],
                    [0.89, 1.55, 0.73, 0.88, 0.55, 0.80, 0.00]])

def read2C(filename):
    """Read 2C input data"""
    with open(filename, "r") as file:
        D = int(file.readline().strip())
        n = np.loadtxt(file)  # Default type is float
        return D, n

def heirarchicalclustering(D, n):
    #n single-element clusters labeled 1, ... , n
    clusters = []
    for i in range(0, D):
        clusters.append([i])
        
    matrix = n
    
    while len(matrix) > 1: #find the two closest clusters Ci and Cj 
        closest_distance = float("inf")
        for i in range(0, len(matrix)):
            for j in range(i+1, len(matrix)):
                distance = matrix[i,j] #average_distance(matrix[i], matrix[j])
                if distance < closest_distance:
                    row1 = i
                    row2 = j
                    closest_distance = distance
        #print(row1)
        #print(row2)
        #index = np.where(clusters == np.amin(clusters))
        #print(index)
        #row1 = index[0]
        #row2 = index[1]
        #print(row1)
        #print(row2)
        top = np.array([len(clusters[row2])*matrix[row2,:]+len(clusters[row1]) *matrix[row1, :]])
        #print("top")
        #print(top)
        bottom = np.array([len(clusters[row2]) + len(clusters[row1])])
        #print("bottom")
        #print(bottom)
        #matrix_rnew = rowwise_averages(matrix[row1], matrix[row2])
        matrix_rnew = top/bottom
        matrix_cnew = np.append(matrix_rnew[0], np.array([0]))
        #print(matrix_rnew)
        cluster_rnew = clusters[row1] + clusters[row2]
        
        matrix = np.append(matrix, matrix_rnew, axis = 0)
        matrix = np.concatenate((matrix, matrix_cnew.reshape(-1,1)), axis = 1)
        #matrix = np.
        #print(matrix)
        matrix = np.delete(matrix, row2, axis = 0)
        matrix = np.delete(matrix, row1, axis = 0)
        matrix = np.delete(matrix, row2, axis = 1)
        matrix = np.delete(matrix, row1, axis = 1)
        
        
        assert row2 > row1
        clusters.pop(row2)
        clusters.pop(row1)
#         for i in range(0, len(clusters)):
#             for j in range(0, len(clusters[i])):
#                 if j == row2:
#                     #print(row2)
#                     clusters[i].pop(j)
#                 if j == row1:
#                     #print(row1)
#                     clusters[i].pop(j)
            # clusters[i].pop(row2)
            # clusters[i].pop(row1)
        clusters.append(cluster_rnew)
        
        #print(clusters)
        print(*[i+1 for i in cluster_rnew])

    return clusters

# def average_distance(C1, C2):
#     davg = ((C1-C0)+(C2-C0))/2
#     return davg

# def rowwise_averages(Crow1, Crow2):
#     Crow_new = []
#     print(Crow1)
#     print(Crow2)
#     for i in range(0, len(Crow1)):
#         Crow_new.append(((i * Crow1[i]) + (i * Crow2[i]))/(i + i))
#         #print("Crow1")
#         print(Crow1[i])
#         print("Crow2")
#         print(Crow2[i])
#         print(((i * Crow1[i]) + (i * Crow2[i]))/(2*i))
#     return Crow_new
        

## TEST
D, n = read2C("rosalind_ba8e.txt")
sample_output = heirarchicalclustering(D, n)
#print(sample_output)

## OUTPUT
#D, n = read2B("rosalind_ba8e.txt")
#output = heirarchicalclustering(D, n)
#print(output)