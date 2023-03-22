import numpy as np
import math
import sys

sample_k = 2
sample_m = 2
sample_stiffness = 2.7
sample_data = np.array([[1.3, 1.1],
                       [1.3, 0.2],
                       [0.6, 2.8],
                       [3.0, 3.2],
                       [1.2, 0.7],
                       [1.4, 1.6],
                       [1.2, 1.0],
                       [1.2, 1.1],
                       [0.6, 1.5],
                       [1.8, 2.6],
                       [1.2, 1.3],
                       [1.2, 1.0],
                       [0.0, 1.9]])

def read2C(filename):
    """Read 2C input data"""
    with open(filename, "r") as file:
        k, m = file.readline().strip().split()
        stiffness = float(file.readline().strip())
        data = np.genfromtxt(file, dtype=np.longdouble)  # Default type is float
        return int(k), int(m), stiffness, data

def centers_of_gravity(clusters):
    centers = []
    for cluster in clusters:
        center = np.average(cluster)
        centers.append(center)
    
    return centers

def soft_clustering(k, m, stiffness, data):
    centers = data[:k]
    clusters = [[]]
    
    for x in range(0, 100):
        hidden_matrix = np.zeros((k, len(data)), dtype = np.longdouble) 
        
        #centers to clusters
        for x in range(len(data)): #data:
            for y in range(len(centers)): #centers:
                euclidD = np.linalg.norm(data[x]-centers[y]) #euclidean distance
                hidden_matrix[y,x] = np.exp(euclidD*(stiffness*-1))
        hidden_matrix /= np.sum(hidden_matrix, axis=0)[np.newaxis,:]
        new_centers = np.matmul(hidden_matrix, data)
        centers = new_centers/np.sum(hidden_matrix, axis=1)[:,np.newaxis]
        
    return centers
            
            
## TEST
#sample_output = soft_clustering(sample_k, sample_m, sample_stiffness, sample_data)
#np.savetxt(sys.stdout, sample_output, fmt="%0.3f")

## OUTPUT
k, m, stiffness, data = read2C("rosalind_ba8d.txt")
output = soft_clustering(k, m, stiffness, data)
np.savetxt(sys.stdout, output, fmt="%0.3f")