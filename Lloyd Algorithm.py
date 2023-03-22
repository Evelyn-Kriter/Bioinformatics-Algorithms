import numpy as np
import sys

sample_k = 2
sample_m = 2
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
        data = np.genfromtxt(file, dtype=np.longdouble)  # Default type is float
        return int(k), int(m), data

# def centers_of_gravity(clusters):
#     centers = np.array()
#     for cluster in clusters:
#         center = np.average(cluster)
#         centers.append(center)
#     
#     return centers

def lloyd(k, m, data):
    centers = np.full((k,m), fill_value = 0.0)
    
    #initial centers
    for i in range(0, k):
        centers[i] = data[i]
    
    changing = True
    
    while changing == True:
        clusters = np.full((k,m), fill_value = 0.0)
        count = np.full((k,), fill_value = 0)
        
        #centers to clusters
        for point in data:
            smallest_distance = 9999999
            closest_center = -1
            for i, center in enumerate(centers):
                distance = np.linalg.norm(point-center) #euclidean distance
                if distance < smallest_distance:
                    smallest_distance = distance
                    closest_center = i
            clusters[closest_center] += point
            count[closest_center] += 1
        
        old_centers = centers
        #print(count)
        
        centers = clusters / count[:,np.newaxis]

        
        #clusters to centers
        #centers = centers_of_gravity(clusters)
        
        #track change
#         if len(centers) == len(old_centers):
#             counter = 0
#             for i in range(0, len(centers)):
#                 for j in range(len(centers[i])):
#                 percent_change = (centers[i,j] - old_centers[i,j])/old_centers[i,j]
#                 if percent_change < 0.001:
#                     counter += 1
#             if counter == len(centers) * m:
#                 changing = False
        if np.allclose(centers, old_centers) == True:
            changing = False
        
    return centers
            
## TEST
# sample_output = lloyd(sample_k, sample_m, sample_data)
# sample_processed_output= np.savetxt(sys.stdout, np.transpose(sample_output),  fmt="%.3f")
# print(sample_processed_output)

## OUTPUT
k, m, data = read2C("rosalind_ba8c.txt")
print(k)
print(m)
print(data)
output = lloyd(k, m, data)
processed_output= np.savetxt(sys.stdout, output,  fmt="%.3f")
print(processed_output)