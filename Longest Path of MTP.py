import itertools
import numpy as np

def read2B(filename):
    """Read rosalind.info 2B input data"""
    with open(filename, "r") as file:
        # Recall that the file object maintains an internal pointer so you can use
        # multiple readline calls to read individual lines.
        n, m = file.readline().split()
        
        num_n = int(n)
        num_m = int(m)
        
        # Similarly readlines will read the remaining lines in the file.
        down = np.genfromtxt(file, max_rows= num_n)
        right = np.genfromtxt(file, skip_header = 1)
        return num_n, num_m, down, right


def longest_path_length(n, m, down, right):
    dag_array = np.full((n+1, m+1), 0) #keep track of array of values
    
    #initialize values for first row in dag_array
    for i in range(1, m+1):
        dag_array[0, i] = right[0, i-1] + dag_array[0, i-1]
    
    #initialize values for first column in dag_array
    for j in range(1, n+1):
        dag_array[j, 0] = down[j-1, 0] + dag_array[j-1, 0]
    
    #fill in values for rest of array
    for y in range(1, n+1): #iterate through row
        for x in range(1, m+1): #iterate through column
            dag_array[y, x] = max(dag_array[y-1,x]+down[y-1,x], dag_array[y,x-1]+right[y,x-1])
    
    return dag_array[n,m]
        


n, m, down, right = read2B("rosalind_ba5b.txt")
length = longest_path_length(n, m, down, right)
print(length)
