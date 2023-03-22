import numpy as np
import sys

STRING = "zyxxxxyxzz"
SIGMA = ["x",
         "y",
         "z"]
STATES = ["A", "B"]
TRANSITION = np.array([[0.911, 0.089],[0.228, 0.772]])
EMISSION = np.array([[0.356, 0.191, 0.453],[0.040, 0.467, 0.493]]) 

def read2B(filename):
    """Read rosalind.info 2B input data"""
    with open(filename, "r") as file:
        string = file.readline().strip()
        file.readline()
        sigma = file.readline().strip().split()
        file.readline()
        states = file.readline().strip().split()
        trans = np.genfromtxt(file, dtype=np.longdouble, skip_header=2, max_rows=len(states), usecols=range(1, 1+len(states)))
        emit = np.genfromtxt(file, dtype=np.longdouble, skip_header=2, max_rows=len(states), usecols=range(1, 1+len(sigma)))
        
        return string, sigma, states, trans, emit

def soft_decode(string, sigma, states, trans, emit):
    
    forward = np.full((len(states), len(string)), 0, dtype = np.longdouble) #)[[len(states)]*len(string)]
    backward = np.full((len(states), len(string)), 0, dtype = np.longdouble) #[[len(states)]*len(string)]
    
    for i in range(0, len(states)):
        forward[i,0] = (1/len(states)) * emit[i,sigma.index(string[0])]
     
    #compute the rest of the columns
    for i in range(1, len(string)):
        for k in range(0, len(states)):
            forward[k,i] = np.sum(forward[:,i-1] * trans[:,k] * emit[k,sigma.index(string[i])])
    
    sink = np.sum(forward[:,-1])
    
    backward[:,-1]= 1.0
    for p in range(len(string)-2, -1, -1):
        for q in range(len(states)):
            #for r in range(0, len(string)):
            backward[q, p] = np.sum(backward[:, p+1]*trans[q, :]*emit[:,sigma.index(string[p+1])])

    probability = forward*backward/sink
    
    return probability


### TEST
# sample_probability = soft_decode(STRING, SIGMA, STATES, TRANSITION, EMISSION)
# print(*STATES)
# output= np.savetxt(sys.stdout, np.transpose(sample_probability),  fmt="%.4f")

    

### OUTPUT
string, sigma, states, trans, emit = read2B("rosalind_ba10j.txt")
probability = soft_decode(string, sigma, states, trans, emit)
print(*states)
output= np.savetxt(sys.stdout, np.transpose(probability),  fmt="%.4f")
