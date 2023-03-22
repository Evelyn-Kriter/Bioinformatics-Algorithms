import numpy as np

STRING = "xzyyzzyzyy"
SIGMA = ["x",
         "y",
         "z"]
STATES = ["A", "B"]
TRANSITION = np.array([[0.303, 0.697],[0.831, 0.169]])
EMISSION = np.array([[0.533, 0.065, 0.402],[0.342, 0.334, 0.324]]) 

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

def hmm_prob(string, sigma, states, trans, emit):
    score = np.full((len(states),len(string)), 0, dtype=np.longdouble)
    
    #compute first column
    for k in range(0, len(states)):
        score[k,0] = (1/len(states)) * emit[k, sigma.index(string[0])]
    
    for i in range(1, len(string)):
        for j in range(0, len(states)):
            for l in range(0, len(states)):
                score[j,i] += score[l][i-1] * trans[l][j] * emit[j, sigma.index(string[i])]
                                
    probability = 0
    for m in range(0, len(states)):
        probability += score[m][-1]
    
    return probability



### TEST
#sample_probability = hmm_prob(STRING, SIGMA, STATES, TRANSITION, EMISSION) 
#print(sample_probability)
    

### OUTPUT
string, sigma, states, trans, emit = read2B("rosalind_ba10d.txt")
probability = hmm_prob(string, sigma, states, trans, emit)
print(probability)