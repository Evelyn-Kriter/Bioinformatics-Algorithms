import numpy as np
import math
from scipy.stats import norm

def read2B(filename):
    """Read rosalind.info 2B input data"""
    with open(filename, "r") as file:
        samples = np.genfromtxt(filename, skip_header=1, usecols=0, dtype=str)
        exome_targets = np.genfromtxt(filename, skip_footer = len(samples), dtype=str)[1:]
        matrix = np.genfromtxt(filename, skip_header=1, usecols=range(1,1+len(exome_targets)), dtype=np.longdouble)

        
        return exome_targets, matrix, samples

def xhmm(exome_target, matrix, samples):
    
    p = 0.00000001 #CNV rate
    q = 0.166 #1/mean targets per CNV
    M = [0, 1, 2]
    trans =np.array([[1-q, q, 0],
             [p, 1-(2*p), p],
             [0, q, 1-q]])

    for i in range(0, len(samples)):
        score = np.full((3, len(matrix[i])), 0, dtype=np.longdouble)
        score[0,:] = norm.pdf(matrix[i], -3, 1)
        score[1,:] = norm.pdf(matrix[i], 0, 1)
        score[2,:] = norm.pdf(matrix[i], 3, 1)
        
        path = viterbi_path(trans, score)
        if set(path) == {1}:
            continue
        
        print("SAMPLE CNV START END Q_EXACT")
        start = -1
        end = -1
        for j in range(0, len(path)):
            if start == -1:
                if path[j] != 1:
                    start = j
            else:
                if path[j] == 1:
                    end = j - 1
                    break

        probability = soft_decode(trans, score, start, end, path[start:end+1])

        
        #Phred score
        phred = -10*math.log10(1-probability)
        print(samples[i], "DUP", exome_target[start], exome_target[end], phred)
    
    return probability

def viterbi_path(trans, emit):
    states, targets = emit.shape
    score = np.full((states,targets), 0, dtype=np.longdouble)
    backtrace = np.full((states,targets), 0, dtype=np.int_)
    
    #compute first column (not uniform)
    for i in range(0, states):
        score[i,0] = trans[1,i] * emit[i,0]
     
    #compute the rest of the columns
    for i in range(1, targets):
        for k in range(0, states):
            score[k,i] = np.max(score[:,i-1] * trans[:,k] * emit[k,i])
            
            backtrace[k,i] = np.argmax(score[:,i-1] * trans[:,k] * emit[k,i])
    
    k = np.argmax(score[:,-1])
    l = targets-1
    output = []
    output.append(k) # = k + output
    while l > 0:
        k = backtrace[k,l]
        l += -1
        output.append(k) #= states[k] + output
    
    return output[::-1]

def soft_decode(trans, emit, start, end, path):
    states, targets = emit.shape
    
    forward = np.full((states, targets), 0, dtype = np.longdouble) 
    backward = np.full((states, targets), 0, dtype = np.longdouble)
    
    for i in range(0, states):
        forward[i,0] = trans[1,i] * emit[i,0]
     
    #compute the rest of the columns
    for i in range(1, targets):
        for k in range(0, states):
            forward[k,i] = np.sum(forward[:,i-1] * trans[:,k] * emit[k,i])
    
    sink = np.sum(forward[:,-1])
    backward[:,-1]= 1.0
    
    for i in range(targets-2, -1, -1):
        for j in range(states):
            backward[j, i] = np.sum(backward[:, i+1]*trans[j, :]*emit[:,i+1])
    
    middle = 1
    for i in range(1, end-start+1):
        middle = middle * trans[path[i-1], path[i]]*emit[path[i], start+i]
    
    probability = (forward[path[0],start]*middle*backward[path[-1], end])/sink
    
    return probability

### OUTPUT
exome_targets, matrix, samples = read2B("XHMM.in.txt")
xhmm(exome_targets, matrix, samples)
#with open("output_XHMM.txt", "w") as file:
   #print(CNVs, file=file)