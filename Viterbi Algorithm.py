import numpy as np

sample = "HHTTHT"
SIGMA = ["H",
         "T"]
STATES = ["F", "B"]
TRANSITION = np.array([[0.641, 0.359],[0.729, 0.271]])
EMISSION = np.array([[0.117, 0.691, 0.192],[0.097, 0.42, 0.483]]) 

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

def viterbi_path(string, sigma, states, trans, emit):
    score = np.full((len(states),len(string)), 0, dtype=np.longdouble)
    backtrace = np.full((len(states),len(string)), 0, dtype=np.int_)
    
    #convert string to numbers
    
    #compute first column
    for i in range(0, len(states)):
        score[i,0] = (1/len(states)) * emit[i,sigma.index(string[0])]
     
    #compute the rest of the columns
    for i in range(1, len(string)):
        for k in range(0, len(states)):
            score[k,i] = np.max(score[:,i-1] * trans[:,k] * emit[k,sigma.index(string[i])])
            
            backtrace[k,i] = np.argmax(score[:,i-1] * trans[:,k] * emit[k,sigma.index(string[i])])
    
    k = np.argmax(score[:,-1])
    l = len(string)-1
    output = ""
    output = states[k] + output
    while l > 0:
        k = backtrace[k,l]
        l += -1
        output = states[k] + output
    
    return output
            
###TEST
# sample_output = viterbi_path(sample, SIGMA, STATES, TRANSITION, EMISSION)
# print(sample_output)

### OUTPUT
string, sigma, states, trans, emit = read2B("rosalind_ba10c.txt")
output = viterbi_path(string, sigma, states, trans, emit)
print(output)
with open("output_ba10c.txt", "w") as file:
   print(output, file=file)