import itertools
import numpy as np
from Bio.Align import substitution_matrices

sample_string1 = "GTAGGCTTAAGGTTA"
sample_string2 = "TAGATA"
#string1 = "TACAGTA"
#string2 = "GGATACG"
penalty = 1
dna_mat = substitution_matrices.Array(alphabet="ACGT", data=np.array([[1, -1, -1, -1],[-1, 1, -1, -1], [-1, -1, 1, -1],[-1, -1, -1, 1]]))


def highest_scoring_fitting(string1, string2, penalty):
    alignment1 = ""
    alignment2 = ""
    dag_array = np.full((len(string2)+1, len(string1)+1), 0) #keep track of array of values
    bt_array = np.full((len(string2)+1, len(string1)+1), 0)
    
    #initialize first value
    dag_array[0,0] = 0
    
    #initialize values for first row in dag_array
    for i in range(1, len(string1)+1):
        dag_array[0, i] = max(0, dag_array[0,i-1]-penalty)
        bt_array[0, i] = 2
        
    #initialize values for first column in dag_array
    for j in range(1, len(string2)+1):
        dag_array[j, 0] = dag_array[j-1, 0]-penalty
        bt_array[j, 0] = 1
    
    #populate backtracking array
    for i in range(1, len(string1)+1):
        for j in range(1, len(string2)+1):
            option1 = dag_array[j-1,i]-penalty
            option2 = dag_array[j,i-1]-penalty
            option3 = dag_array[j-1,i-1]+dna_mat[(string2[j-1], string1[i-1])]
            best_option = max(option1, option2, option3)
            dag_array[j,i] = best_option
            
            if best_option == option1:
                bt_array[j, i] = 1
            
            if best_option == option2:
                bt_array[j, i] = 2
                
            if best_option == option3:
                bt_array[j, i] = 3
    
    #print out the backtracking
    i = np.argmax(dag_array[-1,:])
    j = len(string2)
    orig_i = i
    orig_j = j
    while j > 0:  
        if bt_array[j,i] == 1:
            alignment1 = "-" + alignment1
            alignment2 = string2[j-1] + alignment2
            #score += best_option
            j -= 1
        if bt_array[j,i] == 2:
            alignment1 = string1[i-1] + alignment1
            alignment2 = "-" + alignment2
            #score += best_option
            i -= 1
        if bt_array[j,i] == 3:
            alignment1 = string1[i-1] + alignment1
            alignment2 = string2[j-1] + alignment2
            #score += best_option
            j -= 1
            i -= 1

    return (dag_array[orig_j, orig_i]+(dag_array[j,i]*-1)), alignment1, alignment2
        

score, alignment1, alignment2 = highest_scoring_fitting(sample_string1, sample_string2, penalty)
print(score)
print(alignment1)
print(alignment2)