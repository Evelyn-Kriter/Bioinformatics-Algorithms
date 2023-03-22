import itertools
import numpy as np

sample_string1 = "SMNDWIMLEMATSTYQWPKFHGQEHPRKQMECHLLTLVPEKYWCRRRARVMFRTVINQTWEFHLHKTNKVCPHQWF"
sample_string2 = "SMNDDIMLETATSTQGDADWHEHPRKQMECCLLTLVPEKYWCLQGPGCAKVMYNHLHPWYCRGIICPHHAI"
#sample_string1 = "PRTEINS"
#sample_string2 = "PRTWPSEIN"

sample_opening = 11
sample_extension = 1

from Bio.Align import substitution_matrices
sub_mat = substitution_matrices.load("BLOSUM62")


def affine_gap(string1, string2, opening, extension):
    lower = np.full((len(string2)+1, len(string1)+1), 0)
    middle = np.full((len(string2)+1, len(string1)+1), 0)
    upper = np.full((len(string2)+1, len(string1)+1), 0)
    #dag_array = np.full((len(string2)+1, len(string1)+1), 0) #keep track of array of values
    bt_lower = np.full((len(string2)+1, len(string1)+1), 0)
    bt_middle = np.full((len(string2)+1, len(string1)+1), 0)
    bt_upper = np.full((len(string2)+1, len(string1)+1), 0)
    
    #initialize first value
    lower[0,0] = 0
    middle[0,0] = 0
    upper[0,0] = 0
    
    #initialize values for first row in array
    for i in range(1, len(string1)+1):
        #dag_array[0, i] = dag_array[0,i-1]-opening
        lower[0, i] = opening - (extension*i)
        middle[0, i] = 0
        bt_middle[0, i] = 2
    
    #initialize values for first column in array
    for j in range(1, len(string2)+1):
        upper[j, 0] = opening - (extension*j)
        middle[j, 0] = 0
        bt_middle[j, 0] = 1

    
    for i in range(1, len(string1)+1):
        for j in range(1, len(string2)+1):
            lower[j, i] = max(lower[j, i-1] - extension, middle[j, i-1] - opening)
            upper[j, i] = max(upper[j-1, i] - extension, middle[j-1, i] - opening)
            middle[j, i] = max(lower[j, i], middle[j-1, i-1] + sub_mat[(string2[j-1], string1[i-1])], upper[j, i])
            
            #deletion
            if lower[j, i] == lower[j, i-1] - extension:
                bt_lower[j, i] = 1
                #i-1, no matrix switching, 
            elif lower[j, i] == middle[j, i-1] - opening:
                bt_lower[j, i] = 2
                #i-1, bt = middle
                
            if upper[j, i] == upper[j-1, i] - extension:
                bt_upper[j, i] = 3
                #j-1, no matrix switching
            elif upper[j, i] == middle[j-1, i] - opening:
                bt_upper[j, i] = 4
                #j-1 , switch to middle
            
            if middle[j, i] == lower[j, i]:
                bt_middle[j, i] = 5
                #bt = bt_lower
            elif middle[j, i] == middle[j-1, i-1] + sub_mat[(string2[j-1], string1[i-1])]:
                bt_middle[j, i] = 6
                #i-1, j-1
            elif middle[j, i] == upper[j, i]:
                bt_middle[j, i] = 7
                #bt_upper
                
#             if bt_array[j-1,i-1] == 2 or bt_array[j-1,i-1] == 3:
#                 option1 = dag_array[j-1,i]-extension
#             else:
#                 option1 = dag_array[j-1,i]-opening
#                 
#             if bt_array[j-1,i-1] == 2 or bt_array[j-1,i-1] == 3:
#                 option2 = dag_array[j,i-1]-extension
#             else:
#                 option2 = dag_array[j,i-1]-opening
#                 
#             option3 = dag_array[j-1,i-1]+sub_mat[(string2[j-1], string1[i-1])]
#             
#             best_option = max(option1, option2, option3)
#             dag_array[j,i] = best_option
# 
#             
#             if best_option == option1:
#                 bt_array[j, i] = 1
#             
#             if best_option == option2:
#                 bt_array[j, i] = 2
#                 
#             if best_option == option3:
#                 bt_array[j, i] = 3
    alignment1 = ""
    alignment2 = ""
    i = len(string1)
    j = len(string2)
    bt = bt_middle
    while i > 0 or j > 0:
        back = bt[j,i]
        #print(i)
        #print(j)
        if back == 1:
            alignment1 = string1[i-1] + alignment1
            alignment2 = "-" + alignment2
            i -= 1
            
        elif back == 2:
            alignment1 = string1[i-1] + alignment1
            alignment2 = "-" + alignment2
            i -= 1
            bt = bt_middle
            
        elif back == 3:
            alignment1 = "-" + alignment1
            alignment2 = string2[j-1] + alignment2
            j -= 1
            
        elif back == 4:
            alignment1 = "-" + alignment1
            alignment2 = string2[j-1] + alignment2
            j -= 1
            bt = bt_middle
                
        elif back == 5:
            bt = bt_lower
        
        elif back == 6:
            alignment1 = string1[i-1] + alignment1
            alignment2 = string2[j-1] + alignment2
            i -= 1
            j -= 1
            
        else: #back = 7
            bt = bt_upper

    return middle[-1,-1], alignment1, alignment2
        
### TEST
sample_score, sample_alignment1, sample_alignment2 = affine_gap(sample_string1, sample_string2, sample_opening, sample_extension)
print(sample_score)
print(sample_alignment1)
print(sample_alignment2)

### OUTPUT
