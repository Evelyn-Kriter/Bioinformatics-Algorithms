import itertools
import numpy as np
from Bio.Align import substitution_matrices

string1 = "AATCTGCCTTGTGGTCGGCTCTTACCTTCAGGCTGCTCTGAGCCCAGAGCAGAATGGTCATCACAGCTCTCCTCAACTTGGCATTGCCTGAGATCAGGATGGCTGCATGCCCAGAGGGACAAGCTGCCATTATCCCAACACAAACCATCACCCCTATTTTGTCGCGCCACAGAATCAGTAGGGGCACAGAGATGAAGGCAGC"
#sample_string2 = "TAGATA"
dna_mat = substitution_matrices.Array(alphabet="ACGT", data=np.array([[1, -1, -1, -1],[-1, 1, -1, -1], [-1, -1, 1, -1],[-1, -1, -1, 1]]))


# import pysam
# with pysam.AlignmentFile("NA12878-ready.sample.bam", mode="rb") as bam_reader:
#     for read in bam_reader:
#         print(read.query_sequence)

def read2B(filename):
    """Read text input data"""
    with open(filename, "r") as file:
        sequences = [line.strip() for line in file.readlines()]
        return sequences

def align_read(string1, string2, penalty):
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
                
    #print(dag_array)
    #print(bt_array)
    
    #print out the backtracking
    i = np.argmax(dag_array[-1,:])
    j = len(string2)
    orig_i = i
    orig_j = j
    while j > 0:  
        if bt_array[j,i] == 1: #insertion
            #do nothing
            j -= 1
        if bt_array[j,i] == 2: #deletion
            alignment1 = string1[i-1] + alignment1
            alignment2 = "-" + alignment2
            i -= 1
        if bt_array[j,i] == 3: #match
            alignment1 = string1[i-1] + alignment1
            alignment2 = string2[j-1] + alignment2
            j -= 1
            i -= 1
    
    prefix = string1[:orig_i-len(alignment1)]
    suffix = string1[orig_i:]
    alignment1 = prefix + alignment1 + suffix
    alignment2 = ("-" * len(prefix)) + alignment2 + ("-" * len(suffix))

    return (alignment1, alignment2)

### TEST
# sample_alignment1, sample_alignment2 = align_read(sample_string1, sample_string2, 1)
# print(sample_alignment1)
# print(sample_alignment2)
    
### OUTPUT
sequences = read2B("NA12878.ready.sample.txt")

with open("output_Read_Alignment.txt", "w") as file:
    print(string1, file=file)
    for f in range(0, len(sequences)):
#         string1 = sequences[0]
#         string2 = sequences[1]
#         if len(sequences[0]) > len(sequences[f]):
#             string1 = sequences[0]
#             string2 = sequences[f]
#         else:
#             string1 = sequences[f]
#             string2 = sequences[0]
        alignment1, alignment2 = align_read(string1, sequences[f], 1)
        print(alignment2, file=file)