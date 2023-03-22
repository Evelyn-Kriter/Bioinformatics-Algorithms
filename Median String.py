

##############################################find median string
import itertools
import numpy as np

sample_DNA = [
'AAATTGACGCAT',
'GACGACCACGTT',
'CGTCAGCGCCTG',
'GCTGAGCACCGG',
'AGTACGGGACAG'
]

test_DNA = [
    'CTGTCTAGGTCGTTGAAACCCCGAAACTGAGTTGAGACACGA',
    'GGTCTGACCGGCGACCTGGACTGATCTCGATTGAAATAGAAA',
    'ATCCGGCAACGTATGAAAATGAGTCTATTAATCTGAATCAAG',
    'CATTGTCAATGTTGCCCCTGCAAAGAAAGTAGCACCCTGAAA',
    'GTGAAAGCTTCGGGTCTGAGCTACGAGATTAGATCCTGCTTC',
    'TCGCCCCTGAAAGTAATGAGGCTTTAGTGGTTTACGAAGGGA',
    'CGCTACCCATCCATGCCATTGAAACCTCTGTGTACTTGAGCC',
    'CTGAAAATTGGTTCGAACCAATGCTCGCCTGAGTGCCATTAC',
    'TAGTAACGAGCATTAATCGTTAGTAGGACGACTAACGTGAAA',
    'CTGAAAAGTAGTTACATTCGCTATGTCGGCCAGACAACGGGT'
]
                         

def all_kmers(k):
    #Generate list of possible DNA k-mers of length k
    return ["".join(kmer) for kmer in itertools.product("ACGT",repeat=k)]

def unique_kmer_composition(sequence, k):
    #Return set of unique k-mers in sequence
    #sequence is a string
    #k is an int
    unique_kmers = []
    for i in range(0, len(sequence)-k+1):
        unique_kmers.append(sequence[i:i+k])
    return unique_kmers

def hamming_distance(pattern1, pattern2):
    #Compute the hamming distance for patterns of identical length
    #pattern1 is a string
    #pattern2 is a string
    h_d = 0
    for i in range(0, len(pattern1)):
        if pattern1[i] != pattern2[i]:
            h_d = h_d + 1
    return h_d
    

def min_distance_to_sequence(pattern, sequence, k):
    #Minimum distance between pattern of length k and all k-mers in sequence
    #pattern is a string
    #sequence is a string
    #k is an int
    list_of_kmers = unique_kmer_composition(sequence, k)
    #winning_kmer = ""
    lowest_score = len(pattern) + 1
    
    for i in range(0, len(list_of_kmers)):
        distance = hamming_distance(pattern, list_of_kmers[i])
        if distance < lowest_score:
            #winning_kmer = list_of_kmers[i]
            lowest_score = distance
    return lowest_score
            
        

def distance_to_all_sequences(pattern, sequences, k):
    #Sum of all distances between pattern and strings in sequences
    #pattern is a string
    #sequences is a list
    all_row_sum_distance = 0
    
    for i in range(0, len(sequences)):
        row_min = min_distance_to_sequence(pattern, sequences[i], k)
        all_row_sum_distance = all_row_sum_distance + row_min
    
    return all_row_sum_distance

def median_string(k, DNA):
    #k is an int
    #DNA is a list of strings
    possible_kmers = all_kmers(k)
    winner = ""
    lowest_score = k*len(DNA)
    
    for i in range(0, len(possible_kmers)):
        score = distance_to_all_sequences(possible_kmers[i], DNA, k)
        #print(possible_kmers[i])
        if score < lowest_score:
            lowest_score = score
            winner = possible_kmers[i]
            
    return winner

print("The winner is: "+ median_string(6, test_DNA) + "!")


    
    



    