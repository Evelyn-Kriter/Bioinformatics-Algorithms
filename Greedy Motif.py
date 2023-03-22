########################################Greedy motif search with/out pseudocounts
import itertools
import numpy as np

DNA="ACGT"

sample_DNA = ([
    'GGCGTTCAGGCA',
    'AAGAATCAGTCA',
    'CAAGGAGTTCGC',
    'CACGTCAATCAC',
    'CAATAATATTCG',
    ])

def read2B(filename):
    """Read rosalind.info 2B input data"""
    with open(filename, "r") as file:
        # Recall that the file object maintains an internal pointer so you can use
        # multiple readline calls to read individual lines.
        k, t = file.readline().split()
        
        # Similarly readlines will read the remaining lines in the file.
        # Here I am using a list comprehension to create a list of strings with any leading
        # or trailing whitespace removed.
        kmers = [line.strip() for line in file.readlines()]
        return kmers, int(k), int(t)

def unique_kmer_composition(sequence, k):
    #Return set of unique k-mers in sequence
    #sequence is a string
    #k is an int
    unique_kmers = []
    for i in range(0, len(sequence)-k+1):
        unique_kmers.append(sequence[i:i+k])
    return unique_kmers

def construct_profile(motifs, k, pseudocount):
    profile = np.full((4,k), pseudocount, dtype=np.int_) # Create a 4xk array of zeros/ones
    #denominator = 0 #keeps track of bottom integer
    #profile = pseudocount_array
    
    for pattern in motifs:
        for j in range(0, len(pattern)): #create a profile for the kmer
            profile["ACGT".find(pattern[j]), j] += 1 #increment each spot by one
                
    #denominator += 1
    #prob_profile = profile/(len(motifs)+(4*pseudocount)) #turn profile tally into a profile of probabilities
    prob_profile = profile/(len(motifs)+(4*pseudocount))
    return prob_profile
    

def hamming_distance(pattern1, pattern2):
    #Compute the hamming distance for patterns of identical length
    #pattern1 is a string
    #pattern2 is a string
    h_d = 0
    for i in range(0, len(pattern1)):
        if pattern1[i] != pattern2[i]:
            h_d = h_d + 1
    return h_d

def compute_kmer_probability(kmer, profile):
    coefficient = 1
    for n in range(0, len(kmer)): #for each character in the kmer string
        #Access the probability for character at index m in the k-mer
        index_prob = profile["ACGT".find(kmer[n]), n]
        coefficient = coefficient*index_prob
    return coefficient #probability of kmer based on profile given

def highest_probable_kmer(list_of_kmers, profile, pseudocount):
    current_highest_score = -1
    current_winning_kmer = "" #return first kmer if none
    for m in range(0, len(list_of_kmers)): #for each possible kmer in sequence
        kmer_prob = compute_kmer_probability(list_of_kmers[m], profile)
        if kmer_prob > current_highest_score:
            current_highest_score = kmer_prob
            current_winning_kmer = list_of_kmers[m]
    return current_winning_kmer

def consensus(motifs, profile):
    consensus = ""
    for i in np.argmax(counts, axis=0):
        consensus += "ACGT"[i]

    return consensus

def greedy_motif_search(DNA, k, t, pseudocount):
    #DNA is a list of strings
    #k is the kmer
    #pseudocount takes 0 or 1
    
    possible_patterns = unique_kmer_composition(DNA[0], k) #list of kmers to test
    lowest_pattern_score = float("infinity") #number of mismatches for current lowest score pattern
    lowest_pattern_score_index = 0 #index of current best
    kmer_sequence_list = []
    
    for i in range(0, len(possible_patterns)): #for each possible pattern run the algorithm
        target_pattern = possible_patterns[i] 
        #target_pattern_score = 0 #number of total mismatches for this pattern
        winning_kmer_list = [target_pattern] #list of winning kmers for each sequence
        
        
        profile = construct_profile(winning_kmer_list, k, pseudocount)
        
        for l in range(1, len(DNA)): #for each sequence
            list_of_kmers = unique_kmer_composition(DNA[l], k) #list of kmers in sequence
            winning_kmer = highest_probable_kmer(list_of_kmers, profile, pseudocount) #kmer with highest probability in sequence
            #mismatches = hamming_distance(winning_kmer, target_pattern)
            #target_pattern_score += mismatches
            winning_kmer_list.append(winning_kmer) #append winning_kmer to list
            profile = construct_profile(winning_kmer_list, k, pseudocount)
            
            #kmer_sequence_list.append(winning_kmer_list) #append list of winning kmers to list of lists
                    
        consensus = ""
        for z in np.argmax(profile, axis=0):
            consensus += "ACGT"[z]
        target_pattern_score = 0
        for winning_kmer in winning_kmer_list:
            target_pattern_score += hamming_distance(winning_kmer, consensus)
        kmer_sequence_list.append(winning_kmer_list) #append list of winning kmers to list of lists
            #target_pattern_score += mismatches
            #winning_kmer_list.append(winning_kmer) #append winning_kmer to list
        
        if target_pattern_score < lowest_pattern_score:
            lowest_pattern_score = target_pattern_score
            lowest_pattern_score_index = i    
            
    return kmer_sequence_list[lowest_pattern_score_index]
            

#with pseudocounts
DNA1, k, t = read2B("rosalind_ba2d.txt")
#DNA1, k, t = read2B("input_5.txt")
#motifs = greedy_motif_search(DNA1, k, t, 0)
#for motif in motifs:
    #print(motif)

#without pseudocounts
motifs = greedy_motif_search(DNA1, k, t, 0)
for motif in motifs:
    print(motif)
