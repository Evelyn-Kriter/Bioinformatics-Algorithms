import itertools
import numpy as np
import random

sample_DNA = ([
    'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
    'CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA',
    'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
    'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
    'AATCCACCAGCTCCACGTGCAATGTTGGCCTA'
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
    
def random_motif_picker(DNA, k, t):
    Motif_indexes = [random.randint(1, len(DNA[0])-k) for l in range(0, t)]
    Motifs = []
    for m in range(0, t):
        sequence = DNA[m]
        n = Motif_indexes[m]
        Motifs.append(sequence[n:n+k])
    return Motifs

def construct_profile(motifs, k):
    profile = np.full((4,k), 1, dtype=np.int_) # Create a 4xk array of ones
    for pattern in motifs:
        for j in range(0, len(pattern)): #create a profile for the kmer
            profile["ACGT".find(pattern[j]), j] += 1 #increment each spot by one
                
    prob_profile = profile/(len(motifs)+4*1) #turn profile tally into a profile of probabilities
    return prob_profile

def unique_kmer_composition(sequence, k):
    #Return set of unique k-mers in sequence
    #sequence is a string
    #k is an int
    unique_kmers = []
    for i in range(0, len(sequence)-k+1):
        unique_kmers.append(sequence[i:i+k])
    return unique_kmers

def profile_most_probable_kmer(DNA, k, profile):
    list_of_kmers = unique_kmer_composition(DNA, k)
    
    high_score = -1
    winning_kmer = ""
    
    for i in range(0, len(list_of_kmers)):
        pattern = list_of_kmers[i] #pattern should be a string
        pattern_score = 1
        
        for j in range(0, len(pattern)):
            # Access the probability for pattern[j] at index j in the k-mer
            index_prob = profile["ACGT".find(pattern[j]), j]
            pattern_score = pattern_score*index_prob
            
        if pattern_score > high_score:
            high_score = pattern_score
            winning_kmer = pattern
    
    return winning_kmer

def Score(Motifs, profile):    
    tot_score = 0
    
    for i in range(0, len(Motifs)):
        pattern = Motifs[i] #pattern should be a string
        pattern_score = 1
        
        for j in range(0, len(pattern)):
            # Access the probability for pattern[j] at j in the k-mer
            index_prob = profile["ACGT".find(pattern[j]), j]
            pattern_score = pattern_score*index_prob
        
        tot_score += pattern_score
    
    return tot_score
        

def randomized_motif_search(DNA, k, t):
    Motifs = random_motif_picker(DNA, k, t) #list of randomly selected motifs
    BestMotifs = Motifs
    temp_profile = construct_profile(BestMotifs, k) 
    BestMotifScore = Score(BestMotifs, temp_profile)
    
    while True:
        profile = construct_profile(Motifs, k)
        
        probable_kmers = []*t
        for x in range(0, t):
            probable_kmer = profile_most_probable_kmer(DNA[x], k, profile)
            probable_kmers.append(probable_kmer)
        
        Motifs = probable_kmers
        new_profile = construct_profile(Motifs, k)
        MotifScore = Score(Motifs, new_profile)
        
        if BestMotifScore < MotifScore:
            BestMotifs = Motifs[:]
            BestMotifScore = MotifScore
        else:
            return BestMotifs, BestMotifScore

def run_algorithm(number):
    overall_score = 0 #record score
    best_answer = []
    
    for times in range(0,number):
        answers, answers_score = randomized_motif_search(sample_DNA, 8, 5)
        if answers_score > overall_score:
            overall_score = answers_score
            best_answer = answers
    return best_answer

#sample_DNA, k, t = read2B("rosalind_ba2f.txt")
motifs = run_algorithm(1001)

for g in range(0, len(motifs)):
    print(motifs[g])
