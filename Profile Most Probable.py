############################################Profile most probable k-mer
import itertools
import numpy as np

sample_DNA = "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT"
sample_profile = np.array([
    [0.2, 0.2, 0.3, 0.2, 0.3],
    [0.4, 0.3, 0.1, 0.5, 0.1],
    [0.3, 0.3, 0.5, 0.2, 0.4],
    [0.1, 0.2, 0.1, 0.1, 0.2]
])

def read2C(filename):
    """Read 2C input data"""
    with open(filename, "r") as file:
        sequence = file.readline().strip()
        k = int(file.readline().strip())
        profile = np.loadtxt(file)  # Default type is float
        return sequence, k, profile

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

sequence, k, profile = read2C("rosalind_ba2c.txt")
print(profile_most_probable_kmer(sequence, k, profile))