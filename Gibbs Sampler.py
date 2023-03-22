import itertools
import numpy as np
import random

sample_DNA = ([
    'CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA',
    'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
    'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
    'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
    'AATCCACCAGCTCCACGTGCAATGTTGGCCTA',
    ])

def read2B(filename):
    """Read rosalind.info 2B input data"""
    with open(filename, "r") as file:
        k, t, N = file.readline().split()
        
        kmers = [line.strip() for line in file.readlines()]
        return kmers, int(k), int(t), int(N)
    
def random_motif_picker(DNA, k, t):
    Motif_indexes = [random.randint(1, len(DNA[0])-k) for l in range(0, t)]
    Motifs = []
    for m in range(0, t):
        sequence = DNA[m]
        n = Motif_indexes[m]
        Motifs.append(sequence[n:n+k])
    return Motifs

def construct_profile(motifs, k):
    profile = np.full((4,k), 1, dtype=np.int_) # Create a 4xk array with pseudocounts
    for pattern in motifs:
        for j in range(0, len(pattern)): #create a profile for the kmer
            profile["ACGT".find(pattern[j]), j] += 1 #increment each spot by one
                
    prob_profile = profile/(len(motifs)+4*1) #turn profile tally into a profile of probabilities
    return prob_profile

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

def unique_kmer_composition(sequence, k):
    #Return set of unique k-mers in sequence
    #sequence is a string
    #k is an int
    unique_kmers = []
    for i in range(0, len(sequence)-k+1):
        unique_kmers.append(sequence[i:i+k])
    return unique_kmers

def Gibbs_Sampler(DNA, k, t, N):
    Motifs = random_motif_picker(DNA, k, t) #list of randomly selected motifs
    BestMotifs = Motifs
    temp_profile = construct_profile(BestMotifs, k) 
    BestMotifScore = Score(BestMotifs, temp_profile)
    
    for j in range(0, N):
        #delete random motif in DNA
        i = random.randint(0, t-1)
        new_Motifs = Motifs[:]
        del new_Motifs[i]
        
        #construct a profile without that motif
        profile = construct_profile(new_Motifs, k)
        
        #find scores of kmers in deleted sequence based on the profile
        kmers_in_del_seq = unique_kmer_composition(DNA[i], k)
        list_of_kmer_scores = []
        for p in range(0, len(kmers_in_del_seq)):
            #score_kmer = Score(kmers_in_del_seq[p], profile)
            score_kmer = Score([kmers_in_del_seq[p]], profile)

            list_of_kmer_scores.append(score_kmer)
        
        #randomly choose kmer based on probabilities 
        kmer_probs_array = np.array(list_of_kmer_scores)
        kmer_probs_array /= np.sum(kmer_probs_array)
        random_kmer = np.random.choice(range(len(kmer_probs_array)), p=kmer_probs_array)
        
        #replace deleted kmer with new random kmer
        Motifs[i] = DNA[i][random_kmer:random_kmer+k]
        
        new_profile = construct_profile(Motifs, k)
        MotifScore = Score(Motifs, new_profile)
        
        if BestMotifScore < MotifScore:
            BestMotifs = Motifs[:]
            BestMotifScore = MotifScore
        
    return BestMotifs, BestMotifScore

def run_algorithm(number):
    #DNA, k, t, N= read2B("rosalind_ba2g.txt")
    overall_score = 0 #record score
    best_answer = []
    
    for times in range(0, number):
        #answers, answers_score = Gibbs_Sampler(DNA, k, t, N)
        answers, answers_score = Gibbs_Sampler(sample_DNA, 8, 5, 100)

        if answers_score > overall_score:
            overall_score = answers_score
            best_answer = answers
    return best_answer

motifs = run_algorithm(20)
for motif in range(0, len(motifs)):
    print(motifs[motif])
