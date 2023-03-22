sample_genome = "CAATCCAAC"
sample_k = 5


def read2B(filename):
    """Read rosalind.info 2B input data"""
    with open(filename, "r") as file:
        k = file.readline().split()
        
        sequence = [line.strip() for line in file.readlines()]
        return sequence, int(k)

def kmer_composition(sequence, k):
    #Return set of unique k-mers in sequence
    #sequence is a string
    #k is an int
    unique_kmers = []
    for i in range(0, len(sequence)-k+1):
        unique_kmers.append(sequence[i:i+k])
        
    unique_kmers.sort()
    return unique_kmers


genome, k = read2B("rosalind_ba3a.txt")
kmers = kmer_composition(genome, k)
for x in range(0, len(kmers)):
    print(kmers[x])