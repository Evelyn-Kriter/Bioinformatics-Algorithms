sample_genome = (["ACCGA",
                  "CCGAA",
                  "CGAAG",
                  "GAAGC",
                  "AAGCT"])


def read2B(filename):
    """Read rosalind.info 2B input data"""
    with open(filename, "r") as file:
        kmers = [line.strip() for line in file.readlines()]
        return kmers

def reconstruct_genome(kmers):
    genome = kmers[0]
    i = 0
    for strips in range(0, len(kmers)):
        genome = genome[:i] + kmers[strips]
        i += 1
    return genome

print(reconstruct_genome(sample_genome))

#list_of_kmers = read2B("rosalind_ba3b.txt")
#print(reconstruct_genome(list_of_kmers))