# DNA Sequence Statistics-1

## Function for Reading FASTA/FASQ File

import screed # A Python library for reading FASTA and FASQ file format.
def readFastaFile(inputfile):
    """
    Reads and returns file as FASTA format with special characters removed.
    """
    with screed.open(inputfile) as seqfile:
        for read in seqfile:
            seq = read.sequence
    return seq

## Read FASTA/FASQ File

den1 = readFastaFile("../dengue/den1.fasta")
den2 = readFastaFile("../dengue/den2.fasta")
den3 = readFastaFile("../dengue/den3.fasta")
den4 = readFastaFile("../dengue/den4.fasta")

## Statistical Aanalysis

### Length of a DNA sequence

# length of DEN1 
len(den1)

# length of DEN2
len(den2)

# length of DEN3
len(den3)

len(den4)

### Base composition of a DNA sequence¶


#### By creating a function

def count_base(seq): 
    base_counts = {} 
    for base in seq: 
        if base in base_counts: 
            base_counts[base] += 1 
        else: 
            base_counts[base] = 1 
    return base_counts

# Let's call the function count_base
counts = count_base(den1)
print(counts)

#### By collections modules

from collections import Counter
def count_base(seq): 
    base_counts = Counter(seq)
    return base_counts

bases = count_base(den1)
print(bases)

### Frequency with Counter

from collections import Counter
# Base composition of DEN1
freq1 = Counter(den1)
print(freq1)

# Base composition of DEN2
freq2 = Counter(den2)
print(freq2)

# Base composition of DEN3
freq3 = Counter(den3)
print(freq3)

# Base composition of DEN4
freq4 = Counter(den4)
print(freq2)

### GC Content of DNA

import screed # A Python library for reading FASTA and FASQ file format.
def readFastaFile(inputfile):
    """
    Reads and returns file as FASTA format with special characters removed.
    """
    with screed.open(inputfile) as seqfile:
        for read in seqfile:
            seq = read.sequence
    return seq

def calculate_gc_content(seq): 
    """
    Take DNA sequence as input and calculate the GC content.
    """
    no_of_g = seq.count("G")
    no_of_c = seq.count("C")
    total = no_of_g + no_of_c 
    gc = total/len(seq) * 100
    
    return gc 

# Calculate GC Content of DEN1
den1 = readFastaFile("../dengue/den1.fasta")
result1 = calculate_gc_content(den1)
print(result1)

# Calculate GC Content of DEN2
den2 = readFastaFile("../dengue/den2.fasta")
result2 = calculate_gc_content(den2)
print(result2)

# Calculate GC Content of DEN3
den3 = readFastaFile("../dengue/den3.fasta")
result3 = calculate_gc_content(den3)
print(result3)

# Calculate GC Content of DEN4
den4 = readFastaFile("../dengue/den4.fasta")
result4 = calculate_gc_content(den4)
print(result4)

### AT Content of DNA

import screed # A Python library for reading FASTA and FASQ file format.
def readFastaFile(inputfile):
    """
    Reads and returns file as FASTA format with special characters removed.
    """
    with screed.open(inputfile) as seqfile:
        for read in seqfile:
            seq = read.sequence
    return seq

def calculate_at_content(seq): 
    """
    Take DNA sequence as input and calculate the AT content.
    """
    no_of_a = seq.count("A")
    no_of_t = seq.count("T")
    total = no_of_a + no_of_t
    at = total/len(seq) * 100
    
    return at 

# Calculate AT Content of DEN1
den1 = readFastaFile("../dengue/den1.fasta")
result1 = calculate_at_content(den1)
print(result1)

# Calculate AT Content of DEN2
den2 = readFastaFile("../dengue/den2.fasta")
result2 = calculate_at_content(den2)
print(result2)

# Calculate AT Content of DEN3
den3 = readFastaFile("../dengue/den3.fasta")
result3 = calculate_at_content(den3)
print(result3)

# Calculate AT Content of DEN4
den4 = readFastaFile("../dengue/den4.fasta")
result4 = calculate_at_content(den4)
print(result4)

### A sliding window analysis of GC content

import screed # A Python library for reading FASTA and FASQ file format.
def readFastaFile(inputfile):
    """
    Reads and returns file as FASTA format with special characters removed.
    """
    with screed.open(inputfile) as seqfile:
        for read in seqfile:
            seq = read.sequence
    return seq

def calculate_gc_content(seq): 
    """
    Take DNA sequence as input and calculate the GC content.
    """
    no_of_g = seq.count("G")
    no_of_c = seq.count("C")
    total = no_of_g + no_of_c 
    gc = total/len(seq) * 100
    
    return gc 

# Calculate the GC content of nucleotides 1-2000 of the Dengue genome(DEN-1)
seq = readFastaFile("../dengue/den1.fasta")
calculate_gc_content(seq[1:2001])

# Calculate the GC content of nucleotides 2001-4000 of the Dengue genome(DEN-1)
seq = readFastaFile("../dengue/den1.fasta")
calculate_gc_content(seq[2001:4001])

# Calculate the GC content of nucleotides 4001-6000 of the Dengue genome(DEN-1)
seq = readFastaFile("../dengue/den1.fasta")
calculate_gc_content(seq[4001:6001])

# Calculate the GC content of nucleotides 6001-8000 of the Dengue genome(DEN-1)
seq = readFastaFile("../dengue/den1.fasta")
calculate_gc_content(seq[6001:8001])

# Calculate the GC content of nucleotides 8001-10000 of the Dengue genome(DEN-1)
seq = readFastaFile("../dengue/den1.fasta")
calculate_gc_content(seq[8001:10001])

# Calculate the GC content of nucleotides 10001-10735 of the Dengue genome(DEN-1)
seq = readFastaFile("../dengue/den1.fasta")
calculate_gc_content(seq[10001:10735])

## GC content of the non-overlapping sub-sequences of size k 

def calculate_gc(seq):
    """ 
    Returns percentage of G and C nucleotides in a DNA sequence.
    """
    gc = 0 
    for base in seq:
        if base in "GCgc":
            gc += 1 
        else:
            gc = 1 
    return gc/len(seq) * 100

def gc_subseq(seq, k=2000):
    """
    Returns GC content of non − overlapping sub− sequences of size k.
    The result is a list.
    """
    res = [] 
    for i in range(0, len(seq)-k+1, k):
        subseq = seq[i:i+k]
        gc = calculate_gc(subseq)
        res.append(gc)
    return gc 
    

seq = readFastaFile("../dengue/den1.fasta")

gc_subseq(seq, 2000)

## K-mer Analysis

def build_kmers(sequence, ksize):
    kmers = []
    n_kmers = len(sequence) - ksize + 1

    for i in range(n_kmers):
        kmer = sequence[i:i + ksize]
        kmers.append(kmer)

    return kmers

seq = readFastaFile("../dengue/den1.fasta")

# Dimer
km2 = build_kmers(seq, 2)

# Count dimer
from collections import Counter
dimer = Counter(km2)
print(dimer)

# Trimer 
km3 = build_kmers(seq, 3)

# Count trimer
trimer = Counter(km3)
print(trimer)

## References 
- https://computationalgenomics.blogs.bristol.ac.uk/case_studies/haemophilus_demo