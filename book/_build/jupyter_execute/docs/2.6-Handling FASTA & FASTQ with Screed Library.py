# Introduction to Sequence Data Analysis with Python

import screed
with screed.open("../data/Haemophilus_influenzae.fasta") as seqfile: 
    for read in seqfile:
        seq = read.sequence    
        name = read.name 

print(name)

seq[:500]

len(seq)

## Template for Handling FASTA and FASTQ with Screed 

import screed # A Python library for reading FASTA and FASQ file format.
def readFastaFile(inputfile):
    """
    Reads and returns file as FASTA format with special characters removed.
    """
    with screed.open(inputfile) as seqfile:
        for read in seqfile:
            seq = read.sequence
    return seq

# read data 
seqs = readFastaFile("../data/Haemophilus_influenzae.fasta")

seqs[:200]

import screed # A Python library for reading FASTA and FASQ file format.
def readFastqFile(inputfile):
    """
    Reads and returns file as FASTA format with special characters removed.
    """
    with screed.open(inputfile) as seqfile:
        for read in seqfile:
            seq = read.sequence
    return seq

seqs = readFastqFile("../data/SRR835775_1.first1000.fastq")

seqs 