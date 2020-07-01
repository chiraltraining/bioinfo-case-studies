# screed a library for reading in FASTA/FASTQ
import screed
with screed.open("../data/sequence1.fasta") as seqfile:
    for read in seqfile:
        print(read.sequence)

seq = read.sequence 

print(seq)

len(seq)

def readFastaFile(seqfile):
    import screed
    with screed.open("../data/sequence1.fasta") as seqfile:
        for read in seqfile:
            seq = read.sequence
    return seq


seq1 = readFastaFile("../data/sequence1.fasta")

seq2 = readFastaFile("../data/sequence2.fasta")

seq1

len(seq1)

seq2

len(seq2)

from Bio import pairwise2
from Bio.pairwise2 import format_alignment

seq1 = readFastaFile("../data/sequence1.fasta")
seq3 = readFastaFile("../data/sequence3.fasta")

alignments = pairwise2.align.globalxx(seq1[10:21], seq3[30:41])

for a in alignments:
    print(format_alignment(*a))

seq3 = readFastaFile("../data/sequence3.fasta")

len(seq3)

alignments = pairwise2.align.localxx(seq1[10:30], seq3[30:41])
for a in alignments:
    print(format_alignment(*a))

alignments = pairwise2.align.globalms(seq1[10:30], seq3[30:41], 2, -1, -0.5, -0.1)
for a in alignments:
    print(format_alignment(*a))