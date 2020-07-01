# BioPython-I
## Agenda 
- Sequence processing
- sequence operations
- Reading and writing different file formats
- Sequence annotations

# Load Biopython library & Functions
import Bio
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from Bio.Seq import transcribe, back_transcribe, translate, complement, reverse_complement 

# Check Biopython version 
Bio.__version__

## Sequence Operations

# Sequence 
seq = Seq("GGACCTGGAACAGGCTGAACCCTTTATCCACCTCTCTCCAATTATACCTATCATCCTAACTTCTCAGTGGACCTAACAATCTTCTCCCTTCATCTAGCAGGAGTC")

# Alphabet
seq.alphabet

# Check type 
type(seq.alphabet)

# Find sub-sequence: if TRUE <- SubSeq Position, else <- return -1 
seq.find("ATC")

seq.find("ATGC")

# Number of `A`
seq.count("A")

# Number of `C`
seq.count("C")

# Number of `T`
seq.count("T")

# Number of `G`
seq.count("G")

# K-mer analysis, K = 2(AA)<--dimer
seq.count("AA")

# K-mer analysis, K = 3(AAA)<--trimer
seq.count("AAA")

# Count frequency of nucleotides 
from collections import Counter
freq = Counter(seq)
print(freq)

# Reverse 
print(f'RefSeq: {seq}')
rev = str(seq[::-1])
print(f'RevSeq: {rev}')

# Complement
print(f'RefSeq: {seq}')
com = seq.complement()
print(f'ComSeq: {com}')

# Reverse complement 
print(f'RefSeq: {seq}')
rev_com = seq.reverse_complement()
print(f'RevCom: {rev_com}')

# Transcription(DNA ==> RNA)
print(f'DNA: {seq}')
rna = seq.transcribe()
print(f'RNA: {rna}')

# Back Transcription(RNA ==> DNA)
print(f'RNA: {rna}')
dna = rna.back_transcribe()
print(f'DNA: {dna}')

# Translation(DNA ==> Protein)
print(f'DNA: {seq}')
prt = seq.translate()
print(f'Protein: {prt}')

# Let's varify the protein with length property
len(seq)

# Make codons 
len(seq) % 3 

# Number of codons 
len(seq) / 3 

# Now varify the protein length 
len(prt)

# Translation(DNA ==> Protein) Stop translation when found stop codon 
print(f'DNA: {seq}')
prt = seq.translate(to_stop=True)
print(f'Protein: {prt}')

# Translation(DNA ==> Protein) for Mitochondrial DNA 
print(f'DNA: {seq}')
prt = seq.translate(to_stop=True, table=2)
print(f'Protein: {prt}')

## Handling Files

for seq_record in SeqIO.parse("./dengue/den1.fasta", "fasta"):
    ID = seq_record.id 
    seqs = seq_record.seq[:100]
    rep = repr(seq_record)
    length = len(seq_record)

# ID
print(ID)

# Sequence 
print(seqs)

# Representation
print(rep)

# Length 
print(length)

# Print the first nucleotide of each codon 
seqs[0::3]

# Print the first codon position
seqs[1::3]

# Print the second codon position
seqs[2::3]

# Sequence Length Comparison
seq1 = Seq("TTGTGGCCGCTCAGATCAGGCAGTTTAGGCTTA")
seq2 = Seq("ATTTATAGAAATGTGGTTATTTCTTAAGCATGGC")
seq1 == seq2

# Mutable sequence 
mut_seq = MutableSeq("TTGTGGCCGCTCAGATCAGGCAGTTTAGGCTTA")
print(f'MutSeq: {mut_seq}')
mut_seq[5] == "C"
print(mut_seq)
mut_seq.remove("T") 
print(mut_seq)
mut_seq.reverse()
print(mut_seq)

!wget http://d28rh4a8wq0iu5.cloudfront.net/ads1/data/SRR835775_1.first1000.fastq

# Working with Fastq files 
for record in SeqIO.parse("SRR835775_1.first1000.fastq", "fastq"):
    print(record)
    
    print(record.seq)
    print(record.letter_annotations['phred_quality'])
    
    

quals = [record.letter_annotations['phred_quality'] for record in SeqIO.parse("SRR835775_1.first1000.fastq", "fastq")]

import matplotlib.pyplot as plt 
plt.hist(quals, bins=10)
plt.title("Distribution of Phred Quality Score")
plt.xlabel("Base Position")
plt.ylabel("Phred Score")
plt.show()

sequences = [record.seq for record in SeqIO.parse("SRR835775_1.first1000.fastq", "fastq")]

sequences[:100]