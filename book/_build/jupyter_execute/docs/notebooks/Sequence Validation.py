def validate_seq(seq): 
    """
    Checks if the DNA sequence is valid. Returns True if sequence is valid or False otherwise 
    """
    if seq.upper():
        total_counts = seq.count("A") + seq.count("T") + seq.count("G") + seq.count("C")
    else:
        seq = seq.upper()
        total_counts = seq.count("A") + seq.count("T") + seq.count("G") + seq.count("C")
        
    if total_counts == len(seq):
        return True 
    else:
        return False


with open("../data/dna.txt") as f:
    dna = f.read()
    

dna

validate_seq(dna)

seq1 = dna.replace("\n", "")

validate_seq(seq1)

with open("../output/out.txt", "w") as f:
    f.write(seq1) 

with open("../output/out.txt", "r") as f: 
    seq2 = f.read()

seq2.lower()

seq3 = validate_seq(seq2)

seq3

seq2

validate_seq("ATATC")

def validate_seq(dna_seq):
    seqm = dna_seq.upper() 
    
    total_counts = dna_seq.count("A") + dna_seq.count("T") + dna_seq.count("G") + dna_seq.count("C")
        
    if total_counts == len(dna_seq):
        return True 
    else:
        return False
    

validate_seq("aaatttctctcttcttctttgggct")

validate_seq("atcg")

## Protein visualization

Store_all = []
with open("../data/protein.txt") as protein:
    for lines in protein:
        if "ATOM   " in lines:
            lines = lines.split()
            #'ATOM', '1', 'N', 'LEU', 'A', '125', '4.329', '-12.012', '2.376', '1.00', '0.00', 'N'
            Store_all.append(map(float, lines[3:6]))

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

x,y,z = zip(*Store_all)

fig = plt.figure()
ax = Axes3D(fig)

ax.plot(x,y,z, "o")
ax.axis("off")

plt.show()


from Bio.SeqUtils import GC
from Bio import SeqIO 

for seq in SeqIO.parse("../dengue/den1.fasta", "fasta"): 
    seq = seq.seq



starts = list(range(1, len(seq)-2000, 2000))

starts

def sliding_window(sequence, winSize, step):
        numOfChunks = ((len(sequence)-winSize)/step)+1
        for i in range(0,numOfChunks*step,step):
                yield sequence[i:i+winSize]

print(sliding_window(seq, 2000, 2000))

