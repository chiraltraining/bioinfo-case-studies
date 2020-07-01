 # `FASTQ` Format Handling

## Tasks 
- Exploring Data 
- Reading Data 
- Cleaning Data 
- `FASTQ` Format Handling Template

# get data 
!wget http://d28rh4a8wq0iu5.cloudfront.net/ads1/data/SRR835775_1.first1000.fastq

## Exploring Data 

# examine first few lines 
!head ../data/SRR835775_1.first1000.fastq

# examine first few lines 
!tail ../data/SRR835775_1.first1000.fastq

# convert quality score character to numeric 
ord("?")

ord("B")

ord("#")

## Reading Data 

# open file 
with open("../data/SRR835775_1.first1000.fastq", "r") as f: 
    seq = f.read()

seq[:500]

# remove name line 
with open("../data/SRR835775_1.first1000.fastq", "r") as f: 
    seq = f.readline()
    seq = f.read()

seq[:500]

# read base sequence 
with open("../data/SRR835775_1.first1000.fastq", "r") as f: 
    seq = f.readline()
    seq = f.read()

seq[:500]

seq = seq.replace("\n", "")

seq[:500]

seq = seq.replace("\t", "")

seq[:500]

## FASTQ Format Handling Template

def readFastq(filename):
    """Reads FASTQ file and remove the special characters!"""
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline() # skip name line
            seq = fh.readline().rstrip() # read base sequence
            fh.readline() # skip placeholder line
            qual = fh.readline().rstrip() #base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

seqs, quals = readFastq('../data/SRR835775_1.first1000.fastq')

seqs[:20]

len(seqs)

quals[:20]

ord('#') 

ord('?')