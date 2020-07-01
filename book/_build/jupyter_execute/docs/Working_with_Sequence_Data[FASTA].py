# Introduction to Sequence Data Analysis with Python


## FASTA Format Handling
- Exploring Data 
- Reading Data 
- Cleaning Data 
- `FASTA` Format Handling Template 

## Exploring Data 

# examine first few lines  
!head ../data/Haemophilus_influenzae.fasta

# examine last few lines 
!tail ../data/Haemophilus_influenzae.fasta

## Reading Data

# open file in reading mood 
with open("../data/Haemophilus_influenzae.fasta", "r") as f: 
    # reading data 
    seq = f.read() 

seq 

# inspect specific range of sequence 
seq[1:200] 

# inspect specific range of sequence 
seq[200:500] 

# check length of sequence 
len(seq)

__Note__
- 1882752 is not actual length of this sequence 
- 1856176 is the actual sequence of this sequence

# find the no. of extra characters in sequence 
# present length - actual length[NCBI]
1882752 - 1856176

## Cleaning Data 

# remove name line / info line  
with open("../data/Haemophilus_influenzae.fasta", "r") as f: 
    seq = f.readline() # skip name line 
    seq = f.read() # read data 

seq[1:200]

# check length 
len(seq) 

__Note__
- Still 1882694 is not actual length. it means the dataset is not perfect! 
- We have to remove special characters from this sequence 

# find the no. of extra characters in sequence again 
# present length - actual length[NCBI]
1882694 - 1856176

# remove newline 
seq = seq.replace("\n", "")

# see 1 to 200
seq[1:200]

# remove tab 
seq = seq.replace("\t", "")

# see 1 to 200bp 
seq[1:200]

# check leng
len(seq)

# clean dataset? 
1856176 - 1856176

## FASTA Format Handling Template

def readFASTA(inputfile): 
    """Reads a FASTA file and removes name line and special characters""" 
    # open file 
    with open(inputfile, "r") as f: 
        # remove name line / info line 
        seq = f.readline()
        # read data 
        seq = f.read() 
        # remove special character 
        seq = seq.replace("\n", "")
        seq = seq.replace("\t", "") 
        
    return seq 

# read FASTA 
s = readFASTA("../data/Haemophilus_influenzae.fasta")

len(s) 

s[1:200]