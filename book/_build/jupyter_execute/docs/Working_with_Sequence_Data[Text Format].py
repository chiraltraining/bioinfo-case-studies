# Text Data Handling in Python

## Steps in Text Data Handling
- Exploring Data 
- Reading Data 
- Cleaning Data 
- Text Data Handling Template

## Data Exploring

# Examine first few lines  
!head ../data/dna.txt

# Exmaine last few lines  
!tail ../data/dna.txt

## Reading Data 

# open file in reading mood 
with open("../data/dna.txt", "r") as seqfile: 
    # read data and store as seq 
    seq = seqfile.read() 

# Show dataset 
seq

# Part of dataset 
seq[1:100]

# length 
len(seq)

# find the special characters in sequence 
1175-1157 

## Cleaning Data 
- Remove `\n`(newlines) 
- Remove `\t`(tabs) 

# remove newline 
seq = seq.replace("\n", "")

# take a look after removing newline 
seq 

# now remove tab character 
seq = seq.replace("\t", "") 

# print entire sequence 
seq 

# we get the result: the actual length of sequence 
len(seq)

## Template for Text Data Handling

def readSeq(inputfile): 
    """Reads the DNA sequence and returns as string. Removes all the special characters.""" 
    # opne data 
    with open(inputfile, "r") as seqfile: 
        # read data 
        seq = seqfile.read()
        # remove special characters \n and \t 
        seq = seq.replace("\n", "")
        seq = seq.replace("\t", "") 
    return seq 

# read data with readSeq function
clean_seq = readSeq("../data/dna.txt")

clean_seq

# slicing of sequence 
clean_seq[1:100]

# check length 
len(clean_seq) 