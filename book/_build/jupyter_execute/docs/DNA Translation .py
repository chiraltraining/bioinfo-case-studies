# DNA Translation  ⇒ DNA to Protein
**Translation Theory : DNA ⇒ RNA ⇒ Protein**

## Tasks
1. Manually download DNA and Protein(to check our solution) sequence from NCBI.
2. Import DNA data into Python 
3. Create an algorithm for DNA translation
4. Check if translation matches your dowloaded protein sequence

## Task-1: Manually Download Sequence File
- Go to NCBI and select Nucleotide database 
- Enter the accession code, **NM_207618.2**
- Click on **FASTA**
- Copy and paste(without sequence info) the sequence into a text editor.
- Save both sequence as text file(.txt) into a folder called **data**

## Task-2: Import DNA Sequence Data

inputfile = "../data/dna.txt"

# Read the DNA sequence file 
f = open(inputfile, "r")
seq = f.read()

# To see the DNA sequence 
seq

# Check length of DNA seq.
len(seq)

> The sequence length is not correct! the original size of this sequence is 1157. The 18 characters are the special characters like **\n**(new line), **\r**(tab). The special characters affect the sequence quality and increase the size of this sequence. So, we have to process/clean this sequence file.

# Print the DNA sequence 
print(seq)

# New lines replace with empty space to remove the extra characters from sequence file  
seq = seq.replace("\n", "")

seq

# Check sequence length again! 
len(seq)

That's correct! same as NCBI.

# To see sequence again
print(seq)

# Tabs replace with empty space to remove extra character from sequence file 
seq = seq.replace("\r", "")

len(seq)

print(seq)

# Codon table for translation
table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
    }

# To see the table 
table

# Lookup / Extract infromation from the codon table 
table[CAA]

Error occurs, because "CCA" is a string but we tried to access without quote. That's why Python interpreter give an error.

# Extract "CCA"
table["CAA"]

table["CCT"]

## Task-3: Create an Algorithm
- Step-1: Check the length of sequence is divisible by 3.
- Step-2: Look each 3 letter strings in table and store result.
- Step-3: Continue lookups untill reading end of sequence.

protein = ""
if len(seq) % 3 == 0: 
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        protein += table[codon]


### Create a function to translate DNA sequence into protein.

def translate(seq): 
    """
    Reads sequence file as input. Transcribe input sequence file into RNA, then translate into Protein.
    """
    
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
    }
    
    protein = ""
    if len(seq) % 3 == 0: 
        for i in range(0, len(seq), 3):
            codon = seq[i:i+3]
            protein += table[codon]
    return protein


# Test translate function
translate("ATA")

translate("TCA")

def translate(seq): 
    """
    Reads sequence file as input. Transcribe input sequence file into RNA, then translate into Protein.
    """
    
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
    }
    
    protein = ""
    if len(seq) % 3 == 0: 
        for i in range(0, len(seq), 3):
            codon = seq[i:i+3]
            protein += table[codon]
    return protein


help(translate)

translate("GCC")

# Slicing sequence 
seq[40:50]

### Create a function to read sequence file 

def read_seq(inputfile): 
    """Reads and returns the input sequences with special characters removed"""
    with open(inputfile, "r") as f: 
        seq = f.read() 
    seq = seq.replace("\n", "")
    seq = seq.replace("\r", "") 
    
    return seq 

# Read DNA sequence by using read_seq() function
dna = read_seq("../data/dna.txt")

# Call translate function
translate(dna)

 >Python gives us empty string! because **sequence is not divisible by 3**. NCBI translation starts with 21 and ends with 938. but in Python translation should start at 20 position and end at 938

# Empty string cause of length problem
len(seq) % 3 

# Call translation from 20 to 938
translate(dna[20:938])

## Task-4: Comparison

# Read protein sequence
prt = read_seq("../data/protein.txt")


# Not identical cause of stop codon
prt

translate(dna[20:935])

prt 

# Comparison
prt == translate(dna[20:935])

# Exclude the last character
prt == translate(dna[20:938])[:-1]