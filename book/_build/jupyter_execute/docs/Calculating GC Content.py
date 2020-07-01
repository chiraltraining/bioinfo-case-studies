# Calculating GC Content

## GC Content
In molecular biology and genetics, GC-content (or guanine-cytosine content) is the percentage of nitrogenous bases in a DNA or RNA molecule that are either guanine (G) or cytosine (C). GC-content may be given for a certain fragment of DNA or RNA or for an entire genome.

GC-content may be given for a certain fragment of DNA or RNA or for an entire genome. When it refers to a fragment, it may denote the GC-content of an individual gene or section of a gene (domain), a group of genes or gene clusters, a non-coding region, or a synthetic oligonucleotide such as a primer.
[See More on Wikipedia](https://en.wikipedia.org/wiki/GC-content)

## Calculation
GC-content is usually expressed as a **percentage value**, but sometimes as a **ratio** (called G+C ratio or GC-ratio). GC-content percentage is calculated as

$ GC = \frac {G+C} {Len(DNA Sequence)} \times 100% $

## Applications


### Molecular biology
In polymerase chain reaction (PCR) experiments, the GC-content of short oligonucleotides known as **primers** is often used to predict their annealing temperature to the template DNA. A higher GC-content level indicates a relatively higher melting temperature.

### Systematics
The species problem in **prokaryotic** taxonomy has led to various suggestions in **classifying** bacteria, and the adhoc committee on reconciliation of approaches to bacterial systematics has recommended use of GC-ratios in higher-level hierarchical classification.For example, the Actinobacteria are characterised as "high GC-content bacteria".In Streptomyces coelicolor A3(2), GC-content is 72%

# Calculating GC Content
1. Download a DNA sequence from [NCBI](https://www.ncbi.nlm.nih.gov/)
2. Reading and cleaning the DNA sequence.
3. Create an algorithms for calculating **GC** content.

## Download a DNA Sequence from NCBI
You can download the dataset from [NCBI](https://www.ncbi.nlm.nih.gov/). Accession ID: NM_207618.2, Enter the accession id and click on **FASTA**, copy all the sequence without *>NM_207618.2 Mus musculus vomeronasal 1 receptor, D18 (V1rd18), mRNA*.then save it as **.txt** file.

## Reading and cleaning the DNA sequence

# Read the DNA sequence 
with open("../data/dna.txt", "r") as f: 
    seq = f.read()

seq

The output is not a clean dataset because it contains a special character("\n"). It also may contain a special character("\n"). The special characters affect our analysis and may cause errors. That's why we have to clean the dataset by removing the special characters.

# Check the DNA sequence length 
len(seq)

It is incorrect! because the original length of this DNA sequence is **1157**. The special characters 18. Let's clean the dataset.

# Replace "\n" character with empty space.
seq = seq.replace("\n", "")

seq

# Now check the length
len(seq)

YES! it is the original length of the DNA sequence. But we have to remove "\r"

seq = seq.replace("\r", "")

seq

len(seq)

## Create a Function To Read Sequence File 

def read_seq(inputfile):
    """
    Reads the DNA sequence file and returns as strings. Removes all the special characters.
    """
    with open(inputfile, "r") as f:
        seq = f.read()
        seq = seq.replace("\n", "")
        seq = seq.replace("r", "")
    return seq 

# Read the seq file with rea_seq function
dna = read_seq("../data/dna.txt")

dna

# Check length 
len(dna)

## Calculate GC Content

def calculate_gc_content(seq): 
    """
    Take DNA sequence as input and calculate the GC content.
    """
    no_of_g = seq.count("G")
    no_of_c = seq.count("C")
    total = no_of_g + no_of_c 
    gc = total/len(seq) * 100
    
    return gc 

calculate_gc_content(dna)

## Final Code 

def read_seq(inputfile):
    """
    Reads the DNA sequence file and returns as strings. Removes all the special characters.
    """
    with open(inputfile, "r") as f:
        seq = f.read()
        seq = seq.replace("\n", "")
        seq = seq.replace("r", "")
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

# Read the DNA sequence file 
dna = read_seq("../data/dna.txt")
# Calculate GC Content 
result = calculate_gc_content(dna)
print(result)

So the GC content is 39.5852% <br>