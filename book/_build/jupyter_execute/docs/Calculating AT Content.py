def read_seq(inputfile):
    """
    Reads the DNA sequence file and returns as strings. Removes all the special characters.
    """
    with open(inputfile, "r") as f:
        seq = f.read()
        seq = seq.replace("\n", "")
        seq = seq.replace("r", "")
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

# Read the DNA sequence file 
dna = read_seq("../data/dna.txt")
# Calculate GC Content 
result = calculate_at_content(dna)
print(result)