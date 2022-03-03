# Motif-Mark

Given a a list of gene sequences and motifs, this program will produce an image to visualize the locations of each motif across the gene sequences.

### Input
1. FASTA file with gene sequences 
    - Maximum of 10 sequences (≤1000 bases each)
    - Exons should be indicated by uppercase letters and introns by lowercase letters
2. Text file with a list of motifs to search 
    - One motif per line
    - Maximum of 5 motifs (≤10 bases each)

#### Argparse options:
    -f, --fasta: required arg, file path to FASTA file
    -m, --motifs: required arg, file path to motif text file

### How to format input files:

FASTA:
>\>INSR chr19:7150261-7150808
>ctctgtcctcaaaggcgttggttttgtttccacagAAAAACCTCTTCAGGCACTGGTGCCGAGGACCCTAGgtatgactcacctgtgcgacccctgg\
>\>MBNL chr3:152446461-152447003
>atgttaatgcgcttgaaccccactggcccattgccatcatgtgctcgctgcctgctaattaagACTCAGTCGGCTGTCAAATCACTGAAGCGACCCCTCGAGGCAACCTTTGACCTGgtactatgacctttcaccttttagcttggcatgtagctttattgtagatacaagttttttttt

Motifs:
>ygcy\
>GCAUG\
>catag\
>YYYYYYYYYY

### Example of how to run the program:
```
python motif-mark-oop.py -f Figure_1.dnas -m Fig_1_motifs.txt 
```

### Output
The program will output a single image (.png) containing each of the inputted sequences with motifs, introns, and exons to scale depending on their lengths. An example of the output can be seen below where exons are indicated by the black boxes and introns are indicated by the grey line.

![Example Output](Figure_1.png)
