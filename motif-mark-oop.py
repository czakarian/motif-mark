#!/usr/bin/env python

"""Given a a list of gene sequences and motifs, this program will produce figures to visualize 
the locations of each motif throughout a gene sequence. """

import argparse
from curses.ascii import islower
import Bioinfo
import cairo
import itertools
import re 

iupac_symbols = {'w':'at', 's':'cg', 'm':'ac', 'k':'gt', 'r':'ag', 'y':'ct', 'b':'cgt', 'd':'agt', 'h':'act', 'v':'acg', 'n':'acgt'}

def find_motif(seq, motif):
    """This function takes a sequence string and a motif string and returns a list of the positions of that motif in the sequence"""
    regex_str = ""
    positions = []
    for c in motif:
        if c in iupac_symbols:
            regex_str += "[" + iupac_symbols[c] + "]"
        else:
            regex_str += c
    for match in re.finditer(regex_str, seq):
        positions.append(match.start())
    return positions 

#print(find_motif('acgtcgdkdfkjdfktcgtsdfkjjsacgcfkdjfajjt', 'wcgy'))

class Exon:
    def __init__(self, start, end):
        '''This is how an Exon is made.'''
        
        self.start = start
        self.end = end 

    def get_start(self):
        return self.start
    
    def get_end(self):
        return self.end

class Motif:
    def __init__(self, motif):
        '''This is how a Motif is made.'''

        self.motif = motif.lower()
        self.motif_regex = self.__get_regex()
        self.color = 'red'

    def __get_regex(self):
        """This function takes takes the motif string and generates the regex expression to search for it."""
        regex_str = ""
        for c in self.motif:
            if c in iupac_symbols:
                regex_str += "[" + iupac_symbols[c] + "]"
            else:
                regex_str += c
        return regex_str

    def find_motif(self, Sequence):
        """This function takes a Sequence object and returns a list of the positions of the motif in the given sequence."""
        positions = []
        for match in re.finditer(self.motif_regex, Sequence.sequence):
            positions.append(match.start())
        return positions 

class Sequence:
    def __init__(self, header, sequence):
        '''This is how a Sequence is made.'''

        self.header = header
        self.sequence = sequence
        self.exon_coords = self.find_exon()
    
    def find_exon(self):
        start = 0
        end = 0 
        start_found = False
        end_found = False
        for i,c in enumerate(self.sequence):
            if(start_found and end_found):
                break
            elif(c.isupper() and not start_found):
                start = i
                start_found = True
            elif(c.islower() and start_found and not end_found):
                end = i - 1
                end_found = True       
        self.exon = Exon(start, end)

class Draw:
    def __init__(self):
        ''''''
        self.seqs = []

    def generate_image(self):
        print("hi")

    
def get_args():
    """This function returns the parser arguments entered in command line"""
    parser = argparse.ArgumentParser(description="A program to visualize motifs on sequences")
    parser.add_argument("-f", "--fasta", help="Input FASTQ file", required=True)
    parser.add_argument("-m", "--motifs", help="Input motifs file", required=True)
    return parser.parse_args()

# store the command line args in variables
args = get_args()
fasta= args.fasta
motifs = args.motifs

# read in list of motifs
motif_objs = []
with open(motifs, "r") as fr:
    for line in fr:
        line = line.strip()
        motif_objs.append(Motif(line))

for i in motif_objs:
    print(i.motif_regex)


# turn multiline fasta sequences in one-liners
ol_output = "OL_" + fasta
Bioinfo.oneLineFasta(fasta, ol_output)




# will contain sequence objects, what should key be?
seq_dict = {}
header = ""
with open(ol_output, "r") as fr:
    for line in fr:
        line = line.strip()
        if line[0] == ">":
            header = line
            seq_dict[header] = ""
        else:
            seq_dict[header] = line
