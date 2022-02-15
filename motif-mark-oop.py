#!/usr/bin/env python

""" This program  """

import argparse
from curses.ascii import islower
import Bioinfo
import cairo
import itertools

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
    def __init__(self, seq):
        '''This is how a Motif is made.'''

        self.seq = seq.lower()
        self.possible_motifs = set()
        self.__generate_motifs()

    def __generate_motifs(self):
        y_indices = [] # will store the indices in motif string where pyrimidines (y) are located
        for i,char in enumerate(self.seq):
            if(char == "y"):
                y_indices.append(i)
        # generate a list of the possible 'c' and 't' combinations to use in place of 'y'
        combos = list(itertools.product("ct", repeat=len(y_indices))) 

        for c in combos:
            seq_list = list(self.seq)
            for i,n in enumerate(y_indices):
                seq_list[n] = c[i]
            self.possible_motifs.add("".join(seq_list))

    def get_motifs(self):
        return self.possible_motifs

class Sequence:
    def __init__(self, header, sequence):
        '''This is how a Sequence is made.'''

        self.header = header
        self.sequence = sequence
        self.exon = ""
    
    def __find_exon(self):
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

    def get_exon(self):
        return self.exon

    def get_sequence(self):
        return self.sequence

    def get_header(self):
        return self.sequence

class MotifImage:
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
    print(len(i.get_motifs()))


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
