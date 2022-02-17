#!/usr/bin/env python

"""Given a a list of gene sequences and motifs, this program will produce figures to visualize 
the locations of each motif throughout a gene sequence. """

import argparse
from asyncio import start_server
from curses.ascii import islower
import Bioinfo
import cairo
import itertools
import re 

iupac_symbols = {'w':'at', 's':'cg', 'm':'ac', 'k':'gt', 'r':'ag', 'y':'ct', 'b':'cgt', 'd':'agt', 'h':'act', 'v':'acg', 'n':'acgt'}

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
        self.exon_start = ""
        self.exon_length = ""
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
        self.exon_start = start
        self.exon_length = end - start 

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

# turn multiline fasta sequences in one-liners
ol_output = "OL_" + fasta
Bioinfo.oneLineFasta(fasta, ol_output)

# will contain sequence objects, what should key be?
seq_objs = []
header = ""
with open(ol_output, "r") as fr:
    for line in fr:
        line = line.strip()
        if line[0] == ">":
            header = line
        else:
            seq_objs.append(Sequence(header, line))

for i in seq_objs:
    print(i.header)

# read in list of motifs
motif_objs = []
with open(motifs, "r") as fr:
    for line in fr:
        line = line.strip()
        motif_objs.append(Motif(line))

for i in motif_objs:
    print(i.motif_regex)
    print(i.find_motif(seq_objs[0]))



# draw 
surface = cairo.SVGSurface("seq.svg", 1000, 1000)
context = cairo.Context(surface)


line_start_x = 50
line_start_y = 100 

for i in seq_objs:
    seq_length = len(i.sequence)
    exon_start = i.exon_start
    exon_length = i.exon_length

    line_end_x = seq_length + line_start_x 
    line_end_y = line_start_y 

    # headers
    context.move_to(line_start_x, line_start_y - 15)
    context.show_text(i.header)

    # draw a line
    context.set_line_width(1)
    context.move_to(line_start_x, line_start_y) 
    context.line_to(line_end_x, line_end_y)
    context.stroke()

    rect_start_x = exon_start  
    rect_start_y = line_start_y - 10
    rect_length = exon_length 
    rect_width = 20

    # draw a rectangle
    context.rectangle(rect_start_x, rect_start_y, exon_length, rect_width)
    context.fill()

    line_start_y += 50 

surface.write_to_png('seq.png')
surface.finish()