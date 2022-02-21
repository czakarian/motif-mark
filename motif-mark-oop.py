#!/usr/bin/env python

"""Given a a list of gene sequences and motifs, this program will produce figures to visualize 
the locations of each motif throughout a gene sequence. """

import argparse
from asyncio import start_server
from curses.ascii import islower
import Bioinfo
import cairo
import re 
import numpy as np

iupac_symbols = {'w':'at', 's':'cg', 'm':'ac', 'k':'gt', 'r':'ag', 'y':'ct', 'b':'cgt', 'd':'agt', 'h':'act', 'v':'acg', 'n':'acgt'}

class Motif:
    def __init__(self, motif):
        '''This is how a Motif is made.'''

        self.motif = motif.lower()
        self.motif_regex = self.get_regex()
        self.color = ""

    def get_regex(self):
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
        self.exon_start = 0
        self.exon_length = 0
        self.find_exon()
    
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
        self.exon_start = start
        self.exon_length = end - start 

class Draw:
    def __init__(self, seqs):
        ''''''
        self.seqs = seqs

        self.surface = cairo.SVGSurface("seq.svg", 1000, 1000)
        self.context = cairo.Context(surface)

        self.line_start_x = 50
        self.line_start_y = 100 

    def draw_header(self):
        pass

    def draw_lines(self):
        print("hi")

    def draw_exons(self):
        print("hi")

    def draw_motifs(self):
        print("hi")

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

# read in list of Sequence objects
seq_objs = []
header = ""
with open(ol_output, "r") as fr:
    for line in fr:
        line = line.strip()
        if line[0] == ">":
            header = line
        else:
            seq_objs.append(Sequence(header, line))

# read in list of motifs
motif_objs = []
with open(motifs, "r") as fr:
    for line in fr:
        line = line.strip()
        motif_objs.append(Motif(line))



out_filename = fasta.split('.')[0]



# draw 
surface_height = len(seq_objs)*100
surface = cairo.SVGSurface(out_filename + ".svg", 1100, surface_height)
context = cairo.Context(surface)

line_start_x = 15
line_start_y = 100 

for seq in seq_objs:
    
    # add header text
    context.set_source_rgb(0, 0, 0)
    context.set_font_size(13)
    context.move_to(line_start_x, line_start_y - 20)
    context.show_text(seq.header)

    # draw sequence line
    seq_length = len(seq.sequence)
    line_end_x = seq_length + line_start_x 
    line_end_y = line_start_y 

    context.set_source_rgb(0, 0, 0)
    context.set_line_width(1)
    context.move_to(line_start_x, line_start_y) 
    context.line_to(line_end_x, line_end_y)
    context.stroke()

    # draw exon
    exon_start = seq.exon_start
    exon_length = seq.exon_length
    rect_start_x = exon_start  
    rect_start_y = line_start_y - 10
    rect_length = exon_length 
    rect_width = 20

    context.set_source_rgb(0, 0, 0)
    context.rectangle(rect_start_x, rect_start_y, exon_length, rect_width)
    context.fill()

    # draw motifs 
    colors = np.array([(189,69,71), (71,122,198), (69,150,62), (214,145,33), (151,79,176)])/255

    for c,m in enumerate(motif_objs):
        positions = m.find_motif(seq)
        m.color = colors[c]
        print(seq.header)
        print(m.motif)
        print(positions)
        context.set_source_rgb(m.color[0], m.color[1], m.color[2])
        for p in positions:     
            context.rectangle(line_start_x + p, rect_start_y, len(m.motif), rect_width)
            context.fill()
    line_start_y += 60 

# make legend
leg_pos_x = 15 
leg_pos_y = 40
context.move_to(leg_pos_x, leg_pos_y)   
context.set_source_rgb(0, 0, 0)
context.set_font_size(13)
context.select_font_face("Courier", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
context.show_text("Motifs:")
leg_pos_x += 60
for m in motif_objs:
    context.rectangle(leg_pos_x, leg_pos_y - 11, 15, 15) 
    context.set_source_rgb(m.color[0], m.color[1], m.color[2])
    context.fill()
    context.move_to(leg_pos_x + 20, leg_pos_y)   
    context.show_text(m.motif)
    leg_pos_x += 100

surface.write_to_png(out_filename + '.png')
surface.finish()