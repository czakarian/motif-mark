#!/usr/bin/env python

"""Given a a list of gene sequences and motifs, this program will produce figures to visualize 
the locations of each motif throughout a gene sequence. """

import argparse
from ast import Pass
from asyncio import start_server
from curses.ascii import islower
import Bioinfo
import cairo
import re 
import numpy as np

iupac_symbols = {'w':'at', 's':'cg', 'm':'ac', 'k':'gt', 'r':'ag', 'y':'ct', 'b':'cgt', 'd':'agt', 'h':'act', 'v':'acg', 'n':'acgt', 'u':'t'}
colors = np.array([(189,69,71), (71,122,198), (69,150,62), (214,145,33), (151,79,176)])/255

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
        pattern = re.compile(self.motif_regex)
        match = pattern.search(Sequence.sequence.lower())
        while match:
            positions.append(match.start())
            match = pattern.search(Sequence.sequence.lower(), match.start() + 1) # takes care of overlaps 
        return positions

    def set_color(self, color):
        """This function sets the color of the motif for drawingg."""
        self.color = color


class Sequence:
    def __init__(self, header, sequence):
        '''This is how a Sequence is made.'''
        self.header = header
        self.sequence = sequence
        self.exon_start = 0
        self.exon_length = 0
        self.find_exon()
        self.format_header()

    def format_header(self):
        """This function takes in the fasta header and reformats it [genename chromosome x (xxx - xxx)]"""
        header_parts = self.header.split()
        gene = header_parts[0][1:]
        chr = header_parts[1].split(":")[0].split("chr")[1]
        pos = header_parts[1].split(":")[1]
        self.header = gene + " chromosome " + chr + " (" + pos + ")"
        
    def find_exon(self):
        """This function uses the fasta sequence to set the exon start location and length of the exon."""
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

class MotifMark:
    def __init__(self, fasta, motifs):
        '''This class parses the fasta and motif files and uses pycairo to generate the image.'''

        self.fasta = fasta
        self.motifs = motifs
        self.seq_objs = self.get_seq_objs()
        self.motif_objs = self.get_motif_objs()       

        self.surface = cairo.SVGSurface(fasta.split('.')[0] + ".svg", self.max_seq_len() + 25, len(self.seq_objs)*80)
        self.context = cairo.Context(self.surface)
        self.draw_position_x = 15
        self.draw_position_y = 85

    def max_seq_len(self):
        """This function returns the length of the longest Sequence object."""
        max_len = len(self.seq_objs[0].sequence)
        for s in self.seq_objs:
            if max_len < len(s.sequence):
                max_len = len(s.sequence)
        return max_len         

    def get_seq_objs(self):
        """This function parses the fasta file and generates a list of Sequence objects"""
        
        # turn multiline fasta seqs into one-liners 
        ol_output = "OL_" + fasta
        Bioinfo.oneLineFasta(fasta, ol_output)  

        seq_objs = []
        header = ""
        with open(ol_output, "r") as fr:
            for line in fr:
                line = line.strip()
                if line[0] == ">":
                    header = line
                else:
                    seq_objs.append(Sequence(header, line))
        return seq_objs

    def get_motif_objs(self):
        """This function parses the motif file and generates a list of Motif objects"""
        motif_objs = []
        with open(motifs, "r") as fr:
            for line in fr:
                line = line.strip()
                motif_objs.append(Motif(line))
        return motif_objs
        
    def draw_header(self, seq):
        """This function draws the header for a Sequence object"""
        self.context.set_source_rgb(0, 0, 0)
        self.context.set_font_size(13)
        self.context.move_to(self.draw_position_x, self.draw_position_y - 20)
        self.context.show_text(seq.header)

    def draw_line(self, seq):
        """This function draws the sequence line for a Sequence object"""
        seq_length = len(seq.sequence)
        line_end_x = seq_length + self.draw_position_x 
        line_end_y = self.draw_position_y 

        self.context.set_source_rgb(0, 0, 0)
        self.context.set_line_width(1)
        self.context.move_to(self.draw_position_x, self.draw_position_y) 
        self.context.line_to(line_end_x, line_end_y)
        self.context.stroke()

    def draw_exon(self, seq):
        """This function draws the exon rectangle for a Sequence objects"""
        rect_start_x = self.draw_position_x + seq.exon_start 
        rect_start_y = self.draw_position_y - 7
        rect_length = seq.exon_length 
        rect_width = 14

        self.context.set_source_rgb(0, 0, 0)
        self.context.rectangle(rect_start_x, rect_start_y, rect_length, rect_width)
        self.context.fill()

    def draw_motifs(self, seq):
        """This function draws out the motifs for a Sequence objects"""
        for c,m in enumerate(self.motif_objs):
            positions = m.find_motif(seq)
            m.set_color(colors[c])
            self.context.set_source_rgba(m.color[0], m.color[1], m.color[2], 0.6)
            for p in positions:     
                self.context.rectangle(self.draw_position_x + p, self.draw_position_y - 10, len(m.motif), 20)
                self.context.fill()
        self.draw_position_y += 60 

    def draw_legend(self):
        """This function draws the legend with labeled motifs"""
        start_leg_pos_x = 20
        start_leg_pos_y = 30
        leg_pos_x = 20
        leg_pos_y = 30
        self.context.move_to(leg_pos_x, leg_pos_y)   
        self.context.set_source_rgb(0, 0, 0)
        self.context.set_font_size(13)
        self.context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
        self.context.show_text("Motifs:")
        leg_pos_x += 50
        for m in self.motif_objs:
            self.context.rectangle(leg_pos_x, leg_pos_y - 11, 15, 15) 
            self.context.set_source_rgba(m.color[0], m.color[1], m.color[2], 0.6)
            self.context.fill()
            self.context.move_to(leg_pos_x + 20, leg_pos_y)   
            self.context.show_text(m.motif)
            leg_pos_x += len(m.motif) * 7 + 30

        # draw box around the legend
        self.context.set_line_width(0.5)
        self.context.set_source_rgb(0, 0, 0)
        self.context.move_to(start_leg_pos_x - 5, start_leg_pos_y - 15) 
        self.context.line_to(leg_pos_x, start_leg_pos_y - 15)
        self.context.stroke()
        self.context.move_to(start_leg_pos_x - 5, start_leg_pos_y + 10) 
        self.context.line_to(leg_pos_x, start_leg_pos_y + 10)
        self.context.stroke()
        self.context.move_to(start_leg_pos_x - 5, start_leg_pos_y - 15) 
        self.context.line_to(start_leg_pos_x - 5, start_leg_pos_y + 10)
        self.context.stroke()
        self.context.move_to(leg_pos_x, start_leg_pos_y - 15) 
        self.context.line_to(leg_pos_x, start_leg_pos_y + 10)
        self.context.stroke()
        

    def generate_image(self):
        """This function iterates through the Sequence objects, draws each of the components (header, line, exon, motifs, legend)
        and ouputs the finished .png image."""  
        for s in self.seq_objs:
            self.draw_header(s)
            self.draw_line(s)
            self.draw_exon(s)
            self.draw_motifs(s)
        self.draw_legend()   
   
        self.surface.write_to_png(fasta.split('.')[0] + '.png')
        self.surface.finish()


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

# create MotifMark object with fasta and motif file as inputs and generate the output image 
mm = MotifMark(fasta, motifs)
mm.generate_image()
