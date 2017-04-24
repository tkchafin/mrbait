#!/usr/bin/python

import sys
from Bio import AlignIO
from Bio.Alphabet import IUPAC

#Less shitty consensus function than BioPython has..
def make_consensus(alignment="", threshold=0.1):
	print("here you go")
	aln_len = alignment.get_alignment_length()
	
	for i in range(aln_len):
		nuc_counts = {} #dictionary to track occurences
		nuc_types = 0 #track number of things we found
		print(i)
		
