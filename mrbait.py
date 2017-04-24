#!/usr/bin/python

import sys
import sqlite3
import getopt
import Bio
from Bio import AlignIO
from Bio.AlignIO import MafIO
from mrbait_menu import display_help
from mrbait_menu import parseArgs

############################### MAIN ###################################

#Parse Command line arguments
params = parseArgs()
#print("Cov is",params.cov," and maf is ",params.maf)

for aln in AlignIO.parse(params.maf, "maf"):
	if len(aln) < params.cov:
		print("Alignment only has coverage of ", len(aln), ", skipping")
		continue
	elif aln.get_alignment_length() < params.minlen or aln.get_alignment_length() < params.bmin:
		print("Alignment only has length of ", aln.get_alignment_length(), ", skipping")
		continue
	#for seq in aln:
		#print("starts at %s on the %s strand of seq %s in length, and runs for %s bp" % \
			#(seq.annotations["start"],
			#seq.annotations["strand"],
			#seq.annotations["srcSize"],
			#seq.annotations["size"]))
