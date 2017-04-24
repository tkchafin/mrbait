#!/usr/bin/python

import sys
import sqlite3
import getopt
import Bio
from Bio import AlignIO
from Bio.AlignIO import MafIO
from mrbait_menu import display_help
from mrbait_menu import parseArgs
import manage_bait_db as m

############################### MAIN ###################################


#  To Do:
#    -Input .loci
#    -Input VCF?
#    -Input raw read map
#       --Could then filter by coverage
#


#Parse Command line arguments
params = parseArgs()

#Intiate database connection
sqlite_file = params.db
conn = sqlite3.connect(sqlite_file)
c = conn.cursor()

#Initialize empty databases
#if conn.empty() or something like that 
m.init_new_db(c, conn)

#Parse MAF file
for aln in AlignIO.parse(params.maf, "maf"):
	#Skip if too few individuals in alignment
	if len(aln) < params.cov:
		print("Alignment only has coverage of ", len(aln), ", skipping")
		continue
	#Skip if alignment length too low
	elif aln.get_alignment_length() < params.minlen or aln.get_alignment_length() < params.bmin:
		print("Alignment only has length of ", aln.get_alignment_length(), ", skipping")
		continue
	
	#Add each locus to database
	
	#for
	#for seq in aln:
		#print("starts at %s on the %s strand of seq %s in length, and runs for %s bp" % \
			#(seq.annotations["start"],
			#seq.annotations["strand"],
			#seq.annotations["srcSize"],
			#seq.annotations["size"]))
			
conn.commit()
conn.close()
