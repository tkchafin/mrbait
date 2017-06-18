#!/usr/bin/python

import sys
import sqlite3
import getopt
import Bio
import re
from Bio import AlignIO
from mrbait_menu import display_help
from mrbait_menu import parseArgs
import manage_bait_db as m
import alignment_tools as a
import sequence_tools as s
import pandas as pd
import numpy as np

############################### MAIN ###################################

#Parse Command line arguments
params = parseArgs()

#Intiate database connection
conn = m.create_connection(params.db)
c = conn.cursor()

#Initialize empty databases
#if conn.empty() or something like that 
m.init_new_db(conn)

#Parse MAF file and create database
for aln in AlignIO.parse(params.alignment, "maf"):
	cov = len(aln)
	alen = aln.get_alignment_length()
	
	#Add each locus to database
	locus = a.consensAlign(aln, threshold=params.thresh)
	#consensus = str(a.make_consensus(aln, threshold=params.thresh)) #Old way
	locid = m.add_locus_record(conn, cov, locus.conSequence, 0)
	
	#Extract variable positions for database
	for var in locus.alnVars:
		m.add_variant_record(conn, locid, var.name, var.position, var.value)

#First-pass bait design on loci passing pre-filters
#Pre-filters: Length, alignment depth 
c.execute("UPDATE loci SET pass=1 WHERE length < %s OR depth < %s"""%(params.minlen,params.cov))
passedLoci = pd.read_sql_query("""SELECT consensus FROM loci WHERE pass=0""", conn)
#returns pandas dataframe


print("Gs allowed: ",params.numG)
print("Ns allowed: ", params.numN)
#Target region discovery according to params set 
for seq in passedLoci.itertuples():
	start = 0
	stop = 0
	print("\nConsensus: ", seq[1], "\n")
	generator = s.slidingWindowGenerator(seq[1], params.win_shift, params.win_width)
	for window_seq in generator():
		#print()
		seq_temp = re.sub('[ACGT]', '', (window_seq[0]).upper())
		seq_norm = seq_temp.translate(str.maketrans("RYSWKMBDHV", "**********"))
		#print(window_seq[0], ": ", seq_norm)
		counts = s.seqCounterSimple(seq_norm)
		
		#If window passes filters, extend current bait region
		#print("Start is ", start, " and stop is ",stop) #debug print
		if counts['*'] <= params.var_max and counts['N'] <= params.numN and counts['-'] <= params.numG:
			stop = window_seq[2]	
		else:
			#If window fails, check if previous bait region passes to submit to DB
			#print (stop-start)
			if (stop - start) > params.blen:
				target = (seq[1])[start:stop]
				print("	Target region: ", target)
				#Submit target region to database
				#If target region selected, set start of next window to end of current TR
				generator.setI(stop)
				
			#If bait fails, set start to start point of next window
			start = generator.getI()+1
				#print("	--Bait failed.", end='')
			#Current window fails, update start to 
			#print("Window fails")
print()

			



#NOTE: parallelize bait discovery in future!!


#Next:
#	Find all possible bait regions: Contiguous bases
#c.execute("SELECT * FROM loci")
#print (pd.read_sql_query("SELECT * FROM loci", conn))
#print (pd.read_sql_query("SELECT * FROM variants", conn))
#print (pd.read_sql_query("SELECT * FROM samples", conn))			
conn.commit()
conn.close()








