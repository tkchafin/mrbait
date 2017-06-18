#!/usr/bin/python

import sys
import sqlite3
import getopt
import Bio
from Bio import AlignIO
from mrbait_menu import display_help
from mrbait_menu import parseArgs
import manage_bait_db as m
import alignment_tools as a
import sequence_tools as s
import pandas as pd
import numpy as np

############################# FUNCTIONS ################################





############################### MAIN ###################################

#BELOW IS WORKFLOW FOR UCE DESIGN, FINISH AND THEN CONVERT TO FUNCTIONS

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


#Target region discovery according to params set 
#looping through passedLoci only
for seq in passedLoci.itertuples():
	start = 0
	stop = 0
	print("\nConsensus: ", seq[1], "\n")
	generator = s.slidingWindowGenerator(seq[1], params.win_shift, params.win_width)
	for window_seq in generator():

		seq_norm = s.simplifySeq(window_seq[0])
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
				tr_counts = s.seqCounterSimple(s.simplifySeq(target))
				#Check that there aren't too many SNPs
				if tr_counts["*"] <= params.vmax_r:
					print("	Target region: ", target)
					#Submit target region to database
					
				#set start of next window to end of current TR
				generator.setI(stop)
				
			#If bait fails, set start to start point of next window
			start = generator.getI()+params.win_shift
				#print("	--Bait failed.", end='')
			#Current window fails, update start to 
			#print("Window fails")
print()

#Filter target regions 
#If multiple regions NOT allowed, need to choose which to keep
if params.mult_reg == 0:	
	print("Multiple regions allowed")	
	#Apply --select_r filters 
#Either way, need to apply --filter_r filters


#NOTE: parallelize bait discovery in future!!


#Next:
#	Find all possible bait regions: Contiguous bases
#c.execute("SELECT * FROM loci")
#print (pd.read_sql_query("SELECT * FROM loci", conn))
#print (pd.read_sql_query("SELECT * FROM variants", conn))
#print (pd.read_sql_query("SELECT * FROM samples", conn))			
conn.commit()
conn.close()








