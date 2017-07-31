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

def loadMAF(conn, params):
	#Parse MAF file and create database
	for aln in AlignIO.parse(params.alignment, "maf"):
		#NOTE: Add error handling, return error code
		cov = len(aln)
		alen = aln.get_alignment_length()
		
		#Add each locus to database
		locus = a.consensAlign(aln, threshold=params.thresh)
		#consensus = str(a.make_consensus(aln, threshold=params.thresh)) #Old way
		locid = m.add_locus_record(conn, cov, locus.conSequence, 0)
		
		#Extract variable positions for database
		for var in locus.alnVars:
			m.add_variant_record(conn, locid, var.name, var.position, var.value)


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

#load alignment to database 
if params.alignment:
	print("Loading MAF file:",params.alignment)
	loadMAF(conn, params)
else:
	#Option to load .loci alignment goes here!
	print("No alignment input found. .loci and .fasta support not added yet!")


#First-pass bait design on loci passing pre-filters
#Pre-filters: Length, alignment depth 
c.execute("UPDATE loci SET pass=1 WHERE length < %s OR depth < %s"""%(params.minlen,params.cov))
passedLoci = pd.read_sql_query("""SELECT id, consensus FROM loci WHERE pass=0""", conn) #returns pandas dataframe


#Target region discovery according to params set 
#looping through passedLoci only
for seq in passedLoci.itertuples():
	start = 0
	stop = 0
	print("\nConsensus: ", seq[2], "ID is: ", seq[1], "\n")
	generator = s.slidingWindowGenerator(seq[2], params.win_shift, params.win_width)
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
				target = (seq[2])[start:stop]
				tr_counts = s.seqCounterSimple(s.simplifySeq(target))
				#Check that there aren't too many SNPs
				#if tr_counts["*"] <= params.vmax_r:
				print("	Target region: ", target)
				#Submit target region to database
				m.add_region_record(conn, int(seq[1]), start, stop, target, tr_counts)	
				#set start of next window to end of current TR
				generator.setI(stop)
				
			#If bait fails, set start to start point of next window
			start = generator.getI()+params.win_shift
print()

#Filter target regions 
#If multiple regions NOT allowed, need to choose which to keep
if params.mult_reg == 0:	
	print("Multiple regions NOT allowed")	
	#Apply --select_r filters 
#Either way, need to apply --filter_r filters
for option in params.filter_r_objects: 
	print("Select Region Option: ", option.o1)
	rand = 0 #false
	rand_num = 0
	if option.o1 is "r":
		rand = 1
		rand_num = option.o2
	elif option.o1 is "g": 
		c.execute("UPDATE regions SET pass=1 WHERE gap > %s"%int(option.o2))
	elif option.o1 is "n":
		c.execute("UPDATE regions SET pass=1 WHERE bad > %s"%int(option.o2))
	elif option.o1 is "m": 
		m.regionFilterMinVar(conn, option.o3, option.o2)
	elif option.o1 is "M":
		m.regionFilterMaxVar(conn, option.o3, option.o2)
	else: 
		assert False, "Unhandled option %r"%option 
		

#Pre-filters: Length, alignment depth 
#c.execute("UPDATE loci SET pass=1 WHERE length < %s OR depth < %s"""%(params.minlen,params.cov))
#passedLoci = pd.read_sql_query("""SELECT consensus FROM loci WHERE pass=0""", conn)

#NOTE: parallelize bait discovery in future!!
#NOTE: Add option for end-trimming off of baits/regions to remove Ns, gaps, etc

#Next:
#	Find all possible bait regions: Contiguous bases
#c.execute("SELECT * FROM loci")
print (pd.read_sql_query("SELECT * FROM regions", conn))
#print (pd.read_sql_query("SELECT * FROM loci", conn))
#print (pd.read_sql_query("SELECT * FROM variants", conn))
#print (pd.read_sql_query("SELECT * FROM samples", conn))			
conn.commit()
conn.close()








