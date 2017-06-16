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
print(passedLoci)
	#print(cons)

#NOTE: parallelize bait discovery in future!!



#Next:
#	Find all possible bait regions: Contiguous bases
#c.execute("SELECT * FROM loci")
print (pd.read_sql_query("SELECT * FROM loci", conn))
print (pd.read_sql_query("SELECT * FROM variants", conn))
print (pd.read_sql_query("SELECT * FROM samples", conn))			
conn.commit()
conn.close()








