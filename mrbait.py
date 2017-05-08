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
import alignment_tools as a
import pandas as p

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
conn = m.create_connection(params.db)
c = conn.cursor()

#Initialize empty databases
#if conn.empty() or something like that 
m.init_new_db(conn)

#Parse MAF file and create database
for aln in AlignIO.parse(params.maf, "maf"):
	cov = len(aln)
	alen = aln.get_alignment_length()
	#Skip if too few individuals in alignment
	#if cov < params.cov:
		#print("Alignment only has coverage of ", cov, ", skipping")
		#continue
	#Skip if alignment length too low
	#elif alen < params.minlen or alen < params.blen:
		#print("Alignment only has length of ", alen, ", skipping")
		#continue
	
	#Add each locus to database
	locus = a.consensAlign(aln, threshold=params.thresh)
	#consensus = str(a.make_consensus(aln, threshold=params.thresh)) #Old way
	locid = m.add_locus_record(conn, cov, locus.conSequence, 0)
	
	#Extract variable positions for database
	for var in locus.alnVars:
		m.add_variant_record(conn, locid, var.name, var.position, var.value)

#First-pass bait design on loci passing pre-filters
	#Pre-filters: Length, alignment depth 
	
#c.execute("SELECT * FROM loci")
print (p.read_sql_query("SELECT * FROM loci", conn))
print (p.read_sql_query("SELECT * FROM positions", conn))
print (p.read_sql_query("SELECT * FROM variants", conn))

#print (c.fetchall())
			
conn.commit()
conn.close()
