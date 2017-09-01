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

#Function to load a MAF file into database
def loadMAF(conn, params):
	#Parse MAF file and create database
	for aln in AlignIO.parse(params.alignment, "maf"):
		#NOTE: Add error handling, return error code
		cov = len(aln)
		alen = aln.get_alignment_length()

		#Add each locus to database
		locus = a.consensAlign(aln, threshold=params.thresh)
		#consensus = str(a.make_consensus(aln, threshold=params.thresh)) #Old way
		locid = m.add_locus_record(conn, cov, locus.conSequence, 1)

		print("Loading Locus #:",locid)

		#Extract variable positions for database
		for var in locus.alnVars:
			m.add_variant_record(conn, locid, var.position, var.value)

#Function to load .loci file into database.
def loadLOCI(conn, params):
	#Parse LOCI file and create database
	for aln in read_loci(params.loci):
		#NOTE: Add error handling, return error code
		cov = len(aln)
		alen = aln.get_alignment_length()

		#Add each locus to database
		locus = a.consensAlign(aln, threshold=params.thresh)
		#consensus = str(a.make_consensus(aln, threshold=params.thresh)) #Old way
		locid = m.add_locus_record(conn, cov, locus.conSequence, 1)

		print("Loading Locus #:",locid)

		#Extract variable positions for database
		for var in locus.alnVars:
			m.add_variant_record(conn, locid, var.position, var.value)

#Generator function by ZDZ to parse a .loci file
def read_loci(infile):
	# make emptyp dictionary
	loci = Bio.Align.MultipleSeqAlignment([])

	# read file from command line
	try:
		f = open(infile)
	except IOError as err:
		print("I/O error({0}): {1}".format(err.errno, err.strerror))
	except:
		print("Unexpected error:", sys.exec_info()[0])

	with f as file_object:

		for line in file_object:

			if line[0] == ">":
				identifier = line.split()[0]
				sequence = line.split()[1]
				loci.add_sequence(identifier, sequence)

			else:
				yield(loci)
				loci = Bio.Align.MultipleSeqAlignment([])


#Function to filter target regions by --filter_R arguments
def filterTargetRegions(conn, params):
	rand = 0 #false
	rand_num = 0
	for option in params.filter_r_objects:
		print("Select Region Option: ", option.o1)

		if option.o1 == "r":
			#Set 'rand' to TRUE for random selection AFTER other filters
			rand_num = int(option.o2)
			assert rand_num > 0, "Number for random TR selection must be greater than zero!"
		elif option.o1 == "g":
			c.execute("UPDATE regions SET pass=0 WHERE gap > %s"%int(option.o2))
		elif option.o1 == "n":
			c.execute("UPDATE regions SET pass=0 WHERE bad > %s"%int(option.o2))
		elif option.o1 == "m":
			m.regionFilterMinVar(conn, option.o2, option.o3)
		elif option.o1 == "M":
			m.regionFilterMaxVar(conn, option.o2, option.o3)
			#m.printVarCounts(conn, option.o3)
		else:
			assert False, "Unhandled option %r"%option

	#If 'random' select is turned on, then apply AFTER resolving conflicts
	if rand_num:
		return(rand_num)


#Function to filter target regions by --filter_R arguments
def selectTargetRegions(conn, params):
	#For all alignments over --min_mult:
		#If TRs are within --dist
	print("Select TR criterion is: ",params.select_r)
	print("Select dist is: ", params.select_r_dist)
	print("Minimum mult_reg dist is: ",params.min_mult)

	#Build conflict tables
	#print(pd.read_sql_query("SELECT * FROM regions", conn))
	if params.mult_reg == 0:
		print("Multiple regions NOT allowed, apply --select_r within whole loci")
		#TODO: Need function call to buid conflict_blocks by whole loci
		m.fetchConflictTRs_NoMult(conn)
	else:
		print("Multiple TRs allowed, apply --select_r within conflict_blocks <--dist_r>")
		#Build conflict_blocks according to --dist_r and --min_mult parameters
		m.fetchConflictTRs(conn, params.min_mult, params.dist_r)

	#NEXT: Need to select TRs within conflict_blocks
	#Apply select_r filters for all conflicting TRs
	if params.select_r == "r":
		#Do it later
		pass
	elif params.select_r == "s":
		#Select based on SNPs flanking in "d" dist
		print("--select_r is SNP, dist is ",params.select_r_dist)
		try:
			m.regionSelect_SNP(conn,params.select_r_dist)
		except ValueError as err:
			sys.exit(err.args)
	#	except:
		#	sys.exit(sys.exc_info()[0])
	elif params.select_r == "m":
		#Select based on minimizing SNPs, Ns and gaps in TR region, otherwise randomly
		print("--select_r is MINVAR_TR")
	elif params.select_r == "b":
		#Select based on least Ns and gaps in "d" flanking bases
		print("--select_r is MINBAD, dist is ", params.select_r_dist)
	elif params.select_r == "c":
		#Select based on minimizing SNPs in flanking region
		print("select_r is MINVAR_FLANK, dist is ", params.select_r_dist)
	else:
		assert False, "Unhandled option %r"%params.select_r

	#randomly resolve any remaining conflicts
	print("--select_r is RANDOM")
	try:
		m.regionSelectRandom(conn)
	except ValueError as err:
		sys.exit(err.args)
	except:
			sys.exit(sys.exc_info()[0])
	print(pd.read_sql_query("SELECT * FROM conflicts", conn))
	#NEXT: Randomly resolve any remaining conflicts
	#NEXT: Push conflicts to change "pass" attribute in regions table

#Function to check that target regions table is valid to continue
def checkTargetRegions(conn):
	#Fetch number of entries to TR table
	total = m.getNumTRs(conn)
	if (int(total) <= 0):
		sys.exit("Program killed: No Target Regions were found.")
	passed = m.getNumPassedTRs(conn)
	if (int(passed) <= 0):
		sys.exit("Program killed: No Target Regions passed selection/filtering.")
	print("Checking target regions")

############################### MAIN ###################################

#BELOW IS WORKFLOW FOR UCE DESIGN, FINISH AND THEN CONVERT TO FUNCTIONS
#ADD GFF FUNCTIONALITY LATER
#Add multithreading support later... Each thread will need its own db conn
#If TR too short, can add option to include variable flanking regions?
#TODO: Option for first and second pass over database (e.g. first conservative, second of only failed TRs??)
#TODO: Add better checking to make sure database isn't empty before proceeding (e.g. if filters too stringent)
#TODO: Add "flow control" options, e.g. only make db, load previous db, only TR, etc

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
elif params.loci:
	print("Loading LOCI file:",params.loci)
	loadLOCI(conn, params)
else:
	#Option to load .loci alignment goes here!
	print("No alignment input found. .fasta, .gff, and .phylip support not added yet!")


#First-pass bait design on loci passing pre-filters
#PASS=1 is PASS=FALSE
#Pre-filters: Length, alignment depth
m.filterLoci(conn, params.minlen, params.cov)
passedLoci = m.getPassedLoci(conn) #returns pandas dataframe
if passedLoci.shape[0] <= 0:
	sys.exit("Program killed: No loci passed filtering.")

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

#Assert that there are TRs chosen, and that not all have been filtered out
checkTargetRegions(conn)

#Filter target regions
#If multiple regions NOT allowed, need to choose which to keep
print("Starting: Target Region Selection...")

#First pass filtering of target regions
rand = filterTargetRegions(conn, params)

#Apply --select_r filters to sort any conflicting TRs
selectTargetRegions(conn, params)

#If RANDOM filter for TRs was applied, apply if AFTER TR conflict resolution
if rand:
	print("randomly filtering all TRs")
	m.regionFilterRandom(conn, rand)

#Pre-filters: Length, alignment depth
#c.execute("UPDATE loci SET pass=1 WHERE length < %s OR depth < %s"""%(params.minlen,params.cov))
#passedLoci = pd.read_sql_query("""SELECT consensus FROM loci WHERE pass=0""", conn)

#NOTE: parallelize bait discovery in future!!
#NOTE: Add option for end-trimming off of baits/regions to remove Ns, gaps, etc

print("\n\nProgram ending...Here are some results\n\n")
#Next:
#	Find all possible bait regions: Contiguous bases
#c.execute("SELECT * FROM loci")
#print (pd.read_sql_query("SELECT * FROM loci", conn))
print (pd.read_sql_query("SELECT * FROM regions", conn))
#print (pd.read_sql_query("SELECT * FROM variants", conn))
conn.commit()
conn.close()
