#!/usr/bin/python

import sys
import sqlite3
import getopt
import Bio
from Bio import AlignIO
from mrbait_menu import display_help
from mrbait_menu import parseArgs
from substring import SubString
import manage_bait_db as m
import alignment_tools as a
import sequence_tools as s
import misc_utils as utils
import pandas as pd
import numpy as np
import emboss_align as emb


############################# FUNCTIONS ################################

#Function to load a MAF file into database
def loadMAF(conn, params):
	#Parse MAF file and create database
	for aln in AlignIO.parse(params.alignment, "maf"):
		#NOTE: Add error handling, return error code
		cov = len(aln)
		alen = aln.get_alignment_length()

		#Add each locus to database
		locus = a.consensAlign(aln, threshold=params.thresh, mask=params.mask)
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
		locus = a.consensAlign(aln, threshold=params.thresh, mask=params.mask)
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

#Function to discover target regions using a sliding windows through passedLoci
def targetDiscoverySlidingWindow(conn, params, loci):

	#looping through passedLoci only
	for seq in loci.itertuples():
		#print(seq)
		start = 0
		stop = 0
		#print("\nConsensus: ", seq[2], "ID is: ", seq[1], "\n")
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
					n_mask = utils.n_lower_chars(target)
					n_gc = s.gc_content(target)
					#Check that there aren't too many SNPs
					#if tr_counts["*"] <= params.vmax_r:
					#print("	Target region: ", target)
					#Submit target region to database
					m.add_region_record(conn, int(seq[1]), start, stop, target, tr_counts, n_mask, n_gc)
					#set start of next window to end of current TR
					generator.setI(stop)

				#If bait fails, set start to start point of next window
				start = generator.getI()+params.win_shift

#Function to filter target regions by --filter_R arguments
def filterTargetRegions(conn, params):
	rand = 0 #false
	rand_num = 0

	#Filter by --vmax_r
	m.varMaxFilterTR(conn, params.vmax_r)

	for option in params.filter_r_objects:
		print("Select Region Option: ", option.o1)
		if option.o1 == "rand":
			#Set 'rand' to TRUE for random selection AFTER other filters
			rand_num = int(option.o2)
			assert rand_num > 0, "Number for random TR selection must be greater than zero!"
		elif option.o1 == "gap":
			c.execute("UPDATE regions SET pass=0 WHERE gap > %s"%int(option.o2))
		elif option.o1 == "bad":
			c.execute("UPDATE regions SET pass=0 WHERE bad > %s"%int(option.o2))
		elif option.o1 == "min":
			m.regionFilterMinVar(conn, val=option.o2, flank=option.o3)
		elif option.o1 == "max":
			m.regionFilterMaxVar(conn, val=option.o2, flank=option.o3)
			#m.printVarCounts(conn, option.o3)
		elif option.o1 == "mask":
			min_mask_prop = option.o2
			max_mask_prop = option.o3
			m.regionFilterMask(conn, minprop=min_mask_prop, maxprop=max_mask_prop)
		elif option.o1 == "gc":
			min_mask_prop = option.o2
			max_mask_prop = option.o3
			print("Filter regions by GC, but option not yet implemented")
			m.regionFilterGC(conn, minprop=min_mask_prop, maxprop=max_mask_prop)
		elif option.o1 == "len":
			minlen = option.o2
			maxlen= option.o3
			assert minlen < maxlen, "<--filter_r> suboption \"len\": Min must be less than max"
			m.lengthFilterTR(conn, maxlen, minlen)
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
	#TODO: First fetch how many conflicts, if there are none, then EXIT FUNCTION

	#Apply select_r filters for all conflicting TRs
	if params.select_r == "rand":
		print("--select_r is RANDOM")
		#Do it later
		pass
	elif params.select_r == "snp":
		#Select based on SNPs flanking in "d" dist
		print("--select_r is SNP, dist is ",params.select_r_dist)
		try:
			m.regionSelect_SNP(conn,params.select_r_dist)
		except ValueError as err:
			sys.exit(err.args)
	#	except:
		#	sys.exit(sys.exc_info()[0])
	elif params.select_r == "conw":
		#Select based on minimizing SNPs, Ns and gaps in TR region, otherwise randomly
		print("--select_r is MINVAR_TR")
		try:
			m.regionSelect_MINVAR_TR(conn)
		except ValueError as err:
			sys.exit(err.args)
	elif params.select_r == "bad":
		#Select based on least Ns and gaps in "d" flanking bases
		print("--select_r is MINBAD, dist is ", params.select_r_dist)
		try:
			m.regionSelect_MINBAD(conn,params.select_r_dist)
		except ValueError as err:
			sys.exit(err.args)
	elif params.select_r == "conf":
		#Select based on minimizing SNPs in flanking region
		print("select_r is MINVAR_FLANK, dist is ", params.select_r_dist)
		try:
			m.regionSelect_MINSNP(conn,params.select_r_dist)
		except ValueError as err:
			sys.exit(err.args)
	else:
		assert False, "Unhandled option %r"%params.select_r

	#randomly resolve any remaining conflicts
	try:
		m.regionSelectRandom(conn)
	except ValueError as err:
		sys.exit(err.args)
	except:
			sys.exit(sys.exc_info()[0])
	print(pd.read_sql_query("SELECT * FROM conflicts", conn))

	#NEXT: Push conflicts to change "pass" attribute in regions table
	m.pushResolvedConflicts(conn)

#Function to check that target regions table is valid to continue
def checkTargetRegions(conn):
	#Fetch number of entries to TR table
	total = m.getNumTRs(conn)
	if (int(total) <= 0):
		sys.exit("Program killed: No Target Regions were found.")
	passed = m.getNumPassedTRs(conn)
	if (int(passed) <= 0):
		sys.exit("Program killed: No Target Regions passed selection/filtering.")

#Function for deduplication of targets by pairwise alignment
def pairwiseAlignDedup(params, seqs):
	for a in seqs.itertuples():
		print(a[2])



#function for sliding window bait generation
def baitSlidingWindow(conn, source, sequence, overlap, length):
	generator = s.slidingWindowGenerator(sequence, overlap, length)
	for window_seq in generator():
		#Don't need to do a bunch of filtering, because all was checked when TRs built
		#print(window_seq)
		if (len(window_seq[0]) == length):
			m.add_bait_record(conn, source, window_seq[0], window_seq[1], window_seq[2])

#function for sliding window bait generation, with custom coordinates
def baitSlidingWindowCoord(conn, source, sequence, overlap, length, start):
	generator = s.slidingWindowGenerator(sequence, overlap, length)
	for window_seq in generator():
		#Don't need to do a bunch of filtering, because all was checked when TRs built
		#print(window_seq)
		if (len(window_seq[0]) == length):
			start_coord = start + window_seq[1]
			stop_coord = start_coord + length
			m.add_bait_record(conn, source, window_seq[0], start_coord, stop_coord)

#Function to discover target regions
def baitDiscovery(conn, params, targets):
	print("Params.overlap is ", params.overlap)
	print("Params.bait_shift is", params.bait_shift)
	#Design baits based on specified selection criterion (default is to tile at 2X)
	if params.select_b == "tile":
		#looping through passedLoci only
		for seq in targets.itertuples():
			#seq[1] is the regid; seq[2] is the target sequence
			baitSlidingWindow(conn, seq[1], seq[2], params.bait_shift, params.blen)
	elif params.select_b == 'center':
		#print("Designing centered baits...")
		#First calculate union length needed, if this is longer than target, just
		#tile all of it
		union = utils.calculateUnionLengthFixed(params.select_b_num, params.blen, params.overlap)
		print("Union length is",union)
		#looping through passedLoci only
		for seq in targets.itertuples():
			length = len(seq[2])
			#If the target is too short, just do a full sliding window
			if union >= length:
				baitSlidingWindow(conn, seq[1], seq[2], params.bait_shift, params.blen)
			else:
				center = len(seq[2]) // 2 #Divide by two and round down
				start = center - (union // 2)
				stop = start+union
				#print("Starting at:",start," and stopping at:",stop)
				subseq = (seq[2])[start:stop]
				#print(subseq)
				baitSlidingWindowCoord(conn, seq[1], subseq, params.bait_shift, params.blen, start)
	elif params.select_b == "flank":
		#First calculate union length needed, if this is longer than target, just
		#tile all of it
		#union x 2 because we do this on both sides
		union = (utils.calculateUnionLengthFixed(params.select_b_num, params.blen, params.overlap))
		for seq in targets.itertuples():
			length = len(seq[2])
			#If the target is too short, just do a full sliding window
			if union*2 >= length:
				baitSlidingWindow(conn, seq[1], seq[2], params.bait_shift, params.blen)
			else:
				#Need to: Substring both ends (start + union and stop - union)
				subseq1 = (seq[2])[0:union] #right
				subseq2 = (seq[2])[length-union:] #left
				#print(subseq1)
				#print(subseq2)
				#Right side
				baitSlidingWindowCoord(conn, seq[1], subseq1, params.bait_shift, params.blen, 0)
				#Left side
				baitSlidingWindowCoord(conn, seq[1], subseq2, params.bait_shift, params.blen, length-union)
	# elif params.select_b == "rand":
	# 	#Here, union is the MINIMUM length required to make the specified number of baits with maximum overlap
	# 	union = (utils.calculateUnionLengthFixed(params.select_b_num, params.blen, params.overlap))
	# 	union_noOverlap = (params.select_b_num * params.blen)
	# 	for seq in targets.itertuples():
	# 		length = len(seq[2])
	# 		print("Union is ",union, ", and seq length is ", length)
	# 		#If the target is too short, just do a full sliding window
	# 		if union >= length:
	# 			print("Too short, make all.")
	# 			baitSlidingWindow(conn, seq[1], seq[2], params.overlap, params.blen)
	# 		else:
	# 			pass
	# 			#TODO: Make substring class, with method to search list of already gathered substring to check for overlap
	# 			#While overlap is to high, generate and add/check another random substring until X passing are gathered.
	# 			#Loop through that list and commit all to table
	# 			subseqs = []
	# 			while len(subseqs) < 3:
	# 				new = SubString()
	# 				new.randomDrawSubstring(seq[2], params.blen)
	# 				#Check if new substring overlaps with already picked substrings
	# 				if new.checkMatch(subseqs, params.overlap) == 0:
	# 					subseqs.append(new)
	# 					print("Length of subseqs is now: ", len(subseqs))
	# 					print("Desired length is ", params.select_b_num)
	# 			#if (len(subseqs) > 1):
	# 				#SubString.sortSubStrings(subseqs)
	# 			for s in subseqs:
	# 				print("From regid:",seq[1], " -- ",s.string, ":", s.start, s.stop)
				#new.printAll()

			#	randomDrawSubstring
	else:
		assert False, "Unhandled option %r"%params.select_b

	#At end, make sure to filter out any baits which are too short

############################### MAIN ###################################

#BELOW IS WORKFLOW FOR UCE DESIGN, FINISH AND THEN CONVERT TO FUNCTIONS
#ADD GFF FUNCTIONALITY LATER
#Add multithreading support later... Each thread will need its own db conn
#If TR too short, can add option to include variable flanking regions?
#TODO: Option for first and second pass over database (e.g. first conservative, second of only failed TRs??)
#TODO: Add better checking to make sure database isn't empty before proceeding (e.g. if filters too stringent)
#TODO: Add "flow control" options, e.g. only make db, load previous db, only TR, etc
#TODO: Add BRANCHING option for Flow Control. e.g. Branch RUN1 from Regions level to try new method of bait design, etc
#TODO: For whole genome option, need to read an mpileup or similar to capture variant information.
#TODO: Could also set minimum coverage thresholds for whole genome?
#TODO: Some form of duplicate screening. Screen targets for dupe or screen baits? Not sure.
#TODO: Way to filter targets by masking in flanking region???
#TODO: Filter targets by distance to GFF element
#TODO: Option to constrain targets to ONLY in GFF elements
#TODO: Parallelize loadLOCI and loadMAF
#TODO: Parallelize TR discovery somehow??
#TODO: FIlter by flanking masked, and flanking GFF elements
#TODO: Dedup Targets using either smith-waterman or needleman-wunsch (from Bio.EMBOSS)
#TODO: Could also have Target filtering by BLAST to a local database (provide as FASTA)???

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
targetDiscoverySlidingWindow(conn, params, passedLoci)
print()

#Assert that there are TRs chosen, and that not all have been filtered out
checkTargetRegions(conn)

#Filter target regions
#If multiple regions NOT allowed, need to choose which to keep
print("Starting: Target Region Selection...")

#First pass filtering of target regions
rand = filterTargetRegions(conn, params)

#Check again that not all have been filtered out
checkTargetRegions(conn)

#Apply --select_r filters to sort any conflicting TRs
selectTargetRegions(conn, params)

#If RANDOM filter for TRs was applied, apply if AFTER TR conflict resolution
if rand:
	print("randomly filtering all TRs")
	m.regionFilterRandom(conn, rand)

#Target region deduplication by pairwise alignment
passedTargets = m.getPassedTRs(conn)
pairwiseAlignDedup(params, passedTargets)
sys.exit()
#Bait discovery
print("Starting probe design...")

#Check again that not all have been filtered out
checkTargetRegions(conn)

#Sliding window bait discovery
passedTargets = m.getPassedTRs(conn)
if passedTargets.shape[0] <= 0:
	sys.exit("Program killed: No targets passed filtering.")
baitDiscovery(conn, params, passedTargets)

print("\n\nProgram ending...Here are some results\n\n")

print(m.getLoci(conn))
print(m.getRegions(conn))
print(m.getBaits(conn))

conn.commit()
conn.close()
