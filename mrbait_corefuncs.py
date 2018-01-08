import sys
import sqlite3
import getopt
import Bio
import os
from Bio import AlignIO
from mrbait_menu import display_help
from mrbait_menu import parseArgs
from substring import SubString
import manage_bait_db as m
import alignment_tools as a
import sequence_tools as s
import misc_utils as utils
import seq_graph as graph
import aln_file_tools
import pandas as pd
import numpy as np
import blast as b
import vsearch
import subprocess
import vcf_tools
import gff3_parser as gff

############################# FUNCTIONS ################################

#Function to load a MAF file into database
def loadMAF(conn, params):
	#Parse MAF file and create database
	num = 1
	for aln in AlignIO.parse(params.alignment, "maf"):
		#NOTE: Add error handling, return error code
		cov = len(aln)
		alen = aln.get_alignment_length()

		#Add each locus to database
		locus = a.consensAlign(aln, threshold=params.thresh, mask=params.mask)
		#consensus = str(a.make_consensus(aln, threshold=params.thresh)) #Old way
		locid = m.add_locus_record(conn, cov, locus.conSequence, 1, num)
		num+=1

		print("Loading Locus #:",locid)

		#Extract variable positions for database
		for var in locus.alnVars:
			m.add_variant_record(conn, locid, var.position, var.value)

#Function to load .loci file into database.
def loadLOCI(conn, params):
	locnum = 1
	#Parse LOCI file and create database
	for aln in aln_file_tools.read_loci(params.loci):
		#NOTE: Add error handling, return error code
		cov = len(aln)
		alen = aln.get_alignment_length()

		#Add each locus to database
		locus = a.consensAlign(aln, threshold=params.thresh, mask=params.mask)
		#consensus = str(a.make_consensus(aln, threshold=params.thresh)) #Old way
		locid = m.add_locus_record(conn, cov, locus.conSequence, 1, locnum)
		locnum += 1
		#print("Loading Locus #:",locid)

		#Extract variable positions for database
		for var in locus.alnVars:
			m.add_variant_record(conn, locid, var.position, var.value)

#Function to load FASTA into database
def loadFASTA(conn, params):
	for contig in aln_file_tools.read_fasta(params.assembly):
		#print("Reading contig:",contig[0])
		#print("Sequence is:",contig[1])
		locid = m.add_locus_record(conn, 1, contig[1], 1, contig[0])

		#Parse consensus for vars, submit those vars to db
		for var in a.get_vars(contig[1]):
			m.add_variant_record(conn, locid, var.position, var.value)

#Function to load GFF file into database
def loadGFF(conn, params):
	#For each GFF record in params.gff
	for record in gff.read_gff(params.gff):
		#Skip any records that are missing the sequence ID, or coordinates
		if record.seqid == "NULL" or record.start == "NULL" or record.end == "NULL":
			continue
		if record.start > record.end:
			temp = record.start
			record.start = record.end
			record.end = temp
		#Get the alias, if it exists
		alias = ""
		if record.getAlias(): #returns false if no alias
			alias = record.getAlias()
		else:
			alias = "NULL"
		#NOTE: This function ONLY inserts GFFRecords where record.seqid matches an existing locus in the loci table
		m.add_gff_record(conn, record.seqid, record.type.lower(), record.start, record.end, alias)

	#Check if all GFF records fall within bounds of
	m.validateGFFRecords(conn)


#Function to load VCF variants file
def loadVCF(conn, params):

	loci = m.getPassedLoci(conn) #get DF of passed loci
	print(loci)
	chrom_lookup = set(loci["chrom"].tolist()) #lookup table of locus IDs
	loci = loci.set_index(['chrom']) #index loci DF by chrom column
	#print(loci)
	passed=0 #To track number of VCF records for which no locus exists
	failed=0

	for reclist in vcf_tools.read_vcf(params.vcf):
		rec_chrom = reclist[0].CHROM
		print("Starting locus",rec_chrom)
		if rec_chrom in chrom_lookup:
			locid = int(loci.loc[rec_chrom]["id"])
			#print("locid is",locid)
			passed+=1
			#for rec in reclist:
			#	print(rec.CHROM, rec.POS, rec.REF, rec.ALT, len(rec.samples), rec.call_rate, rec.aaf)
			#Grab DF record for the matching CHROM
			sub = loci.loc[rec_chrom]
			#Get new consensus sequence given VCF records
			new_cons = vcf_tools.make_consensus_from_vcf(sub['consensus'],reclist, params.thresh)
			#Update new consensus seq in db
			if len(new_cons) != len(sub['consensus']): #Check length first
				print("Warning: New consensus sequence for locus %s (locid=<%s>) is the wrong length! Skipping."%(rec_chrom, locid))
			else:
				m.updateConsensus(conn, locid, new_cons)
				#Delete old vars for locus, and parse new consensus
				m.purgeVars(conn, locid)
				#Get new vars and add records to table
				for var in a.get_vars(new_cons):
					m.add_variant_record(conn, locid, var.position, var.value)

		else:
			#print(rec_chrom, "not found.")
			failed+=1
	if failed > 0:
		print("WARNING:%s/%s records in <%s> referenced sequences not found in <%s> FASTA headers"%(failed, failed+passed, params.vcf, params.assembly))

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
					n_gc = s.gc_counts(target)
					#Check that there aren't too many SNPs
					#if tr_counts["*"] <= params.vmax_r:
					#print("	Target region: ", target)
					#Submit target region to database
					m.add_region_record(conn, int(seq[1]), start, stop, target, tr_counts, n_mask, n_gc)
					#set start of next window to end of current TR
					generator.setI(stop)

				#If bait fails, set start to start point of next window
				start = generator.getI()+params.win_shift
	#Now update regions table to include information for flanking regions if available
	m.flankDistParser(conn, params.flank_dist)

#Function to filter target regions by --filter_R arguments
def filterTargetRegions(conn, params):

	rand = 0 #false
	rand_num = 0
	aln = 0
	blast = 0
	minid = None
	mincov = None

	for option in params.filter_r_objects:
		print("Filter Region Option: ", option.o1)
		if option.o1 == "rand":
			#Set 'rand' to TRUE for random selection AFTER other filters
			rand_num = int(option.o2)
			assert rand_num > 0, "Number for random TR selection must be greater than zero!"
		elif option.o1 == "gap":
			m.simpleFilterTargets_gap(conn, int(option.o2))
		elif option.o1 == "bad":
			m.simpleFilterTargets_bad(conn, int(option.o2))
		elif option.o1 == "snp":
			m.simpleFilterTargets_SNP(conn, int(option.o2), int(option.o3))
		elif option.o1 == "mask":
			max_mask_prop = option.o2
			m.regionFilterMask(conn, maxprop=max_mask_prop)
		elif option.o1 == "gc":
			min_mask_prop = option.o2
			max_mask_prop = option.o3
			m.regionFilterGC(conn, minprop=min_mask_prop, maxprop=max_mask_prop)
		elif option.o1 == "len":
			minlen = option.o2
			maxlen= option.o3
			assert minlen < maxlen, "<--filter_r> suboption \"len\": Min must be less than max"
			m.lengthFilterTR(conn, maxlen, minlen)
		elif option.o1 == "pw":
			aln = 1
			minid = option.o2
			mincov= option.o3
		elif option.o1 in ("blast_x", "blast_i", "blast_a"):
			#if blast database given as fasta, make a blastdb:
			db_path = None
			if (params.blastdb):
				db_path = params.blastdb
			elif (params.fastadb):
				db_path = params.workdir + "/blastdb/" + params.out
				b.makeblastdb(params.makedb, params.fastadb, db_path)
				params.blastdb = db_path
			elif(not params.blastdb and not params.fastadb):
				print("WARNING: No blast database was provided. Skipping <--filter_r> option %s"%option.o1)
				break
			print("BLASTDB PATH IS: ", db_path)
			#Get targets, print to fasta
			seqs = m.getPassedTRs(conn)
			fas = params.workdir + "/.temp.fasta"
			aln_file_tools.writeFasta(seqs, fas)
			outfile = params.workdir + "/.temp.blast"
			if option.o1 == "blast_x":
				blacklist = b.blastExcludeMatch(params, db_path, fas, option.o2, option.o3, outfile)
				m.removeRegionsByList(conn, blacklist)
			elif option.o1 == "blast_i":
				whitelist = b.blastIncludeMatch(params, db_path, fas, option.o2, option.o3, outfile)
				m.removeRegionsByWhitelist(conn, whitelist)
			elif option.o1 == "blast_a":
				sys.exit("blast_a option is not yet implemented.")
				#blacklist = b.blastExcludeAmbig(params, db_path, fas, option.o2, option.o3, outfile)
				#m.removeRegionsByList(conn, blacklist)
			os.remove(fas)
			os.remove(outfile)
			#sys.exit()
		elif option.o1 in ("gff", "gff_a"):
			if params.gff and params.assembly:
				if option.o1 == "gff":
					m.regionFilterGFF(conn, option.o2, params.flank_dist)
				elif option.o1 == "gff_a":
					m.regionFilterGFF_Alias(conn, option.o2, params.flank_dist)
			else:
				sys.exit("ERROR: Filtering targets on proximity to GFF elements requires FASTA <-A> and GFF <-G> inputs!")
		else:
			assert False, "Unhandled option %r"%option

	#Perform pairwise alignment AFTER all other filters because it is analytically more expensive
	#Target region deduplication by pairwise alignment
	if aln:
		passedTargets = m.getPassedTRs(conn)
		minid = option.o2
		mincov= option.o3
		assert (0.0 <= minid <= 1.0), "Minimum ID for pairwise alignment must be between 0.0 and 1.0"
		assert (0.0 <= mincov <= 1.0), "Minimum alignment coverage for pairwise alignment must be between 0.0 and 1.0"
		blacklist_edges = pairwiseAlignDedup(conn, params, passedTargets, minid, mincov)
		if (len(blacklist_edges) > 0):
			revised_blacklist = dupEdgeResolution(conn, params, blacklist_edges)
			if len(revised_blacklist) > 0:
				m.removeRegionsByList(conn, revised_blacklist)

	#If 'random' select is turned on, then apply AFTER resolving conflicts (--select_r)
	if rand_num:
		return(rand_num)

#Function to filter target regions by --filter_R arguments
def selectTargetRegions(conn, params):
	#For all alignments over --min_mult:
		#If TRs are within --dist
	print("Select TR criterion is: ",params.select_r)
	print("Flank dist is: ", params.flank_dist)
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
	if m.getNumConflicts(conn) > 0:
	#Apply select_r filters for all conflicting TRs
		if params.select_r == "rand":
			print("--select_r is RANDOM")
			#Do it later
			pass
		elif params.select_r == "snp":
			#Select based on SNPs flanking in "d" dist
			try:
				#TODO: Change to parse flank first and populate in table
				m.regionSelect_SNP(conn)
			except ValueError as err:
				sys.exit(err.args)
		#	except:
			#	sys.exit(sys.exc_info()[0])
		elif params.select_r == "bad":
			#Select based on least Ns and gaps in "d" flanking bases
			try:
				m.regionSelect_MINBAD(conn)
			except ValueError as err:
				sys.exit(err.args)
		elif params.select_r == "cons":
			#Select based on minimizing SNPs in flanking region
			#TODO: Implement flankDistParser first!!!
			try:
				m.regionSelect_MINSNP(conn)
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
		#print(pd.read_sql_query("SELECT * FROM conflicts", conn))

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
def pairwiseAlignDedup(conn, params, seqs, minid, mincov):
	"""Seqs must be a pandas DF where cols: 0=index, 1=name, 2=sequence"""
	fas = params.workdir + "/.temp_ref.fasta"
	aln_file_tools.writeFasta(seqs, fas)

	#First sort FASTA by size
	sor = params.workdir + "/.temp.sort"
	try:
		vsearch.sortByLength(params.vsearch, fas, sor)
	except KeyboardInterrupt:
		sys.exit("Process aborted: Keyboard Interrupt")
	except subprocess.CalledProcessError as err:
		print("VSEARCH encountered a problem.")
		sys.exit(err.args)
	except NameError as err:
		sys.exit(err.args)
	except OSError as err:
		print("Exception: OSError in VSEARCH call. Check that you are using the correct executable.")
		sys.exit(err.args)
	except:
		sys.exit(sys.exc_info()[0])
	os.remove(fas)

	#Pairwise align sorted FASTA (sorted so the shorter seq is always 'target' amd longer is 'query')
	pw = params.workdir + "/" + params.out + ".pw"
	try:
		vsearch.allpairsGlobal(params.vsearch, params.vthreads, sor, minid, mincov, pw)
	except KeyboardInterrupt:
		sys.exit("Process aborted: Keyboard Interrupt")
	except subprocess.CalledProcessError as err:
		print("VSEARCH encountered a problem.")
		sys.exit(err.args)
	except OSError as err:
		print("Exception: OSError in VSEARCH call. Check that you are using the correct executable.")
		sys.exit(err.args)
	except:
		sys.exit(sys.exc_info()[0])
	os.remove(sor)

	#Finally, parse the output of pairwise alignment, to get 'bad matches'
	#Returns a list of edges
	ret = vsearch.parsePairwiseAlign(pw)
	os.remove(pw)
	return(ret)


#Function to perform conflict resolution based on VSEARCH results
def dupEdgeResolution(conn, params, blacklist_edges):
	#OUTLINE:
	#--Get blacklist
	#--Query table to get weights (number of vars)
	#	If --noGraph:
	#		Remove all blacklist
	#	Else:
	#		if --noWeights:
	#			conflict resolve, keeping right neighbor
	#		else:
	#			query table to get weights
	#			conflict resolve with list of weights (make a hash lookup in-function)
	#		In above:
	#			Keep in mind the --graphApproximate and --graphMax, which determine
	#			which algorithm for finding the maximal independent set we will use.

	num_nodes = len(blacklist_edges)
	if (params._noGraph):
		#If _noGraph: Just delete them all!
		blacklist_nodes = graph.listFromEdges(blacklist_edges) #Should return list of nodes from a list of edges
		return(blacklist_nodes)
	else:
		revised = []
		#If _noWeightGraph, or too many nodes, just run the approximate algorithm
		if (params._noWeightGraph) or (num_nodes > params._weightMax):
			revised = graph.edgeResolveApproximate(blacklist_edges) #returns list
			#print(revised)
		else:
			blacklist_nodes = graph.listFromEdges(blacklist_edges)
			weights_df = m.getRegionWeightsByList_VAR(conn, blacklist_nodes)
			# print("weights:")
			# print(weights_df)
			if params._weightByMin:
				weights_df = m.getRegionWeightsByList_BAD(conn, blacklist_nodes)
			#convert pandas df of weights to a dict
			weights = utils.dictFromDF(weights_df)
			revised = graph.edgeResolveWeighted(blacklist_edges, weights) #returns list
		return(revised)
	#os.remove(pw)



#function for sliding window bait generation
def baitSlidingWindow(conn, source, sequence, overlap, length):
	generator = s.slidingWindowGenerator(sequence, overlap, length)
	for window_seq in generator():
	#Don't need to do a bunch of filtering, because all was checked when TRs built
	#print(window_seq)
		if (len(window_seq[0]) == length):
			n_mask = utils.n_lower_chars(window_seq[0])
			n_gc = s.gc_counts(window_seq[0])
			m.add_bait_record(conn, source, window_seq[0], window_seq[1], window_seq[2], n_mask, n_gc)

#function for sliding window bait generation, with custom coordinates
def baitSlidingWindowCoord(conn, source, sequence, overlap, length, start):
	generator = s.slidingWindowGenerator(sequence, overlap, length)
	for window_seq in generator():
		#Don't need to do a bunch of filtering, because all was checked when TRs built
		#print(window_seq)
		if (len(window_seq[0]) == length):
			start_coord = start + window_seq[1]
			stop_coord = start_coord + length
			n_mask = utils.n_lower_chars(window_seq[0])
			n_gc = s.gc_counts(window_seq[0])
			m.add_bait_record(conn, source, window_seq[0], start_coord, stop_coord, n_mask, n_gc)

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

#Function to filter target regions by --filter_R arguments
def filterBaits(conn, params):
	rand = 0 #false
	rand_num = 0
	aln = 0

	for option in params.filter_b_objects:
		print("Select Region Option: ", option.o1)
		if option.o1 == "rand":
		#Set 'rand' to TRUE for random selection AFTER other filters
			rand_num = int(option.o2)
			assert rand_num > 0, "Number for random bait selection must be greater than zero!"
		elif option.o1 == "mask":
			min_mask_prop = option.o2
			m.baitFilterMask(conn, maxprop=max_mask_prop)
		elif option.o1 == "gc":
			min_mask_prop = option.o2
			max_mask_prop = option.o3
			m.baitFilterGC(conn, minprop=min_mask_prop, maxprop=max_mask_prop)
		elif option.o1 == "pw":
			aln = 1
		else:
			assert False, "Unhandled option %r"%option

	#Perform pairwise alignment AFTER all other filters because it is analytically more expensive
	if aln:
		passedBaits = m.getPassedBaits(conn)
		minid = option.o2
		mincov= option.o3
		assert (0.0 < minid < 1.0), "Minimum ID for pairwise alignment must be between 0.0 and 1.0"
		assert (0.0 < mincov < 1.0), "Minimum alignment coverage for pairwise alignment must be between 0.0 and 1.0"
		pairwiseAlignDedup(params, passedBaits, minid, mincov)

	#If 'random' select is turned on, then apply AFTER all other options
	if rand and not rand_num:
		m.baitFilterRandom(conn, rand_num)


#Function to print baits in final output
def printBaits(conn, params):
	pass
