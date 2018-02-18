import sys
import sqlite3
import getopt
import Bio
import os
import time
from Bio import AlignIO
from mrbait_menu import display_help
from mrbait_menu import parseArgs
from substring import SubString
from functools import partial
import manage_bait_db as m
import alignment_tools as a
import sequence_tools as s
import misc_utils as utils
import seq_graph as graph
import aln_file_tools
import pandas as pd
import numpy as np
import blast as b
import multiprocessing

"""

Parallel versions of some of the MrBait corefuncs.

"""


#Generic function to load whicever alignment file is given by
#first chunking, then arsing in parallel
def loadPARALLEL(conn, params, ftype):
	"""
	Format:
	multiprocessing pool.
	Master:
		splits file into n chunks
		creates multiprocessing pool
	Workers:
		read file chunk
		calculate consensus
		grab lock
		INSERT data to SQL database
		release lock
	"""

	if ftype not in ["loci", "maf"]:
		sys.exit("ERROR:",ftype,"not supported by loadPARALLEL()")

	t = int(params.threads)

	if ftype in ["loci", "maf"]:
		numLoci = getLociCount(params, ftype)
		if numLoci < 10000:
			print("\t\t\tReading",numLoci,"alignments.")
		else:
			print("\t\t\tReading",numLoci,"alignments... This may take a while.")

	#file chunker call
	file_list = loadChunker(params, t, ftype)

	#print("Files are:",file_list)
	#Initialize multiprocessing pool
	#if 'lock' not in globals():
	lock = multiprocessing.Lock()
	with multiprocessing.Pool(t,initializer=init, initargs=(lock,)) as pool:
		func = getWorkerFunc(params, ftype)
		results = pool.map(func, file_list)
	pool.close()
	pool.join()

	#reset_lock()
	#Remove chunkfiles
	aln_file_tools.removeChunks(params.workdir)

#Function to chunk input aln file
def loadChunker(params, t, ftype):
	if ftype == "loci":
		return(aln_file_tools.loci_chunker(params.loci, t, params.workdir))
	elif ftype == "maf":
		print("calling maf chunker")
		return(aln_file_tools.maf_chunker(params.alignment, t, params.workdir))

#Function to build and return partial function for worker calls
def getWorkerFunc(params, ftype):
	if ftype == "loci":
		f = partial(loadLOCI_worker,params.db, params.cov, params.minlen, params.thresh, params.mask)
		return(f)
	elif ftype == "maf":
		f = partial(loadMAF_worker, params.db, params.cov, params.minlen, params.thresh, params.mask)
		return(f)

#Function to get count of loci in aln file
def getLociCount(params, ftype):
	if ftype == "loci": #pyRAD LOCI
		return(aln_file_tools.countLoci(params.loci))
	elif ftype == "maf": #MAF
		return(aln_file_tools.countMAF(params.alignment))

#INitialize a global lock. Doing it this way allows it to be inherited by the child processes properly
#Found on StackOverflow: https://stackoverflow.com/questions/25557686/python-sharing-a-lock-between-processes
#Thanks go to SO user dano
def init(l):
    global lock
    lock = l

#Function to reset lock
def reset_lock():
	global lock
	del lock

#NOTE: 'params' object can't be pickled, so I have to do it this way.
#worker function version of loadMAF
def loadMAF_worker(db, params_cov, params_minlen, params_thresh, params_mask, chunk):
	connection = sqlite3.connect(db)
	#Parse MAF file and create database
	for aln in AlignIO.parse(chunk, "maf"):
		#NOTE: Add error handling, return error code
		cov = len(aln)
		alen = aln.get_alignment_length()

		if cov < params_cov or alen < params_minlen:
			continue
		#Add each locus to database
		locus = a.consensAlign(aln, threshold=params_thresh, mask=params_mask)
		lock.acquire()
		locid = m.add_locus_record(connection, cov, locus.conSequence, 1, "NULL")
		lock.release()
		#Extract variable positions for database
		#for var in locus.alnVars:
			#m.add_variant_record(connection, locid, var.position, var.value)


	connection.close()

#Worker function for loadLOCI_parallel
def loadLOCI_worker(db, params_cov, params_minlen, params_thresh, params_mask, chunk):
	connection = sqlite3.connect(db)
	#Parse LOCI file and create database
	for aln in aln_file_tools.read_loci(chunk):
		#NOTE: Add error handling, return error code
		cov = len(aln)
		alen = aln.get_alignment_length()

		#Skip if coverage or alignment length too short
		if cov < params_cov or alen < params_minlen:
			#print("Locus skipped")
			continue
		else:
			#Add each locus to database
			locus = a.consensAlign(aln, threshold=params_thresh, mask=params_mask)

			#Acquire lock, submit to Database
			lock.acquire()
			locid = m.add_locus_record(connection, cov, locus.conSequence, 1, "NULL")
			#print("Loading Locus #:",locid)

			#Extract variable positions for database
			#for var in locus.alnVars:
				#m.add_variant_record(connection, locid, var.position, var.value)
			lock.release()
	connection.close()


#Function to discover target regions using a sliding windows through passedLoci
def targetDiscoverySlidingWindow_parallel(conn, params, loci):
	"""
	Format:
	1. Write pandas DF to n chunk files
	2. List of chunk file names
	3. Pass 1 chunk file to each worker in a multiprocessing pool.
	Master:
		creates n chunk files
		creates multiprocessing pool
	Workers:
		read file chunk
		calculate consensus
		grab lock
		INSERT data to SQL database
		release lock
	"""
	t = int(params.threads)
	chunk = 1
	loci_num = int(loci.shape[0])
	#print("number of loci:",loci_num)
	#print("number of threads:",t)
	chunks = 0
	if loci_num < t:
		chunks = loci_num
	else:
		chunks = t
	chunk_size = loci_num // chunks
	remainder = loci_num % chunks
	#print("Chunk size is:",chunk_size)
	#print("remainder is:",remainder)

	start = 0
	stop = 0

	files = list()
	#Split loci DataFrame into chunks, and keep list of chunk files
	for df_chunk in np.array_split(loci, chunks):
		size = df_chunk.shape[0]
		#print("size of chunk",chunk,"is:",size)
		chunk_file = params.workdir + "/." + str(chunk) + ".chunk.hdf"
		df_chunk.to_hdf(chunk_file, "NULL", mode="w", format="f")
		files.append(chunk_file)
		chunk += 1

	#Initialize multiprocessing pool
	#if 'lock' not in globals():
	lock = multiprocessing.Lock()
	with multiprocessing.Pool(t,initializer=init, initargs=(lock,)) as pool:
		func = partial(targetDiscoverySlidingWindow_worker, params.db, params.win_shift, params.win_width, params.var_max, params.numN, params.numG, params.blen, params.flank_dist)
		results = pool.map(func, files)
	pool.close()
	pool.join()

	#reset_lock()
	#Remove chunkfiles
	d = os.listdir(params.workdir)
	for item in d:
	    if item.endswith(".chunk.hdf"):
	        os.remove(os.path.join(params.workdir, item))


#Function to discover target regions using a sliding windows through passedLoci
def targetDiscoverySlidingWindow_worker(db, shift, width, var, n, g, blen, flank_dist, chunk):
	connection = sqlite3.connect(db)
	#print("process: reading hdf from",chunk)
	loci = pd.read_hdf(chunk)
	#looping through passedLoci only
	#print("process: starting iteration")
	for seq in loci.itertuples():
		#print(seq)
		start = 0
		stop = 0
		#print("\nConsensus: ", seq[2], "ID is: ", seq[1], "\n")
		generator = s.slidingWindowGenerator(seq[2], shift, width)
		for window_seq in generator():

			seq_norm = s.simplifySeq(window_seq[0])
			counts = s.seqCounterSimple(seq_norm)

			#If window passes filters, extend current bait region
			#print("Start is ", start, " and stop is ",stop) #debug print
			if counts['*'] <= var and counts['N'] <= n and counts['-'] <= g:
				stop = window_seq[2]
			else:
				#If window fails, check if previous bait region passes to submit to DB
				#print (stop-start)
				if (stop - start) > blen:
					target = (seq[2])[start:stop]
					tr_counts = s.seqCounterSimple(s.simplifySeq(target))
					n_mask = utils.n_lower_chars(target)
					n_gc = s.gc_counts(target)
					#Check that there aren't too many SNPs
					#if tr_counts["*"] <= params.vmax_r:
					#print("	Target region: ", target)
					#Submit target region to database
					#print("process: grabbing lock")'
					flank_counts = s.getFlankCounts(seq[2], start, stop, flank_dist)
					lock.acquire()
					m.add_region_record(connection, int(seq[1]), start, stop, target, tr_counts, flank_counts, n_mask, n_gc)
					#print("process: releasing lock")
					lock.release()
					#set start of next window to end of current TR
					generator.setI(stop)

				#If bait fails, set start to start point of next window
				start = generator.getI()+shift
	connection.close()


#Function to get DataFrame of targets + flank regions, and calculate some stuff
def flankDistParser_parallel(conn, dist):
	#Call manage_bait_db function to return DataFrame
	targets = m.getTargetFlanks(conn, dist)
