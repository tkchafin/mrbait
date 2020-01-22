#!/usr/bin/python
import sys
import sqlite3
import getopt
import Bio
import os
import time
from Bio import AlignIO
from mrbait import mrbait_menu
from mrbait import substring
from mrbait.substring import SubString
from functools import partial
from mrbait import manage_bait_db as m
from mrbait import alignment_tools as a
from mrbait import sequence_tools as s
from mrbait import misc_utils as utils
from mrbait import seq_graph as graph
from mrbait import aln_file_tools
from mrbait import vcf_tools
from mrbait import vsearch
from mrbait import gff3_parser as gff
from mrbait import blast as b
import subprocess
import pandas as pd
import numpy as np
import multiprocessing

"""

Parallel versions of some of the MrBait corefuncs.

Much thanks to SO user 'dano' for 2014 post on how to share lock in multiprocessing pool:
https://stackoverflow.com/questions/25557686/python-sharing-a-lock-between-processes

"""

#Function to load a XMFA file into database
def loadXMFA_parallel(conn, params):

	t = int(params.threads)
	numLoci = aln_file_tools.countXMFA(params.xmfa)
	if numLoci < 10000:
		print("\t\t\tReading",numLoci,"alignments.")
	else:
		print("\t\t\tReading",numLoci,"alignments... This may take a while.")

	#file chunker call
	file_list = aln_file_tools.xmfa_chunker(params.xmfa, t, params.workdir)

	#print("Files are:",file_list)
	#Initialize multiprocessing pool
	#if 'lock' not in globals():
	lock = multiprocessing.Lock()
	try:
		with multiprocessing.Pool(t,initializer=init, initargs=(lock,)) as pool:
			func = partial(loadXMFA_worker, params.db, params.cov, params.minlen, params.thresh, params.mask, params.maf)
			results = pool.map(func, file_list)
	except Exception as e:
		pool.close()
	pool.close()
	pool.join()

	#reset_lock()
	#Remove chunkfiles
	aln_file_tools.removeChunks(params.workdir)

#worker function version of loadMAF
def loadXMFA_worker(db, params_cov, params_minlen, params_thresh, params_mask, params_maf, chunk):
	try:
		connection = sqlite3.connect(db)
		#Parse MAF file and create database
		for aln in AlignIO.parse(chunk, "mauve"):
			#NOTE: Add error handling, return error code
			cov = len(aln)
			alen = aln.get_alignment_length()

			if cov < params_cov or alen < params_minlen:
				continue
			#Add each locus to database
			locus = a.consensAlign(aln, threshold=params_thresh, mask=params_mask, maf=params_maf)
			lock.acquire()
			locid = m.add_locus_record(connection, cov, locus.conSequence, 1, "NULL")
			lock.release()
		connection.close()
	except Exception as e:
		raise Exception(e.message)


#Function to load LOCI file in parallel
def loadLOCI_parallel(conn, params):
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
	t = int(params.threads)

	numLoci = aln_file_tools.countLoci(params.loci)
	if numLoci < 10000:
		print("\t\t\tReading",numLoci,"alignments.")
	else:
		print("\t\t\tReading",numLoci,"alignments... This may take a while.")

	#file chunker call
	file_list = aln_file_tools.loci_chunker(params.loci, t, params.workdir)

	#print("Files are:",file_list)
	#Initialize multiprocessing pool
	#if 'lock' not in globals():
	lock = multiprocessing.Lock()
	try:
		with multiprocessing.Pool(t,initializer=init, initargs=(lock,)) as pool:
			func = partial(loadLOCI_worker, params.db, params.cov, params.minlen, params.thresh, params.mask, params.maf)
			results = pool.map(func, file_list)
	except Exception as e:
		pool.close()
	pool.close()
	pool.join()

	#reset_lock()
	#Remove chunkfiles
	aln_file_tools.removeChunks(params.workdir)

#Function to load MAF file in parallel
def loadMAF_parallel(conn, params):

	t = int(params.threads)
	numLoci = aln_file_tools.countMAF(params.alignment)
	if numLoci < 10000:
		print("\t\t\tReading",numLoci,"alignments.")
	else:
		print("\t\t\tReading",numLoci,"alignments... This may take a while.")

	#file chunker call
	file_list = aln_file_tools.maf_chunker(params.alignment, t, params.workdir)

	#print("Files are:",file_list)
	#Initialize multiprocessing pool
	#if 'lock' not in globals():
	lock = multiprocessing.Lock()
	try:
		with multiprocessing.Pool(t,initializer=init, initargs=(lock,)) as pool:
			func = partial(loadMAF_worker, params.db, params.cov, params.minlen, params.thresh, params.mask, params.maf)
			results = pool.map(func, file_list)
	except Exception as e:
		pool.close()
	pool.close()
	pool.join()

	#reset_lock()
	#Remove chunkfiles
	aln_file_tools.removeChunks(params.workdir)

# #first chunking, then arsing in parallel
# def loadVCF_parallel(conn, params):
# 	t = int(params.threads)
# 	#file chunker call
# 	file_list = vcf_tools.vcf_chunker(params.vcf, t, params.workdir)
#
# 	print("Files are:",file_list)
# 	#Initialize multiprocessing pool
# 	#if 'lock' not in globals():
# 	lock = multiprocessing.Lock()
# 	try:
# 		with multiprocessing.Pool(t,initializer=init, initargs=(lock,)) as pool:
# 			func = partial(loadVCF_worker, params.db, params.thresh)
# 			results = pool.map(func, file_list)
# 	except Exception as e:
# 		pool.close()
# 	pool.close()
# 	pool.join()
#
# 	#reset_lock()
# 	#Remove chunkfiles
# 	#aln_file_tools.removeChunks(params.workdir)


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
def loadMAF_worker(db, params_cov, params_minlen, params_thresh, params_mask, params_maf, chunk):
	try:
		connection = sqlite3.connect(db)
		#Parse MAF file and create database
		for aln in AlignIO.parse(chunk, "maf"):
			#NOTE: Add error handling, return error code
			cov = len(aln)
			alen = aln.get_alignment_length()

			if cov < params_cov or alen < params_minlen:
				continue
			#Add each locus to database
			locus = a.consensAlign(aln, threshold=params_thresh, mask=params_mask, maf=params_maf)
			lock.acquire()
			locid = m.add_locus_record(connection, cov, locus.conSequence, 1, "NULL")
			#print(locid)
			lock.release()
			#Extract variable positions for database
			#for var in locus.alnVars:
				#m.add_variant_record(connection, locid, var.position, var.value)
		connection.close()
	except Exception as e:
		raise Exception(e.message)

# #Function to load VCF variants file
# def loadVCF_worker(db, threshold, chunk):
# 	try:
# 		#Each worker opens unique connection to db
# 		connection = sqlite3.connect(db)
# 		#Lock DB and read loci, then release lock
# 		lock.acquire()
# 		loci = m.getPassedLoci(connection) #get DF of passed loci
# 		lock.release()
#
# 		chrom_lookup = loci.set_index('chrom')['id'].to_dict()
# 		loci.set_index('id', inplace=True)
#
# 		passed=0 #To track number of VCF records for which no locus exists
# 		failed=0
# 		for reclist in vcf_tools.read_vcf(chunk):
# 			rec_chrom = reclist[0].CHROM
# 			if rec_chrom in chrom_lookup:
# 				locid = chrom_lookup[rec_chrom]
# 				passed+=1
# 				#Grab DF record for the matching CHROM
# 				seq = loci.loc[locid,'consensus']
# 				#Get new consensus sequence given VCF records
# 				new_cons = vcf_tools.make_consensus_from_vcf(seq,rec_chrom,reclist, threshold)
# 				print(new_cons)
# 				#Update new consensus seq in db
# 				if len(new_cons) != len(seq): #Check length first
# 					print("\t\t\tWarning: New consensus sequence for locus %s (locid=<%s>) is the wrong length! Skipping."%(rec_chrom, locid))
# 				else:
# 					#Lock database for update, then relase lock
# 					lock.acquire()
# 					m.updateConsensus(connection, locid, new_cons)
# 					lock.release()
# 			else:
# 				failed+=1
# 		if failed > 0:
# 			print("\t\t\tWARNING:%s/%s records in <%s> don't match any reference sequences"%(failed, failed+passed, chunk))
# 		#close connection
# 		connection.close()
# 	except Exception as e:
# 		raise Exception(e.message)

#Worker function for loadLOCI_parallel
def loadLOCI_worker(db, params_cov, params_minlen, params_thresh, params_mask, params_maf, chunk):
	try:
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
				locus = a.consensAlign(aln, threshold=params_thresh, mask=params_mask, maf=params_maf)

				#Acquire lock, submit to Database
				lock.acquire()
				locid = m.add_locus_record(connection, cov, locus.conSequence, 1, "NULL")
				lock.release()
				#print("Loading Locus #:",locid)

				#Extract variable positions for database
				#for var in locus.alnVars:
					#m.add_variant_record(connection, locid, var.position, var.value)
		connection.close()
	except Exception as e:
		raise Exception(e.message)

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
		chunk_file = params.workdir + "/." + str(chunk) + ".chunk"
		#print(df_chunk)
		df_chunk.to_csv(chunk_file, mode="w", index=False)
		files.append(chunk_file)
		chunk += 1

	#Initialize multiprocessing pool
	#if 'lock' not in globals():
	lock = multiprocessing.Lock()
	with multiprocessing.Pool(t,initializer=init, initargs=(lock,)) as pool:
		func = partial(targetDiscoverySlidingWindow_worker, params.db, params.win_shift, params.win_width, params.var_max, params.numN, params.numG, params.blen, params.flank_dist, params.target_all)
		results = pool.map(func, files)
	pool.close()
	pool.join()

	#reset_lock()
	#Remove chunkfiles
	d = os.listdir(params.workdir)
	for item in d:
		if item.endswith(".chunk"):
			os.remove(os.path.join(params.workdir, item))


#Function to discover target regions using a sliding windows through passedLoci
def targetDiscoverySlidingWindow_worker(db, shift, width, var, n, g, blen, flank_dist, target_all, chunk):
	connection = sqlite3.connect(db)
	#print("process: reading hdf from",chunk)
	loci = pd.read_csv(chunk)
	#print(loci)
	for seq in loci.itertuples():
		#print(seq)
		start = 0
		stop = 0
		if target_all:
			#print("target_all")
			#submit full locus as target
			seq_norm = s.simplifySeq(seq[2])
			counts = s.seqCounterSimple(seq_norm)
			if counts['*'] <= var and counts['N'] <= n and counts['-'] <= g:
				target = seq[2]
				tr_counts = s.seqCounterSimple(seq_norm)
				n_mask = utils.n_lower_chars(seq[2])
				n_gc = s.gc_counts(seq[2])
				#NOTE: flank count set to number of variable sites in whole locus
				#print(int(seq[1]), 0, len(seq[2]), seq[2], tr_counts, tr_counts, n_mask, n_gc)
				lock.acquire()
				m.add_region_record(connection, int(seq[1]), 0, len(seq[2]), seq[2], tr_counts, tr_counts, n_mask, n_gc)
				lock.release()
		else:
			#print("\nConsensus: ", seq[2], "ID is: ", seq[1], "\n")
			generator = s.slidingWindowGenerator(seq[2], shift, width)
			for window_seq in generator():

				seq_norm = s.simplifySeq(window_seq[0])
				counts = s.seqCounterSimple(seq_norm)

				#If window passes filters, extend current bait region
				#print("Start is ", start, " and stop is ",stop) #debug print
				if counts['*'] <= var and counts['N'] <= n and counts['-'] <= g:
					stop = window_seq[2]
					#if this window passes BUT is the last window, evaluate it
					if stop == len(seq[2]):
						#print("last window")
						if (stop - start) >= blen:
							target = (seq[2])[start:stop]
							tr_counts = s.seqCounterSimple(s.simplifySeq(target))
							#print("candidate:",window_seq[0])
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
				else:
					#If window fails, check if previous bait region passes to submit to DB
					#print (stop-start)
					if (stop - start) >= blen:
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
