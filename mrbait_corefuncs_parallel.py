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
	sys.exit()
	#print("Files are:",file_list)
	#Initialize multiprocessing pool
	#if 'lock' not in globals():
	lock = multiprocessing.Lock()
	pool = multiprocessing.Pool(t,initializer=init, initargs=(lock,))
	func = getWorkerFunc(params, ftype)
	results = pool.map(func, file_list)
	pool.close()
	pool.join()

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
		pass

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
			for var in locus.alnVars:
				m.add_variant_record(connection, locid, var.position, var.value)
			lock.release()
	connection.close()
