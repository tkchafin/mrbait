#!/usr/bin/python

from subprocess import Popen, PIPE, CalledProcessError
import os
import sys
import Bio

"""Includes utilities for calling VSEARCH and parsing output of pairwise alignments"""

#Function to build blastn command and make subprocess call
def blastn(binary, fasta, blastdb, method, gapopen, gapextend, e_value, threads, hits, outfile):
	blastn_cli = [binary,
		"-task", str(method),
		"-query", str(fasta),
		"-db", str(blastdb),
		"-evalue", str(e_value),
		"-outfmt", "6",
		"-num_threads", str(threads),
		"-max_target_seqs", str(hits),
		"-out", str(outfile)]

	command = " ".join(blastn_cli)
	print("Command is:",command )

	#Vsearch subprocess
	print("Running blast...")
	proc = Popen(blastn_cli, stdout=PIPE, stdin=PIPE, env={'PATH': os.getenv('PATH')})

	#WRap to enable keyboard interrupe
	try:
		t = proc.communicate()[0]
	except KeyboardInterrupt:
		proc.kill()
		raise KeyboardInterrupt

	#Get return code from process
	if proc.returncode:
		raise CalledProcessError ("BLASTN exited with non-zero status")

#Function to build makeblastdb command and make subprocess call
def makeblastdb(binary, reference, outfile):
	makedb = [binary,
		"-in", str(reference),
		"-dbtype", "nucl",
		"-out", str(outfile)]

	command = " ".join(makedb)
	print("Command is:",command )

	#Vsearch subprocess
	print("Running blast...")
	proc = Popen(makedb, stdout=PIPE, stdin=PIPE, env={'PATH': os.getenv('PATH')})

	#WRap to enable keyboard interrupe
	try:
		t = proc.communicate()[0]
	except KeyboardInterrupt:
		proc.kill()
		raise KeyboardInterrupt

	#Get return code from process
	if proc.returncode:
		raise CalledProcessError ("BLASTN exited with non-zero status")


#TODO: Set "hits" to 2 if looking for seqs with ambiguous mappings
try:
	makeblastdb("./bin/ncbi-makeblastdb-2.6.0-macos", "ref.fasta", "/Users/tkchafin/test/test")
except:
	sys.exit(sys.exc_info()[0])
blastn("./bin/ncbi-blastn-2.6.0-macos","baits.fasta","/Users/tkchafin/test/test","blastn", 5, 2, 0.00001, 4, 1, "test.out")
