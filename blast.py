#!/usr/bin/python

from subprocess import Popen, PIPE, CalledProcessError
import os
import Bio

"""Includes utilities for calling VSEARCH and parsing output of pairwise alignments"""

def blastn(binary, fasta, blastdb, e_value, threads, outfile):
	blastn_cli = [binary,
		"-task", "blastn",
		"-query", str(fasta),
		"-db", str(blastdb),
		"-evalue", str(e_value),
		"-outfmt", "6",
		"-num_threads", str(threads),
		"-max_target_seqs", "1",
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

#Function to parse output of allpairsGlobal
def parsePairwiseAlign(filename):
	#Function assumes id's are the first 2 positions, and prepended with "id_"
	#Also assumes file is tab delimited
	print("parsing pw file")
	bad_ids = []
	filehandle = open(filename, "r")
	for line in filehandle:
		array = re.split(r'\t+', line)
		bad1 = re.sub('id_', '', array[0])
		bad2 = re.sub('id_', '', array[1])
		bad_ids.append([bad1, bad2])
	return(bad_ids)


blastn("./bin/ncbi-blastn-2.6.0-macos","baits.fasta","ref.fasta", 0.00001, 4, "test.out")
