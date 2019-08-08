#!/usr/bin/python

import subprocess
from subprocess import Popen, PIPE, CalledProcessError
import os
import sys
import Bio
import pandas as pd
from mrbait import misc_utils as utils

"""Includes utilities for calling NCBI BLAST and parsing outputs"""

#Function to parse params object to call blastn, and parse output to remove all matching queries
#Parses a MrBait parseArg object
def blastExcludeMatch(opts, db, fas, pid, qcov, out):
	dust = "yes"
	if opts.nodust:
		dust = "no"
	try:
		blastn(opts.blastn, fas, db, opts.blast_method,opts.gapopen, opts.gapextend, opts.evalue, opts.threads, opts.max_hits, dust, out)
	except KeyboardInterrupt:
		sys.exit("Process aborted: Keyboard Interrupt")
	except subprocess.CalledProcessError as err:
		print("BLASTN encountered an unknown problem in subprocess call: ", end="")
		sys.exit(err.args)
	except OSError as err:
		print("Exception: OSError in BLASTN call. Check that you are using the correct executable.", end="")
		sys.exit(err.args)
	except AttributeError as err:
		print("AttributeError in blastn call (internal problem, contact developer): ", end="")
		sys.exit(err.args)
	except TypeError as err:
		print("TypeError in blastn call (internal problem, contact developer): ", end="")
		sys.exit(err.args)
	except:
		print("BLASTN encountered an unknown runtime problem: ", end="")
		sys.exit(sys.exc_info()[0])

	#Parse output and return blacklisted IDs
	results = getBlastResults(out)
	if (len(results) < 1):
		return(list())
	pid_adj = pid*100
	blacklist = results[(results.pident > pid_adj) & (results.qcov > qcov)]["qseqid"].tolist()
	blacklist = [s.replace('id_', '') for s in blacklist]
	return(list(set(blacklist)))

#Function to parse params object to call blastn, and parse output to remove all matching queries
#Parses a MrBait parseArg object
#removes targets having more than 1 NON-OVERLAPPING alignment to a given genome
def blastExcludeAmbig(opts, db, fas, pid, qcov, out):
	dust = "yes"
	if opts.nodust:
		dust = "no"
	try:
		blastn(opts.blastn, fas, db, opts.blast_method,opts.gapopen, opts.gapextend, opts.evalue, opts.threads, opts.max_hits, dust, out)
	except KeyboardInterrupt:
		sys.exit("Process aborted: Keyboard Interrupt")
	except subprocess.CalledProcessError as err:
		print("BLASTN encountered an unknown problem in subprocess call: ", end="")
		sys.exit(err.args)
	except OSError as err:
		print("Exception: OSError in BLASTN call. Check that you are using the correct executable.", end="")
		sys.exit(err.args)
	except AttributeError as err:
		print("AttributeError in blastn call (internal problem, contact developer): ", end="")
		sys.exit(err.args)
	except TypeError as err:
		print("TypeError in blastn call (internal problem, contact developer): ", end="")
		sys.exit(err.args)
	except:
		print("BLASTN encountered an unknown runtime problem: ", end="")
		sys.exit(sys.exc_info()[0])

	#Parse output and return blacklisted IDs
	results = getBlastResults(out)
	if (len(results) < 1):
		return(list())
	pid_adj = pid*100
	
	#reduce to matches passing thresholds, and only queries with >1 hit
	subset = results[(results.pident > pid_adj) & (results.qcov > qcov)]
	subset = subset[subset.groupby('qseqid').qseqid.transform(len) > 1]
	
	#
	blacklist = list()
	unique_queries = dict(list(subset.groupby(['qseqid','sseqid'])))
	if (len(unique_queries) > 1):
		index = next(iter(unique_queries))
		#print(index[0])
		blacklist.append(index[0])
	else:
	# print("All hits are to same target: Check for overlaps!")
		for qseqid, qseqid_group in unique_queries.items():
			nonoverlapping_hits = 0
			for (indx1,row1),(indx2,row2) in zip(qseqid_group[:-1].iterrows(),qseqid_group[1:].iterrows()):
				if utils.calcOverlap(row1["sstart"], row1["send"], row2["sstart"], row2["send"]):
					#print(qseqid[0])
					blacklist.append(qseqid[0])
					break
					
	blacklist = [s.replace('id_', '') for s in blacklist]
	return(list(set(blacklist)))


#Function to parse params object to call blastn, and parse output to INCLUDE all matching queries
#Parses a MrBait parseArg object
#returns a whitelist of IDs to KEEP
def blastIncludeMatch(opts, db, fas, pid, qcov, out):
	dust = "yes"
	if opts.nodust:
		dust = "no"
	try:
		blastn(opts.blastn, fas, db, opts.blast_method,opts.gapopen, opts.gapextend, opts.evalue, opts.threads, opts.max_hits, dust, out)
	except KeyboardInterrupt:
		sys.exit("Process aborted: Keyboard Interrupt")
	except subprocess.CalledProcessError as err:
		print("BLASTN encountered an unknown problem in subprocess call: ", end="")
		sys.exit(err.args)
	except OSError as err:
		print("Exception: OSError in BLASTN call. Check that you are using the correct executable.", end="")
		sys.exit(err.args)
	except AttributeError as err:
		print("AttributeError in blastn call (internal problem, contact developer): ", end="")
		sys.exit(err.args)
	except TypeError as err:
		print("TypeError in blastn call (internal problem, contact developer): ", end="")
		sys.exit(err.args)
	except:
		print("BLASTN encountered an unknown runtime problem: ", end="")
		sys.exit(sys.exc_info()[0])

	#Parse output and return blacklisted IDs
	results = getBlastResults(out)
	if (len(results) < 1):
		return(list())
	pid_adj = pid*100
	whitelist = list()
	whitelist = results[(results.pident > pid_adj) & (results.qcov > qcov)]["qseqid"].tolist()
	whitelist = [s.replace('id_', '') for s in whitelist]
	return(list(set(whitelist)))

#Function reads a blast oufmt-6 table and returns pandas dataframe
def getBlastResults(out):
	r = pd.read_csv(out, sep="\t",names = ["qseqid", "sseqid", "pident", "qlen", "length", "evalue", "bitscore", "sstart", "send"])
	r["qcov"] = r["length"] / r["qlen"]
	return(r)


#Function to build blastn command and make subprocess call
def blastn(binary, fasta, blastdb, method, gapopen, gapextend, e_value, threads, hits, dust, outfile):
	#print("starting")
	#NOTE: I changed the columns included in outfmt 6.
	mask = "true"
	if dust == "no":
		mask = "false"

	blastn_cli = [binary,
		"-task", str(method),
		"-query", str(fasta),
		"-db", str(blastdb),
		"-evalue", str(e_value),
		"-outfmt", "6 qseqid sseqid pident qlen length evalue bitscore sstart send",
		"-num_threads", str(threads),
		"-max_target_seqs", str(hits),
		"-max_hsps", str(1),
		"-dust", str(dust),
		"-soft_masking", str(mask),
		"-out", str(outfile)]

	command = " ".join(blastn_cli)
	#print("BLASTN command is:",command )

	#Vsearch subprocess
	#print("Running blast...")
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
	#print("Command is:",command )

	#Vsearch subprocess
	#print("Running blast...")
	proc = Popen(makedb, stdout=PIPE, stdin=PIPE, env={'PATH': os.getenv('PATH')})

	#WRap to enable keyboard interrupe
	try:
		t = proc.communicate()[0]
	except KeyboardInterrupt:
		proc.kill()
		raise KeyboardInterrupt

	#Get return code from process
	if proc.returncode:
		raise CalledProcessError ("MAKEBLASTDB exited with non-zero status")
