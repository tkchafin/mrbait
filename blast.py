#!/usr/bin/python

import subprocess
from subprocess import Popen, PIPE, CalledProcessError
import os
import sys
import Bio
import pandas as pd

"""Includes utilities for calling NCBI BLAST and parsing outputs"""

#Function to parse params object to call blastn, and parse output to remove all matching queries
#Parses a MrBait parseArg object
def blastExcludeMatch(opts, db, fas, pid, qcov, out):
	dust = "yes"
	if opts.nodust:
		dust = "no"
	try:
		blastn(opts.blastn, fas, db, opts.blast_method,opts.gapopen, opts.gapextend, opts.evalue, opts.threads, 1, dust, out)
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
	pid_adj = pid*100
	blacklist = results[(results.pident > pid_adj) & (results.qcov > qcov)]["qseqid"].tolist()
	blacklist = [s.replace('id_', '') for s in blacklist]
	return(blacklist)

#Function to parse params object to call blastn, and parse output to remove all matching queries
#Parses a MrBait parseArg object
#removes targets having more than 1 NON-OVERLAPPING alignment to a given genome
def blastExcludeAmbig(opts, db, fas, pid, qcov, out):
	dust = "yes"
	if opts.nodust:
		dust = "no"
	try:
		blastn(opts.blastn, fas, db, opts.blast_method,opts.gapopen, opts.gapextend, opts.evalue, opts.threads, 1, dust, out)
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
	pid_adj = pid*100
	subset = results[(results.pident > pid_adj) & (results.qcov > qcov)]
	subset = findNonOverlappingMatches(subset)
	#blacklist = [s.replace('id_', '') for s in blacklist]
	return(blacklist)


#Function to parse params object to call blastn, and parse output to INCLUDE all matching queries
#Parses a MrBait parseArg object
#returns a whitelist of IDs to KEEP
def blastIncludeMatch(opts, db, fas, pid, qcov, out):
	dust = "yes"
	if opts.nodust:
		dust = "no"
	try:
		blastn(opts.blastn, fas, db, opts.blast_method,opts.gapopen, opts.gapextend, opts.evalue, opts.threads, 1, dust, out)
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
	pid_adj = pid*100
	whitelist = results[(results.pident > pid_adj) & (results.qcov > qcov)]["qseqid"].tolist()
	whitelist = [s.replace('id_', '') for s in whitelist]
	return(whitelist)

#Function reads a blast oufmt-6 table and returns pandas dataframe
def getBlastResults(out):
	r = pd.read_csv(out, sep="\t",names = ["qseqid", "sseqid", "pident", "qlen", "length", "evalue", "bitscore"])
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
		"-outfmt", "6 qseqid sseqid pident qlen length evalue bitscore",
		"-num_threads", str(threads),
		"-max_target_seqs", str(hits),
		"-max_hsps", str(1),
		"-dust", str(dust),
		"-soft_masking", str(mask),
		"-out", str(outfile)]

	command = " ".join(blastn_cli)
	print("BLASTN command is:",command )

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
		raise CalledProcessError ("MAKEBLASTDB exited with non-zero status")
