#!/usr/bin/python
import sys
import sqlite3
import getopt
import Bio
import os
import time
from Bio import AlignIO
from decimal import *
from mrbait import mrbait_menu
from mrbait import substring
from mrbait.substring import SubString
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



############################# FUNCTIONS ################################

#Function to load a XMFA file into database
def loadXMFA(conn, params):
	numLoci = aln_file_tools.countXMFA(params.xmfa)
	if numLoci < 10000:
		print("\t\t\tReading",numLoci,"alignments.")
	else:
		print("\t\t\tReading",numLoci,"alignments... This may take a while.")
	#Parse XMFA file and create database

	num = 1
	for aln in AlignIO.parse(params.xmfa, "mauve"):
		#NOTE: Add error handling, return error code
		#print(aln)
		cov = len(aln)
		alen = aln.get_alignment_length()

		#Add each locus to database
		locus = a.consensAlign(aln, threshold=params.thresh, mask=params.mask, maf=params.maf)

		locid = m.add_locus_record(conn, cov, locus.conSequence, 1, num)
		num+=1


#Function to load a MAF file into database
def loadMAF(conn, params):
	numLoci = aln_file_tools.countMAF(params.alignment)
	if numLoci < 10000:
		print("\t\t\tReading",numLoci,"alignments.")
	else:
		print("\t\t\tReading",numLoci,"alignments... This may take a while.")
	#Parse MAF file and create database
	num = 1
	for aln in AlignIO.parse(params.alignment, "maf"):
		#NOTE: Add error handling, return error code
		cov = len(aln)
		alen = aln.get_alignment_length()

		#Add each locus to database
		locus = a.consensAlign(aln, threshold=params.thresh, mask=params.mask, maf=params.maf)
		#consensus = str(a.make_consensus(aln, threshold=params.thresh)) #Old way
		locid = m.add_locus_record(conn, cov, locus.conSequence, 1, num)
		num+=1

	#	print("Loading Locus #:",locid)

		"""deprecated"""
		#Extract variable positions for database
		#for var in locus.alnVars:
			#m.add_variant_record(conn, locid, var.position, var.value)

#Function to load .loci file into database.
def loadLOCI(conn, params):
	numLoci = aln_file_tools.countLoci(params.loci)
	if numLoci < 10000:
		print("\t\t\tReading",numLoci,"alignments.")
	else:
		print("\t\t\tReading",numLoci,"alignments... This may take a while.")
	#Parse LOCI file and create database
	for aln in aln_file_tools.read_loci(params.loci):
		#NOTE: Add error handling, return error code
		cov = len(aln)
		alen = aln.get_alignment_length()

		#Skip if coverage or alignment length too short
		if cov < params.cov or alen < params.minlen:
			#print("Locus skipped")
			continue
		else:
			#Add each locus to database
			locus = a.consensAlign(aln, threshold=params.thresh, mask=params.mask, maf=params.maf)
			#consensus = str(a.make_consensus(aln, threshold=params.thresh)) #Old way
			locid = m.add_locus_record(conn, cov, locus.conSequence, 1, "NULL")
			#print("Loading Locus #:",locid)

			#Extract variable positions for database
			#or var in locus.alnVars:
				#m.add_variant_record(conn, locid, var.position, var.value)

#Function to load FASTA into database
def loadFASTA(conn, params):
	for contig in aln_file_tools.read_fasta(params.assembly):
		#print("Reading contig:",contig[0])
		#print("Sequence is:",contig[1])
		locid = m.add_locus_record(conn, 1, contig[1], 1, contig[0])

		#Parse consensus for vars, submit those vars to db
		"""Deprecated"""
		#for var in a.get_vars(contig[1]):
			#m.add_variant_record(conn, locid, var.position, var.value)

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

#Function to load BED file
def loadBED(conn, params):

	with open(params.bed)as f:
		count=0
		for line in f:
			line = line.strip()
			if not line:
				continue
			count+=1
			if count <= params.bed_header:
				continue
			content = line.split()

			#NOTE: This function ONLY inserts BEDRecords where record.seqid matches an existing locus in the loci table
			m.add_bed_record(conn, content[0], content[1], content[2])

	#remove BED records not falling within our loci
	#print(m.getBED(conn))
	m.validateBEDRecords(conn)
	#print(m.getBED(conn))


#Function to load VCF variants file
def loadVCF(conn, params):

	loci = m.getPassedLoci(conn) #get DF of passed loci
	#print(loci)
	#chrom_lookup = set(loci["chrom"].tolist()) #lookup table of locus IDs
	chrom_lookup = loci.set_index('chrom')['id'].to_dict()
	loci.set_index('id', inplace=True)
	#print(loci)
	#print(chrom_lookup)
	#Using set_index causes a bug for some reason.
	#loci.set_index('chrom', inplace=True) #index loci DF by chrom column

	#print(loci)
	passed=0 #To track number of VCF records for which no locus exists
	failed=0

	for reclist in vcf_tools.read_vcf(params.vcf):
		rec_chrom = reclist[0].CHROM
		#print("Starting locus",rec_chrom)
		if rec_chrom in chrom_lookup:
			#print("chrom:",rec_chrom, " - id:",chrom_lookup[rec_chrom])
			locid = chrom_lookup[rec_chrom]
			#print("locid is",locid)
			passed+=1
			#for rec in reclist:
			#	print(rec.CHROM, rec.POS, rec.REF, rec.ALT, len(rec.samples), rec.call_rate, rec.aaf)
			#Grab DF record for the matching CHROM
			seq = loci.loc[locid,'consensus']
			#Get new consensus sequence given VCF records
			new_cons = vcf_tools.make_consensus_from_vcf(seq,rec_chrom,reclist, params.thresh, params.vcfALT)

			#Update new consensus seq in db
			if len(new_cons) != len(seq): #Check length first
				print("\t\t\tWarning: New consensus sequence for locus %s (locid=<%s>) is the wrong length! Skipping."%(rec_chrom, locid))
			else:
				m.updateConsensus(conn, locid, new_cons)
				#Delete old vars for locus, and parse new consensus
				"""
				#Deprecated
				#m.purgeVars(conn, locid)
				#Get new vars and add records to table
				#for var in a.get_vars(new_cons):
					#m.add_variant_record(conn, locid, var.position, var.value)
				"""

		else:
			#print(rec_chrom, "not found.")
			failed+=1
	if failed > 0:
		print("\t\t\tWARNING:%s/%s records in <%s> don't match any reference sequences"%(failed, failed+passed, params.vcf))

#Function to discover target regions using a sliding windows through passedLoci
def targetDiscoverySlidingWindow(conn, params, loci):

	#looping through passedLoci only
	for seq in loci.itertuples():
		#print(seq)
		start = 0
		stop = 0
		#print(params.win_shift)
		#print("\nConsensus: ", seq[2], "ID is: ", seq[1], "\n")
		if params.target_all:
			#print("target_all")
			#submit full locus as target
			seq_norm = s.simplifySeq(seq[2])
			counts = s.seqCounterSimple(seq_norm)
			if counts['*'] <= params.var_max and counts['N'] <= params.numN and counts['-'] <= params.numG:
				target = seq[2]
				tr_counts = s.seqCounterSimple(seq_norm)
				n_mask = utils.n_lower_chars(seq[2])
				n_gc = s.gc_counts(seq[2])
				#NOTE: flank count set to number of variable sites in whole locus
				#print(int(seq[1]), 0, len(seq[2]), seq[2], tr_counts, tr_counts, n_mask, n_gc)
				m.add_region_record(conn, int(seq[1]), 0, len(seq[2]), seq[2], tr_counts, tr_counts, n_mask, n_gc)
		else:
			#print("\nConsensus: ", seq[2], "ID is: ", seq[1], "\n")
			generator = s.slidingWindowGenerator(seq[2], params.win_shift, params.win_width)
			for window_seq in generator():

				seq_norm = s.simplifySeq(window_seq[0])
				counts = s.seqCounterSimple(seq_norm)

				#If window passes filters, extend current bait region
				#print("Start is ", start, " and stop is ",stop) #debug print
				if counts['*'] <= params.var_max and counts['N'] <= params.numN and counts['-'] <= params.numG:
					stop = window_seq[2]
					#if this window passes BUT is the last window, evaluate it
					if stop == len(seq[2]):
						#print("last window")
						if (stop - start) >= params.blen:
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
							flank_counts = s.getFlankCounts(seq[2], start, stop, params.flank_dist)
							m.add_region_record(conn, int(seq[1]), start, stop, target, tr_counts, flank_counts, n_mask, n_gc)
							#set start of next window to end of current TR
							generator.setI(stop)
				else:
					#If window fails, check if previous bait region passes to submit to DB
					#print (stop-start)
					if (stop - start) >= params.blen:
						target = (seq[2])[start:stop]
						tr_counts = s.seqCounterSimple(s.simplifySeq(target))
						n_mask = utils.n_lower_chars(target)
						n_gc = s.gc_counts(target)
						#Check that there aren't too many SNPs
						#if tr_counts["*"] <= params.vmax_r:
						#print("	Target region: ", target)
						#Submit target region to database
						#print("process: grabbing lock")'
						flank_counts = s.getFlankCounts(seq[2], start, stop, params.flank_dist)
						m.add_region_record(conn, int(seq[1]), start, stop, target, tr_counts, flank_counts, n_mask, n_gc)
						#set start of next window to end of current TR
						generator.setI(stop)

					#If bait fails, set start to start point of next window
					start = generator.getI()+params.win_shift
	#Now update regions table to include information for flanking regions if available
	#m.flankDistParser(conn, params.flank_dist)



#Function to filter target regions by --filter_R arguments
def filterTargetRegions(conn, params):

	rand = 0 #false
	rand_num = 0
	aln = 0
	blast = 0
	minid = None
	mincov = None
	if (len(params.filter_r_objects) <= 0):
		pass
	else:
		for option in params.filter_r_objects:
			#print("Filter Region Option: ", option.o1)
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
					print("\t\t\tWARNING: No blast database provided. Skipping <--filter_r> option %s"%option.o1)
					break
				#print("BLASTDB PATH IS: ", db_path)
				#Get targets, print to fasta
				seqs = m.getPassedTRs(conn)
				fas = params.workdir + "/.temp.fasta"
				aln_file_tools.writeFasta(seqs, fas)
				outfile = params.workdir + "/.temp.blast"
				if option.o1 == "blast_x":
					blacklist = b.blastExcludeMatch(params, db_path, fas, option.o2, option.o3, outfile)
					if len(blacklist) > 1:
						m.removeRegionsByList(conn, blacklist)
				elif option.o1 == "blast_i":
					whitelist = b.blastIncludeMatch(params, db_path, fas, option.o2, option.o3, outfile)
					if len(blacklist) > 1:
						m.removeRegionsByWhitelist(conn, whitelist)
				elif option.o1 == "blast_a":
					blacklist = b.blastExcludeAmbig(params, db_path, fas, option.o2, option.o3, outfile)
					if len(blacklist) > 1:
						m.removeRegionsByList(conn, blacklist)
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
			elif option.o1 in ("bed_x", "bed_i"):
				if params.bed and params.assembly:
					if option.o1 == "bed_x":
						m.regionFilterBED_exclude(conn, params.flank_dist)
					elif option.o1 == "bed_i":
						m.regionFilterBED_include(conn, params.flank_dist)
				else:
					sys.exit("ERROR: Filtering targets on proximity to BED elements requires FASTA <-A> and BED <-B> inputs!")
			else:
				assert False, "Unhandled option %r"%option

		#Perform pairwise alignment AFTER all other filters because it is analytically more expensive
		#Target region deduplication by pairwise alignment
		if aln:
			passedTargets = m.getPassedTRs(conn)
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
def filterTargetRegions_verbose(conn, params):

	rand = 0 #false
	rand_num = 0
	aln = 0
	blast = 0
	minid = None
	mincov = None
	if (len(params.filter_r_objects) <= 0):
		print("\t\t\tNo filtering criteria provided. Retaining all targets.")
	else:
		for option in params.filter_r_objects:
			#print("Filter Region Option: ", option.o1)
			if option.o1 == "rand":
				#Set 'rand' to TRUE for random selection AFTER other filters
				rand_num = int(option.o2)
				assert rand_num > 0, "Number for random TR selection must be greater than zero!"
			elif option.o1 == "gap":
				print("\t\t\tFiltering criterion: Maximum of",option.o2,"gaps")
				m.simpleFilterTargets_gap(conn, int(option.o2))
			elif option.o1 == "bad":
				print("\t\t\tFiltering criterion: Maximum of",option.o2,"ambiguities")
				m.simpleFilterTargets_bad(conn, int(option.o2))
			elif option.o1 == "snp":
				print("\t\t\tFiltering criterion: Between",option.o2,"and",option.o3,"flanking SNPs")
				m.simpleFilterTargets_SNP(conn, int(option.o2), int(option.o3))
			elif option.o1 == "mask":
				print("\t\t\tFiltering criterion: Maximum",option.o2,"proportion masked")
				max_mask_prop = option.o2
				m.regionFilterMask(conn, maxprop=max_mask_prop)
			elif option.o1 == "gc":
				print("\t\t\tFiltering criterion: Between",option.o2,"and",option.o3,"GC content")
				min_mask_prop = option.o2
				max_mask_prop = option.o3
				m.regionFilterGC(conn, minprop=min_mask_prop, maxprop=max_mask_prop)
			elif option.o1 == "len":
				print("\t\t\tFiltering criterion: Target length between",option.o2,"and",option.o3)
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
					print("\t\t\tWARNING: No blast database provided. Skipping <--filter_r> option %s"%option.o1)
					break
				#print("BLASTDB PATH IS: ", db_path)
				#Get targets, print to fasta
				seqs = m.getPassedTRs(conn)
				fas = params.workdir + "/.temp.fasta"
				aln_file_tools.writeFasta(seqs, fas)
				outfile = params.workdir + "/.temp.blast"
				if option.o1 == "blast_x":
					print("\t\t\tFiltering criterion: BLAST exclusion")
					print("\t\t\t  --blastn path:",params.blastn)
					print("\t\t\t  --Database:",db_path)
					print("\t\t\t  --Percent identity:",option.o2)
					print("\t\t\t  --Query coverage:",option.o3)
					print("\t\t\t  --N threads:",params.threads)
					print("\t\t\t  --Word size",params.word_size)
					print("\t\t\t  --Gap open penalty:",params.gapopen)
					print("\t\t\t  --Gap extend:",params.gapextend)
					print("\t\t\t  --E-value cutoff:",params.evalue)
					print("\t\t\t  --Maximum returned hits:",params.max_hits)
					if (params.nodust):
						print("\t\t\t  --DUST: False")
					else:
						print("\t\t\t  --DUST: True")
					print("\t\t\t  --Method:",params.blast_method)
					print("\t\t\t  --Method:",params.blast_method)
					local_db_path = db_path
					if params._bx_db or params._bx_fdb:
						if params._bx_db:
							print("\t\t\t  --OVERRIDING BLAST DB:",params._bx_db)
							local_db_path = params._bx_db
						elif params._bx_fdb:
							print("\t\t\t  --OVERRIDING BLAST DB:",params._bx_fdb)
							local_db_path = params.workdir + "/blastdb/" + params.out
							b.makeblastdb(params.makedb, params._bx_fdb, local_db_path)
					blacklist = b.blastExcludeMatch(params, local_db_path, fas, option.o2, option.o3, outfile)
					if len(blacklist) > 1:
						m.removeRegionsByList(conn, blacklist)
				elif option.o1 == "blast_i":
					print("\t\t\tFiltering criterion: BLAST inclusion")
					print("\t\t\t  --blastn path:",params.blastn)
					print("\t\t\t  --Database:",db_path)
					print("\t\t\t  --Percent identity:",option.o2)
					print("\t\t\t  --Query coverage:",option.o3)
					print("\t\t\t  --N threads:",params.threads)
					print("\t\t\t  --Word size",params.word_size)
					print("\t\t\t  --Gap open penalty:",params.gapopen)
					print("\t\t\t  --Gap extend:",params.gapextend)
					print("\t\t\t  --E-value cutoff:",params.evalue)
					print("\t\t\t  --Maximum returned hits:",params.max_hits)
					if (params.nodust):
						print("\t\t\t  --DUST: False")
					else:
						print("\t\t\t  --DUST: True")
					print("\t\t\t  --Method:",params.blast_method)
					local_db_path = db_path
					if params._bi_db or params._bi_fdb:
						if params._bi_db:
							print("\t\t\t  --OVERRIDING BLAST DB:",params._bi_db)
							local_db_path = params._bi_db
						elif params._bi_fdb:
							print("\t\t\t  --OVERRIDING BLAST DB:",params._bi_fdb)
							local_db_path = params.workdir + "/blastdb/" + params.out
							b.makeblastdb(params.makedb, params._bi_fdb, local_db_path)
					whitelist = b.blastIncludeMatch(params, local_db_path, fas, option.o2, option.o3, outfile)
					if len(whitelist) > 1:
						m.removeRegionsByWhitelist(conn, whitelist)
				elif option.o1 == "blast_a":
					print("\t\t\tFiltering criterion: BLAST ambiguous map exclusion")
					print("\t\t\t  --blastn path:",params.blastn)
					print("\t\t\t  --Database:",db_path)
					print("\t\t\t  --Percent identity:",option.o2)
					print("\t\t\t  --Query coverage:",option.o3)
					print("\t\t\t  --N threads:",params.threads)
					print("\t\t\t  --Word size",params.word_size)
					print("\t\t\t  --Gap open penalty:",params.gapopen)
					print("\t\t\t  --Gap extend:",params.gapextend)
					print("\t\t\t  --E-value cutoff:",params.evalue)
					print("\t\t\t  --Maximum returned hits:",params.max_hits)
					if (params.nodust):
						print("\t\t\t  --DUST: False")
					else:
						print("\t\t\t  --DUST: True")
					print("\t\t\t  --Method:",params.blast_method)
					local_db_path = db_path
					if params._ba_db or params._ba_fdb:
						if params._ba_db:
							print("\t\t\t  --OVERRIDING BLAST DB:",params._ba_db)
							local_db_path = params._ba_db
						elif params._ba_fdb:
							print("\t\t\t  --OVERRIDING BLAST DB:",params._ba_fdb)
							local_db_path = params.workdir + "/blastdb/" + params.out
							b.makeblastdb(params.makedb, params._ba_fdb, local_db_path)
					blacklist = b.blastExcludeAmbig(params, local_db_path, fas, option.o2, option.o3, outfile)
					if len(blacklist) > 1:
						m.removeRegionsByList(conn, blacklist)
				os.remove(fas)
				os.remove(outfile)
				#sys.exit()
			elif option.o1 in ("gff", "gff_a"):
				if params.gff and params.assembly:
					if option.o1 == "gff":
						print("\t\t\tFiltering criterion: Proximity to",option.o2,"GFF elements")
						m.regionFilterGFF(conn, option.o2, params.flank_dist)
					elif option.o1 == "gff_a":
						print("\t\t\tFiltering criterion: Proximity to",option.o2,"GFF aliases")
						m.regionFilterGFF_Alias(conn, option.o2, params.flank_dist)
				else:
					sys.exit("ERROR: Filtering targets on proximity to GFF elements requires FASTA <-A> and GFF <-G> inputs!")
			elif option.o1 in ("bed_x", "bed_i"):
				if params.bed and params.assembly:
					if option.o1 == "bed_x":
						print("\t\t\tFiltering criterion: Exclude by proximity to BED elements")
						m.regionFilterBED_exclude(conn, params.flank_dist)
					elif option.o1 == "bed_i":
						print("\t\t\tFiltering criterion: Include by proximity to BED elements")
						m.regionFilterBED_include(conn, params.flank_dist)
				else:
					sys.exit("ERROR: Filtering targets on proximity to BED elements requires FASTA <-A> and BED <-B> inputs!")
			else:
				assert False, "Unhandled option %r"%option

		#Perform pairwise alignment AFTER all other filters because it is analytically more expensive
		#Target region deduplication by pairwise alignment
		if aln:
			print("\t\t\tFiltering criterion: Pairwise alignment")
			passedTargets = m.getPassedTRs(conn)
			assert (0.0 <= minid <= 1.0), "Minimum ID for pairwise alignment must be between 0.0 and 1.0"
			assert (0.0 <= mincov <= 1.0), "Minimum alignment coverage for pairwise alignment must be between 0.0 and 1.0"
			print("\t\t\t  --VSEARCH executable:",params.vsearch)
			print("\t\t\t  --VSEARCH threads:",params.vthreads)
			print("\t\t\t  --Percent identity:",minid)
			print("\t\t\t  --Query coverage:",mincov)
			print("\t\t\t  --Masking method:",params.vsearch_qmask)
			blacklist_edges = pairwiseAlignDedup(conn, params, passedTargets, minid, mincov)
			if (len(blacklist_edges) > 0):
				print("\t\t\tResolving edges...")
				if(params._noGraph):
					print("\t\t\t  --Resolve by graph: False (deleting all conflicts)")
				else:
					print("\t\t\t  --Resolve by graph: True")
					if (params._noWeightGraph):
						print("\t\t\t  --MIS method: Approximate")
					else:
						print("\t\t\t  --MIS method: Weighted")
						if (params._weightByMin):
							print("\t\t\t  --Weights: Minimum ambiguity")
						else:
							print("\t\t\t  --Weights: Maximum variation")
				revised_blacklist = dupEdgeResolution(conn, params, blacklist_edges)
				if len(revised_blacklist) > 0:
					m.removeRegionsByList(conn, revised_blacklist)

		#If 'random' select is turned on, then apply AFTER resolving conflicts (--select_r)
		if rand_num:
			return(rand_num)

#Function to fetch target conflicts
def findTargetConflicts(conn, params):
	if params.mult_reg == 0:
		#print("Multiple regions NOT allowed, apply --select_r within whole loci")
		m.fetchConflictTRs_NoMult(conn)
	else:
		#print("Multiple TRs allowed, apply --select_r within conflict_blocks <--dist_r>")
		#Build conflict_blocks according to --dist_r and --min_mult parameters
		m.fetchConflictTRs(conn, params.min_mult, params.dist_r)
	return(m.getNumConflicts(conn))

#Function to filter target regions by --filter_R arguments
def selectTargetRegions(conn, params):
	#For all alignments over --min_mult:
		#If TRs are within --dist
	#print("Select TR criterion is: ",params.select_r)
	#print("Flank dist is: ", params.flank_dist)
	#print("Minimum mult_reg dist is: ",params.min_mult)

	# #Build conflict tables
	# #print(pd.read_sql_query("SELECT * FROM regions", conn))
	# if params.mult_reg == 0:
	# 	print("Multiple regions NOT allowed, apply --select_r within whole loci")
	# 	#TODO: Need function call to buid conflict_blocks by whole loci
	# 	m.fetchConflictTRs_NoMult(conn)
	# else:
	# 	print("Multiple TRs allowed, apply --select_r within conflict_blocks <--dist_r>")
	# 	#Build conflict_blocks according to --dist_r and --min_mult parameters
	# 	m.fetchConflictTRs(conn, params.min_mult, params.dist_r)

	#NEXT: Need to select TRs within conflict_blocks
	if m.getNumConflicts(conn) > 0:
	#Apply select_r filters for all conflicting TRs
		if params.select_r == "rand":
			#print("--select_r is RANDOM")
			#Do it later
			pass
		elif params.select_r == "snp":
			#Select based on SNPs flanking in "d" dist
			try:
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
	aln_file_tools.writeFastaNogap(seqs, fas)

	#First sort FASTA by size
	sor = params.workdir + "/.temp.sort"
	try:
		vsearch.sortByLength(params.vsearch, fas, sor, params.blen)
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
		vsearch.allpairsGlobal(params.vsearch, params.vthreads, sor, minid, mincov, pw, params.blen, params.vsearch_qmask)
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
	#os.remove(pw)
	return(ret)

#Function to call VSEARCH to mask loci
def fastxMaskLoci(conn, params, seqs):
	"""Seqs must be a pandas DF where cols: 0=index, 1=name, 2=sequence"""
	fas = params.workdir + "/.temp_premask.fasta"
	aln_file_tools.writeFastaNogap(seqs, fas)


	#First sort FASTA by size
	mask = params.workdir + "/.temp_masked.fasta"
	try:
		vsearch.fastxMask(params.vsearch, fas, mask)
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
	#
	# Finally, parse masked fasta for new mask values
	newMask = list()
	for contig in aln_file_tools.read_fasta(mask):
		id = contig[0].replace("id_","")
		m = s.mask_content(contig[1])
		tup = (id, m)
		newMask.append(tup)
	os.remove(mask)
	#print(newMask)
	#os.remove(pw)
	return(newMask)


#Function for deduplication of targets by pairwise alignment
def pairwiseAlignReverseComp(conn, params, seqs, minid, mincov):
	"""Seqs must be a pandas DF where cols: 0=index, 1=name, 2=sequence"""
	#write sequences
	fas = params.workdir + "/.temp_ref.fasta"
	aln_file_tools.writeFastaNogap(seqs, fas)

	# #First sort FASTA by size
	sor = params.workdir + "/.temp.sort"
	try:
		vsearch.sortByLength(params.vsearch, fas, sor, params.blen)
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


	#VSEARCH call to reverse complement .temp_ref.fasta
	# revcomp = params.workdir + "/.temp_revcomp.fasta"
	# try:
	# 	vsearch.fastxRevcomp(params.vsearch, sor, revcomp)
	# except KeyboardInterrupt:
	# 	sys.exit("Process aborted: Keyboard Interrupt")
	# except subprocess.CalledProcessError as err:
	# 	print("VSEARCH encountered a problem.")
	# 	sys.exit(err.args)
	# except NameError as err:
	# 	sys.exit(err.args)
	# except OSError as err:
	# 	print("Exception: OSError in VSEARCH call. Check that you are using the correct executable.")
	# 	sys.exit(err.args)
	# except:
	# 	sys.exit(sys.exc_info()[0])

	#custom function call to reverse complement FASTA
	revcomp = params.workdir + "/.temp_revcomp.fasta"
	aln_file_tools.reverseComplementFasta(sor, revcomp)

	# global align of baits against reverse complements
	pw = params.workdir + "/" + params.out + ".rc"
	try:
		vsearch.usearchGlobal(params.vsearch, params.vthreads, sor, revcomp, minid, mincov, pw, params.blen, params.vsearch_qmask)
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
	os.remove(revcomp)

	#Finally, parse the output of pairwise alignment, to get 'bad matches'
	#Returns a list of edges
	ret = vsearch.parsePairwiseAlign(pw)
	#os.remove(pw)
	#print(ret)
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
			#print("add")
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
	#print("Params.overlap is ", params.overlap)
	#print("Params.bait_shift is", params.bait_shift)
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
		#print("Union length is",union)
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
	elif params.select_b == "calc":
		#calculate the minimum target length to meet requirement
		union = utils.calculateUnionLengthFixed(params.select_b_num, params.blen, params.overlap)
		#print("union:",union)
		#print("max",params.overlap)
		#NOTE: Here select_b_num is a desired amount and bait_shift is a overlap is a MAXIMUM allowable overlap
		for seq in targets.itertuples():
			length = len(seq[2])
			#if target shorter than union length, tile whole thing
			if length < union:
				#print(length,"<",union,"- tiling all")
				baitSlidingWindow(conn, seq[1], seq[2], params.bait_shift, params.blen)
			#otherwise if target is LONGER, calculate new overlap
			else:
				#bait_shift = how far to shift each starting point
				getcontext().prec = 6
				pad_shift = Decimal(length - union) / Decimal(params.select_b_num) #this is the extra space to add to max_overlap
				new_shift = int(params.bait_shift + pad_shift)
				#print("old_shift",params.bait_shift)
				#print("pas_shift",pad_shift)
				#print("new_shift",new_shift)
				baitSlidingWindow(conn, seq[1], seq[2], new_shift, params.blen)

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
def filterBaits_verbose(conn, params):
	rand_num = None
	alnP = False
	minidP = None
	mincovP = None
	alnR = False
	minidR = None
	mincovR = None
	if (len(params.filter_b_objects) <= 0):
		print("\t\t\tNo filtering criteria provided. Retaining all baits.")
		return(0)
	else:
		for option in params.filter_b_objects:
		#	print("Bait filtering option: ", option.o1)
			if option.o1 == "rand":
			#Set 'rand' to TRUE for random selection AFTER other filters
				rand_num = int(option.o2)
				assert rand_num > 0, "Number for random bait selection must be greater than zero!"
			elif option.o1 == "mask":
				print("\t\t\tFiltering criterion: Maximum",option.o2,"proportion masked")
				max_mask_prop = option.o2
				m.baitFilterMask(conn, maxprop=max_mask_prop)
			elif option.o1 == "gc":
				print("\t\t\tFiltering criterion: Between",option.o2,"and",option.o3,"GC content")
				min_mask_prop = option.o2
				max_mask_prop = option.o3
				m.baitFilterGC(conn, minprop=min_mask_prop, maxprop=max_mask_prop)
			elif option.o1 == "pw":
				#print("Filtering baits by pairwise alignment")
				alnP = True
				minidP = option.o2
				mincovP = option.o3
			elif option.o1 == "rc":
				alnR = True
				minidR = option.o2
				mincovR = option.o3
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
				#print("BLASTDB PATH IS: ", db_path)
				#Get targets, print to fasta
				seqs = m.getPassedBaits(conn)
				fas = params.workdir + "/.temp.fasta"
				aln_file_tools.writeFasta(seqs, fas)
				outfile = params.workdir + "/.temp.blast"
				if option.o1 == "blast_x":
					print("\t\t\tFiltering criterion: BLAST exclusion")
					print("\t\t\t  --blastn path:",params.blastn)
					print("\t\t\t  --Database:",db_path)
					print("\t\t\t  --Percent identity:",option.o2)
					print("\t\t\t  --Query coverage:",option.o3)
					print("\t\t\t  --N threads:",params.threads)
					print("\t\t\t  --Word size",params.word_size)
					print("\t\t\t  --Gap open penalty:",params.gapopen)
					print("\t\t\t  --Gap extend:",params.gapextend)
					print("\t\t\t  --E-value cutoff:",params.evalue)
					if (params.nodust):
						print("\t\t\t  --DUST: False")
					else:
						print("\t\t\t  --DUST: True")
					print("\t\t\t  --Method:",params.blast_method)
					local_db_path = db_path
					if params._bx_db or params._bx_fdb:
						if params._bx_db:
							print("\t\t\t  --OVERRIDING BLAST DB:",params._bx_db)
							local_db_path = params._bx_db
						elif params._bx_fdb:
							print("\t\t\t  --OVERRIDING BLAST DB:",params._bx_fdb)
							local_db_path = params.workdir + "/blastdb/" + params.out
							b.makeblastdb(params.makedb, params._bx_fdb, local_db_path)
					blacklist = b.blastExcludeMatch(params, local_db_path, fas, option.o2, option.o3, outfile)
					if len(blacklist) > 1:
						m.removeBaitsByList(conn, blacklist)
				elif option.o1 == "blast_i":
					print("\t\t\tFiltering criterion: BLAST inclusion")
					print("\t\t\t  --blastn path:",params.blastn)
					print("\t\t\t  --Database:",db_path)
					print("\t\t\t  --Percent identity:",option.o2)
					print("\t\t\t  --Query coverage:",option.o3)
					print("\t\t\t  --N threads:",params.threads)
					print("\t\t\t  --Word size",params.word_size)
					print("\t\t\t  --Gap open penalty:",params.gapopen)
					print("\t\t\t  --Gap extend:",params.gapextend)
					print("\t\t\t  --E-value cutoff:",params.evalue)
					if (params.nodust):
						print("\t\t\t  --DUST: False")
					else:
						print("\t\t\t  --DUST: True")
					print("\t\t\t  --Method:",params.blast_method)
					local_db_path = db_path
					if params._bi_db or params._bi_fdb:
						if params._bi_db:
							print("\t\t\t  --OVERRIDING BLAST DB:",params._bi_db)
							local_db_path = params._bi_db
						elif params._bi_fdb:
							print("\t\t\t  --OVERRIDING BLAST DB:",params._bi_fdb)
							local_db_path = params.workdir + "/blastdb/" + params.out
							b.makeblastdb(params.makedb, params._bi_fdb, local_db_path)
					whitelist = b.blastIncludeMatch(params, local_db_path, fas, option.o2, option.o3, outfile)
					if len(whitelist) > 1:
						m.removeBaitsByWhitelist(conn, whitelist)
				elif option.o1 == "blast_a":
					print("\t\t\tFiltering criterion: BLAST ambiguous map exclusion")
					print("\t\t\t  --blastn path:",params.blastn)
					print("\t\t\t  --Database:",db_path)
					print("\t\t\t  --Percent identity:",option.o2)
					print("\t\t\t  --Query coverage:",option.o3)
					print("\t\t\t  --N threads:",params.threads)
					print("\t\t\t  --Word size",params.word_size)
					print("\t\t\t  --Gap open penalty:",params.gapopen)
					print("\t\t\t  --Gap extend:",params.gapextend)
					print("\t\t\t  --E-value cutoff:",params.evalue)
					print("\t\t\t  --Maximum returned hits:",params.max_hits)
					if (params.nodust):
						print("\t\t\t  --DUST: False")
					else:
						print("\t\t\t  --DUST: True")
					print("\t\t\t  --Method:",params.blast_method)
					local_db_path = db_path
					if params._ba_db or params._ba_fdb:
						if params._ba_db:
							print("\t\t\t  --OVERRIDING BLAST DB:",params._ba_db)
							local_db_path = params._ba_db
						elif params._ba_fdb:
							print("\t\t\t  --OVERRIDING BLAST DB:",params._ba_fdb)
							local_db_path = params.workdir + "/blastdb/" + params.out
							b.makeblastdb(params.makedb, params._ba_fdb, local_db_path)
					blacklist = b.blastExcludeAmbig(params, local_db_path, fas, option.o2, option.o3, outfile)
					if len(blacklist) > 1:
						m.removeBaitsByList(conn, blacklist)
				os.remove(fas)
				os.remove(outfile)
			else:
				assert False, "Unhandled option %r"%option

		#Perform pairwise alignment AFTER all other filters because it is analytically more expensive
		if alnP:
			print("\t\t\tFiltering criterion: Pairwise alignment")
			passedBaits = m.getPassedBaits(conn)
			assert (0.0 < minidP < 1.0), "Minimum ID for pairwise alignment must be between 0.0 and 1.0"
			assert (0.0 < mincovP < 1.0), "Minimum alignment coverage for pairwise alignment must be between 0.0 and 1.0"
			print("\t\t\t  --VSEARCH executable:",params.vsearch)
			print("\t\t\t  --VSEARCH threads:",params.vthreads)
			print("\t\t\t  --Percent identity:",minidP)
			print("\t\t\t  --Query coverage:",mincovP)
			print("\t\t\t  --Masking method:",params.vsearch_qmask)

			blacklist_edges = pairwiseAlignDedup(conn, params, passedBaits, minidP, mincovP)

			#print("Blacklisted edges:",blacklist_edges)
			if (len(blacklist_edges) > 0):
				print("\t\t\tResolving edges...")
				if(params._noGraph):
					print("\t\t\t  --Resolve by graph: False (deleting all conflicts)")
				else:
					print("\t\t\t  --Resolve by graph: True")
					if (params._noWeightGraph):
						print("\t\t\t  --MIS method: Approximate")
					else:
						print("\t\t\t  --MIS method: Weighted")
						if (params._weightByMin):
							print("\t\t\t  --Weights: Minimum ambiguity")
						else:
							print("\t\t\t  --Weights: Maximum variation")
				params._noWeightGraph = 1
				revised_blacklist = dupEdgeResolution(conn, params, blacklist_edges)
				if len(revised_blacklist) > 0:
					m.removeBaitsByList(conn, revised_blacklist)

		#Perform reverse complement search
		if alnR:
			print("\t\t\tFiltering criterion: Reverse Complement Pairwise Alignment")
			passedBaits = m.getPassedBaits(conn)
			assert (0.0 < minidR < 1.0), "Minimum ID for pairwise alignment must be between 0.0 and 1.0"
			assert (0.0 < mincovR < 1.0), "Minimum alignment coverage for pairwise alignment must be between 0.0 and 1.0"
			print("\t\t\t  --VSEARCH executable:",params.vsearch)
			print("\t\t\t  --VSEARCH threads:",params.vthreads)
			print("\t\t\t  --Percent identity:",minidR)
			print("\t\t\t  --Query coverage:",mincovR)
			print("\t\t\t  --Masking method:",params.vsearch_qmask)

			blacklist_edges = pairwiseAlignReverseComp(conn, params, passedBaits, minidR, mincovR)

			#print("Blacklisted edges:",blacklist_edges)
			if (len(blacklist_edges) > 0):
				print("\t\t\tResolving edges...")
				if(params._noGraph):
					print("\t\t\t  --Resolve by graph: False (deleting all conflicts)")
				else:
					print("\t\t\t  --Resolve by graph: True")
					if (params._noWeightGraph):
						print("\t\t\t  --MIS method: Approximate")
					else:
						print("\t\t\t  --MIS method: Weighted")
						if (params._weightByMin):
							print("\t\t\t  --Weights: Minimum ambiguity")
						else:
							print("\t\t\t  --Weights: Maximum variation")
				params._noWeightGraph = 1
				revised_blacklist = dupEdgeResolution(conn, params, blacklist_edges)
				if len(revised_blacklist) > 0:
					m.removeBaitsByList(conn, revised_blacklist)

		#If 'random' select is turned on, then apply AFTER all other options
		if rand_num:
			print("\t\t\tRandomly selecting",rand_num,"baits.")
			m.baitFilterRandom(conn, rand_num)

#Function to filter target regions by --filter_R arguments
# def filterBaits(conn, params):
# 	rand_num = None
# 	aln = False
# 	minid = None
# 	mincov = None
# 	if (len(params.filter_b_objects) <= 0):
# 		return(0)
# 	else:
# 		for option in params.filter_b_objects:
# 			#print("Bait filtering option: ", option.o1)
# 			if option.o1 == "rand":
# 			#Set 'rand' to TRUE for random selection AFTER other filters
# 				rand_num = int(option.o2)
# 				assert rand_num > 0, "Number for random bait selection must be greater than zero!"
# 			elif option.o1 == "mask":
# 				max_mask_prop = option.o2
# 				m.baitFilterMask(conn, maxprop=max_mask_prop)
# 			elif option.o1 == "gc":
# 				min_mask_prop = option.o2
# 				max_mask_prop = option.o3
# 				m.baitFilterGC(conn, minprop=min_mask_prop, maxprop=max_mask_prop)
# 			elif option.o1 == "pw":
# 				#print("Filtering baits by pairwise alignment")
# 				aln = True
# 				minid = option.o2
# 				mincov = option.o3
# 			elif option.o1 in ("blast_x", "blast_i"):
# 				#if blast database given as fasta, make a blastdb:
# 				db_path = None
# 				if (params.blastdb):
# 					db_path = params.blastdb
# 				elif (params.fastadb):
# 					db_path = params.workdir + "/blastdb/" + params.out
# 					b.makeblastdb(params.makedb, params.fastadb, db_path)
# 					params.blastdb = db_path
# 				elif(not params.blastdb and not params.fastadb):
# 					print("WARNING: No blast database was provided. Skipping <--filter_r> option %s"%option.o1)
# 					break
# 				#print("BLASTDB PATH IS: ", db_path)
# 				#Get targets, print to fasta
# 				seqs = m.getPassedTRs(conn)
# 				fas = params.workdir + "/.temp.fasta"
# 				aln_file_tools.writeFasta(seqs, fas)
# 				outfile = params.workdir + "/.temp.blast"
# 				if option.o1 == "blast_x":
# 					blacklist = b.blastExcludeMatch(params, db_path, fas, option.o2, option.o3, outfile)
# 					m.removeBaitsByList(conn, blacklist)
# 				elif option.o1 == "blast_i":
# 					whitelist = b.blastIncludeMatch(params, db_path, fas, option.o2, option.o3, outfile)
# 					m.removeBaitsByWhitelist(conn, whitelist)
# 				elif option.o1 == "blast_a":
# 					sys.exit("blast_a option is not yet implemented.")
# 					#blacklist = b.blastExcludeAmbig(params, db_path, fas, option.o2, option.o3, outfile)
# 					#m.removeRegionsByList(conn, blacklist)
# 				os.remove(fas)
# 				os.remove(outfile)
# 			else:
# 				assert False, "Unhandled option %r"%option
#
# 		#Perform pairwise alignment AFTER all other filters because it is analytically more expensive
# 		if aln:
# 			passedBaits = m.getPassedBaits(conn)
# 			assert (0.0 < minid < 1.0), "Minimum ID for pairwise alignment must be between 0.0 and 1.0"
# 			assert (0.0 < mincov < 1.0), "Minimum alignment coverage for pairwise alignment must be between 0.0 and 1.0"
# 			blacklist_edges = pairwiseAlignDedup(conn, params, passedBaits, minid, mincov)
# 			#print("Blacklisted edges:",blacklist_edges)
# 			if (len(blacklist_edges) > 0):
# 				params._noWeightGraph = 1
# 				revised_blacklist = dupEdgeResolution(conn, params, blacklist_edges)
# 				if len(revised_blacklist) > 0:
# 					m.removeBaitsByList(conn, revised_blacklist)
#
# 		#If 'random' select is turned on, then apply AFTER all other options
# 		if rand_num:
# 			#print("Randomly selecting",rand_num,"baits.")
# 			m.baitFilterRandom(conn, rand_num)


#Function to print baits in final output
def printBaits(conn, params):
	df = m.getPrintBaits(conn)
	out = params.workdir + "/" + params.out + ".fasta"
	file_object = open(out, "w")

	for i, r in df.iterrows():
		id=1
		if r.chrom not in ["NULL", "NA"]:
			id=r.chrom
		else:
			id=r.locid

		#If user wants expanded sequences:
		if params.expand:
			var = 1
			#Get fully expanded sequence
			for expanded in s.expandAmbiquousDNA(r.sequence):
				expanded = expanded.replace("-","")
				if params.strand in ("+", "both"):
					header = ">" + str(id) + ":" + str(r.rstart) + "-" + str(r.rstop) +str(r.regid) + ":" + str(r.start) + "-" + str(r.stop) + "_Bait=" +str(r.baitid) + "." + str(var) + "\n"
					seq = expanded + "\n"
					file_object.write(header)
					file_object.write(seq)
				if params.strand in ("-", "both"):
					header = ">" + str(id) + ":" + str(r.rstart) + "-" + str(r.rstop) +str(r.regid) + ":" + str(r.start) + "-" + str(r.stop) + "_Bait=" +str(r.baitid) + "." + str(var) + "_revcomp" + "\n"
					seq = s.reverseComplement(expanded) + "\n"
					file_object.write(header)
					file_object.write(seq)
				var += 1
		#Otherwise just print it
		else:
			if params.strand in ("+", "both"):
				header = ">" + str(id) + ":" + str(r.rstart) + "-" + str(r.rstop) + str(r.regid) + ":" + str(r.start) + "-" + str(r.stop) + "_Bait=" + str(r.baitid) + "\n"
				seq = r.sequence + "\n"
				file_object.write(header)
				file_object.write(seq)
			if params.strand in ("-", "both"):
				header = ">" + str(id) + ":" + str(r.rstart) + "-" + str(r.rstop) +str(r.regid) + ":" + str(r.start) + "-" + str(r.stop) + "_Bait=" +str(r.baitid) + "_revcomp" + "\n"
				seq = s.reverseComplement(r.sequence) + "\n"
				file_object.write(header)
				file_object.write(seq)
	file_object.close()

#Function to print targets in final output
def printTargets(conn, params):
	df = m.getPrintRegions(conn)

	out = params.workdir + "/" + params.out + "_targets.fasta"
	file_object = open(out, "w")

	rel_num = dict() #dictionary to keep target number relative to locus num
	for i, r in df.iterrows():
		id=1
		if r.chrom not in ["NULL", "NA"]:
			id=r.chrom
		else:
			id=r.locid

		#build FASTA header
		p = "T"
		if r["pass"]==0:
			p="F"
		header = ">" + str(id) + ":" + str(r.start) + "-" + str(r.stop) +"_Target=" + str(r.regid) + "_Pass=" + str(p) + "\n"
		seq = r.sequence + "\n"
		file_object.write(header)
		file_object.write(seq)

	file_object.close()

#Function to print loci in final output
def printLoci(conn, params):
	df = m.getLoci(conn)

	out = params.workdir + "/" + params.out + "_catalog.fasta"
	file_object = open(out, "w")

	for i, r in df.iterrows():
		#build FASTA header
		p = "T"
		if r["pass"]==0:
			p="F"
		id=1
		if r.chrom not in ["NULL", "NA"]:
			id=r.chrom
		else:
			id=r.id
		header = ">Locus" + str(id) + "_Pass=" + str(p) + "\n"
		seq = r.consensus + "\n"
		file_object.write(header)
		file_object.write(seq)

	file_object.close()
