#!/usr/bin/python
import sys
import sqlite3
import pandas
import time
from timeit import default_timer as timer
from mrbait import mrbait_menu
from mrbait.mrbait_menu import parseArgs
from mrbait import manage_bait_db as m
from mrbait import mrbait_corefuncs as core
from mrbait import mrbait_corefuncs_parallel as pcore


############################### MAIN ###################################

#TODO: Filter targets and baits by low-complexity (Dusting)
#TODO: Python implementation of SDUST algorithm (gonna be a lot of work, hold off for now)
#TODO: Move some hacker options to full visible options
#TODO: When doing PW align with baits, query targets table to get weights!


def main():
	global_start = timer()
	#Parse Command line arguments
	#Print header information
	mrbait_menu.printHeader()

	print("\n\tParsing command line arguments...")
	params = parseArgs()



	#Intiate database connection
	print ("\tLoading SQLite Database...")
	print ("\t\tEstablishing database connection:", params.db)
	conn = m.create_connection(params.db)

	#Initialize empty databases
	#if conn.empty() or something like that
	if params.resume:
		step = params.resume + 1
	else:
		step = 0
	while step < 6:
		if step == 0:
			#Establishing new database
			start = timer()
			print ("\t\tInitializing empty tables.\n")
			m.init_new_db(conn)
			step = 1
			printTime(start,2)
		elif step == 1:
			start = timer()
			#Loading inputs
			print ("\n\tStep 1: Loading Alignments")
			#Clear database
			if m.getNumPassedLoci(conn) > 0:
				print("\t\tClearing existing records from database")
				m.init_new_db(conn)
			loadAlignments(conn, params)
			#PASS=1 is PASS=FALSE
			
			#optional masking of loci using DUST algorithm
			if params.dustMask:
				print("\t\tMasking consensus sequences using the DUST algorithm in VSEARCH...", end="")
				masked = core.fastxMaskLoci(conn, params, m.getPassedLoci(conn))
				m.updateLociMask(conn, masked)
				print(" Done!\n")

			
			#Pre-filters: Length, alignment depth
			print("\t\tFiltering loci...",end="")
			m.filterLoci(conn, params.minlen, params.cov, params.max_ambig, params.max_mask)
			#print(m.getLoci(conn))
			print(" Done!")
			passedLoci = m.getNumPassedLoci(conn)
			if passedLoci <= 0:
				sys.exit("\nProgram killed: No loci passed filtering.\n")
			else:
				print("\t\t### Results: %s loci passed filtering! ###"%passedLoci)
			if params.print_loc:
				print("\t\tPrinting locus catalog to file...")
				core.printLoci(conn, params)
			step = 2
			printTime(start,2)
		elif step == 2:
			start = timer()
			#Target discovery
			print("\n\tStep 2: Target Discovery")
			#Check that database has loci
			passedLoci = m.getNumPassedLoci(conn)#returns pandas dataframe
			if passedLoci <= 0:
				sys.exit("\nProgram killed: No loci in database.\n")
			else:
				#Check if targets exist
				#If yes, clear them
				numTRs = m.getNumTRs(conn)
				if numTRs > 0:
					print("\t\tWarning: Database already contains targets. Clearing existing records.")
					m.clearBaits(conn)
					m.clearTargets(conn)

				#Target discovery call
				targetDiscovery(conn, params)
				passed = m.getNumPassedTRs(conn)
				if passed <= 0:
					sys.exit("\nProgram killed: No viable targets found.\n")
				else:
					print("\n\t\t### Results: %s potential targets identified! ###"%passed)
				step = 3
				printTime(start,2)
			step = 3
			#print(m.getRegions(conn))
		elif step == 3:
			start = timer()
			print("\n\tStep 3: Target Filtering and Selection")
			#Target filtering and conflict resolution
			passedLoci = m.getNumPassedLoci(conn)#returns pandas dataframe
			passedTargets = m.getNumPassedTRs(conn)
			if passedLoci <= 0:
				sys.exit("\nProgram killed: No loci in database.\n")
			elif passedTargets <= 0:
				sys.exit("\nProgram killed: No targets in database.\n")
			else:
				#Clear baits
				m.clearBaits(conn)
				#reset targets
				m.resetTargets(conn)
				#select: resolve conflicts, apply filters
				selectFilterTargets(conn, params)
			passed = m.getNumPassedTRs(conn)
			if params.print_tr:
				print("\t\tPrinting targets to file...")
				core.printTargets(conn, params)
			if passed <= 0:
				sys.exit("\nProgram killed: No targets passed filtering.\n")
			else:
				print("\n\t\t### Results: %s targets passed filtering! ###"%passed)
			printTime(start,2)
			step = 4
		elif step == 4:
			#Bait discovery
			start = timer()
			print("\n\tStep 4: Bait discovery")
			passedLoci = m.getNumPassedLoci(conn)#returns pandas dataframe
			passedTargets = m.getNumPassedTRs(conn)
			if passedLoci <= 0:
				sys.exit("\nProgram killed: No loci in database.\n")
			elif passedTargets <= 0:
				sys.exit("\nProgram killed: No targets in database.\n")
			else:
				#clear baits
				m.clearBaits(conn)
				baitDiscovery(conn, params)
				passed = m.getNumPassedBaits(conn)
				if passed <= 0:
					sys.exit("\nProgram killed: No baits found.\n")
				else:
					print("\n\t\t### Results: %s potential baits identified! ###"%passed)
			printTime(start,2)
			step = 5
		elif step == 5:
			start = timer()
			print("\n\tStep 5: Bait filtering")
			passedLoci = m.getNumPassedLoci(conn)#returns pandas dataframe
			passedTargets = m.getNumPassedTRs(conn)
			if passedLoci <= 0:
				sys.exit("\nProgram killed: No passed loci in database.\n")
			elif passedTargets <= 0:
				sys.exit("\nProgram killed: No passed targets in database.\n")
			else:
				#Bait filtering
				#reset baits 
				m.resetBaits(conn)
				filterBaits(conn,params)
				passedBaits = m.getNumPassedBaits(conn)
				if passedBaits <= 0:
					sys.exit("\nProgram killed: No baits passed filtering.\n")
				else:
					print("\t\t\t",passedBaits,"passed filtering.")
					#Print baits
					print("\t\tFormatting for printing...")
					formatPrintBaits(conn, params)
					passed = m.getNumPassedBaits(conn)
					if passed <= 0:
						sys.exit("\nProgram killed: No baits passed filtering.\n")
					else:
						print("\n\t\t### Results: %s baits output to file! ###"%passed)
			printTime(start,2)
			step = 6

	print("\n\t=======================================================================")
	out = params.workdir + "/" + params.out + ".fasta"
	print("\t### Final output: %s ###"%out)
	printTime(global_start, 1)
	print()
	conn.close()
	sys.exit()

#Loading inputs
def loadAlignments(conn, params):
	#load alignment to database
	print("\t\tStep 1 parameters:")
	if params.alignment or params.loci or params.xmfa:
		print("\t\t\tMinimum alignment depth (-c):", params.cov)
		print("\t\t\tMinimum alignment length (-l):", params.minlen)
		print("\t\t\tThreshold to call \"N\" or gap consensus base (-q):",params.thresh)
		print("\t\t\tThreshold to retain masked base (-k):",params.mask)
		print("\t\t\tMaximum allowed gap/N bases (-Q):",params.max_ambig)
		print("\t\t\tMaximum allowed masked bases (-K):",params.max_mask)
		#Load alignments
		if params.alignment:
			print("\t\tLoading MAF file:",params.alignment)
			if int(params.threads) > 1:
				print("\t\t\tLoading alignments using",str(params.threads),"parallel processes.")
				pcore.loadMAF_parallel(conn, params)
			else:
				core.loadMAF(conn, params)
		elif params.xmfa:
			print("\t\tLoading XMFA file:",params.xmfa)
			if int(params.threads)>1:
				print("\t\t\tLoading alignments using",str(params.threads),"parallel processes.")
				pcore.loadXMFA_parallel(conn, params)
			else:
				pass
				core.loadXMFA(conn, params)

		elif params.loci:
			print("\t\tLoading LOCI file:",params.loci)
			if int(params.threads) > 1:
				print("\t\t\tLoading alignments using",str(params.threads),"parallel processes.")
				pcore.loadLOCI_parallel(conn, params)
			else:
				core.loadLOCI(conn, params)

	elif params.assembly:
		print("\t\t\tMinimum contig length (-l,--len):", params.minlen)
		print("\t\t\tMaximum allowed gap/N bases (-Q):",params.max_ambig)
		print("\t\t\tMaximum allowed masked bases (-K):",params.max_mask)
		#Load assembly file
		print("\t\tLoading FASTA file:",params.assembly)
		core.loadFASTA(conn, params)

		#If VCF file
		if params.vcf:
			print("\t\tLoading VCF file:",params.vcf)
			# if int(params.threads) > 1:
			# 	print("\t\t\tLoading VCF records using",str(params.threads),"parallel processes.")
			# 	pcore.loadVCF_parallel(conn, params)
			# else:
			if (params.vcfALT):
				print("\t\t\tAttempting to override N/gap alleles using ALT from VCF (--vcfALT=True)\n")
			else:
				print("\t\t\tRetaining N/gap alleles from FASTA reference (--vcfALT=False)\n")
			core.loadVCF(conn, params)

		#if GFF file
		if params.gff:
			print("\t\tLoading GFF file:",params.gff)
			core.loadGFF(conn, params)
			#print(m.getGFF(conn))
	else:
		#Option to load .loci alignment goes here!
		sys.exit("No input files provided.")


#Function to call targetDiscoverySlidingWindow and print relevant params
def targetDiscovery(conn, params):
	print("\t\tStep 2 parameters:")
	print("\t\t\tSliding window width (-b, --bait):", params.win_width)
	print("\t\t\tSliding window shift distance (-w, --win_shift):", params.win_shift)
	print("\t\t\tVariable columns allowed per target (-v, --var_max):", params.var_max)
	print("\t\t\tAmbiguities (N's) allowed per target (-n, --numN):", params.numN)
	print("\t\t\tGap characters allowed per target (-g, --numG):", params.numG)
	print("\t\t\tFlanking distance to parse (-d,--flank_dist):",params.flank_dist)

	#Fetch passed loci
	passedLoci = m.getPassedLoci(conn)
	numPassedLoci = m.getNumPassedLoci(conn)
	#sliding window call
	if params.target_all:
		#NOTE: Target criteria still apply (e.g. number of allowed Ns, gaps, vars)
		print("\t\tSkipping target discovery and applying target filters to",numPassedLoci,"full-length loci...")
	else:
		print("\t\tStarting sliding window target discovery of",numPassedLoci,"loci...")

	if int(params.threads) > 1:
		print("\t\t\tFinding targets using",str(params.threads),"parallel processes...")
		pcore.targetDiscoverySlidingWindow_parallel(conn, params, passedLoci)
	else:
		core.targetDiscoverySlidingWindow(conn, params, passedLoci)

#Function to print params and make calls for target region selection
def selectFilterTargets(conn, params):
	print("\t\tStep 3 parameters:")
	if(params.mult_reg):
		print("\t\t\tMultiple targets allowed per locus (-R, --mult_reg): True")
		print("\t\t\tMinimum locus length for multiple targets (-m, --min_mult):",params.min_mult)
		print("\t\t\tMinimum distance between targets within locus (-D, --dist_r):",params.dist_r)
	else:
		print("\t\t\tMultiple targets allowed per locus (-R, --mult_reg): False")
	print("\t\t\tFlanking distance for target selection (-d,--flank_dist):",params.flank_dist)

	rand = 0
	passed=m.getNumPassedTRs(conn)
	if (passed > 0):
		print("\t\tApplying filters to",passed,"targets.")
		rand = core.filterTargetRegions_verbose(conn, params)
		core.checkTargetRegions(conn)
	else:
		sys.exit("Program killed: No targets in database.")


	#Build conflict tables
	conflicts = core.findTargetConflicts(conn, params)
	if (conflicts > 0):
		print("\t\tResolving conflicting targets... Found",conflicts,"needing resolution.")
		print("\t\t\tSelection criterion (-S, --select_r): ",end="")
		if params.select_r == "rand":
			print("Random selection")
		else:
			if params.select_r == "snp":
				print("Most SNPs w/in",params.flank_dist,"bases")
			elif params.select_r == "bad":
				print("Least ambiguities w/in",params.flank_dist,"bases")
			elif params.select_r == "cons":
				print("Most conserved w/in",params.flank_dist,"bases")
			print("\t\t\tNote: Remaining conflicts will be randomly resolved.")
		core.selectTargetRegions(conn, params)
	else:
		print("\t\tNo conflicting targets found. Skipping target selection.")

	#If RANDOM filter for TRs was applied, apply AFTER TR conflict resolution
	passed=m.getNumPassedTRs(conn)
	if rand:
		if rand > passed:
			print("\t\tYou chose to randomly keep",rand,"targets but only",passed,"remain.")
		else:
			print("\t\tRandomly selecting",rand,"targets from",passed,"passing filtering")
			m.regionFilterRandom(conn, rand)

#function to print and call core functions for probe development
def baitDiscovery(conn, params):
	print("\t\tStep 4 parameters:")
	print("\t\t\tBait length (-b, --blen):",params.blen)

	if (params.select_b == "tile"):
		print("\t\t\tTiling baits with",params.overlap,"base overlap")
	elif (params.select_b == "center"):
		print("\t\t\tCentering",params.select_b_num,"baits with",params.overlap,"overlap")
	elif(params.select_b == "flank"):
		print("\t\t\tDesigning",params.select_b_num,"terminal baits with",params.overlap,"overlap")
	elif(params.select_b=="calc"):
		print("\t\t\tDesigning",params.select_b_num,"baits per target with",params.overlap,"maximum overlap")
	#core function call
	core.baitDiscovery(conn, params, m.getPassedTRs(conn))

#Function to filter baits
def filterBaits(conn, params):
	passed=m.getNumPassedBaits(conn)
	if (passed > 0):
		print("\t\tFound",passed,"potential baits!")
		print("\t\tApplying filters to",passed,"baits...")
		core.filterBaits_verbose(conn, params)
	else:
		sys.exit("Program killed: No baits in database.")

#Print options and core call for printBaits
def formatPrintBaits(conn, params):
	if(params.expand):
		print("\t\t\tExpanding all ambiguities (-X, --expand)")
	if (params.strand == "-"):
		print("\t\t\tOnly reporting reverse complement (--strand -)")
	elif (params.strand == "both"):
		print("\t\t\tReporting bait and reverse complement (--strand both)")
	out = params.workdir + "/" + params.out + ".fasta"
	core.printBaits(conn, params)
	print("\t\tFinal baits output to:",out)

#Function to print runtime given a start time
def printTime(start, tabs):
	t = (timer() - start)
	m, s = divmod(t, 60)
	h, m = divmod(m, 60)
	out = ""
	out = "".join(["\t" for i in range(0,tabs)])
	print("%s### Runtime: %d:%02d:%02d (%2f seconds) ###" % (out,h, m, s, t))

#Function to print runtime given a start time
def printTimeClean(start, tabs):
	t = (timer() - start)
	m, s = divmod(t, 60)
	h, m = divmod(m, 60)
	out = ""
	out = "".join(["\t" for i in range(0,tabs)])
	print("%sRuntime: %d:%02d:%02d (%2f seconds)" % (out,h, m, s, t))


#Call main function
if __name__ == '__main__':
	try:
		main()
	except KeyboardInterrupt:
		sys.exit(1)
