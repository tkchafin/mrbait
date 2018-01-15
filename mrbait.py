import sys
import sqlite3
import pandas
import time
from mrbait_menu import parseArgs
import manage_bait_db as m
import mrbait_corefuncs as core


############################### MAIN ###################################

#BELOW IS WORKFLOW FOR UCE DESIGN, FINISH AND THEN CONVERT TO FUNCTIONS
#ADD GFF FUNCTIONALITY LATER
#Add multithreading support later... Each thread will need its own db conn
#If TR too short, can add option to include variable flanking regions?
#TODO: Add better checking to make sure database isn't empty before proceeding (e.g. if filters too stringent)
#TODO: Add "flow control" options, e.g. only make db, load previous db, only TR, etc
#TODO: Add BRANCHING option for Flow Control. e.g. Branch RUN1 from Regions level to try new method of bait design, etc
#TODO: Parallelize loadLOCI and loadMAF, implement ZDZ's chunking functions
#TODO: Parallelize TR discovery somehow??
#TODO: Filter by flanking masked, and flanking GFF elements
#TODO: Functions to read and write params as JSON? Or maybe just stick in database?
#TODO: Filter targets and baits by low-complexity (Dusting)
#TODO: Option to apply SDUST to all loci, and pass/fail loci/alignments on complexity
#TODO: Python implementation of SDUST algorithm (gonna be a lot of work, hold off for now)
#TODO: Option to output baits as +, -, or both strands
#TODO: Runtime bottleneck profiling: cProfile + pstats or prun (??)
#TODO: Memory usage profiling: Check out mprof, looks easy, maybe look at guppy
#TODO: Progress bar for loadXXXX() functions
#TODO: Implement file chunkers and parallel loading
#NOTE: Not all database drivers support the "?" syntax I use for passing params to SQLite. Be careful.


def main():
	global_start = time.clock()
	#Parse Command line arguments
	print("\nParsing command line arguments...")
	params = parseArgs()

	#Print header information
	core.printHeader(params)

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
	while step < 5:
		if step == 0:
			#Establishing new database
			start = time.clock()
			print ("\t\tInitializing empty tables.")
			m.init_new_db(conn)
			step = 1
			printTime(start,2)
		elif step == 1:
			start = time.clock()
			#Loading inputs
			print ("\n\tStep 1: Loading Alignments")
			#Clear database
			if m.getNumPassedLoci(conn) > 0:
				print("\t\tClearing existing records from database")
				m.init_new_db(conn)
			loadAlignments(conn, params)
			#PASS=1 is PASS=FALSE
			#Pre-filters: Length, alignment depth
			print("\t\tFiltering loci...")
			m.filterLoci(conn, params.minlen, params.cov, params.max_ambig, params.max_mask)
			passedLoci = m.getNumPassedLoci(conn)#returns pandas dataframe
			if passedLoci <= 0:
				sys.exit("\nProgram killed: No loci passed filtering.\n")
			else:
				print("\t\t## %s loci passed filtering! ##"%passedLoci)
			step = 2
			printTime(start,2)
		elif step == 2:
			#Target discovery
			print("\n\tStep 2: Target Discovery")
			#Check that database has loci
			#Clear targets and baits
			step = 3
		elif step == 3:
			#Target filtering and conflict resolution
			#Clear baits
			step = 4
		elif step == 4:
			#Bait discovery
			#clear baits
			step = 5
		elif step == 5:
			#Bait filtering
			#Print baits
			#DONE
			step = 6

	print("\n\t=======================================================================")
	print("\tDone!")
	printTime(global_start, 1)
	print()
	conn.close()
	sys.exit()

#Loading inputs
def loadAlignments(conn, params):
	#load alignment to database
	if params.alignment or params.loci:
		print("\t\tStep 1 parameters:")
		print("\t\t\t--Minimum alignment depth (-c):", params.cov)
		print("\t\t\t--Minimum alignment length (-l):", params.minlen)
		print("\t\t\t--Threshold to call \"N\" or gap consensus base (-q):",params.thresh)
		print("\t\t\t--Threshold to retain masked base (-k):",params.mask)
		print("\t\t\t--Maximum allowed gap/N bases (-Q):",params.max_ambig)
		print("\t\t\t--Maximum allowed masked bases (-K):",params.max_mask)
		#Load alignments
		if params.alignment:
			print("\t\tLoading MAF file:",params.loci)
			core.loadMAF(conn, params)
		elif params.loci:
			print("\t\tLoading LOCI file:",params.loci)
			core.loadLOCI_parallel(conn, params)
	elif params.assembly:
		print("\t\tStep 1 params:")
		print("\t\t\t--Minimum contig length (-l,--len):", params.minlen)
		print("\t\t\t--Maximum allowed gap/N bases (-Q):",params.max_ambig)
		print("\t\t\t--Maximum allowed masked bases (-K):",params.max_mask)
		#Load assembly file
		print("\t\tLoading FASTA file:",params.assembly)
		core.loadFASTA(conn, params)

		#If VCF file
		if params.vcf:
			print("\t\tLoading VCF file:",params.vcf)
			core.loadVCF(conn, params)

		#if GFF file
		if params.gff:
			print("\t\tLoading GFF file:",params.gff)
			core.loadGFF(conn, params)
			print(m.getGFF(conn))
	else:
		#Option to load .loci alignment goes here!
		sys.exit("No input files provided.")


#Function to print runtime given a start time
def printTime(start, tabs):
	t = (time.clock() - start)
	m, s = divmod(t, 60)
	h, m = divmod(m, 60)
	out = ""
	out = "".join(["\t" for i in range(0,tabs)])
	print("%sRuntime: %d:%02d:%02d (%2f seconds)" % (out,h, m, s, t))


#Call main function
if __name__ == '__main__':
    main()






#First-pass bait design on loci passing pre-filters


#TODO: Could add an option here to extract GFF elements, and ONLY design baits for those?

#Target region discovery according to params set
print("Starting Target Region discovery...")
core.targetDiscoverySlidingWindow(conn, params, passedLoci)

#Assert that there are TRs chosen, and that not all have been filtered out
core.checkTargetRegions(conn)

#Filter target regions
#If multiple regions NOT allowed, need to choose which to keep
print("Starting Target Region selection...")

#First pass filtering of target regions
rand = core.filterTargetRegions(conn, params)

#Check again that not all have been filtered out
core.checkTargetRegions(conn)

#Apply --select_r filters to sort any conflicting TRs
core.selectTargetRegions(conn, params)

#If RANDOM filter for TRs was applied, apply if AFTER TR conflict resolution
if rand:
	print("Randomly filtering all TRs")
	m.regionFilterRandom(conn, rand)

numTR = m.getNumPassedTRs(conn)
print("Finished Target Region selection... %r targets passed filtering."%numTR)

#Bait discovery
print("Starting probe design...")

#Check again that not all have been filtered out
core.checkTargetRegions(conn)

#Sliding window bait discovery
if m.getNumPassedTRs(conn) <= 0:
	sys.exit("Program killed: No targets passed filtering.")
core.baitDiscovery(conn, params, m.getPassedTRs(conn))

#Bait Filtering
if m.getNumPassedBaits(conn) > 0:
	print("Filtering probe sequences")
	core.filterBaits(conn, params)

#Create final outputs
if m.getNumPassedBaits(conn) > 0:
	core.printBaits(conn, params)
else:
	sys.exit("Program killed: No baits passed filtering.")

print("\n\nProgram ending...Here are some results\n\n")

print(m.getLoci(conn))
#print(m.getVariants(conn))
print(m.getRegions(conn))
print(m.getBaits(conn))

conn.commit()
conn.close()
