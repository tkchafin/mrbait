#!/usr/bin/python

import sys
import sqlite3
import pandas
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
#TODO: For whole genome option, need to read an mpileup or similar to capture variant information.
#TODO: Could also set minimum coverage thresholds for whole genome?
#TODO: Some form of duplicate screening. Screen targets for dupe or screen baits? Not sure.
#TODO: Filter targets by distance to GFF element
#TODO: Option to constrain targets to ONLY in GFF elements
#TODO: Parallelize loadLOCI and loadMAF, implement ZDZ's chunking functions
#TODO: Parallelize TR discovery somehow??
#TODO: Filter by flanking masked, and flanking GFF elements
#TODO: Functions to read and write params as JSON? Or maybe just stick in database?
#TODO: Option to print LOCI, TARGETS to fasta files?
#TODO: Keep vsearch and blast outfiles, remain to desired paths.
#TODO: Filter targets and baits by low-complexity (Dusting)
#TODO: Python implementation of SDUST algorithm

#Parse Command line arguments
params = parseArgs()

#Intiate database connection
conn = m.create_connection(params.db)
c = conn.cursor()

#Initialize empty databases
#if conn.empty() or something like that
m.init_new_db(conn)

#load alignment to database
if params.alignment:
	print("Loading MAF file:",params.alignment)
	core.loadMAF(conn, params)
elif params.loci:
	print("Loading LOCI file:",params.loci)
	core.loadLOCI(conn, params)
else:
	#Option to load .loci alignment goes here!
	print("No alignment input found. .fasta, .gff, and .phylip support not added yet!")

#First-pass bait design on loci passing pre-filters
#PASS=1 is PASS=FALSE
#Pre-filters: Length, alignment depth
m.filterLoci(conn, params.minlen, params.cov)
passedLoci = m.getPassedLoci(conn) #returns pandas dataframe
if passedLoci.shape[0] <= 0:
	sys.exit("Program killed: No loci passed filtering.")

#Target region discovery according to params set
core.targetDiscoverySlidingWindow(conn, params, passedLoci)

#Assert that there are TRs chosen, and that not all have been filtered out
#core.checkTargetRegions(conn)

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
passedTargets = m.getPassedTRs(conn)
if passedTargets.shape[0] <= 0:
	sys.exit("Program killed: No targets passed filtering.")
core.baitDiscovery(conn, params, passedTargets)

#print

print("\n\nProgram ending...Here are some results\n\n")

print(m.getLoci(conn))
print(m.getRegions(conn))
print(m.getBaits(conn))

conn.commit()
conn.close()
