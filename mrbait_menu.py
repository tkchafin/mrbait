#!/usr/bin/python

import getopt
import sys

def display_help(message=None):
    if message is not None:
        print (message)
    print ("\nMrBait is a program for developing baits for targeted enrichment NGS methods\n")
    print ("Author:		Tyler K. Chafin")
    print ("Contact:	tkchafin@uark.edu")
    print ("Version:	0.01")
    print ("\nUsage: ", sys.argv[0], "-m=</path/to/MAF>")
    print ("\nNOTE: SQLite databases are currently remade from scratch with each", \
    "instance of the program. Support for opening and editing older database files will come soon.")
    print ("\nNOTE: Loci failing the <--cov> and <--len> filters are not included in ",\
    "the database and thus cannot be salvaged, even if subsequent processing of the", \
    "database applies less stringent filtering. Changes in these parameters will require", \
    "re-running the program")
    print ("\nMandatory arguments:")
    print ("	-i, --in	: Input format- currently only MAF supported [More input methods to come later]")
    print ("	-m, --maf	: Input multiple alignment [More input methods to come later]")
    print ("\nOptional arguments:")
    print ("	-p,--ploidy	: PLoidy of your organism [def=2]")
    print ("	-t,--thresh	: Min. proportion for a variant to override the consensus sequence [def=0.1]")
    print ("	-b,--bait	: Desired bait length [default=80]")
    print ("	-c,--cov	: Minimum number of sequences in each alignment [def=1]")
    print ("	-l,--len	: Minimum alignment length [default=1]")
    print ("	-d,--db		: Output database name and path [default=./mrbait.sqlite]")
    print ("	-h,--help	: Displays this help menu [boolean]")
    print()

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try: 
			options, remainder = getopt.getopt(sys.argv[1:], 'p:t:d:l:b:c:m:h', \
			['ploidy=','thresh=','db=','len=','bait=','cov=','maf=','help'])
		except getopt.GetoptError as err:
			print(err)
			display_help("Exiting because getopt returned non-zero exit status.")
			sys.exit(2)
		#Default values for vars
		call_help=0
		self.maf=None
		self.cov=1
		self.blen=80
		self.minlen=80
		self.thresh=0.1
		self.ploidy=2
		self.db="./mrbait.sqlite"
		for opt, arg in options:
			if opt in ('-m', '--maf'):
				self.maf = arg
			elif opt in ('-h', '--help'):
				call_help=1
			elif opt in ('-c', '--cov'):
				self.cov = int(arg)
			elif opt in ('-l', '--len'):
				self.minlen = int(arg)
			elif opt in ('-d', '--db'):
				self.db = arg
			elif opt in ('-t', '--thresh'):
				self.thresh = float(arg)
			elif opt in ('-p', '--ploidy'):
				self.ploidy = int(arg)
			elif opt in ('-b', '--bait'):
				self.blen = int(arg)
			else: 
				assert False, "unhandled option"
		if self.blen > self.minlen:
			self.minlen = self.blen
		#Call help menu if prompted
		if call_help is 1:
			display_help("Exiting because help menu was called.")
			sys.exit(0)
		
