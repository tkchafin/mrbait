#!/usr/bin/python

import getopt
import sys

def display_help(message=None):
    if message is not None:
        print (message)
    print ("\nMrBait is a program for developing baits for targeted enrichment NGS methods\n")
    print ("Author:		Tyler K. Chafin")
    print ("Contact:	tkchafin@uark.edu")
    print ("Version:	alpha 0.01")
    print ("\nUsage: ", sys.argv[0], "-m=</path/to/MAF> <-c><-l><-t>")
    print ("\nNOTE: SQLite databases are currently remade from scratch with each", \
    "instance of the program. Support for opening and editing older database files will come soon.")

    print ("\nMandatory arguments:")
    print ("	-i, --in	: Input format- currently only MAF supported [More to come later]")
    print ("	-a, --alignment	: Input multiple alignment [More input methods to come later]")
    
    print ("\nLocus filtering/ consensus options:")
    print ("	-c,--cov	: Minimum number of sequences in each alignment [def=1]")
    print ("	-l,--len	: Minimum alignment length [default=1]")
    print ("	-t,--thresh	: Min. proportion for a gap/N for inclusion in consensus sequence [def=0.1]")
    
    print ("\nBait design options:")
    print ("	-b,--bait	: Desired bait length [default=80]")
    print ("	-m,--mult	: Allow multiple baits per locus [default=false]")
    print ("		 	  NOTE: Use <-M> to use conditional on fragment length") 
    print ("	-r,--maxreg	: Maximum length of target region to design baits [default=1000]")
    print ("	-n,--numN	: Number of consensus N bases allowed in a bait [default=0]")
    print ("		 	  NOTE: These will be inserted as an \"N\" in bait sequence unless <-N>")
    print ("	-N,--callN	: For any Ns in bait, call as majority nucleotide if possible")
    print ("	-g,--numG	: Number of gaps allowed in a bait")
    print ("	-G,--callG	: For gaps in bait, call as majority nucleotide if possible")
    print ("		 	  NOTE: Default is to insert gaps as \"N\" unless <-G> is used")
    print ("	-w,--win_shift	: Shift distance for sliding windows")
    print ("		 	  NOTE: Default is 1/4 of bait length")  
    print ("	-R,--reg_type	: Constrain bait design to within element type of provided GFF")
    print ("		 	  NOTE: COMING SOON, NOT WORKING YET") 
      
    print ("\nBait optimization options (related to <-m>):")
    print ("	-M,--min_mult	: Minimum alignment length to allow multiple baits")
    print ("		 	  NOTE: Default is 4 * bait length when <-m=TRUE>")  
    print ("	-D,--reg_dist	: Minimum distance between disjunct target regions")
    print ("		 	  NOTE: See manual. [Default=100]") 
    print ("	-O,--overlap	: Overlap for tiled baits within target regions [default=40; asserts -T]")
    print ("	-T,--tile	: Bait tiling for regions >N bases [default=160]")
    print ("		 	  NOTE: Recommended for long alignments")
    print ("	-s,--select	: Which criterion is used to select when multiple bait regions w/in <-D>")
    print ("		 	  Options:")
    print ("				-\"s/snp\": Most SNPs within <-d> bases")
    print ("				-\"r/range\": Between <u> and <l> SNPs within <-d> bases")
    print ("				-\"b/bad\": Least Ns and gaps within <-d> bases")
    print ("				-\"c/cons\": Most conserved within <-d> bases")
    print ("				-\"m/min\": Minimum non-conserved (snp, gap, N) in bait region")
    print ("		 	  NOTE: Default behavior is to select first passing bait region per locus")
    print ("	-T,--tile	: Bait tiling for regions >N bases [default=160]")
    
    print ("\nBait region selection (for use with <-s>):")
    print ("	-C,--force_s	: Enforce <-s> criterion as a filter for all bait regions")
    print ("		 	  NOTE: When False, <-s> only chooses best of conflicting bait regions")
    print ("	-d,--opt_dist	: Number of flanking bases to consider when choosing among TRs")
    print ("	-u,--upper_snp	: Maximum SNPs when <-s=range>")
    print ("	-l,--lower_snp	: Minimum SNPs when <-s=range>")
   
    
    print ("\nOutput options:")
    print ("	-e,--expand	: In output bait table, expand all ambiguities")
    print ("	-o,--out	: Prefix for outputs [default=\"./mrbait\"]")
    
    print ("\nGeneral options:")
    print ("	-h,--help	: Displays this help menu [boolean]")
    print()

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try: 
			options, remainder = getopt.getopt(sys.argv[1:], 't:l:b:c:a:h', \
			['thresh=','len=','bait=','cov=','alignment=','help'])
		except getopt.GetoptError as err:
			print(err)
			display_help("Exiting because getopt returned non-zero exit status.")
			sys.exit(2)
		#Default values for params
		call_help=0
		self.alignment=None
		self.cov=1
		self.blen=80
		self.minlen=80
		self.thresh=0.1
		self.ploidy=2
		self.db="./mrbait.sqlite"
				#Parse through arguments and set params 
		for opt, arg in options:
			if opt in ('-a', '--alignment'):
				self.alignment = arg
			elif opt in ('-h', '--help'):
				call_help=1
			elif opt in ('-c', '--cov'):
				self.cov = int(arg)
			elif opt in ('-l', '--len'):
				self.minlen = int(arg)
			elif opt in ('-t', '--thresh'):
				self.thresh = float(arg)
			elif opt in ('-b', '--bait'):
				self.blen = int(arg)
			else: 
				assert False, "unhandled option"
		#Assertions and conditional changes to params
		#Assert that win_shift cannot be larger than blen
		if self.blen > self.minlen:
			self.minlen = self.blen
		#Call help menu if prompted
		if call_help is 1:
			display_help("Exiting because help menu was called.")
			sys.exit(0)
		
