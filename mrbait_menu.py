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

	print("""
Input options:
	-a,--maf	: Input multiple alignment MAF file [more input options to come later]
	-e,--gff	: GFF file containing annotation information [NOT WORKING YET]
	-L,--loci	: For RAD-data, input locus alignment formatted as the \".loci\" output of pyRAD
	-A,--assembly	: Input whole genome assembly as FASTA [NOT WORKING YET]""")
	
	print("""
Locus filtering/ consensus options:
	-c,--cov	: Minimum number of sequences per alignment, for MAF or LOCI input [1] 
	-l,--len	: Minimum alignment length to attempt bait design [80]
	-t,--thresh	: Minimum proportion for gap/N to include in consensus [0.1]""")

	print("""
General Bait Design options:
	-b, --bait	: Bait length [default=80]
	-w, --win_shift	: Shift distance for sliding window discovery of bait regions/baits [1]
			  --Increasing will speed up program but could lower accuracy
			  --Cannot be longer than bait length <-b,--bait>
	-R,--mult_reg	: Allow multiple target regions per locus [false]
			  --See \"Target Region options\" below
	-m,--min_mult	: Minimum alignment length to allow multiple target region [1000]
			  --By default is set to the value of <-l,--len>. Only applies when <-R>
			  --When -R and -T are turned off, only one bait is selected per locus
	-v,--var_max	: Maximum number of SNP columns to be included in a bait [0]
			  --These will be coded as the appropriate IUPAC ambiguity code
			  --Can be expanded in final output using the <--expand> flag
	-n,--numN	: Number of consensus Ns allowed in a bait [0]
			  --These will be inserted as \"N\" unless <-N> is used
	-N,--callN	: For Ns in bait sequence, call as majority nucleotide if possible
	-g,--numG	: Number of consensus gaps (\"-\") allowed in bait [0]
			  --These will be inserted as \"N\" unless <-G> is used
	-G,--callG	: For gaps in bait, call as majority nucleotide if possible
	-E,--gff_type	: Constrain bait design to within element type of provided GFF
			  --COMING SOON, NOT WORKING YET""")

	print("""
Target Region options:
	-T,--tiling	: Turn on tiling to design multiple baits per bait region [false]
			  --When tiling is off, see Bait Selection options below for how baits are sampled
	-O,--overlap	: Overlap for tiled baits within target regions [40]
			  --By default will be set to 1/2 of the bait length <-b>
			  --Asserts <-T> is turned on	  
	-x,--max_reg	: Maximum length of target region to retain baits [500]
	-y,--min_reg	: Minimum length of target region to retain baits [80]
			  --Default will be set to bait length <-b,--bait>
	-V,--vmax_reg	: Maximumum number of SNPs in a target region to throw out all baits [0]
			  --Individual baits are constrained by <-s>
	-D,--dist_reg	: Minimum distance between disjunct target regions [100]
	-t,--tile_min	: Minimum bait region size to allow tiling 
			  --By default will be set to bait length <-b> + overlap <-O>
			  --Asserts <-T> is turned on
	-S, --select	: Which criterion to select which target region to use when w/in <-D> 
			  --Options (e.g. -S s=100 to choose region with most SNPs within 100 bases)
				s=/snp=[d]  : Most SNPs within \"d\" bases (replace [d] with value)
				b=/bad=[d]  : Least Ns and gaps within \"d\" bases
				c=/cons=[d] : Most conserved \"d\" bases
				m/min       : Least non-conserved (SNP,N,gap) bases in bait region
				r/rand      : Randomly choose a bait region to retain [default]""")

	#Need -s option for picking baits within bait regions and for picking bait regions within loci
	print("""
Bait Selection/ Optimization options:

""")
	print ("\nBait optimization/ selection options (for use with <-s>):")
	print ("	-C,--force_s	: Enforce <-s> criterion as a filter for all bait regions")
	print ("		 	  NOTE: When False, <-s> only chooses best of conflicting bait regions")
	print ("	-d,--opt_dist	: Number of flanking bases to consider when choosing among TRs")
	print ("	-u,--upper_snp	: Maximum SNPs when <-s=range>")
	print ("	-l,--lower_snp	: Minimum SNPs when <-s=range>")


	print ("\nOutput options:")
	print ("	-X,--expand	: In output bait table, expand all ambiguities")
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
		
