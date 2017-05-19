#!/usr/bin/python

import getopt
import sys

def display_help(message=None):
	if message is not None:
		print (message)
	print ("\nMrBait version alpha 0.01\n")
	print ("Contact:\n\n\tTyler K. Chafin\n\tUniversity of Arkansas\n\ttkchafin@uark.edu\n")
	print ("\nUsage:\n\n\t", sys.argv[0], "-a </path/to/MAF> <-c><-l><-t>\n")
	print ("Description:\n")
	print("\tCURRENTLY IN DEVELOPMENT\n")
	print("\tMrBait is a program for designing baits for targeted enrichment studies\n",\
	"\tincluding Ultraconverved Elements (UCEs), Anchored enrichment, and RAD-capture.\n",\
	"\tIt can also be used to customize bait design for any of these methods based on\n",\
	"\ta variety of built-in methods and constraints.\n")
	
	print("\tNOTE: SQLite databases are currently remade from scratch with each\n",\
	"\tinstance of the program. Support for opening and editing older database files\n",\
	"\twill come soon.\n")

	print("""
Input options:

	-M,--maf	: Input multiple alignment MAF file [more to come later]
	-e,--gff	: GFF file containing annotation information [NOT WORKING YET]
	-L,--loci	: For RAD-data, as the \".loci\" output of pyRAD
	-A,--assembly	: Input whole genome assembly as FASTA [NOT WORKING YET]""")
	
	print("""
Locus filtering/ consensus options:

	-c,--cov	: Minimum number of sequences per alignment, for MAF or LOCI input [1] 
	-l,--len	: Minimum alignment length to attempt bait design [80]
	-t,--thresh	: Minimum proportion for gap/N to include in consensus [0.1]""")

	print("""
General Bait Design options:

	-b, --bait	: Bait length [default=80]
	-w, --win_shift	: Sliding window shift distance [1]
			  --Increasing will speed up program but could lower accuracy
			  --Cannot be longer than bait length <-b,--bait>
	-R,--mult_reg	: Allow multiple target regions per locus [false]
			  --See \"Target Region options\" below
	-m,--min_mult	: Minimum alignment length to allow multiple regions [1000]
			  --By default is set to the value of <-l>. Only applies when <-R>
			  --When -R and -T are off, only one bait is selected per locus
	-v,--var_max	: Maximum number of SNP columns to be included in a bait [0]
			  --These will be coded as the appropriate IUPAC ambiguity code
			  --Can be expanded in final output using the <--expand> flag
	-n,--numN	: Number of consensus Ns allowed in a bait [0]
			  --These will be inserted as \"N\" unless <-N> is used
	-N,--callN	: For Ns in bait sequence, call majority nucleotide if possible
	-g,--numG	: Number of consensus gaps (\"-\") allowed in bait [0]
			  --These will be inserted as \"N\" unless <-G> is used
	-G,--callG	: For gaps in bait, call as majority nucleotide if possible
	-E,--gff_type	: Constrain bait design to within element type of GFF
			  --COMING SOON, NOT WORKING YET""")

	print("""
Target Region options:

	-T,--tiling	: Turn on tiling to design multiple baits per bait region 
			  --When tiling is off, see Bait Selection options below 
	-O,--overlap	: Overlap for tiled baits within target regions [40]
			  --By default will be set to 1/2 of the bait length <-b>
			  --Asserts <-T> is turned on	  
	-x,--max_r	: Maximum length of target region to retain baits [500]
	-y,--min_r	: Minimum length of target region to retain baits [80]
			  --Default will be set to bait length <-b,--bait>
	-V,--vmax_r	: Maximumum SNPs in a target region to throw out baits [0]
			  --Individual baits are constrained by <-s>
	-D,--dist_r	: Minimum distance between disjunct target regions [100]
	-p,--tile_min	: Minimum bait region size to allow tiling 
			  --By default will be set to bait length <-b> + overlap <-O>
			  --Asserts <-T> is turned on
	-S, --select_r	: Which criterion to select target regions w/in <-D> 
			  --Options 
				s=/snp=[d]   : Most SNPs within \"d\" bases
				b=/bad=[d]   : Least Ns and gaps within \"d\" bases
				c=/cons=[d]  : Most conserved \"d\" bases
				m/min        : Least SNP, N, or gap bases in bait region
				r/rand       : Randomly choose a bait region [default]
				Ex: -S s=100 to choose region with most SNPs w/in 100 bases
	-F,--filter_r	: Include any criteria used to filter ALL bait regions
			  --Warning: May mask selections made using <-S> or <-f>	
			  --Options
				m=/min=[d,x] : Minimum of \"x\" SNPs within \"d\" bases
				M=/max=[d,x] : Maximum of \"x\" SNPs within \"d\" bases
				r=/rand=[x]  : Randomly retain \"x\" target regions w/ baits
				Ex: m=100,1/M=100,10 to sample when 1-10 SNPs w/in 100 bases""")
	#Need -s option for picking baits within bait regions and for picking bait regions within loci
	print("""
Bait Selection/ Optimization options:
	-a,--balign	: Maximum allowable alignment length between baits [40]
			  --By default is set to 1/2 bait length
	-s,--select_b	: Which criterion to select a bait for target region 
			  --NOTE: This option is ignored when <-T> or <-X>
			  --Options 
				s=/snp=[d]   : Most SNPs within \"d\" bases 
				b=/bad=[d]   : Least Ns and gaps within \"d\" bases
				c=/cons=[d]  : Most conserved \"d\" bases
				m/min        : Least SNP, N, or gap bases in bait region
				r/rand       : Randomly choose a bait [default]
				Ex: -S s=100 to choose region with most SNPs within 100 bases
	-f,--filter_b	: Include any criteria used to filter ALL baits
			  --Warning: May mask selections made using <-S> or <-F>	
			  --Options 
				m=/min=[d,x] : Minimum of \"x\" SNPs within \"d\" bases
				M=/max=[d,x] : Maximum of \"x\" SNPs within \"d\" bases
				r=/rand=[x]  : Randomly retain \"x\" baits OVERALL
				Ex: m=100,1/M=100,10 to sample when 1-10 SNPs w/in 100 bases""")


	print("""
Running options/ shortcuts:
	-W,--tile_all	: Ignore target regions and tile baits across all loci
			  --Relevant arguments to set: -b, -O, -v, -n, -g
	-Q,--quiet	: Shut up and run - don't output ANYTHING to stdout 
			  --Errors and assertions are not affected""")
	print("""
Output options:
	-X,--expand	: In output bait table, expand all ambiguities
	-o,--out	: Prefix for outputs [\"mrbait\"]""")
	
	#SPLIT option not figured out yet
	print("""
Output options:
	-Z,--split	: Split large file and read into memory piece-meal
	-h,--help	: Displays this help menu 
	-P,--plot_all	: Plot variable positions for all baits/bait regions""")
	print()

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try: 
			options, remainder = getopt.getopt(sys.argv[1:], 'M:e:L:A:hc:l:t:b:w:Rm:v:n:Ng:GE:TO:x:y:V:D:p:S:F:a:s:f:WQXo:P', \
			["maf=","gff=","loci=","assembly=",'help',"cov=","len=","thresh=",
			"bait=","win_shift=","mult_reg","min_mult=","var_max=","numN=",
			"callN","numG=","callG","gff_type=","tiling","overlap=","max_r",
			"min_r=","vmax_r=","dist_r=","tile_min=","select_r=","filter_r=",
			"balign=","select_b=","filter_b=","tile_all","queit","expand","out=",
			"plot_all"])
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
		
