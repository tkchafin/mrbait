#!/usr/bin/python

import getopt
import sys
import re
import ntpath
import misc_utils as utils

def string_containsAny(str, set):
	for c in set:
		if c in str: return 0;
	return 1;

def bad_opts(message=""):
	print(message)
	print("Invalid options: Exiting program. Please see manual or use <--help>\n")
	sys.exit(1)

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

	-M,--maf	: Input multiple alignment MAF file
	-e,--gff	: GFF file containing annotation information [NOT WORKING YET]
	-L,--loci	: For RAD-data, as the \".loci\" output of pyRAD
	-A,--assembly	: Input whole genome assembly as FASTA [NOT WORKING YET]""")

	print("""
Locus filtering/ consensus options:

	-c,--cov	: Minimum number of sequences per alignment, for MAF or LOCI input [1]
	-l,--len	: Minimum alignment length to attempt bait design [80]
	-t,--thresh	: Minimum proportion for gap/N to include in consensus [0.1]
	-k,--mask	: Minimum proportion for masking (lower case) to include in consensus [0.1]""")

	print("""
General Bait Design options:

	-b, --bait	: Bait length [default=80]
	-w, --win_shift	: Sliding window shift distance [1]
		  	--Increasing will speed up program but could lower accuracy
		  	--Cannot be longer than bait length <-b,--bait>
	-R,--mult_reg	: Allow multiple target regions per locus [false]
			--See \"Target Region options\" below
	-m,--min_mult	: Minimum alignment length to allow multiple regions
			--By default is set to the value of <-l>. Only applies when <-R>
			--When -R and -T are off, only one bait is selected per locus
	-v,--var_max	: Maximum number of SNP columns to be included in a bait [0]
			--These will be coded as the appropriate IUPAC ambiguity code
			--Can be expanded in final output using the <--expand> flag
	-n,--numN	: Number of consensus Ns allowed in a bait [0]
			--These will be inserted as \"N\" unless <-N> is used
	-g,--numG	: Number of consensus indels (\"-\") allowed in bait [0]
	-E,--gff_type	: Constrain bait design to within element type of GFF
			--COMING SOON, NOT WORKING YET""")

	print("""
Target Region options:

	-D,--dist_r	: Minimum distance between target regions [100]
			--Conflicts will be resolved according to --select_r
	-d,--flank_dist	: Distance in flanking regions selection and filtering targets
			--See relevant options in <--select_r> and <--filter_r>
			--By default will be set to [500]
	-S, --select_r	: Which criterion to select target regions w/in <-D>
			--Options
			snp          : Most SNPs w/in \"d\" bases
			bad          : Least Ns and gaps w/n \"d\" bases
			cons         : Most conserved w/in \"d\" bases
			rand         : Randomly choose a bait region [default]
			Ex: -S s=100 to choose region with most SNPs w/in 100 bases
	-F,--filter_r	: Include any criteria used to filter ALL bait regions
			--Options
			len=[x,y]    : Target length between \"x\" (min) and \"y\" (max)
			gap=[x]      : Maximum of \"x\" indels in target region
			bad=[x]      : Maximum of \"x\" Ns in target region
			min=[x]      : Minimum of \"x\" SNPs within \"d\" bases
			max=[x]      : Maximum of \"x\" SNPs within \"d\" bases
			mask=[x,y]   : Proportion masked bases between \"x\" (min) and \"y\" (max)
			gc=[x,y]     : Proportion of G/C bases between \"x\" (min) and \"y\" (max)
			rand=[x]     : Randomly retain \"x\" target regions w/ baits
			aln=[i,q]    : Pairwise alignment, removing when \"i\" identity in \"q\" proportion
			blast_i=[i,q]: Only retain hits over \"i\" identity and \"q\" query coverage to provided db
			blast_x=[i,q]: Remove hits over \"i\" identity and \"q\" query coverage to provided db
			blast_a=[i,q]: Remove targets having ambiguous mapping (more than one hit) to a provided db
			gff=         : ADD LATER!!! Only keep targets within X distance of GFF element
			Ex: -F min=100,1 -F max=100,10 to sample when 1-10 SNPs w/in 100 bases""")

	print("""
Bait Design / Selection options:

	-s,--select_b	: Which criterion to select a bait for target region
			--NOTE: This option is ignored when <-W>
			--Options
		  	tile=[x]     : Tile baits with \"x\" overlapping bases [default; x= bait_length/2]
		  	center=[n,x] : Design \"n\" centered baits with \"x\" overlap
			flank=[n,x]  : Design \"x\" terminal baits (each side) with \"x\" overlap
			Ex: -s tile=40 to tile baits with an overlap of 40 bases
	-f,--filter_b	: Include any criteria used to filter ALL baits
			--Options
			mask=[x,y]   : Proportion masked bases between \"x\" (min) and \"y\" (max)
			gc=[x,y]     : Proportion of G/C bases between \"x\" (min) and \"y\" (max)
			aln=[i,q]    : Pairwise alignment, removing when \"i\" identity in \"q\" proportion
			blast_i=[i,q]: ADD LATER!!! Keep targets hitting BLAST db
			blast_x=[i,q]: ADD LATER!!! Remove targets hitting BLAST db (e.g. contamination)
			rand=[x]     : Randomly retain \"x\" baits""")

	print("""
Pairwise-Alignment Deduplication (use when --select_b or --select_r = \"aln\"):

	--vsearch	: Path to VSEARCH executable if other than provided
			--MrBait will try to detect OS and appropriate exectable to use
	--vthreads	: Number of threads for VSEARCH to run in parallel [default=4]""")

	print("""
BLAST Paramters (use when --select_b or --select_r = \"blast\"):

	--blastdb	: Path and prefix to existing Blast database to use
	--fastadb	: Path and prefix to FASTA formatted sequences to make blast db
	--e_value	: Minimum e-value cutoff for reporting BLAST hits
	--gapopen	: Gap opening penalty for blastn
	--gapextend	: Gap extension penatly for blastn
	--word_size	: Word size for blastn
	--megablast	: Use megablast rather than blastn
	--bthreads	: Number of threads for BLAST to run in parallel [default=4]""")


	print("""
Output options:

	-X,--expand	: In output bait table, expand all ambiguities
			--Gaps are expanded as [ACGT] and absent
			--\"N\"s are expanded as [ACGT]
	-o,--out	: Output directory and prefix [e.g. ./baits/run1]""")

	#SPLIT option not figured out yet
	print("""
General options:

	-W,--tile_all	: Tile baits across all target regions
			--Skips target filtering and selection, and just tiles all
	-K, --no_mask	: Ignore all masking information [boolean]
	-Q,--quiet	: Shut up and run - don't output ANYTHING to stdout
			--Errors and assertions are not affected
	-h,--help	: Displays this help menu
	""")
	print()

#Sub-object for holding filtering options
class subArg():
	def __init__(self, o1, o2, o3=None):
		self.o1 = o1
		self.o2 = o2
		self.o3 = o3

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'M:e:L:A:hc:l:t:b:w:Rm:v:n:Ng:GE:D:p:S:F:s:f:QXo:Pk:Kd:', \
			["maf=","gff=","loci=","assembly=",'help',"cov=","len=","thresh=",
			"bait=","win_shift=","mult_reg","min_mult=","var_max=","numN=",
			"callN","numG=","callG","gff_type=","dist_r=","tile_min=",
			"select_r=","filter_r=",
			"select_b=","filter_b=","quiet","expand","out=",
			"plot_all","mask=","no_mask", "flank_dist=","vsearch=",
			"vthreads=","hacker="])
		except getopt.GetoptError as err:
			print(err)
			display_help("\nExiting because getopt returned non-zero exit status.")
			sys.exit(2)
		#Default values for params
		call_help=0 #boolean

		#Input params
		self.alignment=None
		self.gff=None
		self.loci=None
		self.assembly=None

		#Locus filtering params
		self.cov=1
		self.minlen=None
		self.thresh=0.1
		self.mask=0.1

		#Bait params
		self.blen=80
		self.win_width=None
		self.win_shift=1
		self.mult_reg=0 #boolean
		self.min_mult=None
		self.var_max=0
		self.numN=0
		self.callN=0 #boolean
		self.numG=0
		self.callG=0 #boolean
		self.anchor=None

		#target region options
		self.dist_r=None
		self.flank_dist=500
		self.select_r="rand"
		self.filter_r=0 #bool
		self.filter_t_whole=None
		self.filter_r_objects=[]

		#VSEARCH options - deduplication
		self.vsearch = None
		self.vthreads = 4

		#Bait selection options
		self.overlap=None
		self.bait_shift=None
		self.select_b="tile"
		self.select_b_num=None
		self.filter_b=0 #bool
		self.filter_b_whole=None
		self.filter_b_objects=[]

		#Running options/ shortcuts
		self.no_mask = 0
		self.stfu = 0

		#Output options
		self.expand = 0
		self.out = None
		self.workdir = None

		self.ploidy=2
		self.db="./mrbait.sqlite"

		#HACKER ONLY OPTIONS
		self._noGraph = 0
		self._noWeightGraph = 0
		self._weightMax = 50000 #maximum size to attempt weighted edge resolution
		self._weightByMin = 0

				#Parse through arguments and set params
		for opt, arg_raw in options:
			arg = arg_raw.replace(" ","")
			arg = arg.strip()
			if opt in ('-M', '--maf'):
				self.alignment = arg
			elif opt in ('-h', '--help'):
				call_help = 1

			#Input params
			elif opt in ('-e', '--gff'):
				self.gff = arg
			elif opt in ('-L', '--loci'):
				self.loci = arg
			elif opt in ('-A', '--assembly'):
				self.assembly = arg

			#Locus filtering params
			elif opt in ('-c', '--cov'):
				self.cov = int(arg)
			elif opt in ('-l', '--len'):
				self.minlen = int(arg)
			elif opt in ('-t', '--thresh'):
				self.thresh = float(arg)
			elif opt in ('-k', '--mask'):
				self.mask = float(arg)

			#Bait general params
			elif opt in ('-b', '--bait'):
				self.blen = int(arg)
			elif opt in ('-w', '--win_shift'):
				self.win_shift = int(arg)
			elif opt in ('-R', '--mult_reg'):
				self.mult_reg = 1
			elif opt in ('-m', '--min_mult'):
				self.min_mult = int(arg)
			elif opt in ('-v', '--var_max'):
				self.var_max = int(arg)
			elif opt in ('-n', '--numN'):
				self.numN = int(arg)
			elif opt in ('-N', '--callN'):
				self.callN = 1
			elif opt in ('-g', '--numG'):
				self.numG = int(arg)
			elif opt in ('-G', '--callG'):
				self.callG = 1
			elif opt in ('-E', '--gff_type'):
				self.anchor = arg

			#target region opts
			elif opt in ('-D', '--dist_r'):
				self.dist_r = int(arg)
			elif opt in ('-p', '--tile_min'):
				self.tile_min = int(arg)
				self.tiling = 1
			elif opt in ('-d', '--flank_dist'):
				self.flank_dist = int(arg)
				assert isinstance(self.flank_dist, int), "<--flank_dist> must be an integer"
				assert self.flank_dist >= 0, "<--flank_dist> must be an integer greater than zero!"
			elif opt in ('-S', '--select_r'):
				temp = arg.split('=')
				assert len(temp) == 1, "Invalid specification for <--select_r>: %s"%arg
				self.select_r = (temp[0]).lower()
				chars = (['snp','bad','cons','rand'])
				if self.select_r not in chars:
					raise ValueError("Invalid option \"%r\" for <--select_r>" % self.select_r)
			elif opt in ('-F', '--filter_r'):
				self.filter_r = 1 #turn on region filtering
				#temp = arg.split('/') #parse region filtering options
				self.filter_r_whole = arg
				#for sub in temp:
				subopts = re.split('=|,',arg)
				if subopts[0] in ('snp','mask','gc','len', "aln"): #TODO: Add blast options
					assert len(subopts) == 3, "Incorrect specification of option %r for <--filter_r>" %subopts[0]
					if subopts[0] in ('gc', 'mask','len'):
						assert subopts[1] < subopts[2], "In <--filter_r> suboption \"%s\": Min must be less than max"%subopts[0]
						self.filter_r_objects.append(subArg(subopts[0],float(subopts[1]),float(subopts[2])))
					elif subopts[0] == "aln":
						self.filter_r_objects.append(subArg(subopts[0],float(subopts[1]),float(subopts[2])))
					elif subopts[0] == "snp":
						minS = int(subopts[1])
						maxS = int(subopts[2])
						assert 0 <= minS < maxS, "In <--filter_r> suboption \"%s\": Min must be less than max"%subopts[0]
						assert isinstance(minS, int), "--filter_r suboption snp: Minimum must be an integer"
						assert isinstance(maxS, int), "--filter_r suboption snp: Maximum must be an integer"
						self.filter_r_objects.append(subArg(subopts[0],minS,maxS))
				elif subopts[0] in ('rand','gap','bad'):
					assert len(subopts) == 2, "Incorrect specification of option %r for <--filter_r>" %subopts[0]
					self.filter_r_objects.append(subArg(subopts[0],int(subopts[1])))
				else:
					bad_opts("Invalid option %r for <--filter_r>!" %subopts[0])

			#Bait selection options
			elif opt in ('-s', '--select_b'):
				subopts = re.split('=|,',arg)
				self.select_b = (subopts[0]).lower()
				chars = (['tile', 'center', 'flank'])
				if self.select_b not in chars:
					raise ValueError("Invalid option \"%r\" for <--select_b>" % self.select_b)
				subchars = (['center','flank'])
				if self.select_b in subchars:
					assert len(subopts) == 3, "Incorrect specification of option %r for <--select_b>" %subopts[0]
					self.overlap = int(subopts[2])
					self.select_b_num = int(subopts[1])
					assert self.overlap < self.blen, "Overlap distance cannot be greater than bait length"
				elif self.select_b == "tile":
					self.select_b_num = None
					self.overlap = int(subopts[1])
				#print("select_b is %r" %self.select_b)
				#print("select_b_dist is %r"%self.select_b_dist)
			elif opt in ('-f', '--filter_b'):
				self.filter_b = 1 #turn on region filtering
				#temp = arg.split('/') #parse region filtering options
				self.filter_b_whole = arg
				#for sub in temp:
				subopts = re.split('=|,',arg)
				if subopts[0] in ('min','max','mask','gc'):
					assert len(subopts) == 3, "Incorrect specification of option %r for <--filter_b>" %subopts[0]
					if subopts[0] in ('gc', 'mask'):
						assert subopts[1] < subopts[2], "In <--filter_b> for suboptions \"mask\" and \"gc\": Min must be less than max"
						self.filter_r_objects.append(subArg(subopts[0],float(subopts[1]),float(subopts[2])))
					else:
						self.filter_r_objects.append(subArg(subopts[0],int(subopts[1]),int(subopts[2])))
				elif (subopts[0] is 'rand'):
					assert len(subopts) == 2, "Incorrect specification of option %r for <--filter_b>" %subopts[0]
					self.filter_b_objects.append(subArg(subopts[0],subopts[1]))
				else:
					bad_opts("Invalid option %r for <--filter_b>!" %subopts[0])

			#Running options
			elif opt in ('-Q', '--quiet'):
				self.stfu = 1
			elif opt in ('-K', '--no_mask'):
				self.no_mask = 1

			#vsearch options
			elif opt in ("--vsearch"):
				self.vsearch = str(arg)
			elif opt in ("--vthreads"):
				self.vthreads = int(arg)

			#output options
			elif opt in ('-X', '--expand'):
				self.expand = 1
			elif opt in ('-o', '--out'):
				self.out = arg

			#HACKER ONLY OPTIONS
			elif opt in ('--hacker'):
				print(opt, arg)
				subopts = re.split('=|, ',arg)
				main = subopts[0]
				print("Warning: Using \"hacker only\" option <%s>. Be careful."%main)
				if main == "noGraph":
					self._noGraph = 1
				elif main == "win_width":
					assert len(subopts) == 2, "Warning: HACKER option <win_width> must have two arguments separated by \"=\""
					win_width = int(subopts[1])
				elif main == "noWeightGraph":
					self._noWeightGraph = 1
				elif main == "weightByMin":
					self._weightByMin = 1
				elif main == "weightMax":
					assert len(subopts) == 2, "Warning: HACKER option <graphMax> must have two arguments separated by \"=\""
					self._weightMax = int(subopts[1])
				else:
					assert False, "Unhandled option %r"%main
			else:
				assert False, "Unhandled option %r"%opt

		#DEBUG PRINTS
		for subopt in self.filter_r_objects:
			print("filter_r: Suboption %s has parameters: %s %s" %(subopt.o1,subopt.o2,subopt.o3))

		#for subopt in self.filter_b_objects:
			#print("filter_b: Suboption %s has parameters: %s %s" %(subopt.o1,subopt.o2,subopt.o3))

		#Assertions and conditional changes to params
		if (self.alignment is None) and (self.loci is None) and (self.assembly is None):
			display_help("Input not specified!")
			sys.exit(0)

		assert self.blen > 0, "Bait length cannot be less than or equal to zero!"

		#Set default win_width
		if self.win_width is None:
			self.win_width = self.blen

		#Warn of argument masking
		if self.mult_reg is 0:
			if self.dist_r is not None:
				print("Warning: You have set <--dist_r/-D>, but it is ignored when <--mult_reg/-R> is off. ")

		#Default of dist_r
		if self.dist_r is None:
			self.dist_r = 100

		#Default of overlap
		if self.overlap is None:
			self.overlap = (self.blen // 2)
		self.bait_shift = self.blen - self.overlap

		#if --no_mask, set mask thresh to 1.0
		if self.no_mask:
			self.mask = 1.0

		#Get working dir path and output prefix
		if self.out is None:
			self.out = "mrbait"
			self.workdir = utils.getWorkingDir()
		else:
			self.workdir, self.out = ntpath.split(self.out)
			if self.out == "":
				self.out = "mrbait"
			if self.workdir == "":
				self.workdir = os.getcwd()
		print("Path is: ", self.workdir)
		print("Prefix is: ", self.out)


		#Assert that win_shift cannot be larger than blen
		if self.minlen is None:
			self.minlen = self.blen
		elif self.blen > self.minlen:
			self.minlen = self.blen

		#set default of min_mult
		if self.min_mult is None:
			self.min_mult = self.minlen

		#If vsearch path not given, try to figure it out
		if self.vsearch is None:
			os_platform = utils.getOS()
			print("Found OS platform:", os_platform)
			if os_platform == "linux" or os_platform == "unknown":
				print("Automatically detected LINUX or UNKNOWN platform: Using LINUX VSEARCH executable.")
				self.vsearch = utils.getScriptPath() + "/bin/vsearch-2.4.4-linux"
			elif os_platform == "darwin": #mac os
				print("Automatically detected MACOS platform: Using MACOS VSEARCH executable.")
				self.vsearch = utils.getScriptPath() + "/bin/vsearch-2.4.4-macos"

		#Default bait design behavior
		if self.select_b == "tile" and self.overlap is None:
			self.overlap = self.blen // 2

		#Call help menu if prompted
		if call_help is 1:
			display_help("Exiting because help menu was called.")
			sys.exit(0)
