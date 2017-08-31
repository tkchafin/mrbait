#!/usr/bin/python

import getopt
import sys
import re

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
	-x,--max_r	: Maximum length of target region to retain [500]
	-y,--min_r	: Minimum length of target region to retain [80]
			  --Default will be set to bait length <-b,--bait>
	-V,--vmax_r	: Maximumum SNPs allowed in a target region [0]
			  --Individual baits are constrained by <-v>
			  --By default set to the value of <-v, --var_max>
			  --Also controllable with -F M=0,[x] where x is maximum value
	-D,--dist_r	: Minimum distance between disjunct target regions [100]
	-p,--tile_min	: Minimum bait region size to allow tiling
			  --By default will be set to bait length <-b> + overlap <-O>
			  --Asserts <-T> is turned on
	-S, --select_r	: Which criterion to select target regions w/in <-D>
			  --Options
				s=[d]   : Most SNPs w/in \"d\" bases
				b=[d]   : Least Ns and gaps w/n \"d\" bases
				c=[d]   : Most conserved w/in \"d\" bases
				m       : Least SNP, N, or gap bases in bait region
				r       : Randomly choose a bait region [default]
				Ex: -S s=100 to choose region with most SNPs w/in 100 bases
	-F,--filter_r	: Include any criteria used to filter ALL bait regions
			  --Warning: May mask selections made using <-S> or <-f>
			  --Options
				g=[x]   : Maximum of \"x\" gaps in target region
				b=[x]   : Maximum of \"x\" Ns in target region
				m=[d,x] : Minimum of \"x\" SNPs within \"d\" bases
				M=[d,x] : Maximum of \"x\" SNPs within \"d\" bases
				r=[x]   : Randomly retain \"x\" target regions w/ baits
				Ex: -F m=100,1 -F M=100,10 to sample when 1-10 SNPs w/in 100 bases""")
	#Need -s option for picking baits within bait regions and for picking bait regions within loci
	print("""
Bait Selection/ Optimization options:
	-a,--balign	: Maximum allowable alignment length between baits [40]
			  --By default is set to 1/2 bait length
	-s,--select_b	: Which criterion to select a bait for target region
			  --NOTE: This option is ignored when <-T> or <-W>
			  --Options
				s=[d]   : Most SNPs w/in \"d\" bases
				b=[d]   : Least Ns and gaps w/in \"d\" bases
				c=[d]   : Most conserved w/in \"d\" bases
				m       : Least SNP, N, or gap bases in bait region
				r       : Randomly choose a bait [default]
				Ex: -S s=100 to choose region with most SNPs within 100 bases
	-f,--filter_b	: Include any criteria used to filter ALL baits
			  --Warning: May mask selections made using <-S> or <-F>
			  --Options
				m=[d,x] : Minimum of \"x\" SNPs within \"d\" bases
				M=[d,x] : Maximum of \"x\" SNPs within \"d\" bases
				r=[x]   : Randomly retain \"x\" baits OVERALL
				Ex: -f m=100,1 -f M=100,10 to sample when 1-10 SNPs w/in 100 bases""")


	print("""
Running options/ shortcuts:
	-W,--tile_all	: Tile baits across all target regions
			  --Relevant arguments to set: -b, -O, -v, -n, -g, -f
	-Q,--quiet	: Shut up and run - don't output ANYTHING to stdout
			  --Errors and assertions are not affected""")
	print("""
Output options:
	-X,--expand	: In output bait table, expand all ambiguities
			  --Gaps are expanded as [ACGT] and absent
			  --\"N\"s are expanded as [ACGT]
	-o,--out	: Prefix for outputs [\"mrbait\"]""")

	#SPLIT option not figured out yet
	print("""
General options:
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
			options, remainder = getopt.getopt(sys.argv[1:], 'M:e:L:A:hc:l:t:b:w:Rm:v:n:Ng:GE:TO:x:y:V:D:p:S:F:a:s:f:WQXo:Pd:', \
			["maf=","gff=","loci=","assembly=",'help',"cov=","len=","thresh=",
			"bait=","win_shift=","mult_reg","min_mult=","var_max=","numN=",
			"callN","numG=","callG","gff_type=","tiling","overlap=","max_r",
			"min_r=","vmax_r=","dist_r=","tile_min=","select_r=","filter_r=",
			"balign=","select_b=","filter_b=","tile_all","queit","expand","out=",
			"plot_all", "win_width="])
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
		self.tiling=0 #bool
		self.overlap=40
		self.max_r=500
		self.min_r=None
		self.vmax_r=None
		self.dist_r=None
		self.tile_min=None
		self.select_r="r"
		self.select_r_dist=None
		self.filter_r=0 #bool
		self.filter_t_whole=None
		self.filter_r_objects=[]

		#Bait selection options
		self.balign=None
		self.select_b="r"
		self.select_b_dist=None
		self.filter_b=0 #bool
		self.filter_b_whole=None
		self.filter_b_objects=[]

		#Running options/ shortcuts
		self.tile_all = 0
		self.stfu = 0

		#Output options
		self.expand = 0
		self.out = "mrbait"


		self.ploidy=2
		self.db="./mrbait.sqlite"
				#Parse through arguments and set params
		for opt, arg in options:
			arg = arg.replace(" ","")
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

			#Bait general params
			elif opt in ('-b', '--bait'):
				self.blen = int(arg)
			elif opt in ('-w', '--win_shift'):
				self.win_shift = int(arg)
			elif opt in ('-R', '--mult_reg'):
				self.mult_reg = 1
			elif opt in ('-d', '--win_width'):
				print("Warning: You are setting a hidden option <-d/--win_width>")
				self.win_width = int(arg)
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
			elif opt in ('-T', '--tiling'):
				self.tiling = 1
			elif opt in ('-O', '--overlap'):
				self.overlap = int(arg)
			elif opt in ('-x', '--max_r'):
				self.max_r = int(arg)
			elif opt in ('-y', '--min_r'):
				self.min_r = int(arg)
			elif opt in ('-V', '--vmax_r'):
				self.vmax_r = int(arg)
			elif opt in ('-D', '--dist_r'):
				self.dist_r = int(arg)
			elif opt in ('-p', '--tile_min'):
				self.tile_min = int(arg)
				self.tiling = 1
			elif opt in ('-S', '--select_r'):
				temp = arg.split('=')
				self.select_r = (temp[0]).lower()
				chars = (['s','b','c','m','r'])
				if self.select_r not in chars:
					raise ValueError("Invalid option \"%r\" for <--select_r>" % self.select_r)
				subchars = (['s','b','c'])
				if string_containsAny(self.select_r, subchars) == 0:
					if (len(temp) > 1):
						assert len(temp) is 2, "invalid use of <--select_r>"
						self.select_r_dist = int(temp[1])
						assert isinstance(self.select_r_dist, int), "select_r_dist must be an integer"
						assert self.select_r_dist >= 0, "select_r_dist must be an integer greater than zero!"
					else:
						self.select_r_dist = 100
				else:
					self.select_r_dist = None
				#print("select_r is %r" %self.select_r)
				#print("select_r_dist is %r"%self.select_r_dist)
			elif opt in ('-F', '--filter_r'):
				self.filter_r = 1 #turn on region filtering
				#temp = arg.split('/') #parse region filtering options
				self.filter_r_whole = arg
				#for sub in temp:
				subopts = re.split('=|,',arg)
				if subopts[0] in ('m','M'):
					assert len(subopts) == 3, "Incorrect specification of option %r for <--filter_r>" %subopts[0]
					self.filter_r_objects.append(subArg(subopts[0],subopts[1],subopts[2]))
				elif subopts[0] in ('r','g','n'):
					assert len(subopts) == 2, "Incorrect specification of option %r for <--filter_r>" %subopts[0]
					self.filter_r_objects.append(subArg(subopts[0],subopts[1]))
				else:
					bad_opts("Invalid option %r for <--filter_r>!" %subopts[0])

			#Bait selection options
			elif opt in ('-a', '--balign'):
				self.balign = int(arg)
			elif opt in ('-s', '--select_b'):
				temp = arg.split('=')
				self.select_b = (temp[0]).lower()
				chars = (['s','b','c','m','r'])
				if self.select_b not in chars:
					raise ValueError("Invalid option \"%r\" for <--select_r>" % self.select_r)
				subchars = (['s','b','c'])
				if string_containsAny(self.select_b, subchars) == 0:
					if (len(temp) > 1):
						self.select_b_dist = int(temp[1])
						assert self.select_b_dist >= 0, "select_b_dist must be an integer greater than zero!"
					else:
						self.select_b_dist = 100
				else:
					self.select_b_dist = None
				print("select_b is %r" %self.select_b)
				print("select_b_dist is %r"%self.select_b_dist)
			elif opt in ('-f', '--filter_b'):
				self.filter_b = 1 #turn on region filtering
				#temp = arg.split('/') #parse region filtering options
				self.filter_b_whole = arg
				#for sub in temp:
				subopts = re.split('=|,',arg)
				if subopts[0] in ('m','M'):
					assert len(subopts) == 3, "Incorrect specification of option %r for <--filter_b>" %subopts[0]
					self.filter_b_objects.append(subArg(subopts[0],subopts[1],subopts[2]))
				elif (subopts[0] is 'r'):
					assert len(subopts) == 2, "Incorrect specification of option %r for <--filter_b>" %subopts[0]
					self.filter_b_objects.append(subArg(subopts[0],subopts[1]))
				else:
					bad_opts("Invalid option %r for <--filter_b>!" %subopts[0])

			#Running options
			elif opt in ('-W', '--tile_all'):
				self.tile_all = 1
			elif opt in ('-Q', '--quiet'):
				self.stfu = 1

			#output options
			elif opt in ('-X', '--expand'):
				self.expand = 1
			elif opt in ('-o', '--out'):
				self.out = arg

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

		#Default vmax_r
		if self.vmax_r is None:
			self.vmax_r = self.var_max

		#Warn of argument masking
		if self.mult_reg is 0:
			if self.dist_r is not None:
				print("Warning: You have set <--dist_r/-D>, but it is ignored when <--mult_reg/-R> is off. ")
			if self.min_mult is not None:
				print("Warning: You have set <--min_mult/-m>, but it is ignored when <--mult_reg/-R> is off. ")

		#Default of dist_r
		if self.dist_r is None:
			self.dist_r = 100

		#set default of min_mult
		if self.min_mult is None:
			self.min_mult = 1000

		#Assert that win_shift cannot be larger than blen
		if self.minlen is None:
			self.minlen = self.blen
		elif self.blen > self.minlen:
			self.minlen = self.blen

		#Set default value for balign
		if self.balign is None:
			self.balign = self.blen / 2

		#Set minimum target region size
		if self.min_r is None:
			self.min_r = self.blen

		#Set tile_min default if not given
		if (self.tiling == 1) and (self.tile_min is None):
			self.tile_min = self.blen + self.overlap

		#Call help menu if prompted
		if call_help is 1:
			display_help("Exiting because help menu was called.")
			sys.exit(0)
