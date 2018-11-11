#!/usr/bin/python

import getopt
import sys
import os
import re
import ntpath
from mrbait import misc_utils as utils

def string_containsAny(str, set):
	for c in set:
		if c in str: return 0;
	return 1;

def bad_opts(message=""):
	print(message)
	print("Invalid options: Exiting program. Please see manual or use <--help>\n")
	sys.exit(1)

#Function to print header information
def printHeader():
	print("""

	=======================================================================
	    MrBait: Universal Probe Design for Targeted-Enrichment Methods
	=======================================================================

	Version: 1.1.5
	Author: Tyler K. Chafin
	Contact: tkchafin@uark.edu
	License: GNU Public License v3.0

	Releases: https://github.com/tkchafin/mrbait/releases

	Citation: Chafin TK, Douglas MR, Douglas ME. 2018. MrBait: Universal
		identification and design of targeted-enrichment capture probes.
		Bioinformatics. https://doi.org/10.1093/bioinformatics/bty548

	Source code: https://github.com/tkchafin/mrbait
	Documentation: https://mrbait.readthedocs.io/
	Conda package: https://anaconda.org/tylerkchafin/mrbait

	=======================================================================
	""")

def display_help(message=None):
	if message is not None:
		print (message)
	print ("\nMrBait Help Menu\n\nFor more detailed information, please see the manual.\n")
	print ("\nUsage:\n\n\t", sys.argv[0], "-A </path/to/MAF> <-c><-l><-t><...>\n")
	print ("Description:\n")
	print("\tMrBait is a program for designing baits for targeted enrichment studies\n",\
	"\tincluding Ultraconverved Elements (UCEs), Anchored enrichment, and RAD-capture.\n",\
	"\tIt can also be used to customize bait design for any of these methods based on\n",\
	"\ta variety of built-in methods and constraints.\n")

	print("\tNOTE: SQLite databases are currently remade from scratch with each\n",\
	"\tinstance of the program. Support for opening and editing older database files\n",\
	"\twill come soon.\n")
	#SPLIT option not figured out yet
	print("""
General options:

	-r,--resume	: Resumes from a database pre-loaded with alignments (integer)
		--Options
			-r 1 : Continues pipeline beginning after Step 1 (Loading alignments)
			-r 2 : Continues after Step 2 (Target discovery)
			-r 3 : Continues after Step 3 (Target filtering/ selection)
			-r 4 : Continues after step 4 (Bait discovery)
	--db		: .sqlite file containing pre-existing database. For use with --resume
	-T,--threads	: Number of threads to use for processes that can run in parallel [1]
	-h,--help	: Displays this help menu
	""")
	print("""
Input options:

	-M,--maf	: Input multiple alignment MAF file
	-L,--loci	: For RAD-data, as the \".loci\" output of pyRAD
	-A,--assembly	: Input whole genome assembly as FASTA
	-X,--xmfa	: Input whole genome alignments as XMFA""")
	print("""
Assembly input options (for use only with -A <genome.fasta>):

	-V,--vcf	: VCF file containing variant information
	-G,--gff	: GFF file containing annotation information
	--vcfALT	: If using VCF and reference allele is gap/N, attempt
		 	   to call new consensus using VCF ALT alleles [default=OFF]""")
	print("""
Alignment filtering/ consensus options (use with -M or -L inputs):

	-c,--cov	: Minimum number of sequences per alignment, for MAF or LOCI input [1]
	-l,--len	: Minimum alignment length to attempt bait design [80]
	-q,--thresh	: Threshold proportion for gap/N to include in consensus [0.1]
	-Q,--max_ambig	: Maximum proportion of gap/N bases allowed in a consensus sequence [0.5]
	-k,--mask	: Threshold proportion for masking (lower case) to include in consensus [0.1]
	-K,--max_mask	: Maximum proportion of masked bases allowed in a consensus sequence [0.5]""")

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
	-v,--var_max	: Maximum number of SNP columns to be allowed in a bait [0]
		--These will be coded as the appropriate IUPAC ambiguity code
		--Can be expanded in final output using the <--expand> flag
	-n,--numN	: Number of consensus Ns allowed in a bait [0]
		--These will be inserted as \"N\" unless <-N> is used
	-g,--numG	: Number of consensus indels (\"-\") allowed in bait [0]""")

	print("""
Target Region options:

	-D,--dist_r	: Minimum distance between target regions [100]
		--Conflicts will be resolved according to --select_r
	-d,--flank_dist	: Distance in flanking regions selection and filtering targets
		--See relevant options in <--select_r> and <--filter_r>
		--By default will be set to [500]
		--Appopriate value likely will correspond to your sequencing insert size
	-S, --select_r	: Which criterion to select target regions w/in <-D>
		--Options
			snp          : Most SNPs w/in \"d\" bases
			bad          : Least Ns and gaps w/n \"d\" bases
			cons         : Most conserved w/in \"d\" bases
			rand         : Randomly choose a bait region [default]
			Ex: -S snp -d 100 to choose region with most SNPs w/in 100 bases
	-F,--filter_r	: Include any criteria used to filter ALL bait regions
		--Options
			len=[x,y]    : Target length between \"x\" (min) and \"y\" (max)
			gap=[x]      : Maximum of \"x\" indels in target region
			bad=[x]      : Maximum of \"x\" Ns in target region
			snp=[x,y]    : Min. of \"x\" and max \"y\" SNPs within --flank_dist
			mask=[x]     : Maximum of \"x\" proportion of masked bases allowed
			gc=[x,y]     : Proportion of G/C bases between \"x\" (min) and \"y\" (max)
			rand=[x]     : Randomly retain \"x\" target regions w/ baits
			pw=[i,q]     : Pairwise alignment, removing when \"i\" identity in \"q\" proportion
			blast_i=[i,q]: Only retain hits over \"i\" identity and \"q\" query coverage to provided db
			blast_x=[i,q]: Remove hits over \"i\" identity and \"q\" query coverage to provided db
			gff=[type]   : Only retain targets within \"d\" distance of \"type\" GFF elements
			               Use \"all\" to target all types, or see GFF3 specifications for SO feature types:
			               https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
			gff_a=[alias]: Only retain targets within \"d\" distance of GFF records with attributes of \"Alias\"
			               Note that \"gff\" and \"gffa\" options are case-insensitive
		Ex1: -F snp=1,10 -d 100 to sample when 1-10 SNPs w/in 100 bases
		Ex2: -F gc=0.2,0.8 -F rand=100 to randomly sample 100 targets with GC between 20-80%
		Ex3: -F mask=0.0,0.1 to remove targets with >10% \masked bases
		Ex4: -F gff=exon -d 1000 to constrain targets to within 1000 bases of exons
		NOTE: Values for gc, mask, pw, and blast* methods are given as proportions (e.g. 90% = 0.9)""")

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
			mask=[x]     : Maximum of \"x\" proportion of masked bases allowed
			gc=[x,y]     : Proportion of G/C bases between \"x\" (min) and \"y\" (max)
			pw=[i,q]     : Pairwise alignment, removing when \"i\" identity in \"q\" proportion
			blast_i=[i,q]: Only retain hits over \"i\" identity and \"q\" query coverage to provided db
			blast_x=[i,q]: Remove hits over \"i\" identity and \"q\" query coverage to provided db
			rand=[x]     : Randomly retain \"x\" baits""")

	print("""
VSEARCH Parameters (use when --select_b or --select_r = \"aln\"):

	--vsearch	: Path to VSEARCH executable if other than provided
		--MrBait will try to detect OS and appropriate exectable to use
	--vthreads	: Number of threads for VSEARCH
		--Will set to the value of <--threads> if not defined""")

	print("""
BLAST Parameters (use when --select_b or --select_r = \"blast\"):

	--blastdb	: Path and prefix to existing Blast database to use
	--fastadb	: Path and prefix to FASTA formatted sequences to make blast db
	--e_value	: Minimum e-value cutoff for reporting BLAST hits [0.000001]
	--gapopen	: Gap opening penalty for blastn [5]
	--gapextend	: Gap extension penatly for blastn [2]
	--word_size	: Word size for blastn
	--nodust	: Turn off low-complexity filter for blastn
	--megablast	: Use megablast rather than blastn
		--Default word size will be set to 28 if not provided
	--blastn	: Path to blastn binary if different than default
	--makedb	: Path to makeblastdb binary if different than default
		--MrBait will try to detect OS and appropriate exectable to use""")


	print("""
Output options:

	-x,--expand	: In output bait table, expand all ambiguities
		--Gaps are expanded as [ACGT] and absent
		--\"N\"s are expanded as [ACGT]
	--strand	: Output strand. Options: "+", "-", "both"
		--"-" print antiparallel (reverse complement)
	-t,--print_tr	: Boolean. Print target regions to FASTA file
	--print_loc	: Boolean. Print consensus locus catalog to FASTA file
	-o,--out	: Prefix for output files""")

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
			options, remainder = getopt.getopt(sys.argv[1:], 'M:G:V:L:A:hc:l:q:Q:b:w:Rm:v:n:Ng:E:D:p:S:F:s:f:xo:Pk:K:d:r:T:tX:', \
			["maf=","gff=","vcf=","loci=","assembly=",'help',"cov=","len=","thresh=", "max_ambig="
			"bait=","win_shift=","mult_reg","min_mult=","var_max=","numN=",
			"callN","numG=","callG","gff_type=","dist_r=","tile_min=",
			"select_r=","filter_r=", "threads=", "blastdb=", "fastadb=",
			"select_b=","filter_b=","quiet","expand","out=",
			"plot_all","mask=","max_mask=","flank_dist=","vsearch=",
			"vthreads=","hacker=", "evalue=", "e_value=", "gapopen=", "gapextend=",
			"word_size=", "megablast", "blastn=", "makedb=", "gap_extend=",
			"word=", "mega", "gap_open=", "blast_db=", "fasta_db=", "wordsize=", "nodust", "strand=",
			"resume=","db=", "print_tr", "xmfa=", "print_loc", "vcfALT"])
		except getopt.GetoptError as err:
			print(err)
			display_help("\nExiting because getopt returned non-zero exit status.")
			sys.exit(2)
		#Default values for params
		call_help=0 #boolean

		#Input params
		self.alignment=None
		self.gff=None
		self.vcf=None
		self.loci=None
		self.assembly=None
		self.xmfa=None
		self.vcfALT=False

		#Locus filtering params
		self.cov=1
		self.minlen=None
		self.thresh=0.1
		self.max_ambig=0.5
		self.mask=0.1
		self.max_mask=0.5

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
		self.vsearch = "vsearch"
		self.vthreads = None

		#BLAST options - contaminant removal and filtering by specificity
		self.blastdb=None
		self.fastadb=None
		self.evalue = 0.000001
		self.gapopen = 5
		self.gapextend=2
		self.word_size = None
		self.blast_method = "blastn"
		self.blastn = "blastn"
		self.makedb = "makeblastdb"
		self.nodust = None

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
		self.expand = None
		self.strand = "+"
		self.out = ""
		self.print_tr=False
		self.print_loc=False
		self.workdir = ""
		self.threads = 1

		self.ploidy=2
		self.db=""
		self.resume=None

		#HACKER ONLY OPTIONS
		self._noGraph = 0
		self._noWeightGraph = 0
		self._weightMax = 50000 #maximum size to attempt weighted edge resolution
		self._weightByMin = 0
		self._os = None



		#First pass to see if help menu was called
		for o, a in options:
			if o in ("-h", "-help", "--help"):
				display_help("\tExiting because help menu was called.")
				sys.exit(0)

		#Second pass to set all args.
		for opt, arg_raw in options:
			arg = arg_raw.replace(" ","")
			arg = arg.strip()

			if opt=='-M' or opt=='--maf':
				self.alignment = arg
			elif opt=='-X' or opt=='--xmfa':
				self.xmfa=arg
			elif opt=='-h' or opt=='--help':
				pass
			#Input params
			elif opt=='-G' or opt=='--gff':
				self.gff = arg
			elif opt=='-V' or opt== '--vcf':
				self.vcf = arg
			elif opt=='-L' or opt=='--loci':
				self.loci = arg
			elif opt =='-A' or opt== '--assembly':
				self.assembly = arg
			elif opt =="--vcfALT":
				self.vcfALT = True

			#Locus filtering params
			elif opt =='-c' or opt == '--cov':
				self.cov = int(arg)
			elif opt=='-l' or opt=='--len':
				self.minlen = int(arg)
			elif opt=='-q' or opt=='--thresh':
				self.thresh = float(arg)
			elif opt=='-Q' or opt=='--max_ambig':
				self.max_ambig = float(arg)
			elif opt=='-k' or opt=='--mask':
				self.mask = float(arg)
			elif opt=='-K' or opt=='--max_mask':
				self.max_mask = float(arg)

			#Bait general params
			elif opt=='-b' or opt=='--bait':
				self.blen = int(arg)
			elif opt=='-w' or opt=='--win_shift':
				self.win_shift = int(arg)
			elif opt=='-R' or opt=='--mult_reg':
				self.mult_reg = 1
			elif opt=='-m' or opt=='--min_mult':
				self.min_mult = int(arg)
			elif opt =='-v' or opt == '--var_max':
				self.var_max = int(arg)
			elif opt=='-n' or opt=='--numN':
				self.numN = int(arg)
			elif opt=='-g' or opt=='--numG':
				self.numG = int(arg)
			elif opt=='-E' or opt=='--gff_type':
				self.anchor = arg

			#target region opts
			elif opt=='-D' or opt=='--dist_r':
				self.dist_r = int(arg)
			elif opt=='-p' or opt=='--tile_min':
				self.tile_min = int(arg)
				self.tiling = 1
			elif opt=='-d' or opt=='--flank_dist':
				self.flank_dist = int(arg)
				assert isinstance(self.flank_dist, int), "<--flank_dist> must be an integer"
				assert self.flank_dist >= 0, "<--flank_dist> must be an integer greater than zero!"
			elif opt=='-S' or opt=='--select_r':
				temp = arg.split('=')
				assert len(temp) == 1, "Invalid specification for <--select_r>: %s"%arg
				self.select_r = (temp[0]).lower()
				chars = (['snp','bad','cons','rand'])
				if self.select_r not in chars:
					raise ValueError("Invalid option \"%r\" for <--select_r>" % self.select_r)
			elif opt=='-F' or opt=='--filter_r':
				self.filter_r = 1 #turn on region filtering
				#temp = arg.split('/') #parse region filtering options
				self.filter_r_whole = arg
				#for sub in temp:
				subopts = re.split('=|,',arg)
				if subopts[0] in ('snp','gc','len', "pw", "blast_i", "blast_x", "blast_a"):
					assert len(subopts) == 3, "Incorrect specification of option %r for <--filter_r>" %subopts[0]
					if subopts[0] == 'gc':
						assert 0.0 <= float(subopts[1]) < float(subopts[2]) <= 1.0, "In <--filter_r> suboption \"%s\": Min must be less than max"%subopts[0]
						self.filter_r_objects.append(subArg(subopts[0],float(subopts[1]),float(subopts[2])))
					elif subopts[0] == "len":
						assert subopts[1] < subopts[2], "In <--filter_r> suboption \"%s\": Min must be less than max"%subopts[0]
						self.filter_r_objects.append(subArg(subopts[0],int(subopts[1]),int(subopts[2])))
					elif subopts[0] in ("pw", "blast_x", "blast_i", "blast_a"):
						assert 0.0 <= float(subopts[1]) <= 1.0, "In <--filter_r> suboption \"%s\": Values must be given as proportions! "%subopts[0]
						assert 0.0 <= float(subopts[2]) <= 1.0, "In <--filter_r> suboption \"%s\": Values must be given as proportions! "%subopts[0]
						self.filter_r_objects.append(subArg(subopts[0],float(subopts[1]),float(subopts[2])))
					elif subopts[0] == "snp":
						minS = int(subopts[1])
						maxS = int(subopts[2])
						assert 0 <= minS < maxS, "In <--filter_r> suboption \"%s\": Min must be less than max"%subopts[0]
						assert isinstance(minS, int), "--filter_r suboption snp: Minimum must be an integer"
						assert isinstance(maxS, int), "--filter_r suboption snp: Maximum must be an integer"
						self.filter_r_objects.append(subArg(subopts[0],minS,maxS))

				elif subopts[0] in ('rand','gap','bad', 'mask', 'gff', 'gff_a'):
					assert len(subopts) == 2, "Incorrect specification of option %r for <--filter_r>" %subopts[0]
					if subopts[0] in ("gff", "gff_a"):
						self.filter_r_objects.append(subArg(subopts[0],subopts[1]))
					else:
						self.filter_r_objects.append(subArg(subopts[0],int(subopts[1])))
				else:
					bad_opts("Invalid option %r for <--filter_r>!" %subopts[0])

			#Bait selection options
			elif opt=='-s' or opt=='--select_b':
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
			elif opt=='-f' or opt=='--filter_b':
				self.filter_b = 1 #turn on region filtering
				#temp = arg.split('/') #parse region filtering options
				self.filter_b_whole = arg
				#for sub in temp:
				subopts = re.split('=|,',arg)
				if subopts[0] in ('pw','gc'):
					assert len(subopts) == 3, "Incorrect specification of option %r for <--filter_b>" %subopts[0]
					assert 0.0 <= float(subopts[1]) <= 1.0, "In <--filter_b> suboption \"%s\": Value must be between 0.0 and 1.0"%subopts[0]
					assert 0.0 <= float(subopts[2]) <= 1.0, "In <--filter_b> suboption \"%s\": Value must be between 0.0 and 1.0"%subopts[0]
					if subopts[0] == "gc":
						assert subopts[1] < subopts[2], "In <--filter_b> for suboptions \"mask\" and \"gc\": Min must be less than max"
					self.filter_b_objects.append(subArg(subopts[0],float(subopts[1]),float(subopts[2])))
				elif subopts[0] == 'mask':
					assert 0.0 <= float(subopts[1]) <= 1.0, "In <--filter_b> suboption \"%s\": Value must be between 0.0 and 1.0"%subopts[0]
					self.filter_b_objects.append(subArg(subopts[0],float(subopts[1])))
				elif (subopts[0] == 'rand'):
					assert len(subopts) == 2, "Incorrect specification of option %r for <--filter_b>" %subopts[0]
					self.filter_b_objects.append(subArg(subopts[0],subopts[1]))
				else:
					bad_opts("Invalid option %r for <--filter_b>!" %subopts[0])

			#Running options

			#vsearch options
			elif opt == "--vsearch":
				self.vsearch = str(arg)
			elif opt == "--vthreads":
				self.vthreads = int(arg)

			#BLAST options
			elif opt=='--blastdb' or opt=='--blast_db':
				self.blastdb = arg
			elif opt=='--fastadb' or opt=='--fasta_db':
				self.fastadb = arg
			elif opt=='--e_value' or opt=='--evalue':
				self.evalue = float(arg)
			elif opt=='--gapopen' or opt=='--gap_open':
				self.gapopen = int(arg)
			elif opt=='--gapextend' or opt=='--gap_extend':
				self.gapextend = int(arg)
			elif opt=='--word' or opt=='--word_size' or opt=='--wordsize':
				self.word_size = int(arg)
			elif opt=='--mega' or opt=='--megablast':
				self.blast_method = "megablast"
			elif opt == "--blastn":
				self.blastn = arg
			elif opt == "--makedb":
				self.makedb = arg
			elif opt == "--nodust":
				self.nodust = "TRUE"

			#output options
			elif opt=='-x' or opt=='--expand':
				self.expand = 1
			elif opt == "--strand":
				assert arg in ("+", "-", "both"), "Invalid option" + arg + "for <--strand>"
				self.strand = arg
			elif opt=='-r' or opt=='--resume':
				assert int(arg) in (0, 1, 2, 3, 4), "Invalid option" + arg + "for <-r, --resume>. Please specify a step (1-4) to resume pipeline."
				self.resume = int(arg)
			elif opt=='-o' or opt=='--out':
				self.out = arg
			elif opt=='-t' or opt=='--print_tr':
				self.print_tr = True
			elif opt=='--print_loc':
				self.print_loc = True
			elif opt == "--db":
				self.db = str(arg)
			elif opt=='-T' or opt=='--threads':
				self.threads = arg

			#HACKER ONLY OPTIONS
			elif opt=='--hacker':
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
				elif main == "os":
					assert len(subopts) == 2, "Warning: HACKER option <os> must have two arguments separated by \"=\""
					os_in = subopts[1]
					if os_in.lower() in ("linux", "ubuntu", "unknown"):
						self._os = "linux"
					elif os_in.lower() in ("mac", "apple", "macos", "darwin", "unix"):
						self_os = "darwin"
					else:
						assert False, "Unrecognized option %s for <--hacker os=>"%os_in
				else:
					assert False, "Unhandled option %r"%main
			else:
				assert False, "Unhandled option %r"%opt


		for subopt in self.filter_r_objects:
			if subopt.o1 in ("gff", "gff_a"):
				if not self.gff:
					sys.exit("ERROR: You have selected to filter targets on proximity to GFF elements, but have not provided a GFF file!")
			#print("filter_r: Suboption %s has parameters: %s %s" %(subopt.o1,subopt.o2,subopt.o3))

		#for subopt in self.filter_b_objects:
			#print("filter_b: Suboption %s has parameters: %s %s" %(subopt.o1,subopt.o2,subopt.o3))

		#Assertions and conditional changes to params
		if (self.alignment is None) and (self.loci is None) and (self.assembly is None) and (self.xmfa is None):
			display_help("Input not specified!")
			sys.exit(0)

		assert self.blen > 0, "Bait length cannot be less than or equal to zero!"
		assert 0.0 <= self.thresh <= 1.0, "Threshold for including ambiguous bases in consensus <-q, --thresh> must be between 0.0 and 1.0!"
		assert 0.0 <= self.mask <= 1.0, "Threshold for masking bases in consensus <-k, --mask> must be between 0.0 and 1.0!"
		assert 0.0 <= self.max_ambig <= 1.0, "Maximum proportion of ambiguous bases in a consensus <-Q, --max_ambig> must be between 0.0 and 1.0!"
		assert 0.0 <= self.max_mask <= 1.0, "Maximum proportion of masked bases in a consensus <-K, --max_mask> must be between 0.0 and 1.0!"

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
		if self.out == "":
			self.out = "baits"
		self.workdir = utils.getWorkingDir()
		if self.db == "":
			self.db = self.workdir + "/" + self.out + ".sqlite"
		#print("Database is:", self.db)
		#print("Prefix is: ", self.out)


		#Assert that win_shift cannot be larger than blen
		if self.minlen is None:
			self.minlen = self.blen
		elif self.blen > self.minlen:
			self.minlen = self.blen

		#set default of min_mult
		if self.min_mult is None:
			self.min_mult = self.minlen

		self._os = utils.getOS()
		if self._os == "darwin":
			self._os = "darwin (MacOS)"

		#If vsearch path not given, try to figure it out
		# if self.vsearch is None:
		# 	os_platform = utils.getOS()
		# 	#print("Found OS platform:", os_platform)
		# 	if os_platform == "linux" or os_platform == "unknown":
		# 		print("Automatically detected LINUX or UNKNOWN platform: Using LINUX VSEARCH executable.")
		# 		self.vsearch = utils.getScriptPath() + "/bin/vsearch-2.4.4-linux"
		# 	elif os_platform == "darwin": #mac os
		# 		print("Automatically detected MACOS platform: Using MACOS VSEARCH executable.")
		# 		self.vsearch = utils.getScriptPath() + "/bin/vsearch-2.4.4-macos"

		if self.vthreads is None:
			self.vthreads = self.threads


		#BLAST defaults
		if self.word_size is None:
			if self.blast_method == "megablast":
				self.word_size = 28
			else:
				self.word_size = 11

		# if self.blastn is None:
		# 	os_platform = utils.getOS()
		# 	#print("Found OS platform:", os_platform)
		# 	if os_platform == "linux" or os_platform == "unknown":
		# 		print("Automatically detected LINUX or UNKNOWN platform: Using LINUX BLASTN executable.")
		# 		self.blastn = utils.getScriptPath() + "/bin/ncbi-blastn-2.6.0-linux"
		# 	elif os_platform == "darwin": #mac os
		# 		print("Automatically detected MACOS platform: Using MACOS BLASTN executable.")
		# 		self.blastn = utils.getScriptPath() + "/bin/ncbi-blastn-2.6.0-macos"

		# if self.makedb is None:
		# 	os_platform = utils.getOS()
		# 	#print("Found OS platform:", os_platform)
		# 	if os_platform.lower() in ("linux", "ubuntu", "unknown", None):
		# 		print("Automatically detected LINUX or UNKNOWN platform: Using LINUX BLASTN executable.")
		# 		self.makedb = utils.getScriptPath() + "/bin/ncbi-makeblastdb-2.6.0-linux"
		# 	elif os_platform.lower() in ("darwin", "mac", "macos", "unix", "apple"): #mac os
		# 		print("Automatically detected MACOS platform: Using MACOS MAKEBLASTDB executable.")
		# 		self.makedb = utils.getScriptPath() + "/bin/ncbi-makeblastdb-2.6.0-macos"

		#Default bait design behavior
		if self.select_b == "tile" and self.overlap is None:
			self.overlap = self.blen // 2

		#Validate input types
		if self.assembly:
			if self.loci or self.alignment:
				sys.exit("ERROR: FASTA inputs cannot be provided with MAF or LOCI")
			#if self.mask or self.thresh:
				#print("WARNING: <-q,--thresh> and <-k,--mask> inputs have no affect when using a FASTA input file. Ignoring them.")
			self.cov = 0
		if self.gff or self.vcf:
			if not self.assembly:
				sys.exit("ERROR: VCF and GFF inputs require a FASTA assembly file.")
			if self.loci or self.alignment:
				sys.exit("ERROR: VCF and GFF files cannot be used with MAF and LOCI inputs.")
