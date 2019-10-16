#!/usr/bin/python

import sys
import random
from mrbait import misc_utils as utils
from Bio import AlignIO

############################# CLASSES ##################################

class consensAlign():
	'Consensus alignment object'
	#Default constructor
	def __init__(self, alignment, threshold, mask, maf=0.0):
		self.alnVars = []
		self.conSequence = make_consensus(alignment, threshold, mask, maf)
		self.alnVars = get_vars(self.conSequence)

class variablePosition():
	'Object to hold information about a variable position'
	#Default constructor
	def __init__(self, pos=None, val=None):
		self.position = pos
		self.value = val.upper()


	@classmethod
	def from_list(cls, data):
		pos = int(data[0])
		val = str(data[1])
		new = cls(pos, val)
		return new



######################## STATIC FUNCTIONS ##############################

#Less shitty consensus function than BioPython has..
#From an AlignIO alignment object
def make_consensus(alignment, threshold, mask, maf=0.0):
	aln_depth = len(alignment)
	#If only one sequence in alignment, return that seq as consensus
	if aln_depth == 1:
		return alignment[:,0]
	aln_len = alignment.get_alignment_length()
	consensus="" #consensus string to build
	#For each column
	for i in range(aln_len):
		#For each character in column, grab occurences
		col_len = len((alignment[:,i]))
		nuc_count = dict() #dictionary to track occurences
		nuc_types = 0 #track number of things we found
		ismask = 0 #Track if we should mask this column
		#Check number of masked bases
		nlower = utils.n_lower_chars(alignment[:,i])
		prop_mask = float(nlower/aln_depth)
		if prop_mask > mask:
			ismask = 1

		for c in ((alignment[:,i]).upper()):
			for nuc in get_iupac(c):
				nuc_count[nuc] = nuc_count.get(nuc,0)+1
				#print(nuc, end='', flush=True)
		nucs= []
		add = 0
		#print(nuc_count)
		bad_count = 0
		bad_counts = dict()
		for key_raw in nuc_count:
			key = key_raw
			if ismask: #If masked above threshold, set to lower case
				key = key_raw.lower()
			#print(key, " is ", nuc_count[key])
			#If only one nuc type, keep it
			if nuc_types == 1:
				consensus+=str(key)
				add=1
				break
			#If N or gap, call consensus N or gap if above threshold
			elif key in ("N", "n", "-"):
				temp = key
				if key == "n":
					temp = "N"
				bad_count += nuc_count[temp]
				bad_counts[temp] = nuc_count[temp]
				if float((nuc_count[temp])/aln_depth*2) >= threshold:
					#print("Found ", nuc_count[key], key, "'s in alignment")
					#print((nuc_count[key])/aln_depth)
					consensus+=str(key)
					add=1
					break
			else:
				nucs.append(str(key))
		if add == 0:
			if bad_count >= col_len:
				#print("All bad!")
				nuc = utils.getMaxKey(bad_counts)
				#print(nuc)
				consensus+=nuc
			else:
				#call consensus from alleles!
				temp = utils.listToSortUniqueString(nucs)
				if len(temp) > 1 and maf > 0.0:
					#print(temp)
					#first, remove alleles below MAF filter
					nucs_clean = filterListByMAF(nucs, nuc_count, aln_depth*2, maf)
					#recreate unique list
					temp = utils.listToSortUniqueString(nucs_clean)
					#print(temp)
					#print()
				#print(nucs,":",temp)
				consensus+=reverse_iupac_case(temp)
	return(consensus)


def filterListByMAF(l, counts, depth, threshold):
	if threshold <= 0.0:
		return(l)
	elif threshold >= 1.0:
		return(l)
	else:
		ret = []
		for c in utils.listToSortUniqueString(l):
			print(c, depth, counts[c], float((float(counts[c]))/float(depth)))
			if float((float(counts[c]))/float(depth)) >= threshold:
				ret.append(c)
		return(ret)
		#sys.exit(0)


#Function to get a list of variablePositions
def get_vars(con):
	#print("Parsing: ", con)
	var_objects = [] #empty list for holding var objects
	#For each position
	for i in range(len(con)):
		#print(i, " is ", con[i])
		#If not monomorphic
		if con[i].upper() not in {"A", "G", "T", "C", "X"}:
			var_objects.append(variablePosition(i, con[i].upper()))
			#continue
		#else:
		'''
		#For each sequence
		for c in range(len(aln[:,i])):
			#print(aln[c,i], end='', flush=True)
			ref = con[i].upper()
			var = aln[c,i].upper()
			#ref.upper()
			#var.upper()
			#print("Var is ",var, " and Ref is ", ref)
			if var == ref:
			#if var == "-" and ref == "-":
				continue
			#elif var == "N" and ref == "N":
				continue
			else:
				#print(var, end='', flush=True)
				#print(aln[c].id, " has ", aln[c,i], " at pos ", i)
				var_objects.append(variablePosition(aln[c].id, i, aln[c,i]))
		'''
	return var_objects

#Function to split character to IUPAC codes, assuming diploidy
def get_iupac(char):
	iupac = {
		"A"	: ["A", "A"],
		"G"	: ["G", "G"],
		"C"	: ["C", "C"],
		"T"	: ["T", "T"],
		"N"	: ["N", "N"],
		"-"	: ["-", "-"],
		"R"	: ["A","G"],
		"Y"	: ["C","T"],
		"S"	: ["G","C"],
		"W"	: ["A","T"],
		"K"	: ["G","T"],
		"M"	: ["A","C"],
		"B"	: random.sample(["C","G","T"], 2),
		"D"	: random.sample(["A","G","T"], 2),
		"H"	: random.sample(["A","C","T"], 2),
		"V"	: random.sample(["A","C","G"], 2)
	}
	return iupac[char]

#Function to translate a string of bases to an iupac ambiguity code
def reverse_iupac(char):
	iupac = {
		'A':'A',
		'N':'N',
		'-':'-',
		'C':'C',
		'G':'G',
		'T':'T',
		'AG':'R',
		'CT':'Y',
		'AC':'M',
		'GT':'K',
		'AT':'W',
		'CG':'S',
		'CGT':'B',
		'AGT':'D',
		'ACT':'H',
		'ACG':'V',
		'ACGT':'N'
	}
	return iupac[char]

#Function to translate a string of bases to an iupac ambiguity code, retains case
def reverse_iupac_case(char):
	iupac = {
		'A':'A',
		'N':'N',
		'-':'-',
		'C':'C',
		'G':'G',
		'T':'T',
		'AG':'R',
		'CT':'Y',
		'AC':'M',
		'GT':'K',
		'AT':'W',
		'CG':'S',
		'CGT':'B',
		'AGT':'D',
		'ACT':'H',
		'ACG':'V',
		'ACGT':'N',
		'a':'a',
		'n':'n',
		'c':'c',
		'g':'g',
		't':'t',
		'ag':'r',
		'ct':'y',
		'ac':'m',
		'gt':'k',
		'at':'w',
		'cg':'s',
		'cgt':'b',
		'agt':'d',
		'act':'h',
		'acg':'v',
		'acgt':'n'
	}
	return iupac[char]
