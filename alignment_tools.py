#!/usr/bin/python

import sys
import misc_utils as utils
from Bio import AlignIO

############################# CLASSES ##################################

class consensAlign():
	'Consensus alignment object'
	#Default constructor
	def __init__(self, alignment, threshold=0.1, mask=0.1):
		self.alnVars = []
		self.conSequence = make_consensus(alignment, threshold, mask)
		self.alnVars = self.get_vars(self.conSequence, alignment)

	@staticmethod
	def get_vars(con, aln):
		#print("Parsing: ", con)
		var_objects = [] #empty list for holding var objects
		#For each position
		for i in range(len(con)):
			#print(i, " is ", con[i])
			#If not monomorphic
			if con[i] not in {"A", "G", "T", "C"}:
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
def make_consensus(alignment, threshold=0.1, mask=0.1):
	aln_depth = len(alignment)
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
				if float((nuc_count[temp])/aln_depth) >= threshold:
					#print("Found ", nuc_count[key], key, "'s in alignment")
					#print((nuc_count[key])/aln_depth)
					consensus+=str(key)
					add=1
					break
			else:
				nucs.append(str(key))
		if add == 0:
			nucs.sort()
			#print (nucs)
			temp = str(''.join(nucs))
			#print(temp)
			consensus+=reverse_iupac_case(temp)
	return(consensus)

#Function to split character to IUPAC codes, assuing diploidy
def get_iupac(char):
	iupac = {
		"A"	: ["A"],
		"G"	: ["G"],
		"C"	: ["C"],
		"T"	: ["T"],
		"N"	: ["N"],
		"-"	: ["-"],
		"R"	: ["A","G"],
		"Y"	: ["C","T"],
		"S"	: ["G","C"],
		"W"	: ["A","T"],
		"K"	: ["G","T"],
		"M"	: ["A","C"],
		"B"	: ["C","G","T"],
		"D"	: ["A","G","T"],
		"H"	: ["A","C","T"],
		"V"	: ["A","C","G"]
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
