#!/usr/bin/python

import re
import sys
from itertools import product

#Function to split character to IUPAC codes, assuing diploidy
def get_iupac_caseless(char):
	lower = False
	if char.islower():
		lower = True
		char = char.upper()
	iupac = {
		"A"	: ["A"],
		"G"	: ["G"],
		"C"	: ["C"],
		"T"	: ["T"],
		"N"	: ["A", "C", "G", "T"],
		"-"	: ["A", "C", "G", "T", "-"],
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
	ret = iupac[char]
	if lower:
		ret = [c.lower() for c in ret]
	return ret


#Function to expand ambiguous sequences
#Generator function
def expandAmbiquousDNA(sequence):
   for i in product(*[get_iupac_caseless(j) for j in sequence]):
      yield("".join(i))

#Function to return reverse complement of a nucleotide, while preserving case
def get_revComp_caseless(char):
	lower = False
	if char.islower():
		lower = True
		char = char.upper()
	d = {
		"A"	: "T",
		"G"	: "C",
		"C"	: "G",
		"T"	: "A",
		"N"	: "N",
		"-"	: "-",
		"R"	: "Y",
		"Y"	: "R",
		"S"	: "S",
		"W"	: "W",
		"K"	: "M",
		"M"	: "K",
		"B"	: "V",
		"D"	: "H",
		"H"	: "D",
		"V"	: "B"
	}
	ret = d[char]
	if lower:
		ret = ret.lower()
	return ret


#Function to reverse complement a sequence, with case preserved
def reverseComplement(seq):
	comp = []
	for i in (get_revComp_caseless(j) for j in seq):
		comp.append(i)
	return("".join(comp[::-1]))



#Function to simplify a sequence
def simplifySeq(seq):
	temp = re.sub('[ACGT]', '', (seq).upper())
	return temp.translate(str.maketrans("RYSWKMBDHV", "**********"))

#returns dict of character counts
def seqCounter(seq):
	d = {}
	d = {
		'A':0,
		'N':0,
		'-':0,
		'C':0,
		'G':0,
		'T':0,
		"R"	: 0,
		"Y"	: 0,
		"S"	: 0,
		"W"	: 0,
		"K"	: 0,
		"M"	: 0,
		"B"	: 0,
		"D"	: 0,
		"H"	: 0,
		"V"	: 0
	}
	for c in seq:
		if c in d:
			d[c] += 1
	d['VAR'] = d['R'] + d['Y'] + d['S'] + d['W'] \
	+ d['K'] + d['M'] + d['B'] + d['D'] + d['H'] + d['V']
	return d

#Returns dict of character counts from a simplified consensus sequence
def seqCounterSimple(seq):
	d = {}
	d = {
		'.':0,
		'N':0,
		'-':0,
		'*':0
	}
	for c in seq:
		if c in d:
			d[c] += 1
	return d

#Function to get GC content of a provided sequence
def gc_counts(string):
	new = re.sub('[GCgc]','#',string)
	return sum(1 for c in new if c == '#')

#Function to get counts of masked bases
def mask_counts(string):
	return sum(1 for c in string if c.islower())


#Function to get GC content as proportion
def gc_content(string):
	new = re.sub('[GCgc]','#',string)
	count = sum(1 for c in new if c == '#')
	return(count/(len(string)))

#Function to count number of lower case in a string
def mask_content(string):
	count = sum(1 for c in string if c.islower())
	return(count/(len(string)))

#generator to create sliding windows by slicing out substrings
def seqSlidingWindowString(seq, shift, width):
	seqlen = len(seq)
	for i in range(0,seqlen,shift):
		if i+width > seqlen:
			j = seqlen
		else:
			j = i + width
		yield seq[i:j]
		if j==seqlen: break

#generator to create sliding windows by slicing out substrings, returns substring indices
def seqSlidingWindow(seq, shift, width):
	seqlen = len(seq)
	for i in range(0,seqlen,shift):
		if i+width > seqlen:
			j = seqlen
		else:
			j = i + width
		yield [seq[i:j], i, j]
		if j==seqlen: break


#Function to simplify a sequence to SNP, gaps, and Ns and get counts of sliding windows
def countSlidingWindow(seq, shift, width):
	seq_temp = re.sub('[ACGT]', '', seq.upper())
	seq_norm = seq_temp.translate(str.maketrans("RYSWKMBDHV", "**********"))
	for i in windowSub(seq_norm, shift, width):
		#print(i)
		window_seq = "".join(i)
		seqCounterSimple(window_seq)

#Object for creating an iterable slidinw window sampling
class slidingWindowGenerator():
	#Need to come back and comment better...
	def __init__(self, seq, shift, width):
		self.__seq = seq
		self.__seqlen = len(self.__seq)
		self.__shift = shift
		self.__width = width
		self.__i = 0

	def __call__(self):
		self.__seqlen
		while self.__i < self.__seqlen:
			#print("i is ", self.__i, " : Base is ", self.__seq[self.__i]) #debug print
			if self.__i+self.__width > self.__seqlen:
				j = self.__seqlen
			else:
				j = self.__i + self.__width
			yield [self.__seq[self.__i:j], self.__i, j]
			if j==self.__seqlen: break
			self.__i += self.__shift

	@classmethod
	def setI(self, value):
		self.__i = value

	@classmethod
	def getI(self):
		return self.__i

	def __getattr__(self, name):
		return getattr(self.__class__, name)
	def __setattr__(self, name, val):
		return setattr(self.__class__, name, val)
