#!/usr/bin/python

import re

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
def gc_content(string):
	new = re.sub('[GCgc]','#',string)
	return sum(1 for c in new if c == '#')

#generator to create sliding windows by slicing out substrings
def seqSlidingWindowString(seq, shift, width=2):
	seqlen = len(seq)
	for i in range(0,seqlen,shift):
		if i+width > seqlen:
			j = seqlen
		else:
			j = i + width
		yield seq[i:j]
		if j==seqlen: break

#generator to create sliding windows by slicing out substrings, returns substring indices
def seqSlidingWindow(seq, shift, width=2):
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
