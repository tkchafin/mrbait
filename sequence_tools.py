#!/usr/bin/python

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

#generator to create sliding windows by slicing out substrings
def seqSlidingWindow(seq, shift, width=2):
	seqlen = len(seq)
	for i in range(0,seqlen,shift):
		if i+width > seqlen:
			j = seqlen
		else:
			j = i + width
		yield seq[i:j]
		if j==seqlen: break	

#Function to simplify a sequence to SNP, gaps, and Ns and 
def countSlidingWindow(seq, num, shift, width):
	for j in range(num):
		seq_temp = re.sub('[ACGT]', '', seq.upper())
		seq_norm = seq_temp.translate(str.maketrans("RYSWKMBDHV", "**********"))
		for i in windowSub(seq_norm, shift, width):
			#print(i)
			window_seq = "".join(i)
			seqCounterSimple(window_seq)
			
			
			
