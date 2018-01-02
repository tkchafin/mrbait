#!/usr/bin/python

from collections import Counter
from collections import defaultdict
import time
import unicodedata

"""
Conclusions:
- Just converting to upper() or lower() is faster than using unicode normalization function
- For loop and dictionaries are faster than using the collections.Counter class
- defaultdict does not substantially improve performance over the seqCounterLoop function
"""


def normalize_caseless(text):
    return unicodedata.normalize("NFKD", text.casefold())

def time_me(method):
    def wrapper(*args, **kw):
        startTime = int(round(time.time() * 1000))
        result = method(*args, **kw)
        endTime = int(round(time.time() * 1000))

        print(endTime - startTime,'ms')
        return result

    return wrapper

@time_me
def seqCounter(seq, num):
	seq_norm = normalize_caseless(seq)
	print("counter style:")
	for i in range(num):
		c = Counter(seq_norm)
		c['var'] = c['r'] + c['y'] + c['s'] + c['w'] \
		+ c['k'] + c['m'] + c['b'] + c['d'] + c['h'] + c['v']
	print(c)


@time_me
def seqCounter2(seq, num):
	seq_norm = seq.lower()
	print("counter 2 style:")
	for i in range(num):
		c = Counter(seq_norm)
		c['var'] = c['r'] + c['y'] + c['s'] + c['w'] \
		+ c['k'] + c['m'] + c['b'] + c['d'] + c['h'] + c['v']
	print(c)

@time_me
def seqCounterLoop(seq, num):
	d = {}
	seq_norm = seq.upper()
	print("loop style:")
	for i in range(num):
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
		for c in seq_norm:
			if c in d:
				d[c] += 1
		d['VAR'] = d['R'] + d['Y'] + d['S'] + d['W'] \
		+ d['K'] + d['M'] + d['B'] + d['D'] + d['H'] + d['V']
	print(d)

@time_me
def seqCounterDict(seq, num):
	print("dict style")
	seq_norm = seq.upper()
	for i in range(num):
		mydict = defaultdict(int)
		for c in seq_norm:
			mydict[c] += 1

seq = "ATGTGTAA-NRATTYRR-NNNAtattgygygwrttsgstttyn--agagg--gwtrrcacacccncacncgcc-ay-accdhaaca-vhVaaccannNNN"

seqCounter(seq, 100000)
seqCounter2(seq, 100000)
seqCounterLoop(seq, 100000)
seqCounterDict(seq, 100000)
