#!/usr/bin/python

from collections import Counter
from collections import defaultdict
import time
import unicodedata
import os
import sys
import pandas as pd
from collections import namedtuple

"""
Benchmarking two different structures for yielding GFF Records

Results:

-At small numbers, namedtuples and native lists are much faster.
-namedtuple and pandas df do have some more up-front cost in defining them
	But this is inconsequential if the definition won't be repeated
-Adding in 1 access (of seqid) per iteration didn't change much.
-At high number of replicates, they all seem to take about the same time

Testing with namedtuples: 100 iterations
1 ms
Testing with namedtuples: 1000 iterations
2 ms
Testing with namedtuples: 10000 iterations
21 ms
Testing with namedtuples: 100000 iterations
194 ms
Testing with namedtuples: 1000000 iterations
1962 ms

Testing with pandas df: 100 iterations
12 ms
Testing with pandas df: 1000 iterations
13 ms
Testing with pandas df: 10000 iterations
31 ms
Testing with pandas df: 100000 iterations
207 ms
Testing with pandas df: 1000000 iterations
1939 ms

Testing with dict df: 100 iterations
0 ms
Testing with dict df: 1000 iterations
2 ms
Testing with dict df: 10000 iterations
20 ms
Testing with dict df: 100000 iterations
196 ms
Testing with dict df: 1000000 iterations
1942 ms

Conclusion:
-probably will go with a native dict since namedtuple doesn't seem to speed
	things up much. Pandas DF seems overboard for this simple container case
-namedtuple and dict should also theoretically take same space in mem,
	although in this case we should only have 1 record in mem per thread, so
	this is negligible anyways 


"""

def time_me(method):
    def wrapper(*args, **kw):
        startTime = int(round(time.time() * 1000))
        result = method(*args, **kw)
        endTime = int(round(time.time() * 1000))

        print(endTime - startTime,'ms')
        return result

    return wrapper

#Function to split GFF attributes
def splitAttributes(a):
	ret = {}
	for thing in a.split(";"):
		stuff = thing.split("=")
		if len(stuff) != 2: continue #Eats error silently, YOLO
		key = stuff[0]
		value = stuff[1]
		ret[key] = value
	return ret

#function to read each GFF element into a namedtuple Collections object
@time_me
def read_gff_TUPLE(g, num):
	fields = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
	GFFTuple = namedtuple("GFFRecord", fields)

	gf = open(g)
	with gf as file_object:
		for i in range(num):
			for line in file_object:
				if line.startswith("#"): continue
				line = line.strip()
				things = line.split("\t")
				if len(things) != len(fields):
					sys.exit("Fatal error: GFF file is not standard-compatible")
				#line = utils.removeURL(line) #Sanitize any URLs out
				normalizedInfo = {
					"seqid": None if things[0] == "." else things[0],
					"source": None if things[1] == "." else things[1],
					"type": None if things[2] == "." else things[2],
					"start": None if things[3] == "." else int(things[3]),
					"end": None if things[4] == "." else int(things[4]),
					"score": None if things[5] == "." else float(things[5]),
					"strand": None if things[6] == "." else things[6],
					"phase": None if things[7] == "." else things[7],
					"attributes": None if things[8] == "." else splitAttributes(things[8])
				}
				ret = GFFTuple(**normalizedInfo)
				test = GFFTuple.seqid
				test2=GFFTuple.phase
				ret = ""
	gf.close()


#Function to read each GFF element into a pandas df
@time_me
def read_gff_PANDAS(g, num):
	fields = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
	GFFTuple = pd.DataFrame(columns=fields)

	g2 = open(g)
	with g2 as file_object:
		for i in range(num):
			for line in file_object:
				if line.startswith("#"): continue
				line = line.strip()
				things = line.split("\t")
				if len(things) != len(fields):
					sys.exit("Fatal error: GFF file is not standard-compatible")
				#line = utils.removeURL(line) #Sanitize any URLs out

				GFFTuple["seqid"] = None if things[0] == "." else things[0]
				GFFTuple["source"] = None if things[1] == "." else things[1]
				GFFTuple["type"] =  None if things[2] == "." else things[2]
				GFFTuple["start"] = None if things[3] == "." else int(things[3])
				GFFTuple["end"] = None if things[4] == "." else int(things[4])
				GFFTuple["score"] = None if things[5] == "." else float(things[5])
				GFFTuple["strand"] = None if things[6] == "." else things[6]
				GFFTuple["phase"] = None if things[7] == "." else things[7]
				GFFTuple["attributes"] = None if things[8] == "." else things[8]
				things = splitAttributes(str(GFFTuple["attributes"]))

				ret = GFFTuple
				test = GFFTuple["seqid"]
				test2 = GFFTuple["phase"]
				ret = "" #This is where the yield would go
	g2.close()

#Function to read each GFF element into a native dict
@time_me
def read_gff_DICT(g, num):
	GFFTuple = {}
	g3 = open(g)
	with g3 as file_object:
		for i in range(num):
			for line in file_object:
				if line.startswith("#"): continue
				line = line.strip()
				things = line.split("\t")
				if len(things) != 9:
					sys.exit("Fatal error: GFF file is not standard-compatible")
				#line = utils.removeURL(line) #Sanitize any URLs out

				GFFTuple["seqid"] = None if things[0] == "." else things[0]
				GFFTuple["source"] = None if things[1] == "." else things[1]
				GFFTuple["type"] =  None if things[2] == "." else things[2]
				GFFTuple["start"] = None if things[3] == "." else int(things[3])
				GFFTuple["end"] = None if things[4] == "." else int(things[4])
				GFFTuple["score"] = None if things[5] == "." else float(things[5])
				GFFTuple["strand"] = None if things[6] == "." else things[6]
				GFFTuple["phase"] = None if things[7] == "." else things[7]
				GFFTuple["attributes"] = None if things[8] == "." else splitAttributes(things[8])

				ret = GFFTuple
				test = GFFTuple["seqid"]
				test2 = GFFTuple["phase"]
				ret = "" #This is where the yield would go
	g3.close()

#MAIN
gff = "../examples/example.gff"

print("Testing with namedtuples: 100 iterations")
read_gff_TUPLE(gff, 100)
print("Testing with namedtuples: 1000 iterations")
read_gff_TUPLE(gff, 1000)
print("Testing with namedtuples: 10000 iterations")
read_gff_TUPLE(gff, 10000)
print("Testing with namedtuples: 100000 iterations")
read_gff_TUPLE(gff, 100000)
print("Testing with namedtuples: 1000000 iterations")
read_gff_TUPLE(gff, 1000000)

print()

print("Testing with pandas df: 100 iterations")
read_gff_PANDAS(gff, 100)
print("Testing with pandas df: 1000 iterations")
read_gff_PANDAS(gff, 1000)
print("Testing with pandas df: 10000 iterations")
read_gff_PANDAS(gff, 10000)
print("Testing with pandas df: 100000 iterations")
read_gff_PANDAS(gff, 100000)
print("Testing with pandas df: 1000000 iterations")
read_gff_PANDAS(gff, 1000000)

print()

print("Testing with dict df: 100 iterations")
read_gff_DICT(gff, 100)
print("Testing with dict df: 1000 iterations")
read_gff_DICT(gff, 1000)
print("Testing with dict df: 10000 iterations")
read_gff_DICT(gff, 10000)
print("Testing with dict df: 100000 iterations")
read_gff_DICT(gff, 100000)
print("Testing with dict df: 1000000 iterations")
read_gff_DICT(gff, 1000000)
