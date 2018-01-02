#!/usr/bin/env python3

import os
import sys
import urllib.parse

#Function to split GFF attributes
def splitAttributes(a):
	ret = {}
	for thing in a.split(";"):
		stuff = thing.split("=")
		if len(stuff) != 2: continue #Eats error silently, YOLO
		key = stuff[0].lower()
		value = stuff[1].lower()
		ret[key] = value
	return ret

#Function to return a GFF record as a dict
def GFFRecordAsDict(things):
	rec = {}
	#Load up dict, and sanitize inputs
	rec["seqid"] = None if things[0] == "." else urllib.parse.unquote(things[0])
	rec["source"] = None if things[1] == "." else urllib.parse.unquote(things[1])
	rec["type"] =  None if things[2] == "." else urllib.parse.unquote(things[2])
	rec["start"] = None if things[3] == "." else int(things[3])
	rec["end"] = None if things[4] == "." else int(things[4])
	rec["score"] = None if things[5] == "." else float(things[5])
	rec["strand"] = None if things[6] == "." else urllib.parse.unquote(things[6])
	rec["phase"] = None if things[7] == "." else urllib.parse.unquote(things[7])
	rec["attributes"] = None if things[8] == "." else splitAttributes(urllib.parse.unquote(things[8]))
	return rec

#function to read a GFF file
#Generator function, yields individual elements
def read_gff(g):
	bad = 0 #tracker for if we have bad lines
	gf = open(g)
	try:
		with gf as file_object:
			for line in file_object:
				if line.startswith("#"): continue
				line = line.strip() #strip leading/trailing whitespace
				if not line: #skip empty lines
					continue
				things = line.split("\t") #split lines
				if len(things) != 9:
					if bad == 0:
						print("Warning: GFF file does not appear to be standard-compatible. See https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md")
						bad = 1
					elif bad == 1:
						sys.exit("Fatal error: GFF file does not appear to be standard-compatible. See https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md")
				#line = utils.removeURL(line) #Sanitize any URLs out
				GFFRecord = GFFRecordAsDict(things)

				yield(GFFRecord)
	finally:
		gf.close()
