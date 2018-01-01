#!/usr/bin/env python3

import os
import sys
import misc_utils as utils

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

#function to read a GFF file
#Generator function, yields individual elements
def read_gff(g):
	fields = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
	GFFRecord = namedtuple("GFFRecord", gffInfoFields)
	pass
