#!/usr/bin/python

from collections import Counter
from collections import defaultdict
import time
import unicodedata
import os
import sys
import gzip
from collections import namedtuple
import misc_utils as utils
import alignment_tools as aln

#TODO: Decided I don't like any of the available GFF3 parsers so I'm going to
#write my own.

#function to read a GFF file
#Generator function, yields individual elements
def read_gff(g):
	fields = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
	GFFRecord = namedtuple("GFFRecord", gffInfoFields)

	gf = open(g)
	with gf as file_object:
		for line in file_object:
			if line.startswith("#"): continue
			line = line.strip()
			pieces = line.split("\t")
			if len(pieces) != len(gffInfoFields):
				sys.exit("Fatal error: GFF file is not standard-compatible")
			line = utils.removeURL(line) #Sanitize any URLs out
			
