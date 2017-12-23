#!/usr/bin/env python3

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
	pass
