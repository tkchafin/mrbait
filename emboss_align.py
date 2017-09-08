#!/usr/bin/python

from Bio.Emboss.Applications import NeedleCommandline
from Bio.Emboss.Applications import WaterCommandline

#Function for calling needle from EMBOSS
def needleAlign(seq1, seq2, gapopen, gapextend):
    needle = NeedleCommandline()
    needle.asequence = seq1
    needle.bsequence = seq2
    needle.gapopen=gapopen
    needle.gapextend = gapextend
    needle.outfile = "needle.txt"

    stdout, stderr = needle()
    print(stdout)


#Function for calling water from EMBOSS
def waterAlign(seq1, seq2, gapopen, gapextend):
    water = WaterCommandline()
    water.asequence = seq1
    water.bsequence = seq2
    water.gapopen=gapopen
    water.gapextend = gapextend
    water.outfile = "needle.txt"

    stdout, stderr = water()
    print(stdout)
