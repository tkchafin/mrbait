#!/usr/bin/python

import sys
import os
import pandas as pd

#Function calculates union length of overlapping fixed length lines
def calculateUnionLengthFixed(n, l, o):
    #N is number of lines
    #L is length of lines
    #O is overlap of lines
    assert isinstance(n, int)
    assert isinstance(l, int)
    assert isinstance(o, int)
    return ((n*l)-((n-1)*o))


#Function to count number of lower case in a string
def n_lower_chars(string):
    return sum(1 for c in string if c.islower())

#Function to return what OS platform we are using
def getOS():
	platform = "unknown"
	try:
		platform = sys.platform
	except OSError as e:
            # try to catch windows-specific issues... MrBait not currently supported in Windows!!
            if sys.platform.startswith('win'):
                if isinstance(e, WindowsError) and e.winerror == 267:
                    raise InvalidFile('uses Windows special name (%s)' % e)
            raise
	return(sys.platform)

#Function to get directory of running python executable
def getScriptPath():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

#Function to get call path
def getWorkingDir():
    return os.getcwd()


#Function to make dict from a 2-column pandas DF (col1=key, col2=val)
def dictFromDF(df):
	d = {}
	for row in df.iterrows():
		key = int(row[1][0])
		d[key] = int(row[1][1])
	return(d)
