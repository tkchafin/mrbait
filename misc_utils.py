#!/usr/bin/python

import sys
import os
import os.path
import pandas as pd
import re

#Function to check if a file path is valid
def fileCheck(f):
	return (os.path.isfile(f))

#Function calculates union length of overlapping fixed length lines
def calculateUnionLengthFixed(n, l, o):
    #N is number of lines
    #L is length of lines
    #O is overlap of lines
    assert isinstance(n, int)
    assert isinstance(l, int)
    assert isinstance(o, int)
    return ((n*l)-((n-1)*o))

#Function to sanitize strings of any URLs, and replace them with "<REMOVED URL>"
def removeURL(st):
	return(re.sub(r'^https?:\/\/.*[\r\n]*', '', st, flags=re.MULTILINE))

#Function implementing fast way to replace single char in string
#This way is a lot faster than doing it by making a list and subst in list
def stringSubstitute(s, pos, c):
	return(s[:pos] + c + s[pos+1:])

#Function to return sorted unique string from list of chars
def listToSortUniqueString(l):
	l.sort()
	return(str(''.join(l)))

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
