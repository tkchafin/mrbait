#!/usr/bin/python

import sys
import os
import os.path
import pandas as pd
import re

#Function to calculate overlap of two 1D line segments
#Found on StackOverflow
#https://stackoverflow.com/questions/16691524/calculating-the-overlap-distance-of-two-1d-line-segments
def calcOverlap(min1, max1, min2, max2):
	return max(0, min(max1, max2) - max(min1, min2))


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

#Function for checking if two things with start and end coordinates overlap on single axis
#Line segment overlap
def checkOverlap(row1, row2, dist):
	#This could be made much more concise
	x2 = row2["start"]
	y2 = row2["stop"]
	x1 = (row1["start"]-dist)
	y1 = (row1["stop"]+dist)

	if x1 < 0:
		x1 = 0

	# print("x1=",x1)
	# print("y1=",y1)
	# print("x2=",x2)
	# print("y2=",y2)
	#Left edge of row1 overlaps right edge of row2
	if (x2 < x1) and (y2 >= x1):
		return 1
	#Right edge of row1 overlaps left edge of row2
	elif (x2 <= y1) and (y2 > y1):
		return 1
	#Row 2 completely overlapped by row1
	elif (x2 >= x1 and y2 > x1) and (x2 < y1 and y2 <= y1):
		return 1
	else:
		return 0

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
