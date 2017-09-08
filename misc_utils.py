#!/usr/bin/python

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
