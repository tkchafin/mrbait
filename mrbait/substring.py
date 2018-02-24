#!/usr/bin/python
from random import randint
import operator

#Object for creating an iterable slidinw window sampling
class SubString():
    def __init__(self):
        self.string = "NULL"
        self.start = 0
        self.stop = 0

    def getString(self):
        return self.string

    def checkMatch(self, stuff, overlap):
        if len(stuff) <= 0:
            return 0
        #For each thing in the list
        for other in stuff:
            #Check if they overlap:
            if self.calcOverlap(other) >= overlap:
                #If they DO overlap, return 1/TRUE
                #print("They overlap by more than ", overlap)
                #Return 1 if an overlap is found
                return 1
        #Return 0 if no overlaps found
        return 0

    def calcOverlap(self, other):
        #Returns TRUE if two lines do not overlap by given amount
        a = other.start
        b = other.stop
        x = self.start
        y = self.stop
        over = max(0, (min(b, y) - max(a, x)))
        #print("Calc overlap for: ", x, y, " -- ", a, b, " : ", over)
        return(over)

    def randomDrawSubstring(self, source, size):
        size = int(size)
        if len(source) < size:
            raise ValueError("SubString.randomDrawSubstring: Length of source string cannot be greater than length of desired substring")
        elif len(source) == size:
            self.string = source
            self.start = 0
            self.stop = len(source)-1
        else:
            limit = int(len(source) - size)
            first = int(randint(0,limit))
            last = int(first + size)
            seq = source[first:last]
            self.string = seq
            self.start = first
            self.stop = last
            #print(seq)

    @staticmethod
    #Function to sort a given list of SubString objects
    def sortSubStrings(stuff):
        print(stuff.sort(key = lambda x: x.start))
