#!/usr/bin/python

import getopt
import sys

def init_new_db(cursor, connection):
	
	cursor.execute('''DROP TABLE loci''')
	cursor.execute('''DROP TABLE annotations''')
	cursor.execute('''
		CREATE TABLE loci(id INTEGER PRIMARY KEY, depth INTEGER, 
			length INTEGER, sequence TEXT)
	''')
	cursor.execute('''
		CREATE TABLE annotations(id INTEGER PRIMARY KEY, ref, 
			length INTEGER, sequence TEXT)
	''')
	connection.commit()
