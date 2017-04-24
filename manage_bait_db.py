#!/usr/bin/python

import sys

def init_new_db(cursor, connection):
	
	cursor.execute('''DROP TABLE IF EXISTS loci''')
	cursor.execute('''DROP TABLE IF EXISTS annotations''')
	cursor.execute('''
		CREATE TABLE loci(id INTEGER PRIMARY KEY, depth INTEGER, 
			length INTEGER, consensus TEXT, refseq TEXT, pass INTEGER)
	''')
	cursor.execute('''
		CREATE TABLE annotations(id INTEGER PRIMARY KEY, ref, 
			length INTEGER, sequence TEXT)
	''')
	connection.commit()
