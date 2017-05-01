#!/usr/bin/python

import sys
import sqlite3

#Function to create database connection
#Add code later to enable re-running from existing database
def create_connection(db):
	try:
		conn = sqlite3.connect(db)
		return conn
	except Error as e:
		print(e)
	return None

#Initialize empty databases
def init_new_db(connection):
	cursor = connection.cursor()
	cursor.execute('''DROP TABLE IF EXISTS loci''')
	cursor.execute('''DROP TABLE IF EXISTS annotations''')
	cursor.execute('''
		CREATE TABLE loci(id INTEGER PRIMARY KEY, depth INTEGER NOT NULL, 
			length INTEGER NOT NULL, consensus TEXT NOT NULL, pass INTEGER NOT NULL)
	''')
	cursor.execute('''
		CREATE TABLE annotations(id INTEGER PRIMARY KEY, ref, 
			length INTEGER, sequence TEXT)
	''')
	connection.commit()

#Code to add record to 'loci' table
def add_locus_record(conn, depth, consensus, passed=0):
	stuff = []*4
	sql = ''' INSERT INTO loci(depth, length, consensus, pass) 
				VALUES(?,?,?,?) '''
	stuff.insert(0, depth)
	stuff.insert(1, int(len(consensus)))
	stuff.insert(2, str(consensus))
	stuff.insert(3, int(passed))
	cur = conn.cursor()
	cur.execute(sql, stuff)
	return cur.lastrowid

#Code to add to 'annotations' table
def add_annotation_record(conn, name, pos, val):
	stuff = [name,pos,val]
	sql = ''' INSERT INTO annotations(depth, length, consensus, pass) 
			VALUES(?,?,?,?) '''
