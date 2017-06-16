#!/usr/bin/python

import sys
import sqlite3
import pandas as pd

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
	cursor.execute('''DROP TABLE IF EXISTS variants''')
	cursor.execute('''DROP TABLE IF EXISTS samples''')
	#Table holding records for each locus
	cursor.execute('''
		CREATE TABLE loci(id INTEGER PRIMARY KEY, depth INTEGER NOT NULL, 
			length INTEGER NOT NULL, consensus TEXT NOT NULL, pass INTEGER NOT NULL,
			UNIQUE(id))
	''')
	#Table holding variant information
	cursor.execute('''
		CREATE TABLE variants(varid INTEGER NOT NULL, locid INTEGER NOT NULL, 
			sampid INTEGER, column INTGER NOT NULL, value TEXT NOT NULL,
			FOREIGN KEY (locid) REFERENCES loci(id), 
			FOREIGN KEY (sampid) REFERENCES samples(sampid),
			PRIMARY KEY(varid), 
			UNIQUE(varid))
	''')
	
	#Table holding records for each locus
	cursor.execute('''
		CREATE TABLE samples(sampid INTEGER NOT NULL, name TEXT NOT NULL, 
			PRIMARY KEY (sampid),
			UNIQUE(name))
	''')
	
	connection.commit()

#Code to add record to 'loci' table
def add_locus_record(conn, depth, consensus, passed=0):
	stuff = [depth, int(len(consensus)), str(consensus), int(passed)]
	sql = ''' INSERT INTO loci(depth, length, consensus, pass) 
				VALUES(?,?,?,?) '''
	cur = conn.cursor()
	cur.execute(sql, stuff)
	return cur.lastrowid

#Code to add to 'variants' table
def add_variant_record(conn, loc, name, pos, val):
	#Establish cursor 
	cur = conn.cursor()
	
	#Check if sample has a sampid, fetch if it exists
	sql= "INSERT OR IGNORE INTO samples(name) VALUES (%r);"%name
	cur.execute(sql)
	fetch = "SELECT sampid FROM samples WHERE name = %r"%name
	cur.execute(fetch)
	sampid = cur.fetchone()[0] #fecth the sample id
	
	#Check if position has a posid, fetch if it exists
	sql2 = '''INSERT OR IGNORE INTO variants(locid, sampid, column, value) VALUES(?,?,?,?)'''
	stuff = [loc, sampid, pos, val]
	cur.execute(sql2,stuff)

	










