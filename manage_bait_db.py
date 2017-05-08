#!/usr/bin/python

import sys
import sqlite3
import pandas as p

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
	cursor.execute('''DROP TABLE IF EXISTS positions''')
	#Table holding records for each locus
	cursor.execute('''
		CREATE TABLE loci(id INTEGER PRIMARY KEY, depth INTEGER NOT NULL, 
			length INTEGER NOT NULL, consensus TEXT NOT NULL, pass INTEGER NOT NULL)
	''')
	#Table holding position identifiers given to each variant column in a locus
	cursor.execute('''
		CREATE TABLE positions(posid INTEGER NOT NULL, locid INTEGER NOT NULL, 
			column INTGER NOT NULL, 
			FOREIGN KEY (locid) REFERENCES loci(id), 
			PRIMARY KEY(posid))
	''')
	
	#cursor
	#Table holding extracted variants from each alignment 
	cursor.execute('''
		CREATE TABLE variants(posid INTEGER NOT NULL, name TEXT NOT NULL, 
			value TEXT NOT NULL, 
			FOREIGN KEY (posid) REFERENCES positions(posid), 
			PRIMARY KEY(posid, name))
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
	#Check if position has a posid
	check = "SELECT posid, COUNT(posid) FROM positions WHERE locid=%s AND column= %s"%(loc,pos)
	res = p.read_sql_query(check, conn)
	posid = None
	#If not, submit into positions table
	if (res["COUNT(posid)"].values) == 0:
		stuff = [loc, pos] 
		sql = '''INSERT INTO positions(locid, column) VALUES(?,?)'''
		cur.execute(sql, stuff)
		posid = cur.lastrowid
	#If present, use that posid
	else:
		posid = int(res["posid"].values)
	#Using posid, submit to annotations
	stuff = [posid, name, val]
	sql = ''' INSERT INTO variants(posid, name, value) 
			VALUES(?,?,?) '''
	cur = conn.cursor()
	cur.execute(sql, stuff)
	










