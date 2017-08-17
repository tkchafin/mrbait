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
	cursor.execute('''DROP TABLE IF EXISTS regions''')
	cursor.execute('''DROP TABLE IF EXISTS samples''')
	
	#Table holding records for each locus
	cursor.execute('''
		CREATE TABLE loci(id INTEGER PRIMARY KEY, depth INTEGER NOT NULL, 
			length INTEGER NOT NULL, consensus TEXT NOT NULL, pass INTEGER NOT NULL,
			UNIQUE(id))
	''')
	
	#Table holding records for each locus
	cursor.execute('''
		CREATE TABLE regions(regid INTEGER PRIMARY KEY, locid INTEGER NOT NULL, 
			length INTEGER NOT NULL, sequence TEXT NOT NULL, vars INTEGER, 
			bad INTEGER, gap INTEGER, start INTEGER NOT NULL, stop INTEGER NOT NULL, 
			pass INTEGER NOT NULL, 
			FOREIGN KEY (locid) REFERENCES loci(id))
	''')
	
	#Table holding variant information
	cursor.execute('''
		CREATE TABLE variants(varid INTEGER NOT NULL, locid INTEGER NOT NULL, 
			column INTGER NOT NULL, value TEXT NOT NULL,
			FOREIGN KEY (locid) REFERENCES loci(id), 
			PRIMARY KEY(varid), 
			UNIQUE(varid))
	''')
	
	#Table holding records for each locus
	#cursor.execute('''
	#	CREATE TABLE samples(sampid INTEGER NOT NULL, name TEXT NOT NULL, 
	#		PRIMARY KEY (sampid),
	#		UNIQUE(name))
	#''')
	
	connection.commit()

#Code to add record to 'loci' table
def add_locus_record(conn, depth, consensus, passed=0):
	stuff = [depth, int(len(consensus)), str(consensus), int(passed)]
	sql = ''' INSERT INTO loci(depth, length, consensus, pass) 
				VALUES(?,?,?,?) '''
	cur = conn.cursor()
	cur.execute(sql, stuff)
	conn.commit()
	return cur.lastrowid

#Code to add to 'variants' table
def add_variant_record(conn, loc, pos, val):
	#Establish cursor 
	cur = conn.cursor()
	
	#Check if sample has a sampid, fetch if it exists
	#sql= "INSERT OR IGNORE INTO samples(name) VALUES (%r);"%name
	#cur.execute(sql)
	#fetch = "SELECT sampid FROM samples WHERE name = %r"%name
	#cur.execute(fetch)
	#sampid = cur.fetchone()[0] #fecth the sample id
	
	#Check if position has a posid, fetch if it exists
	sql2 = '''INSERT OR IGNORE INTO variants(locid, column, value) VALUES(?,?,?)'''
	stuff = [loc, pos, val]
	cur.execute(sql2,stuff)
	conn.commit()

#Function to add region to regions table
def add_region_record(conn, locid, start, stop, seq, counts):
	#Establish cursor 
	cur = conn.cursor()
	
	#build sql and pack values to insert
	sql = '''INSERT INTO regions(locid, length, sequence, vars, bad, gap,
		start, stop, pass) VALUES (?,?,?,?,?,?,?,?,0)'''
	stuff = [locid, len(seq), seq, counts["*"], counts["N"], counts["-"], start, stop] 
	
	#insert
	cur.execute(sql, stuff)
	conn.commit()

#Function to filter regions relation by minimum flanking SNPs 
def regionFilterMinVar(conn, val, flank):
	cur = conn.cursor()
	#Update pass to "1/fail" where COUNT(vars) in flanking region + TR is less than minimum
	sql = ''' 
		UPDATE regions 
		SET pass = 1 
		WHERE regid in 
			(SELECT regid FROM
			(SELECT 
				regid, 
				COUNT(DISTINCT column) as counts 
			FROM 
				regions INNER JOIN variants ON regions.locid = variants.locid 
			WHERE 
				value != "N" AND value != "-"
			AND 
				((column < (stop+%s)) AND (column > start-%s))
			GROUP BY regid)
			WHERE counts < %s);

	'''%(flank, flank, val)
	cur.execute(sql)
	conn.commit()

#Function to filter regions relation by maximum flanking SNPs 
def regionFilterMaxVar(conn, val, flank):
	cur = conn.cursor()
	#Update pass to "1/fail" where COUNT(vars) in flanking region + TR is greater than maximum
	sql = ''' 
		UPDATE regions 
		SET pass = 1 
		WHERE regid in 
			(SELECT regid FROM
			(SELECT 
				regid, 
				COUNT(DISTINCT column) as counts 
			FROM 
				regions INNER JOIN variants ON regions.locid = variants.locid 
			WHERE 
				value != "N" AND value != "-"
			AND 
				((column < (stop+%s)) AND (column > start-%s))
			GROUP BY regid)
			WHERE counts > %s);

	'''%(flank, flank, val)
	cur.execute(sql)
	conn.commit()

#Function for debug printing of variant counts flanking TRs
def printVarCounts(conn, flank):
	cur = conn.cursor()
	print("Variants within %r bases of TRs:"%flank)
	#FOR TESTING/DEBUGGING
	sql2 = ''' 
	SELECT 
		regid, 
		COUNT(DISTINCT column) as counts 
	FROM 
		regions INNER JOIN variants ON regions.locid = variants.locid 
	WHERE 
		value != "N" AND value != "-"
	AND 
		((column < (stop+%s)) AND (column > start-%s))
	GROUP BY regid
	'''%(flank, flank)
	print (pd.read_sql_query(sql2, conn))

#Function for random selection of TRs
def regionFilterRandom(conn, num):
	cur = conn.cursor()
	num = int(num) #number to keep
	
	#Fetch number of total
	cur.execute("SELECT COUNT(*) FROM REGIONS")
	rows = int(cur.fetchone()[0])
	
	#fetch number already failed
	cur.execute("SELECT COUNT(*) FROM regions WHERE pass=1")
	fails = int(cur.fetchone()[0])	
	
	print("Number of rows:",rows)
	if rows is 0 or rows is None: 
		raise ValueError("There are no rows in <regions>!")
	if num < rows-fails:
		sql = ''' 
			UPDATE regions 
			SET pass = 1 
			WHERE regid in 
				(SELECT 
					regid 
				FROM
					regions 
				WHERE 
					pass=0
				ORDER BY RANDOM() LIMIT(%s - %s - %s)
				)
		'''%(rows,fails,num)
		cur.execute(sql)
		conn.commit()

