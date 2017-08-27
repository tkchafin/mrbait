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
	
#Function to select one TR per alignment based on 
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

#Function to fetch all target regions requiring conflict resolution
def fetchConflictTRs(conn, min_len, dist):
	cur = conn.cursor()
	#Create temporary table to store conflict block information
	sql = ''' 
	CREATE TEMPORARY TABLE conflicts (
		regid INTEGER PRIMARY KEY,
		length INTEGER, 
		conflict_block INTEGER,
		choose INTEGER,
		FOREIGN KEY (regid) REFERENCES regions(regid)
	); 

	'''
	cur.execute(sql)
	conn.commit()
	
	#Populate table using passing TRs
	sql2 = '''
		INSERT INTO conflicts
		SELECT 
			regions.regid, 
			loci.length,
			"NULL" AS conflict_block,
			"NULL" AS choose
		FROM 
			regions INNER JOIN loci ON regions.locid = loci.id 
		WHERE 
			regions.pass = "0"
	'''
	cur.execute(sql2)
	conn.commit()
	
	#Query conflicts temp table to make a pandas dataframe for parsing
	sql_test = '''
	SELECT 
		conflicts.regid,
		conflict_block,
		conflicts.length,
		choose, 
		locid, 
		start, 
		stop
	FROM 
		conflicts INNER JOIN regions ON conflicts.regid = regions.regid
	'''
	df = pd.read_sql_query(sql_test, conn)
	block = max(df["locid"]) + 1
	#for name, row in df.iterrows():
		#If row isn't in a block:
		#if row["conflict_block"] == "NULL":
			#Compare to all rows in same locus 
			
				#If overlap, assign block to both
		#print(df["locid"] == row["locid"])
	groups = df.groupby("locid")
	for group, group_df in groups:
		print("Group is: ",group)
		#If only one TR for alignment, set choose to 1:
		if group_df.count == 1:
			for name, row in group_df.iterrows():
				df.loc[name, "conflict_block"] = row["locid"]
				df.loc[name, "choose"] = 1
		else:
			for name, row in group_df.iterrows():
				#If alignment length below minlen, set conflict_group to locid
				if row["length"] <= min_len:
					df.loc[name, "conflict_block"] = row["locid"]
				else:
					for _name, _row in group_df.iterrows():
						if name == _name:
							continue
						else:
							if _row["stop"] > (row["start"]- dist) or _row["start"] < (row["stop"]+ dist):
								if _row["conflict_block"] == "NULL":
									df.loc[_name, "conflict_block"] = block
									df.loc[name, "conflict_block"] = block
									block+=1
								else:
									df.loc[_name, "conflict_block"] = _row["conflict_block"]
	print(df)
	
def randomChooseRegionMINLEN(conn, min_len):
	#Assign blocks as full alignments below min_len
	sql_minlen = '''
	UPDATE conflicts
	SET choose = 1
	WHERE
		regid IN
			(SELECT 
				regid
			FROM
				(SELECT 
					regid,
					locid
				FROM
					conflicts INNER JOIN loci ON conflicts.locid = loci.id
				WHERE 
					length <= %r
				ORDER BY
					RANDOM() 
				)
			GROUP BY
				locid
			)
	'''%min_len
	cur.execute(sql_minlen)





