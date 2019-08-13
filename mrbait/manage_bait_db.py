#!/usr/bin/python

import sys
import sqlite3
from sqlite3 import OperationalError
import pandas as pd
from mrbait import misc_utils as utils
from mrbait import sequence_tools as s

#Function to create database connection
#Add code later to enable re-running from existing database
def create_connection(db):
	conn = sqlite3.connect(db)
	return conn


#Initialize empty databases
def init_new_db(connection):
	clearLoci(connection)
	#clearVariants(connection)
	clearGFF(connection)
	clearTargets(connection)
	clearBaits(connection)

#Function to clear baits table
def clearBaits(conn):
	cursor = conn.cursor()
	cursor.execute('''DROP TABLE IF EXISTS baits''')
	cursor.execute('''
		CREATE TABLE baits(baitid INTEGER NOT NULL, regid INTEGER NOT NULL,
			sequence TEXT NOT NULL, start INTEGER NOT NULL, stop INTEGER NOT NULL,
			mask REAL, gc REAL, pass INTEGER NOT NULL,
			PRIMARY KEY(baitid),
			FOREIGN KEY (regid) REFERENCES regions(regid)
		)
	''')
	conn.commit()

#Function to clear targets
def clearTargets(conn):
	cursor = conn.cursor()
	cursor.execute('''DROP TABLE IF EXISTS regions''')
	#Table holding records for each locus
	cursor.execute('''
		CREATE TABLE regions(regid INTEGER PRIMARY KEY, locid INTEGER NOT NULL,
			length INTEGER NOT NULL, sequence TEXT NOT NULL, vars INTEGER,
			bad INTEGER, gap INTEGER, mask REAL, gc REAL, vars_flank INTEGER,
			bad_flank INTEGER, gap_flank INTEGER, start INTEGER NOT NULL,
			stop INTEGER NOT NULL,pass INTEGER NOT NULL,
			FOREIGN KEY (locid) REFERENCES loci(id))
	''')
	conn.commit()

#function to initialize clean loci
def clearLoci(conn):
	cursor = conn.cursor()
	cursor.execute('''DROP TABLE IF EXISTS loci''')
	#Table holding records for each locus
	cursor.execute('''
		CREATE TABLE loci(id INTEGER PRIMARY KEY, depth INTEGER NOT NULL,
			length INTEGER NOT NULL, consensus TEXT NOT NULL, pass INTEGER NOT NULL,
			chrom TEXT, ambig REAL, gap REAL, mask REAL, gc REAL)
	''') #chrom only relevant if --assembly
	#UNIQUE(chrom) ON CONFLICT ROLLBACK
	conn.commit()

#Function to clear variants
"""DEPRECATED"""
# def clearVariants(conn):
# 	cursor = conn.cursor()
# 	cursor.execute('''DROP TABLE IF EXISTS variants''')
# 	#Table holding variant information
# 	cursor.execute('''
# 		CREATE TABLE variants(varid INTEGER NOT NULL, locid INTEGER NOT NULL,
# 			column INTEGER NOT NULL, value TEXT NOT NULL,
# 			FOREIGN KEY (locid) REFERENCES loci(id),
# 			PRIMARY KEY(varid),
# 			UNIQUE(varid))
# 	''')
# 	conn.commit()


#Function to clear GFF
def clearGFF(conn):
	cursor = conn.cursor()
	cursor.execute('''DROP TABLE IF EXISTS gff''')
	#Table holding variant information
	#Table holding GFF element information
	cursor.execute('''
		CREATE TABLE gff(gffid INTEGER PRIMARY KEY, locid INTEGER NOT NULL,
			type TEXT NOT NULL, start INTEGER NOT NULL,
			stop INTEGER NOT NULL, alias TEXT, pass INTEGER NOT NULL,
			FOREIGN KEY (locid) REFERENCES loci(id))
	''')
	conn.commit()



################################################################################

#function to reset all targets to passing
def resetTargets(conn):
	cur = conn.cursor()
	sql = '''
	UPDATE
		regions
	SET
		pass = 1
	'''
	cur.execute(sql)
	conn.commit()

#function to reset all baits to passing
def resetBaits(conn):
	cur = conn.cursor()
	sql = '''
	UPDATE
		baits
	SET
		pass = 1
	'''
	cur.execute(sql)
	conn.commit()

#Function to filter loci by coverage and length
def filterLoci(conn, minlen, mincov, max_ambig, max_mask):
	cur = conn.cursor()
	sql = '''
	UPDATE
		loci
	SET
		pass=0
	WHERE
		length < ? OR depth < ? OR (ambig + gap) > ? OR mask > ?
	'''
	stuff = [minlen, mincov, max_ambig, max_mask]
	cur.execute(sql, stuff)
	conn.commit()

#Function returns pandas dataframe of passedLoci
def getPassedLoci(conn):
	return(pd.read_sql_query("""SELECT id, consensus, chrom FROM loci WHERE pass=1""", conn))

#Function returns a Pandas DataFrame of passing target regions
def getPassedTRs(conn):
	return(pd.read_sql_query("""SELECT regid, sequence FROM regions WHERE pass=1""", conn))

#Function to parse fetchone() results (internal)
def parseFetchNum(fet):
	if fet is None:
		return 0
	else:
		return fet[0]

"""DEPRECATED"""
# #returns number of variants
# def getNumVars(conn):
# 	cur = conn.cursor()
# 	cur.execute("""SELECT count(varid) FROM variants""")
# 	return(parseFetchNum(cur.fetchone()))

#returns number of GFF elements
def getNumGFF(conn):
	cur = conn.cursor()
	cur.execute("""SELECT count(gffid) FROM gff""")
	return(parseFetchNum(cur.fetchone()))

#returns number of GFF elements
def getNumPassedGFF(conn):
	cur = conn.cursor()
	cur.execute("""SELECT count(gffid) FROM gff WHERE pass = 1""")
	return(parseFetchNum(cur.fetchone()))


#Function returns count of passed loci
def getNumPassedLoci(conn):
	cur = conn.cursor()
	cur.execute("""SELECT count(id) FROM loci WHERE pass=1""")
	return(parseFetchNum(cur.fetchone()))

#returns number of GFF elements
def getNumPassedBaits(conn):
	cur = conn.cursor()
	cur.execute("""SELECT count(baitid) FROM baits WHERE pass = 1""")
	return(parseFetchNum(cur.fetchone()))

#Function returns pandas dataframe of passed targets
def getNumConflicts(conn):
	cur = conn.cursor()
	try:
		cur.execute("""SELECT count(*) FROM conflicts WHERE choose != 1""")
		return(parseFetchNum(cur.fetchone()))
	except OperationalError:
		return(False)

#Function returns pandas dataframe of passedLoci
def getNumTRs(conn):
	cur = conn.cursor()
	check = '''
		SELECT
			count(*)
		FROM
		 	regions
	'''
	cur.execute(check)
	return(parseFetchNum(cur.fetchone()))

#Function returns pandas dataframe of passedLoci
def getNumPassedTRs(conn):
	cur = conn.cursor()
	check = '''
		SELECT
			count(*)
		FROM
		 	regions
		WHERE
			pass=1
	'''
	cur.execute(check)
	return(parseFetchNum(cur.fetchone()))

#Function returns pandas dataframe of baits to print
def getPrintBaits(conn):
	sql = '''
		SELECT
			locid, baits.regid, baitid, baits.sequence
		FROM
		 	baits INNER JOIN regions ON baits.regid = regions.regid
		WHERE
			baits.pass=1
	'''
	return(pd.read_sql_query(sql, conn))

#Function to return loci table
def getLoci(conn):
	return(pd.read_sql_query("""SELECT * FROM loci """, conn))

#Function to return regions table
def getRegions(conn):
	return(pd.read_sql_query("""SELECT * FROM regions """, conn))

#Function to return baits table
def getBaits(conn):
	return(pd.read_sql_query("""SELECT * FROM baits """, conn))

#Function to return GFF table as pandas DataFrame
def getGFF(conn):
	return(pd.read_sql_query("""SELECT * FROM gff """, conn))

#Function to return baits table
def getPassedBaits(conn):
	return(pd.read_sql_query("""SELECT baitid, sequence FROM baits WHERE pass=1""", conn))

"""DEPRECATED"""
# #Function to return variants table
# def getVariants(conn):
# 	return(pd.read_sql_query("""SELECT * FROM variants """, conn))

#Code to add record to 'loci' table
def add_locus_record(conn, depth, consensus, passed, name):
	if name == None:
		name = "NA"
	seq_norm = s.simplifySeq(consensus)
	counts = s.seqCounterSimple(seq_norm)
	ambig = counts["N"]/len(consensus)
	gap = counts["-"]/len(consensus)
	mask = s.mask_content(consensus)
	gc = s.gc_content(consensus)
	stuff = [depth, int(len(consensus)), str(consensus), int(passed), str(name), float(ambig), float(gap), float(mask), float(gc)]
	try:
		sql = ''' INSERT INTO loci(depth, length, consensus, pass, chrom, ambig, gap, mask, gc)
					VALUES(?,?,?,?,?,?,?,?,?) '''
		cur = conn.cursor()
		cur.execute(sql, stuff)
		conn.commit()
	except sqlite3.IntegrityError as err:
		print("Constraint failed: Skipping locus \"%s\" because it already exists, this is usually caused by duplicate headers when parsing a FASTA file."%name)
	return cur.lastrowid


#Code to add record to 'bait' table
def add_bait_record(conn, reg, seq, start, stop, mask, gc):
	mask_p = float(mask/len(seq))
	gc_p = float(gc/len(seq))
	stuff = [int(reg), seq, int(start), int(stop), float(mask_p), float(gc_p)]
	sql = ''' INSERT INTO baits(regid, sequence, start, stop, mask, gc, pass)
				VALUES(?,?,?,?,?,?,1) '''
	cur = conn.cursor()
	cur.execute(sql, stuff)
	conn.commit()
	return cur.lastrowid

#Code to add record to GFF table
def add_gff_record(conn,seqid, gff_type, start, stop, alias):
	cur = conn.cursor()

	#Query seqid for given chrom
	cur.execute("SELECT id FROM loci WHERE chrom = ?;",(seqid,))
	res = cur.fetchone()
	if res is not None:
		#If any match exists, fetch locid
		locid = res[0]
		#print("Locid for",seqid,"is:",locid)
		#BUILT INSERT SQL
		sql = ''' INSERT INTO gff(locid, type, start, stop, alias, pass)
					VALUES(?,?,?,?,?,1);'''
		stuff = [int(locid), str(gff_type), int(start), int(stop), str(alias)]
		cur.execute(sql, stuff)
		conn.commit()

"""DEPRECATED"""
# #Code to add to 'variants' table
# def add_variant_record(conn, loc, pos, val):
# 	#Establish cursor
# 	cur = conn.cursor()
#
# 	#Check if sample has a sampid, fetch if it exists
# 	#sql= "INSERT OR IGNORE INTO samples(name) VALUES (%r);"%name
# 	#cur.execute(sql)
# 	#fetch = "SELECT sampid FROM samples WHERE name = %r"%name
# 	#cur.execute(fetch)
# 	#sampid = cur.fetchone()[0] #fecth the sample id
#
# 	#Check if position has a posid, fetch if it exists
# 	sql2 = '''INSERT OR IGNORE INTO variants(locid, column, value) VALUES(?,?,?)'''
# 	stuff = [loc, pos, val]
# 	cur.execute(sql2,stuff)
# 	conn.commit()

#Function to add region to regions table
def add_region_record(conn, locid, start, stop, seq, counts, fcounts, mask, gc):
	#Establish cursor
	cur = conn.cursor()

	mask_p = float(mask/len(seq))
	gc_p = float(gc/len(seq))
	#build sql and pack values to insert
	sql = '''INSERT INTO regions(locid, length, sequence, vars, bad, gap, mask, gc,
		vars_flank, bad_flank, gap_flank, start, stop, pass) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,1)'''
	stuff = [locid, len(seq), seq, counts["*"], counts["N"], counts["-"], mask_p, gc_p, fcounts["*"], fcounts["N"], fcounts["-"],start, stop]

	#insert
	cur.execute(sql, stuff)
	conn.commit()


#Function to FAIL any GFF elements that do not overlap with our sequence for the given locus/region
def validateGFFRecords(conn):
	cur = conn.cursor()

	sql = '''
		UPDATE
			gff
		SET
			pass = 0
		WHERE gffid IN
			(SELECT gffid FROM
				gff INNER JOIN loci ON gff.locid = loci.id
			WHERE
				(gff.start > loci.length) AND (gff.stop > loci.length)
			);
	'''
	cur.execute(sql)
	conn.commit()

"""DEPRECATED"""
# #Function to filter regions relation by minimum flanking SNPs
# def regionFilterMinVar(conn, val, flank):
# 	cur = conn.cursor()
# 	#Update pass to "0/FALSE" where COUNT(vars) in flanking region + TR is less than minimum
# 	sql = '''
# 		UPDATE regions
# 		SET pass = 0
# 		WHERE regid in
# 			(SELECT regid FROM
# 			(SELECT
# 				regid,
# 				COUNT(DISTINCT column) as counts
# 			FROM
# 				regions INNER JOIN variants ON regions.locid = variants.locid
# 			WHERE
# 				value != "N" AND value != "-"
# 			AND
# 				((column < (stop+?)) AND (column > start-?))
# 			GROUP BY regid)
# 			WHERE counts < ?);
#
# 	'''
# 	stuff = [flank, flank, val]
# 	cur.execute(sql, stuff)
# 	conn.commit()

"""DEPRECATED"""
# #Function to filter regions relation by maximum flanking SNPs
# def regionFilterMaxVar(conn, val, flank):
# 	cur = conn.cursor()
# 	#Update pass to "0/FALSE" where COUNT(vars) in flanking region + TR is greater than maximum
# 	sql = '''
# 		UPDATE regions
# 		SET pass = 0
# 		WHERE regid in
# 			(SELECT regid FROM
# 			(SELECT
# 				regid,
# 				COUNT(DISTINCT column) as counts
# 			FROM
# 				regions INNER JOIN variants ON regions.locid = variants.locid
# 			WHERE
# 				value != "N" AND value != "-"
# 			AND
# 				((column < (stop+?)) AND (column > start-?))
# 			GROUP BY regid)
# 			WHERE counts > ?);
#
# 	'''
# 	stuff = [flank, flank, val]
# 	cur.execute(sql, stuff)
# 	conn.commit()

"""DEPRECATED"""
# #Function for debug printing of variant counts flanking TRs
# def printVarCounts(conn, flank):
# 	cur = conn.cursor()
# 	print("Variants within %r bases of TRs:"%flank)
# 	#FOR TESTING/DEBUGGING
# 	sql2 = '''
# 	SELECT
# 		regid,
# 		COUNT(DISTINCT column) as counts
# 	FROM
# 		regions INNER JOIN variants ON regions.locid = variants.locid
# 	WHERE
# 		value != "N" AND value != "-"
# 	AND
# 		((column <= (stop+%s)) AND (column >= start-%s))
# 	GROUP BY regid
# 	'''%(flank, flank)
# 	print (pd.read_sql_query(sql2, conn))

#Function for random selection of TRs
def regionFilterRandom(conn, num):
	cur = conn.cursor()
	num = int(num) #number to keep

	#Fetch number of total
	cur.execute("SELECT COUNT(*) FROM REGIONS")
	rows = int(cur.fetchone()[0])

	#fetch number already failed
	cur.execute("SELECT COUNT(*) FROM regions WHERE pass=0")
	fails = int(cur.fetchone()[0])

	#print("Number of rows:",rows)
	if rows is 0 or rows is None:
		raise ValueError("There are no rows in <regions>!")
	if num < rows-fails:
		sql = '''
			UPDATE regions
			SET pass = 0
			WHERE regid in
				(SELECT
					regid
				FROM
					regions
				WHERE
					pass=1
				ORDER BY RANDOM() LIMIT(? - ? - ?)
				)
		'''
		stuff = [rows, fails, num]
		cur.execute(sql, stuff)
		conn.commit()

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
						length <= ?
					ORDER BY
						RANDOM()
					)
				GROUP BY
					locid
				)
	'''
	cur.execute(sql_minlen, (minlen,))

"""DEPRECATED"""
# #Function to delete from variable table on locid
# def purgeVars(conn, key):
# 	cur = conn.cursor()
#
# 	sql = '''
# 	DELETE FROM variants
# 	WHERE locid = ?
# 	'''
#
# 	cur.execute(sql,(key,))
# 	conn.commit()

#Function to update consensus sequence of a locus
def updateConsensus(conn, key, seq):
	cur = conn.cursor()
	sql = '''
		UPDATE loci
		SET
			consensus = ?
		WHERE
			id = ?
	'''
	stuff = [str(seq), int(key)]
	cur.execute(sql, stuff)
	conn.commit()

#Internal function for checking if TRs overlap within distance buffer
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


#Function to initialize TEMPORARY TABLE 'conflicts'
def initializeConflicts(conn):
	cur = conn.cursor()

	#Clear conflicts table if it exists already
	cur.execute("DROP TABLE IF EXISTS conflicts")

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

	check = '''
		SELECT
			count(*)
		FROM
		 	sqlite_temp_master
		WHERE
			type='table' AND name='conflicts'
	'''
	res_check = cur.execute(check)

	#Raise exception if conflicts failed to initialize
	if (cur.fetchone()[0]) is 0:
		raise RuntimeError("SQL Error: Failed to CREATE TEMPORARY TABLE \"conflicts\"")

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
			regions.pass=1
	'''
	cur.execute(sql2)
	conn.commit()

#Function to take a pointer to a pandas dataframe, send to SQL, and update conflicts table
def updateConflictsFromPandas(conn, df):
	cur = conn.cursor()

	#Fetch number of entries in conflict tables
	rows = getConflictNumRows(conn)

	#Make sure there is some data to work on
	if rows is 0 or rows is None:
		raise ValueError("There are no rows in <conflicts>!")

	df.to_sql('t', conn, if_exists='replace')

	#Hacky way to do it, but SQlite doesn't support FROM clause in UPDATEs...
	sql_update = '''
		UPDATE
			conflicts
		SET
			conflict_block = (SELECT t.conflict_block FROM t WHERE t.regid = conflicts.regid)
		WHERE
			EXISTS(SELECT * FROM t WHERE t.regid = conflicts.regid)
	'''
	cur.execute(sql_update)

	#Clear up the temp table t
	cur.execute("DROP TABLE IF EXISTS t")
	conn.commit()

#Function to take a pointer to a pandas dataframe, send to SQL, and update conflicts table
def updateChosenFromPandas(conn, df):
	cur = conn.cursor()

	#Fetch number of entries in conflict tables
	rows = getConflictNumRows(conn)

	#Make sure there is some data to work on
	if rows is 0 or rows is None:
		raise ValueError("There are no rows in <conflicts>!")

	df.to_sql('t', conn, if_exists='replace')

	#Hacky way to do it, but SQlite doesn't support FROM clause in UPDATEs...
	sql_update = '''
		UPDATE
			conflicts
		SET
			choose = (SELECT t.choose FROM t WHERE t.regid = conflicts.regid)
		WHERE
			EXISTS(SELECT * FROM t WHERE t.regid = conflicts.regid)
	'''
	cur.execute(sql_update)

	#Clear up the temp table t
	cur.execute("DROP TABLE IF EXISTS t")
	conn.commit()

#Updates
def updateLociMask(conn, newMask):
	cur = conn.cursor()

	#convert to data frame
	df = pd.DataFrame(newMask, columns=("id", "mask"))

	#create temporary table
	df.to_sql('m', conn, if_exists='replace')
	#print(pd.read_sql_query("SELECT * FROM loci", conn))

	#Hacky way to do it, but SQlite doesn't support FROM clause in UPDATEs...
	sql_update = '''
		UPDATE
			loci
		SET
			mask = (SELECT m.mask FROM m WHERE m.id = loci.id)
		WHERE
			EXISTS(SELECT * FROM m WHERE m.id = loci.id)
	'''
	cur.execute(sql_update)

	#Clear up the temp table t
	cur.execute("DROP TABLE IF EXISTS m")
	conn.commit()
	#print(pd.read_sql_query("SELECT * FROM loci", conn))

#Function to build conflicts table when --R is false
def fetchConflictTRs_NoMult(conn):
	cur = conn.cursor()

	#Function call to build conflicts table
	try:
		num = initializeConflicts(conn)
	except RuntimeError as err:
		print(err.args)

	#Query to UPDATE conflicts with each locus as a conflict_block
	update = '''
		UPDATE
			conflicts
		SET
			conflict_block = (SELECT r.locid FROM regions AS r WHERE r.regid = conflicts.regid)
		WHERE
			EXISTS(SELECT * FROM regions AS r WHERE r.regid = conflicts.regid)
	'''
	cur.execute(update)
	conn.commit()

	#DEBUG print
	#print(pd.read_sql_query("SELECT * FROM conflicts", conn))

#Function to fetch number of rows in conflicts table
def getConflictNumRows(conn):
	cur = conn.cursor()
	cur.execute("SELECT COUNT(*) FROM conflicts")
	rows = int(cur.fetchone()[0])
	return(rows)

#Function to fetch all target regions requiring conflict resolution
def fetchConflictTRs(conn, min_len, dist):
	cur = conn.cursor()
	#TODO: Check and first make sure conflicts table exists, or otherwise implement error handling

	#Function call to build conflicts table
	try:
		num = initializeConflicts(conn)
	except RuntimeError as err:
		print(err.args)

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
		WHERE
			regions.pass=1
	'''

	df = pd.read_sql_query(sql_test, conn)
	#Set block num start at 1+ last locid (to prevent overlap in block nums)
	block = int(max(df["locid"]) + 1) #make sure to enforce integer type, or it fucks up in sql later
	#Split df into locid groups (retains INDEX of each entry, but in separate dfs)
	groups = df.groupby("locid")
	for group, _df in groups:
		group_df = _df.sort_values(by="start", ascending=True, kind='mergesort')
		#If only one TR for alignment, set choose to 1:
		if group_df.shape[0] == 1:
			for name, row in group_df.iterrows():
				#modify original dataframe
				df.loc[name, "conflict_block"] = row["locid"]
				df.loc[name, "choose"] = 1
		else:
			#print("\n########Group is: ",group,"\n")
			#print(group_df)
			#For each TR in locus
			#Create iterable
			row_iterator = group_df.iterrows()
			_i, last_row = row_iterator.__next__() #grab first item
			#print("First row:",_i,last_row)
			#Make sure it is from a long enough locus
			if last_row["length"] <= min_len:
				df.loc[_i, "conflict_block"] = last_row["locid"]
			#For each next row (starting at row 2)
			for i, next_row in row_iterator:
				#print("Comparing row ",_i," and ", i)
				if next_row["length"] <= min_len:
					df.loc[i, "conflict_block"] = next_row["locid"]
				#If they overlap, assign them both to same conflict block
				if checkOverlap(last_row, next_row, dist) == 1:
					#print("overlap!")
					#Fetch conflict_block currents from parent df
					cb = df.loc[i, "conflict_block"]
					_cb = df.loc[_i, "conflict_block"]
					#If neither is assigned, assign both to new block
					if (cb == "NULL") and (_cb == "NULL"):
						#print("Assigning block: ",block)
						df.loc[_i, "conflict_block"] = block
						df.loc[i, "conflict_block"] = block
						block+=1
					#Else if row1 is assigned, give row2 the block of row1
					#If last_row is assigned, move next_row to same block
					elif cb == "NULL":
						df.loc[i, "conflict_block"] = df.loc[_i, "conflict_block"]
					#Or if row2 has block, assign it to row1
					elif _cb == "NULL":
						df.loc[_i, "conflict_block"] = df.loc[i, "conflict_block"]
				#If they don't overlap, assign last_row to its own conflict block
				#current next_row will then get checked with its next neighbor
				#in the next iteration.
				else:
					df.loc[_i, "conflict_block"] = block
					block += 1

	#Set choose=1 for any which have no conflicts
	df.loc[~df.duplicated("conflict_block",keep=False),"choose"] = 1
	#print(df)
	#Next step, send df back to SQLite and update conflicts table
	updateConflictsFromPandas(conn, df)
	updateChosenFromPandas(conn, df)
	#DEBUG print
	#print(pd.read_sql_query("SELECT * FROM conflicts", conn))

#Function for random selection of TRs within conflict blocks
def regionSelectRandom(conn):
	cur = conn.cursor()

	#Fetch number of entries in conflict tables
	rows = getConflictNumRows(conn)

	#Make sure there is some data to work on
	if rows is 0 or rows is None:
		raise ValueError("There are no rows in <conflicts>!")

	sql = '''
	UPDATE
		conflicts
	SET
		choose=1
	WHERE
		regid IN
			(SELECT
				regid
			FROM
				(SELECT
					*
				FROM
					conflicts
				WHERE
					choose="NULL"
				ORDER BY
					RANDOM()
				)
			GROUP BY
				conflict_block
			)
	'''
	#print(pd.read_sql_query(sql, conn))
	cur.execute(sql)
	conn.commit()

	#Set "unchosen" regions to 0/FALSE
	cur.execute("UPDATE conflicts SET choose=0 WHERE choose='NULL'")
	conn.commit()


#Function for resolving conflict blocks by number of flanking SNPs
def regionSelect_SNP(conn):
	cur = conn.cursor()

	#Fetch number of entries in conflict tables
	rows = getConflictNumRows(conn)

	#Make sure there is some data to work on
	if rows <= 0 :
		raise ValueError("There are no rows in <conflicts>!")

	# Template:
	# -Fetch joined conflict table (get flanking SNPs)
	# -Split pandas df into groups
	# -Loop through each group, tracking top region(s)
	# -If one region wins, set choose=1 and choose=0 for the rest
	#-Set chosen to 1
	#-If two are chosen, set others to 0 and keep the "tied" as NULL
	#-Any region in a block that had NO vars will be kept as NULL.
	#-----AFTER pushing pandas df back and updating, need to update
	#-----all unchosen regions to 0 ONLY IF the block alreday has a selected region
	#-----final call to RANDOM select by parent function will sort out the rest

	#-Use RANDOM to resolve any remaining NULLs.
	#print("DIST IS ",dist)
	#Query conflicts temp table to make a pandas dataframe for parsing
	sql = '''
		SELECT
			c.regid,
			conflict_block,
			choose,
			(vars+vars_flank) as counts
		FROM
			conflicts AS c INNER JOIN regions AS r USING (regid)
		WHERE
			c.choose="NULL"
	'''

	#printFlankingSNPs_conflicts(conn, dist)
	#printFlankingSNPCounts_conflicts(conn, dist)
	df = pd.read_sql_query(sql, conn)

	#Split df into locid groups (retains INDEX of each entry, but in separate dfs)
	#Loop through groups and select highest
	try:
		df = parseCountsMax(df)
	except:
		raise

	#Push modified DF to SQL as temp table
	updateChosenFromPandas(conn, df)


#Function to prints flanking SNPs for conflicting regions...
#Function for use when debugging
# def printFlankingSNPs_conflicts(conn):
# 	cur = conn.cursor()
#
# 	#Query conflicts temp table to make a pandas dataframe for parsing
# 	sql = '''
# 		SELECT
# 			c.regid,
# 			conflict_block,
# 			choose,
# 			vars
# 		FROM
# 			conflicts AS c INNER JOIN regions AS r USING (regid)
# 		WHERE
# 			c.choose="NULL"
# 	'''
# 	print(pd.read_sql_query(sql, conn))


#Function to loop though pandas DF of conflict counts and choose best by "minimum"
def parseCountsMin(df):
	#Loop through groups and select lowest
	groups = df.groupby("conflict_block")
	for group, group_df in groups:
		#print("\n########Conflict Block is: ",group,"\n")
		#print(group_df)
		#If only one TR for alignment, set choose to 1:
		if group_df.shape[0] == 1:
			for name, row in group_df.iterrows():
				#modify original dataframe
				df.loc[name, "choose"] = 1
		else:
			#For each TR in locus
			chosen_ones = []
			best = 0
			track = 0
			#Loop through and figure out which is best
			for name, row in group_df.iterrows():
				if track == 0:
					best = row["counts"]
					chosen_ones = [name]
					track = 1
				else:
					if row["counts"] > best:
						continue
					elif row["counts"] < best:
						#If this is the best we've seen, replace chosen_ones with new choice
						best = row["counts"]
						chosen_ones = [name]
					elif row["counts"] == best:
						#If its a tie, add to list of chosen_ones
						chosen_ones.append(name)

			#Use chosen list to assign values
			#If one chosen, set chosen to 1/True and the rest to 0/False
			if len(chosen_ones) == 1:
				df.loc[df["conflict_block"] == group, "choose"] = 0
				df.loc[chosen_ones[0], "choose"] = 1
			#Elif, multiple were chosen
			elif len(chosen_ones) > 1:
				df.loc[df["conflict_block"] == group, "choose"] = 0
				for i in chosen_ones:
					df.loc[i, "choose"] = "NULL"
			# else:
			# 	message = "Error: No TR chosen for block " + group + ", something went wrong"
			# 	raise RuntimeError(message)
	return(df)

#Function to go through pandas df of conflict counts and choose best by "maximum"
def parseCountsMax(df):
	groups = df.groupby("conflict_block")
	for group, group_df in groups:
		#print("\n########Conflict Block is: ",group,"\n")
		#print(group_df)
		#If only one TR for alignment, set choose to 1:
		if group_df.shape[0] == 1:
			for name, row in group_df.iterrows():
				#modify original dataframe
				df.loc[name, "choose"] = 1
		else:
			#For each TR in locus
			chosen_ones = []
			best = 0
			#Loop through and figure out which is best
			for name, row in group_df.iterrows():
				if row["counts"] < best:
					continue
				elif row["counts"] > best:
					#If this is the best we've seen, replace chosen_ones with new choice
					best = row["counts"]
					chosen_ones = [name]
				elif row["counts"] == best and best != 0:
					#If its a tie, add to list of chosen_ones
					chosen_ones.append(name)
			#Use chosen list to assign values
			#If one chosen, set chosen to 1/True and the rest to 0/False
			if len(chosen_ones) == 1:
				df.loc[df["conflict_block"] == group, "choose"] = 0
				df.loc[chosen_ones[0], "choose"] = 1
			#Elif, multiple were chosen
			elif len(chosen_ones) > 1:
				df.loc[df["conflict_block"] == group, "choose"] = 0
				for i in chosen_ones:
					df.loc[i, "choose"] = "NULL"
			# else:
			# 	message = "Error: No TR chosen for block %s, something went wrong"%(group)
			# 	raise RuntimeError(message)
	return(df)

#Function to prints flanking SNPs for conflicting regions...
#Function for use when debugging
def printFlankingSNPCounts_conflicts(conn):
	cur = conn.cursor()

	#Query conflicts temp table to make a pandas dataframe for parsing
	sql = '''
		SELECT
			c.regid,
			conflict_block,
			choose,
			(vars+vars_flank) AS counts
		FROM
			conflicts AS c INNER JOIN regions AS r USING (regid)
		WHERE
			c.choose="NULL"

	'''
	print(pd.read_sql_query(sql, conn))

#Function for resolving conflict blocks by minimizing "bad" bases in flanking region
def regionSelect_MINBAD(conn):
	cur = conn.cursor()

	#Fetch number of entries in conflict tables
	rows = getConflictNumRows(conn)

	#Make sure there is some data to work on
	if rows <= 0 :
		raise ValueError("There are no rows in <conflicts>!")

	#-Use RANDOM to resolve any remaining NULLs.
	#print("DIST IS ",dist)
	#Query conflicts temp table to make a pandas dataframe for parsing
	#This version includes records which have a zero from the COUNT aggregate function by
	#outer joining the 'counts' table by the original table.
	#Probably a better way to go about it but i'm tired of dicking with it...
	sql = '''
		SELECT
			c.regid,
			conflict_block,
			choose,
			(gap + bad + bad_flank + gap_flank) as counts
		FROM
			conflicts AS c INNER JOIN regions AS r USING (regid)
		WHERE
			c.choose="NULL"
	'''

	df = pd.read_sql_query(sql, conn)

	#Split df into locid groups (retains INDEX of each entry, but in separate dfs)
	#Loop through groups and select highest
	try:
		df = parseCountsMin(df)
	except:
		raise

	#Push modified DF to SQL as temp table
	updateChosenFromPandas(conn, df)

#Function for resolving conflict blocks by minimizing all variable bases in flanking region
def regionSelect_MINVAR_TR(conn):
	cur = conn.cursor()

	#Fetch number of entries in conflict tables
	rows = getConflictNumRows(conn)
	#Make sure there is some data to work on
	if rows <= 0 :
		raise ValueError("There are no rows in <conflicts>!")

	#Query conflicts temp table to make a pandas dataframe for parsing
	sql = '''
		SELECT
			c.regid,
			conflict_block,
			choose,
			(vars + flank_vars) AS counts
		FROM
			conflicts AS c INNER JOIN regions AS r USING (regid)
		WHERE
			c.choose="NULL"
	'''

	df = pd.read_sql_query(sql, conn)
	#Split df into locid groups (retains INDEX of each entry, but in separate dfs)
	try:
		df = parseCountsMin(df)
	except:
		raise

	#Push modified DF to SQL as temp table
	updateChosenFromPandas(conn, df)

#Function for resolving conflict blocks by minimizing number of flanking SNPs
def regionSelect_MINSNP(conn):
	cur = conn.cursor()

	#Fetch number of entries in conflict tables
	rows = getConflictNumRows(conn)

	#Make sure there is some data to work on
	if rows <= 0 :
		raise ValueError("There are no rows in <conflicts>!")

	#-Use RANDOM to resolve any remaining NULLs.
	#print("DIST IS ",dist)
	#Query conflicts temp table to make a pandas dataframe for parsing
	sql = '''
		SELECT
			c.regid,
			conflict_block,
			choose,
			(vars + vars_flank) as counts
		FROM
			conflicts AS c INNER JOIN regions AS r USING (regid)
		WHERE
			c.choose="NULL"
	'''

	#printFlankingSNPs_conflicts(conn, dist)
	#printFlankingSNPCounts_conflicts(conn, dist)
	df = pd.read_sql_query(sql, conn)

	#Split df into locid groups (retains INDEX of each entry, but in separate dfs)
	#Loop through groups and select highest
	try:
		df = parseCountsMin(df)
	except:
		raise

	#Push modified DF to SQL as temp table
	updateChosenFromPandas(conn, df)

#Function to push resolved TR conflicts to the regions table
def pushResolvedConflicts(conn):

	#Check that all conflicts are resolved
	cur = conn.cursor()
	unres = pd.read_sql_query("SELECT COUNT(*) FROM conflicts WHERE choose='NULL'", conn)
	rows = unres.shape[0]
	if rows <= 0:
		print("Unresolved conflicts:")
		print(unres)
		raise RuntimeError("Error: There are still unresolved conflicts")
	else:
		#Hacky way to do it, but SQlite doesn't support FROM clause in UPDATEs...
		sql_update = '''
			UPDATE
				regions
			SET
				pass = (SELECT c.choose FROM conflicts c WHERE c.regid = regions.regid)
			WHERE
				regions.regid in (SELECT c.regid FROM conflicts c WHERE c.regid = regions.regid)
			AND
				pass = 1
		'''
		cur.execute(sql_update)

		#Clear up the temp table conflicts
		cur.execute("DROP TABLE IF EXISTS conflicts")
		conn.commit()


#Functon to filter targets by length
def lengthFilterTR(conn, maxlen, minlen):
	cur = conn.cursor()

	sql = '''
	UPDATE
		regions
	SET
		pass=0
	WHERE
		length > ? OR length < ?
	'''
	stuff = [int(maxlen), int(minlen)]
	#print(pd.read_sql_query(sql, conn))
	cur.execute(sql, stuff)
	conn.commit()

#Functon to filter targets by --vmax_r
def varMaxFilterTR(conn, varmax):
	cur = conn.cursor()

	sql = '''
	UPDATE
		regions
	SET
		pass=0
	WHERE
		(vars+vars_flank) AS counts > ?
	'''
	#print(pd.read_sql_query(sql, conn))
	cur.execute(sql, (varmax,))
	conn.commit()

def regionFilterMask(conn, maxprop):
	cur = conn.cursor()
	sql = '''
	UPDATE
		regions
	SET
		pass=0
	WHERE
		(mask > ?)
	'''
	#print(pd.read_sql_query(sql, conn))
	cur.execute(sql, (maxprop,))
	conn.commit()

def regionFilterGC(conn, minprop, maxprop):
	cur = conn.cursor()
	sql = '''
	UPDATE
		regions
	SET
		pass=0
	WHERE
		(gc > ?) OR (gc < ?)
	'''
	stuff = [maxprop, minprop]
	#print(pd.read_sql_query(sql, conn))
	cur.execute(sql, stuff)
	conn.commit()

#Function to remove targets NOT included in list
def removeRegionsByWhitelist(conn, whitelist):

	cur = conn.cursor()
	#If nothing in list, no need to do any work:
	if len(whitelist) <= 0:
		return(0)
	df = pd.DataFrame({"regid" : whitelist})
	df.to_sql('tt', conn, if_exists='replace')

	#Hacky way to do it, but SQlite doesn't support FROM clause in UPDATEs...
	sql_update = '''
		UPDATE
			regions
		SET
			pass = 0
		WHERE
			NOT EXISTS(SELECT * FROM tt WHERE tt.regid = regions.regid)
	'''
	cur.execute(sql_update)
	#Clear up the temp table t
	cur.execute("DROP TABLE IF EXISTS tt")
	conn.commit()

#Function to remove Target Regions given a list of blacklisted regids
def removeRegionsByList(conn, blacklist):

	cur = conn.cursor()
	#If nothing in list, no need to do any work:
	if len(blacklist) <= 0:
		return(0)
	df = pd.DataFrame({"regid" : blacklist})
	df.to_sql('tt', conn, if_exists='replace')

	#Hacky way to do it, but SQlite doesn't support FROM clause in UPDATEs...
	sql_update = '''
		UPDATE
			regions
		SET
			pass = 0
		WHERE
			EXISTS(SELECT * FROM tt WHERE tt.regid = regions.regid)
	'''
	cur.execute(sql_update)
	#Clear up the temp table t
	cur.execute("DROP TABLE IF EXISTS tt")
	conn.commit()

#Function to remove baits given a list of blacklisted regids
def removeBaitsByList(conn, blacklist):

	cur = conn.cursor()

	#If nothing in list, no need to do any work:
	if len(blacklist) <= 0:
		return(0)
	df = pd.DataFrame({"baitid" : blacklist})
	#print(df)

	df.to_sql('b', conn, if_exists='replace')

	#Hacky way to do it, but SQlite doesn't support FROM clause in UPDATEs...
	sql_update = '''
		UPDATE
			baits
		SET
			pass = 0
		WHERE
			EXISTS(SELECT * FROM b WHERE b.baitid = baits.baitid)
	'''
	cur.execute(sql_update)

	#Clear up the temp table t
	cur.execute("DROP TABLE IF EXISTS b")
	conn.commit()

#Function to remove baits NOT included in list
def removeBaitsByWhitelist(conn, whitelist):

	cur = conn.cursor()
	#If nothing in list, no need to do any work:
	if len(whitelist) <= 0:
		return(0)
	df = pd.DataFrame({"baitid" : whitelist})
	df.to_sql('ttt', conn, if_exists='replace')

	#Hacky way to do it, but SQlite doesn't support FROM clause in UPDATEs...
	sql_update = '''
		UPDATE
			baits
		SET
			pass = 0
		WHERE
			NOT EXISTS(SELECT * FROM ttt WHERE ttt.baitid = baits.baitid)
	'''
	cur.execute(sql_update)
	#Clear up the temp table t
	cur.execute("DROP TABLE IF EXISTS ttt")
	conn.commit()


def baitFilterMask(conn, maxprop):
	cur = conn.cursor()
	sql = '''
	UPDATE
		baits
	SET
		pass=0
	WHERE
		(mask > ?)
	'''
	#print(pd.read_sql_query(sql, conn))
	cur.execute(sql, (maxprop,))
	conn.commit()

def baitFilterGC(conn, minprop, maxprop):
	cur = conn.cursor()
	sql = '''
	UPDATE
		baits
	SET
		pass=0
	WHERE
		(gc > ?) OR (gc < ?)
	'''
	stuff = [maxprop, minprop]

	cur.execute(sql, stuff)
	conn.commit()

#Function for random selection of baits
def baitFilterRandom(conn, num):
	cur = conn.cursor()
	num = int(num) #number to keep

	#Fetch number of total
	cur.execute("SELECT COUNT(*) FROM baits")
	rows = int(cur.fetchone()[0])

	#fetch number already failed
	cur.execute("SELECT COUNT(*) FROM baits WHERE pass=0")
	fails = int(cur.fetchone()[0])

	if rows is 0 or rows is None:
		raise ValueError("There are no rows in <baits>!")
	if num < rows-fails:
		sql = '''
			UPDATE baits
			SET pass = 0
			WHERE baitid in
				(SELECT
					baitid
				FROM
					baits
				WHERE
					pass=1
				ORDER BY RANDOM() LIMIT(? - ? - ?)
				)
		'''
		stuff = [rows, fails, num]
		cur.execute(sql, stuff)
		conn.commit()

#Function to filter targets by proximity or overlap with GFF records
def regionFilterGFF(conn, gff_type, dist):
	cur = conn.cursor()

	if getNumGFF(conn) > 0:
		if getNumPassedGFF(conn) > 0:
			df = pd.DataFrame() #empty pandas DF
			#If get GFF by type:
			if gff_type == "all":
				sql = """
					SELECT regid, regions.start, regions.stop, regions.pass,
						gffid, gff.start AS gff_start, gff.stop AS gff_stop
					FROM
						regions INNER JOIN gff ON regions.locid = gff.locid
					WHERE
						regions.pass = 1 AND gff.pass = 1
				"""
				df = pd.read_sql_query(sql, conn)
			else:
				#Query database to get targets with passing GFFs
				sql = """
					SELECT regid, regions.start, regions.stop, regions.pass,
						gffid, gff.start AS gff_start, gff.stop AS gff_stop
					FROM
						regions INNER JOIN gff ON regions.locid = gff.locid
					WHERE
						regions.pass = 1 AND gff.pass = 1 AND gff.type = ?
				"""
				df = pd.read_sql_query(sql, conn, params=(gff_type,))

			whitelist = parseJoinGFFTable(df, dist)
			removeRegionsByWhitelist(conn, whitelist)

			"""
			1. Fetch INNER JOINED db matching criterion
			2. Parse overlaps in Python (see function above)
			3. Return list of ones to keep.
			4. For UPDATE- Set pass to 0 if: pass NOT 1 in returned list, OR if already failed
				This should also fail the case where NO GFF RECORDS TO JOIN or gff record was failed
			"""
		else:
			cur.execute("UPDATE regions SET pass = 0")
			print("WARNING: No GFF records passed quality control. Because you chose to filter target regions on proximity to GFF records, no targets will be retained.")
	else:
		print("WARNING: No GFF records present in database. Skipping target region filtering on proximity to GFF records.")
	conn.commit()

#Function to filter targets by proximity or overlap with GFF records
def regionFilterGFF_Alias(conn, gff_type, dist):
	cur = conn.cursor()

	if getNumGFF(conn) > 0:
		if getNumPassedGFF(conn) > 0:
			df = pd.DataFrame() #empty pandas DF
			#If get GFF by type:
			if gff_type == "all":
				sql = """
					SELECT regid, regions.start, regions.stop, regions.pass,
						gffid, gff.start AS gff_start, gff.stop AS gff_stop
					FROM
						regions INNER JOIN gff ON regions.locid = gff.locid
					WHERE
						regions.pass = 1 AND gff.pass = 1 AND gff.alias != "NULL"
				"""
				df = pd.read_sql_query(sql, conn)
			else:
				#Query database to get targets with passing GFFs
				sql = """
					SELECT regid, regions.start, regions.stop, regions.pass,
						gffid, gff.start AS gff_start, gff.stop AS gff_stop
					FROM
						regions INNER JOIN gff ON regions.locid = gff.locid
					WHERE
						regions.pass = 1 AND gff.pass = 1 AND gff.alias = ?
				"""
				df = pd.read_sql_query(sql, conn, params=(gff_type,))

			#Get list of passed targets, pass list to FAIL all non-whitelisted targets
			whitelist = parseJoinGFFTable(df, dist)
			removeRegionsByWhitelist(conn, whitelist)

		else:
			cur.execute("UPDATE regions SET pass = 0")
			print("WARNING: No GFF records passed quality control. Because you chose to filter target regions on proximity to GFF records, no targets will be retained.")
	else:
		print("WARNING: No GFF records present in database. Skipping target region filtering on proximity to GFF records.")
	conn.commit()

#Function to parse an INNER JOIN regions/gff table for overlapping fragments, and return whitelist
def parseJoinGFFTable(df, dist):
	for index, row in df.iterrows():
		min1 = row["start"] - dist or 0
		min2 = row["gff_start"]
		max1 = row["stop"] + dist
		max2 = row["gff_stop"]
		#If no overlap (overlap distance <= 0), set 'pass' to 0/FAIL
		if utils.calcOverlap(min1, max1, min2, max2) <= 0:
			df.loc[index, "pass"] = 0

	#Get list of passed targets, pass list to FAIL all non-whitelisted targets
	whitelist = df[df["pass"]==1].regid.unique()
	return(whitelist)

"""DEPRECATED"""
# #Function to parse variants table to update regions VARS for flanking information
# def parseFlankVars(conn, dist):
# 	cur = conn.cursor()
# 	sql = '''
# 		WITH other AS
# 			(SELECT
# 				r.regid,
# 				COUNT(v.value) as counts
# 			FROM
# 				regions AS r INNER JOIN variants AS v ON r.locid = v.locid
# 			WHERE
# 				v.value !="N" AND v.value != "-"
# 			AND
# 				v.column <= (r.stop + CAST(? as integer)) AND v.column >= (r.start-CAST(? as integer))
# 			GROUP BY
# 				r.regid)
# 		UPDATE
# 			regions
# 		SET
# 			vars = vars + (SELECT counts FROM other WHERE regions.regid = other.regid)
# 		WHERE
# 			regid IN (SELECT regid FROM other WHERE regions.regid=other.regid)
# 	'''
# 	stuff = [dist, dist]
# 	cur.execute(sql, stuff)
# 	conn.commit()
#
# def parseFlankBad(conn, dist):
# 	cur = conn.cursor()
# 	sql2 = '''
# 		WITH other AS
# 			(SELECT
# 				r.regid,
# 				COUNT(v.value) as counts
# 			FROM
# 				regions AS r INNER JOIN variants AS v ON r.locid = v.locid
# 			WHERE
# 				v.value ="N"
# 			AND
# 				v.column <= (r.stop + CAST(? as integer)) AND v.column >= (r.start-CAST(? as integer))
# 			GROUP BY
# 				r.regid)
# 		UPDATE
# 			regions
# 		SET
# 			bad = bad + (SELECT counts FROM other WHERE regions.regid = other.regid)
# 		WHERE
# 			regid IN (SELECT regid FROM other WHERE regions.regid=other.regid)
# 	'''
# 	stuff = [dist, dist]
# 	cur.execute(sql2, stuff)
# 	conn.commit()
#
# def parseFlankGap(conn, dist):
# 	cur = conn.cursor()
# 	sql3 = '''
# 		WITH other AS
# 			(SELECT
# 				r.regid,
# 				COUNT(v.value) as counts
# 			FROM
# 				regions AS r INNER JOIN variants AS v ON r.locid = v.locid
# 			WHERE
# 				v.value = "-"
# 			AND
# 				v.column <= (r.stop + CAST(? as integer)) AND v.column >= (r.start-CAST(? as integer))
# 			GROUP BY
# 				r.regid)
# 		UPDATE
# 			regions
# 		SET
# 			gap = gap + (SELECT counts FROM other WHERE regions.regid = other.regid)
# 		WHERE
# 			regid IN (SELECT regid FROM other WHERE regions.regid=other.regid)
# 	'''
# 	stuff = [dist, dist]
# 	cur.execute(sql3, stuff)
# 	conn.commit()
#
#
# #Function to parse variants table to population flanking character columns for regions table
# def flankDistParser(conn, dist):
# 	parseFlankVars(conn, dist)
# 	parseFlankBad(conn, dist)
# 	parseFlankGap(conn, dist)

"""DEPRECATED"""
# #Function to parse variants table to update regions VARS for flanking information
# def getTargetFlanks(conn, dist):
# 	sql = '''
# 	SELECT
# 		r.regid, r.start AS start, r.stop AS stop, l.length AS len, l.consensus AS seq
# 	FROM
# 		regions AS r INNER JOIN loci AS l ON r.locid = l.id
# 	'''
# 	#Fetch targets with full loci sequences
# 	targets = pd.read_sql_query(sql, conn)
# 	targets.set_index('regid', inplace=True)
#
# 	#Re-pull targets + flank distances
# 	for index, row in targets.iterrows():
# 		start = row["start"] - dist
# 		if start < 0:
# 			start = 0
# 		end = row["stop"] + dist
# 		if end > row["len"]:
# 			end = row["len"]
# 		new_seq = row["seq"][start:end]
# 		targets.loc[index, 'seq'] = new_seq

#Function to return a pandas DF of regions, vars, and 'bad bases'
def getRegionWeights(conn):
	sql = """
	SELECT
		regid,
		(vars+vars_flank) AS vars,
		(bad + gap + bad_flank + gap_flank) AS sum_bad
	FROM
		regions
	"""
	return(pd.read_sql_query(sql ,conn))

#Function to return a pandas DF of regions, vars, and 'bad bases' of ONLY regids in a given list
def getRegionWeightsByList(conn, fetch):
	cur = conn.cursor()
	if len(fetch) <= 0:
		return(None)
	df = pd.DataFrame({"regid" : fetch})
	df.to_sql('ttt', conn, if_exists='replace')

	sql = """
	SELECT
		regions.regid,
		(vars+vars_flank) AS vars,
		(bad + gap + bad_flank + gap_flank) as sum_bad
	FROM
		regions
	WHERE
		regions.regid IN (SELECT regid FROM ttt)
	"""
	new_df = pd.read_sql_query(sql ,conn)

	cur.execute("DROP TABLE IF EXISTS ttt")
	conn.commit()
	return(new_df)

#Function to return a pandas DF of regions, vars, and 'bad bases' of ONLY regids in a given list
def getRegionWeightsByList_BAD(conn, fetch):
	cur = conn.cursor()
	if len(fetch) <= 0:
		return(None)
	df = pd.DataFrame({"regid" : fetch})
	df.to_sql('ttt', conn, if_exists='replace')

	sql = """
	SELECT
		regions.regid,
		(bad + gap + bad_flank + gap_flank) AS weight
	FROM
		regions
	WHERE
		regions.regid IN (SELECT regid FROM ttt)
	"""
	new_df = pd.read_sql_query(sql ,conn)

	cur.execute("DROP TABLE IF EXISTS ttt")
	conn.commit()

	#adjust weights to be: max(weight)-weight:
	max_weight = new_df["weight"].max()
	new_df["weight"] = max_weight - new_df["weight"]

	return(new_df)

#Function to return a pandas DF of regions, vars, and 'bad bases' of ONLY regids in a given list
def getRegionWeightsByList_VAR(conn, fetch):
	cur = conn.cursor()
	if len(fetch) <= 0:
		return(None)
	df = pd.DataFrame({"regid" : fetch})
	df.to_sql('ttt', conn, if_exists='replace')

	sql = """
	SELECT
		regions.regid,
		(vars+vars_flank) AS weight
	FROM
		regions
	WHERE
		regions.regid IN (SELECT regid FROM ttt)
	"""
	new_df = pd.read_sql_query(sql ,conn)

	cur.execute("DROP TABLE IF EXISTS ttt")
	conn.commit()
	return(new_df)

#Function to update REGIONS table based on existing Gap attribute
def simpleFilterTargets_gap(conn, thresh):
	cur = conn.cursor()
	cur.execute("UPDATE regions SET pass=0 WHERE gap > %s"%thresh)
	conn.commit()


#Function to update REGIONS table based on existing Bad attribute
def simpleFilterTargets_bad(conn, thresh):
	cur = conn.cursor()
	cur.execute("UPDATE regions SET pass=0 WHERE bad > ?", (thresh,))
	conn.commit()

#Function to update REGIONS table based on existing vars attribute
def simpleFilterTargets_SNP(conn, minS, maxS):
	cur = conn.cursor()
	cur.execute("UPDATE regions SET pass=0 WHERE (vars+vars_flank) < ? OR (vars+vars_flank) > ?",(minS,maxS))
	conn.commit()
