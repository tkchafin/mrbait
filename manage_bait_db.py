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
	cursor.execute('''DROP TABLE IF EXISTS baits''')
	cursor.execute('''DROP TABLE IF EXISTS gff''')

	#Table holding records for each locus
	cursor.execute('''
		CREATE TABLE loci(id INTEGER PRIMARY KEY, depth INTEGER NOT NULL,
			length INTEGER NOT NULL, consensus TEXT NOT NULL, pass INTEGER NOT NULL,
			chrom TEXT,
			UNIQUE(chrom) ON CONFLICT ROLLBACK)
	''') #chrom only relevant if --assembly

	#Table holding records for each locus
	cursor.execute('''
		CREATE TABLE regions(regid INTEGER PRIMARY KEY, locid INTEGER NOT NULL,
			length INTEGER NOT NULL, sequence TEXT NOT NULL, vars INTEGER,
			bad INTEGER, gap INTEGER, mask REAL, gc REAL, start INTEGER NOT NULL,
			stop INTEGER NOT NULL,pass INTEGER NOT NULL,
			FOREIGN KEY (locid) REFERENCES loci(id))
	''')

	#Table holding variant information
	cursor.execute('''
		CREATE TABLE variants(varid INTEGER NOT NULL, locid INTEGER NOT NULL,
			column INTEGER NOT NULL, value TEXT NOT NULL,
			FOREIGN KEY (locid) REFERENCES loci(id),
			PRIMARY KEY(varid),
			UNIQUE(varid))
	''')

	#Table holding GFF element information
	cursor.execute('''
		CREATE TABLE gff(gffid INTEGER PRIMARY KEY, locid INTEGER NOT NULL,
			type TEXT NOT NULL, start INTEGER NOT NULL,
			stop INTEGER NOT NULL, alias TEXT, pass INTEGER NOT NULL,
			FOREIGN KEY (locid) REFERENCES loci(id))
	''')

	cursor.execute('''
		CREATE TABLE baits(baitid INTEGER NOT NULL, regid INTEGER NOT NULL,
			sequence TEXT NOT NULL, start INTEGER NOT NULL, stop INTEGER NOT NULL,
			mask REAL, gc REAL, pass INTEGER NOT NULL,
			PRIMARY KEY(baitid),
			FOREIGN KEY (regid) REFERENCES regions(regid)
		)
	''')
	#Table holding records for each locus
	#cursor.execute('''
	#	CREATE TABLE samples(sampid INTEGER NOT NULL, name TEXT NOT NULL,
	#		PRIMARY KEY (sampid),
	#		UNIQUE(name))
	#''')

	connection.commit()

#Function to filter loci by coverage and length
def filterLoci(conn, minlen, mincov):
	cur = conn.cursor()
	sql = '''
	UPDATE
		loci
	SET
		pass=0
	WHERE
		length < ? OR depth < ?
	'''
	stuff = [minlen, mincov]
	cur.execute(sql, stuff)
	conn.commit()

#Function returns pandas dataframe of passedLoci
def getPassedLoci(conn):
	return(pd.read_sql_query("""SELECT id, consensus, chrom FROM loci WHERE pass=1""", conn))

#Function returns pandas dataframe of passed targets
def getPassedTRs(conn):
	return(pd.read_sql_query("""SELECT regid, sequence FROM regions WHERE pass=1""", conn))

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

	return((cur.fetchone()[0]))

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

	return((cur.fetchone()[0]))

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
	return(pd.read_sql_query("""SELECT * FROM baits WHERE pass=1""", conn))

#Function to return variants table
def getVariants(conn):
	return(pd.read_sql_query("""SELECT * FROM variants """, conn))

#Code to add record to 'loci' table
def add_locus_record(conn, depth, consensus, passed, name):
	if name == None:
		name = "NA"
	stuff = [depth, int(len(consensus)), str(consensus), int(passed), str(name)]
	try:
		sql = ''' INSERT INTO loci(depth, length, consensus, pass, chrom)
					VALUES(?,?,?,?,?) '''
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
		print("Locid for",seqid,"is:",locid)
		#BUILT INSERT SQL
		sql = ''' INSERT INTO gff(locid, type, start, stop, alias, pass)
					VALUES(?,?,?,?,?,1);'''
		stuff = [int(locid), str(gff_type), int(start), int(stop), str(alias)]
		cur.execute(sql, stuff)
		conn.commit()

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
def add_region_record(conn, locid, start, stop, seq, counts, mask, gc):
	#Establish cursor
	cur = conn.cursor()

	mask_p = float(mask/len(seq))
	gc_p = float(gc/len(seq))
	#build sql and pack values to insert
	sql = '''INSERT INTO regions(locid, length, sequence, vars, bad, gap, mask, gc,
		start, stop, pass) VALUES (?,?,?,?,?,?,?,?,?,?,1)'''
	stuff = [locid, len(seq), seq, counts["*"], counts["N"], counts["-"], mask_p, gc_p, start, stop]

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

#Function to filter regions relation by minimum flanking SNPs
def regionFilterMinVar(conn, val, flank):
	cur = conn.cursor()
	#Update pass to "0/FALSE" where COUNT(vars) in flanking region + TR is less than minimum
	sql = '''
		UPDATE regions
		SET pass = 0
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
				((column < (stop+?)) AND (column > start-?))
			GROUP BY regid)
			WHERE counts < ?);

	'''
	stuff = [flank, flank, val]
	cur.execute(sql, stuff)
	conn.commit()

#Function to filter regions relation by maximum flanking SNPs
def regionFilterMaxVar(conn, val, flank):
	cur = conn.cursor()
	#Update pass to "0/FALSE" where COUNT(vars) in flanking region + TR is greater than maximum
	sql = '''
		UPDATE regions
		SET pass = 0
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
				((column < (stop+?)) AND (column > start-?))
			GROUP BY regid)
			WHERE counts > ?);

	'''
	stuff = [flank, flank, val]
	cur.execute(sql, stuff)
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
		((column <= (stop+%s)) AND (column >= start-%s))
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
	cur.execute("SELECT COUNT(*) FROM regions WHERE pass=0")
	fails = int(cur.fetchone()[0])

	print("Number of rows:",rows)
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

#Function to delete from variable table on locid
def purgeVars(conn, key):
	cur = conn.cursor()

	sql = '''
	DELETE FROM variants
	WHERE locid = ?
	'''

	cur.execute(sql,(key,))
	conn.commit()

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
	for group, group_df in groups:
		#print("\n########Group is: ",group,"\n")
		#If only one TR for alignment, set choose to 1:
		if group_df.shape[0] == 1:
			for name, row in group_df.iterrows():
				#modify original dataframe
				df.loc[name, "conflict_block"] = row["locid"]
		else:
			#For each TR in locus
			for name, row in group_df.iterrows():
				#print("\nName in group_df:",name)
				#print("Last INDEX in group_df:",group_df.index[-1])

				#If alignment length below minlen, set conflict_group to locid
				if row["length"] <= min_len:
					df.loc[name, "conflict_block"] = row["locid"]
				else:
					#Else, compare with RIGHT NEIGHBOR
					#Assumes ordered left to right within locus
					#Maybe make sure to ORDER BY in query??
					_name = int(name) + 1

					#If this is the last TR in locus (i.e. no right neighbor)
					if name == (group_df.index[-1]):
						if df.loc[name, "conflict_block"] == "NULL":
							df.loc[name, "conflict_block"] = block
							block += 1
						break
					else:
						_row = group_df.loc[_name]
						#print(_row)
						#	print(_row)
						#for _name, _row in group_df.iterrows():
						#if name == _name:
						#	continue
						#else:
						#If within dist_r bases, assign same conflict block
						#print("Comparing row ",name," and ", _name)
						#print()
						#If they overlap, assign them both to same conflict block
						if checkOverlap(row, _row, dist) == 1:
							#print("overlap!")
							#Fetch conflict_block currents from parent df
							cb = df.loc[name, "conflict_block"]
							_cb = df.loc[_name, "conflict_block"]
							#If neither is assigned, assign both to new block
							if (cb == "NULL") and (_cb == "NULL"):
								#print("Assigning block: ",block)
								df.loc[_name, "conflict_block"] = block
								df.loc[name, "conflict_block"] = block
								block+=1
							#Else if row1 is assigned, give row2 the block of row1
							elif _cb == "NULL":
								df.loc[_name, "conflict_block"] = df.loc[name, "conflict_block"]
							#Or if row2 has block, assign to row1
							elif cb == "NULL":
								df.loc[name, "conflict_block"] = df.loc[_name, "conflict_block"]
							#If no overlap, assign Row1 to its own conflict_block
						else:
							df.loc[name, "conflict_block"] = block
							block += 1

	#Next step, send df back to SQLite and update conflicts table
	updateConflictsFromPandas(conn, df)
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
			vars as counts
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
def printFlankingSNPs_conflicts(conn):
	cur = conn.cursor()

	#Query conflicts temp table to make a pandas dataframe for parsing
	sql = '''
		SELECT
			c.regid,
			conflict_block,
			choose,
			vars
		FROM
			conflicts AS c INNER JOIN regions AS r USING (regid)
		WHERE
			c.choose="NULL"
	'''
	print(pd.read_sql_query(sql, conn))


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
			else:
				message = "Error: No TR chosen for block " + group + ", something went wrong"
				raise RuntimeError(message)
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
			else:
				message = "Error: No TR chosen for block " + group + ", something went wrong"
				raise RuntimeError(message)
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
			vars AS counts
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
			SUM(gap + bad) as counts
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

#Function for resolving conflict blocks by minimizing "bad" bases in flanking region
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
			(vars + gap + bad) AS counts
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
			vars as counts
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
		vars > ?
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
	print(df)

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

#Function to parse variants table to update regions VARS for flanking information
def parseFlankVars(conn, dist):
	cur = conn.cursor()
	sql = '''
		WITH other AS
			(SELECT
				r.regid,
				COUNT(v.value) as counts
			FROM
				regions AS r INNER JOIN variants AS v ON r.locid = v.locid
			WHERE
				v.value !="N" AND v.value != "-"
			AND
				v.column <= (r.stop + CAST(? as integer)) AND v.column >= (r.start-CAST(? as integer))
			GROUP BY
				r.regid)
		UPDATE
			regions
		SET
			vars = vars + (SELECT counts FROM other WHERE regions.regid = other.regid)
		WHERE
			regid IN (SELECT regid FROM other WHERE regions.regid=other.regid)
	'''
	stuff = [dist, dist]
	cur.execute(sql, stuff)
	conn.commit()

def parseFlankBad(conn, dist):
	cur = conn.cursor()
	sql2 = '''
		WITH other AS
			(SELECT
				r.regid,
				COUNT(v.value) as counts
			FROM
				regions AS r INNER JOIN variants AS v ON r.locid = v.locid
			WHERE
				v.value ="N"
			AND
				v.column <= (r.stop + CAST(? as integer)) AND v.column >= (r.start-CAST(? as integer))
			GROUP BY
				r.regid)
		UPDATE
			regions
		SET
			bad = bad + (SELECT counts FROM other WHERE regions.regid = other.regid)
		WHERE
			regid IN (SELECT regid FROM other WHERE regions.regid=other.regid)
	'''
	stuff = [dist, dist]
	cur.execute(sql2, stuff)
	conn.commit()

def parseFlankGap(conn, dist):
	cur = conn.cursor()
	sql3 = '''
		WITH other AS
			(SELECT
				r.regid,
				COUNT(v.value) as counts
			FROM
				regions AS r INNER JOIN variants AS v ON r.locid = v.locid
			WHERE
				v.value = "-"
			AND
				v.column <= (r.stop + CAST(? as integer)) AND v.column >= (r.start-CAST(? as integer))
			GROUP BY
				r.regid)
		UPDATE
			regions
		SET
			gap = gap + (SELECT counts FROM other WHERE regions.regid = other.regid)
		WHERE
			regid IN (SELECT regid FROM other WHERE regions.regid=other.regid)
	'''
	stuff = [dist, dist]
	cur.execute(sql3, stuff)
	conn.commit()

#Function to parse variants table to population flanking character columns for regions table
def flankDistParser(conn, dist):
	parseFlankVars(conn, dist)
	parseFlankBad(conn, dist)
	parseFlankGap(conn, dist)

#Function to return a pandas DF of regions, vars, and 'bad bases'
def getRegionWeights(conn):
	sql = """
	SELECT
		regid,
		vars,
		(bad + gap) AS sum_bad
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
		vars,
		(bad + gap) as sum_bad
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
		(bad + gap) AS weight
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
		vars AS weight
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
	cur.execute("UPDATE regions SET pass=0 WHERE vars < ? OR vars > ?",(minS,maxS))
	conn.commit()
