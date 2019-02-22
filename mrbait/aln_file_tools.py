#!/usr/bin/python

import os
import sys
import Bio
import vcf
from mrbait import misc_utils as utils
from Bio import AlignIO

"""Functions for parsing and manipulating sequence alignment files
Functions by Zach Zbinden and Tyler Chafin"""

#Write FASTA from pandas df where col1 is index, col2 is sequence
#seqs must be a pandas df
def writeFasta(seqs, fas):
	with open(fas, 'w') as fh:
		try:
			#Write seqs to FASTA first
			#Assumes that a[0] is index, a[1] is id, and a[2] is sequence
			for a in seqs.itertuples():
				name = ">id_" + str(a[1]) + "\n"
				seq = a[2] + "\n"
				fh.write(name)
				fh.write(seq)
		except IOError as e:
			print("Could not read file:",e)
			sys.exit(1)
		except Exception as e:
			print("Unexpected error:",e)
			sys.exit(1)
		finally:
			fh.close()

#Write FASTA from pandas df where col1 is index, col2 is sequence
#seqs must be a pandas df
#this version replaces gaps with N characters
def writeFastaNogap(seqs, fas):
	with open(fas, 'w') as fh:
		try:
			#Write seqs to FASTA first
			#Assumes that a[0] is index, a[1] is id, and a[2] is sequence
			for a in seqs.itertuples():
				name = ">id_" + str(a[1]) + "\n"
				seq = str(a[2]) + "\n"
				seq = seq.replace("-","N")
				fh.write(name)
				fh.write(seq)
		except IOError as e:
			print("Could not read file:",e)
			sys.exit(1)
		except Exception as e:
			print("Unexpected error:",e)
			sys.exit(1)
		finally:
			fh.close()

#Read genome as FASTA. FASTA header will be used
#This is a generator function
#Doesn't matter if sequences are interleaved or not.
def read_fasta(fas):

	if not utils.fileCheck(fas):
		raise FileNotFoundError("Fatal exception, file %s not found."%fas)
	fh = open(fas)
	try:
		with fh as file_object:
			contig = ""
			seq = ""
			for line in file_object:
				line = line.strip()
				if not line:
					continue
				line = line.replace(" ","")
				#print(line)
				if line[0] == ">": #Found a header line
					#If we already loaded a contig, yield that contig and
					#start loading a new one
					if contig:
						yield([contig,seq]) #yield
						contig = "" #reset contig and seq
						seq = ""
					contig = (line.replace(">",""))
				else:
					seq += line
		#Iyield last sequence, if it has both a header and sequence
		if contig and seq:
			yield([contig,seq])
	finally:
		fh.close()


#This is a GENERATOR function to read through a .loci file
#.loci is the RAD alignment output from the promgram pyRAD
#YIELDS: BioPython MultipleSeqAlignment object
def read_loci(infile):

	if not utils.fileCheck(infile):
		raise FileNotFoundError("Fatal exception, file %s not found."%infile)

	# make emptyp dictionary
	loci = Bio.Align.MultipleSeqAlignment([])
	# read file from command line
	try:
		f = open(infile)
	except IOError as err:
		print("I/O error({0}): {1}".format(err.errno, err.strerror))
	except:
		print("Unexpected error:", sys.exec_info()[0])

	with f as file_object:
		for line in file_object:
			line = line.strip()
			if not line:
				continue
			if line[0] == ">":
				identifier = line.split()[0]
				sequence = line.split()[1]
				loci.add_sequence(identifier, sequence)
			else:
				yield(loci)
				loci = Bio.Align.MultipleSeqAlignment([])
	f.close()

#Function to remove existing CHUNK files
def removeChunks(dir_name):
	test = os.listdir(dir_name)

	for item in test:
	    if item.endswith(".chunk"):
	        os.remove(os.path.join(dir_name, item))

#function to count number of loci alignments in file
def countLoci(loci):
	fh  = open(loci, 'r')
	count=0
	for l in fh:
		line = l.strip()
		if not line:
			continue
		if line.startswith("//"):
			count+=1
	return(count)

#function to count number of loci in FASTA file (by headers)
def countMAF(loci):
	fh  = open(str(loci), 'r')
	count=0
	for l in fh:
		line = l.strip()
		if not line:
			continue
		if line.startswith("a"):
			count+=1
	return(count)

#function to count number of loci in FASTA file (by headers)
def countXMFA(loci):
	fh  = open(str(loci), 'r')
	count=0
	for l in fh:
		line = l.strip()
		if not line:
			continue
		if line.startswith("="):
			count+=1
	return(count)


#Function split .loci file into n chunks
def loci_chunker(infile, chunks, wd):

	chunks = int(chunks)
	loci_count = countLoci(infile)
	if loci_count < chunks:
		chunks = loci_count
	chunk_size = loci_count // chunks
	removeChunks(wd)

	files = list()
	#write .loci file into chunk files
	with open(infile) as file_object:
		max_chunks = chunks
		chunks = 1
		loci_number = 1

		chunk_file = wd + "/." + str(chunks) + ".chunk"
		out_object = open(chunk_file, "w")
		files.append(chunk_file)

		for l in file_object:
			line = l.strip()
			if not line:
				continue
			if chunks < max_chunks:
				if loci_number <= chunk_size:
					if line[0] == ">":
						out = line + "\n"
						out_object.write(out)
					else:
						loci_number = loci_number + 1
						out = line + "\n"
						out_object.write(out)
				else:
					loci_number = 1
					chunks = chunks + 1
					out_object.close()
					chunk_file = wd + "/." + str(chunks) + ".chunk"
					out_object = open(chunk_file, "w")
					files.append(chunk_file)
					out = line + "\n"
					out_object.write(out)
			else:
				#If last chunk, keep writing to final chunk file
				out = line + "\n"
				out_object.write(out)
		out_object.close()
		file_object.close()
			# else:
			# 	chunks = max_chunks
			# 	out_object.write(line.strip())
	return(files)

#Function to split maf file into n chunks
def maf_chunker(infile, chunks, wd):

	chunks = int(chunks)
	loci_count = countMAF(infile)
	if loci_count < chunks:
		chunks = loci_count
	chunk_size = loci_count // chunks
	removeChunks(wd) #clear any existing chunkfiles

	files = list()
	#write .loci file into chunk files
	with open(infile) as file_object:
		max_chunks = chunks
		chunks = 1
		loci_number = 0

		chunk_file = wd + "/." + str(chunks) + ".chunk"
		out_object = open(chunk_file, "w")
		files.append(chunk_file)

		header = ""
		hset = 0
		for l in file_object:
			line = l.strip()
			if not line:
				continue
			#First, get header information
			if hset == 0:
				if line[0] == "a":
					loci_number = 1
					hset=1
					header += "\n"
					out_object.write(header)
					out = "\n" + line + "\n"
					out_object.write(out)
				elif line[0] =="#":
					header +=str(line+"\n")
			else:
				#Write chunk_size alignments to each chunk
				if chunks < max_chunks:
					#If starting new alignment
					if line[0] == "a":
						loci_number += 1 #increment locus number
					 	#If current chunk not full, add locus to chunk
						if loci_number <= chunk_size:
							out = "\n" + line + "\n"
							out_object.write(out)
						#Otherwise, start new chunk
						else:
							loci_number = 1
							chunks = chunks + 1
							out_object.close()
							chunk_file = wd + "/." + str(chunks) + ".chunk"
							out_object = open(chunk_file, "w")
							out_object.write(header)
							files.append(chunk_file)
							out = line + "\n"
							out_object.write(out)
					#If not new alignment, write to current chunk
					else:
						out = line + "\n"
						out_object.write(out)

				#If last chunk, keep writing to final chunk file
				else:
					if line[0] == "a":
						out = "\n" + line + "\n"
					else:
						out = line + "\n"
					out_object.write(out)
		out_object.close()
		file_object.close()
			# else:
			# 	chunks = max_chunks
			# 	out_object.write(line.strip())
	return(files)

#Function to split xmfa file into n chunks
def xmfa_chunker(infile, chunks, wd):

	chunks = int(chunks)
	loci_count = countXMFA(infile)
	if loci_count < chunks:
		chunks = loci_count
	chunk_size = loci_count // chunks
	removeChunks(wd)

	files = list()
	#write .loci file into chunk files
	with open(infile) as file_object:
		max_chunks = chunks
		chunks = 1
		loci_number = 1

		chunk_file = wd + "/." + str(chunks) + ".chunk"
		out_object = open(chunk_file, "w")
		files.append(chunk_file)

		header = ""
		hset = 0
		for l in file_object:
			line = l.strip()
			if not line:
				continue
			#First, get header information
			if hset == 0:
				if line[0] == ">":
					loci_number = 1
					hset=1
					out_object.write(header)
					out = line + "\n"
					out_object.write(out)
				elif line[0] =="#":
					header +=str(line+"\n")
			else:
				if chunks < max_chunks:
					if loci_number <= chunk_size:
						#If its the header for a sequence, start seq
						if line[0] == ">":
							out = line + "\n"
							out_object.write(out)
						#If end of alignment, deposit it
						elif line[0] == "=":
							loci_number = loci_number + 1
							out = line + "\n"
							out_object.write(out)
						#otherwise its a sequence!
						else:
							out = line + "\n"
							out_object.write(out)
					else:
						loci_number = 1
						chunks = chunks + 1
						out_object.close()
						chunk_file = wd + "/." + str(chunks) + ".chunk"
						out_object = open(chunk_file, "w")
						out_object.write(header)
						files.append(chunk_file)
						out = line + "\n"
						out_object.write(out)
				else:
					#If last chunk, keep writing to final chunk file
					out = line + "\n"
					out_object.write(out)
		out_object.close()
		file_object.close()
			# else:
			# 	chunks = max_chunks
			# 	out_object.write(line.strip())
	return(files)
