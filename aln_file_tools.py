#!/usr/bin/python

import os
import sys
import Bio
import vcf
import misc_utils as utils
from Bio import AlignIO

"""Functions for parsing and manipulating sequence alignment files
Functions by Zach Zbinden and Tyler Chafin"""

#Write FASTA from pandas df where col1 is index, col2 is sequence
#seqs must be a pandas df
def writeFasta(seqs, fas):
	file_object = open(fas, "w")
	#Write seqs to FASTA first
	#Assumes that a[0] is index, a[1] is id, and a[2] is sequence
	for a in seqs.itertuples():
		name = ">id_" + str(a[1]) + "\n"
		seq = a[2] + "\n"
		file_object.write(name)
		file_object.write(seq)
	file_object.close()

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

#Function by ZVZ to "chunk" a given MAF alignment file into n number of chunks
def maf_chunker(infile, chunks):
	# maf_chunker creates specified number of files containing equal numbers
	# of loci (unless there are remainder loci, which append to the last
	# chunk.
	# 1 to n '.maf_chunk' files will be created

# read file from command line
	with open(infile) as file_object:
#count number of loci, loci_count = -1 so that header is not counted
		loci_count = -1
		chunks = int(sys.argv[2])

		for line in file_object:
			line = line.strip()
			if len(line) > 0:
				pass
			else:
				loci_count = loci_count+1
		chunk_size = loci_count // chunks

#write .maf file into chunk files, with each chunk beginning with header
#first read header
	with open(infile) as file_object:
		max_chunks = int(sys.argv[2])
		chunks = 0
		loci_number = 0
		individual = 1

		for line in file_object:
			line = line.strip()

#isolate header chunk
			if loci_number == 0:
				if len(line) > 0:
					print(line.strip(), file=open(str(chunks) + ".maf_chunk", "a"))
				else:
					loci_number = loci_number + 1
					chunks = chunks + 1
#move to loci chunks
			else:
				if chunks < max_chunks:
					if loci_number <= chunk_size:
#print contents of header before printing loci of individual 1
						if individual == 1 and chunks == 1:
							with open('0.maf_chunk') as header:
								for var in header:
									print(var.strip(), file=open(str(chunks) + ".maf_chunk", "a"))

							print("", file=open(str(chunks) + ".maf_chunk", "a"))
							print(line.strip(), file=open(str(chunks) + ".maf_chunk", "a"))

							individual = individual + 1
						else:
							if len(line) > 0:
								print(line.strip(), file=open(str(chunks) + ".maf_chunk", "a"))
								individaul = individual + 1

							else:
								loci_number = loci_number + 1
								individual = 1
								print("", file=open(str(chunks) + ".maf_chunk", "a"))
					else:
						loci_number = 1
						chunks = chunks + 1
						individual = 1

						with open('0.maf_chunk') as header:
								for var in header:
									print(var.strip(), file=open(str(chunks) + ".maf_chunk", "a"))
						print("", file=open(str(chunks) + ".maf_chunk", "a"))

						print(line.strip(), file=open(str(chunks) + ".maf_chunk", "a"))
				else:
					chunks = max_chunks
					print(line.strip(), file=open(str(chunks) + ".maf_chunk", "a"))

	os.remove("0.maf_chunk")

#Function to remove existing CHUNK files
def removeChunks(dir_name):
	test = os.listdir(dir_name)

	for item in test:
	    if item.endswith(".chunk"):
	        os.remove(os.path.join(dir_name, item))

#function to count number of loci alignments in file
def countLoci(loci):
	file  = open(loci, 'r').read()
	return(file.count("//"))


#Function by ZDZ to split a given .loci file into n chunks
def loci_chunker(infile, chunks, wd):
	# read file from command line
	with open(infile) as file_object:

		#count number of loci
		loci_count = 1
		chunks = int(chunks)

		for line in file_object:
			if line[0] == "/":
				loci_count = loci_count+1
			else:
				pass
		chunk_size = loci_count // chunks
	file_object.close()
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

		for line in file_object:
			if chunks < max_chunks:
				if loci_number <= chunk_size:
					if line[0] == ">":
						out = line.strip() + "\n"
						out_object.write(out)
					else:
						loci_number = loci_number + 1
						out = line.strip() + "\n"
						out_object.write(out)
				else:
					loci_number = 1
					chunks = chunks + 1
					out_object.close()
					chunk_file = wd + "/." + str(chunks) + ".chunk"
					out_object = open(chunk_file, "w")
					files.append(chunk_file)
					out = line.strip() + "\n"
					out_object.write(out)
			else:
				#If last chunk, keep writing to final chunk file
				out = line.strip() + "\n"
				out_object.write(out)
		out_object.close()
		file_object.close()
			# else:
			# 	chunks = max_chunks
			# 	out_object.write(line.strip())
	return(files)
