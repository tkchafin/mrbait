#!/usr/bin/python
import os
import sys
import vcf
import misc_utils as utils
import alignment_tools as aln

#Read VCF variant calls
#Generator function, yields each locus
def read_vcf(v):

	if not utils.fileCheck(v):
		raise FileNotFoundError("Fatal exception, file %s not found."%v)

	try:
		vfh = vcf.Reader(filename=v)
	except IOError as err:
		print("I/O error({0}): {1}".format(err.errno, err.strerror))
	except:
		print("Unexpected error:", sys.exec_info()[0])

	chrom = ""
	recs = []
	added = 0
	for rec in vfh:
		if not rec.FILTER:
			if chrom:
				if chrom == rec.CHROM:
					recs.append(rec)
					added = 1
				else:
					#print("YIELDING")
					yield recs
					chrom = rec.CHROM
					recs = [rec]
					added = 0
			else:
				chrom = rec.CHROM
				recs.append(rec)
				added = 1
	if added == 0 and recs:
		yield recs
#NOTES:
#If reference base from FASTA is N or gap, we try to call new consensus from VCF
#N or - in VCF still have to pass threshold to be incoporated.
#If called as gap, we overwrite FASTA reference at that position.
#MASKING information is retained from FASTA reference and NOT considered in VCF
#Function to return new consensus sequence given REF and VCF records
def make_consensus_from_vcf(ref, chrom, records, thresh):
	consensus = ""
	for rec in records:
		current_ref = ""
		if consensus:
			current_ref = consensus
		else:
			current_ref = ref
		nucs = aln.get_iupac(current_ref[rec.POS-1].upper())
		rec_alt = []
		for x in rec.ALT:
			rec_alt += str(x).upper()

		cons = ""
		chosen = 0
		if rec.REF in nucs or rec.REF in ["N", "n", "-"] or current_ref[rec.POS-1] == "-":
			#Check if an allele is a gap
			if "-" in rec_alt:
				i = (rec_alt).index("-")
				prop = rec.aaf[i]
				if float(prop) >= float(thresh):
					cons = "-"
					chosen = 1
			if "N" in rec_alt:
				i = (rec_alt).index("N")
				prop = rec.aaf[i]
				if float(prop) >= float(thresh):
					cons = "N"
					chosen = 1
			#If no gap is kept, need to make new consensus
			#Add ALT alleles to list of observed nucs at this position
			if chosen == 0:
				nucs += rec_alt
				temp = utils.listToSortUniqueString(nucs)
				cons = aln.reverse_iupac(temp)
				chosen = 1
		else:
			print("\t\t\tWARNING: CHROM %s position %s (%s) doesn't match REF in VCF record (%s). "%(chrom, rec.POS, current_ref[rec.POS-1],rec.REF))

		#Incorporate new consensus base
		if chosen == 1 and cons:
			#Retain masking info from FASTA reference
			if current_ref[rec.POS-1].islower():
				consensus = utils.stringSubstitute(current_ref, (rec.POS-1), cons.lower())
			else:
				consensus = utils.stringSubstitute(current_ref, (rec.POS-1), cons.upper())
	#print("Reference:",ref)
	#print("Consensus:",consensus)
	return(consensus)

#function to count number of loci in FASTA file (by headers)
def countVCF(loci):
	fh  = open(str(loci), 'r')
	count=0
	for l in fh:
		line = l.strip()
		if not line:
			continue
		if not line.startswith("#"):
			count+=1
	return(count)

#Function to remove existing CHUNK files
def removeChunks(dir_name):
	test = os.listdir(dir_name)

	for item in test:
	    if item.endswith(".chunk"):
	        os.remove(os.path.join(dir_name, item))

#Function to chunk a VCF file
def vcf_chunker(infile, chunks, wd):
	chunks = int(chunks)
	loci_count = countVCF(infile)
	if loci_count < chunks:
		chunks = loci_count
	chunk_size = loci_count // chunks
	removeChunks(wd) #clear any existing chunkfiles

	files = list()
	#write .loci file into chunk files
	with open(infile) as file_object:
		max_chunks = chunks
		chunks = 1
		record = 0

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
				if line[0] == "#":
					header +=str(line+"\n")
				else:
					record = 1
					hset=1
					out_object.write(header)
					out = line + "\n"
					out_object.write(out)
			else:
				#Write chunk_size alignments to each chunk
				if chunks < max_chunks:
					record += 1 #increment locus number
				 	#If current chunk not full, add record to chunk
					if record <= chunk_size:
						out = line + "\n"
						out_object.write(out)
					#Otherwise, start new chunk
					else:
						record= 1
						chunks = chunks + 1
						out_object.close()
						chunk_file = wd + "/." + str(chunks) + ".chunk"
						out_object = open(chunk_file, "w")
						out_object.write(header)
						files.append(chunk_file)
						out = line + "\n"
						out_object.write(out)

				#If last chunk, keep writing to final chunk file
				else:
					out = line + "\n"
					out_object.write(out)
		out_object.close()
		file_object.close()
	return(files)
