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
def make_consensus_from_vcf(ref, records, thresh):
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
			print("WARNING: Nucleotide (%s) at position %s in reference sequence does not match REF allele from VCF file (%s). This is most commonly caused by incorrect indexing in your VCF file."%(current_ref[rec.POS-1],rec.POS,rec.REF))

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
