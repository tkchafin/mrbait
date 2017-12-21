import os
import sys
import vcf
import alignment_tools as aln

#Read VCF variant calls
#Generator function, yields each locus
def read_vcf(v):
	try:
		vfh = vcf.Reader(filename=v)
	except IOError as err:
		print("I/O error({0}): {1}".format(err.errno, err.strerror))
	except:
		print("Unexpected error:", sys.exec_info()[0])

	chrom = ""
	recs = []
	for rec in vfh:
		if not rec.FILTER:
			print(rec)
			if chrom:
				if chrom == rec.CHROM:
					recs.append(rec)
					continue
				else:
					yield recs
				chrom = ""
				recs = []
			else:
				chrom = rec.CHROM
				recs.append(rec)
				continue

#Function to return new consensus sequence given REF and VCF records
def make_consensus_from_vcf(ref, records, thresh):
	for rec in records:
		if rec.REF == ref[rec.POS-1]:
			nucs = aln.get_iupac(ref[rec.POS-1])
			#nucs +=
			pass
		else:
			raise ValueError("Nucleotide (%s) at position %s in reference sequence does not match REF allele from VCF file (%s). This is most commonly caused by incorrect indexing in your VCF file."%(ref[rec.POS-1],rec.POS,rec.REF))
