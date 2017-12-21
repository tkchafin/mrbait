import os
import sys
import vcf

#Read VCF variant calls
#Generator function, yields each locus
def read_vcf(v):
	try:
		vfh = vcf.Reader(filename=v)
	except IOError as err:
		print("I/O error({0}): {1}".format(err.errno, err.strerror))
	except:
		print("Unexpected error:", sys.exec_info()[0])

	num = 1
	for rec in vfh:
		# print("Record",num)
		# print("Type:",rec.var_type)
		# print("Subtype:",rec.var_subtype)
		# samples = rec.samples
		# print("Samples in record:")
		# for samp in samples:
		# 	print(samp.called, samp.gt_alleles)
		# num+=1
		# print()
		if not rec.FILTER:
			print(rec.CHROM, rec.REF, rec.ALT, len(rec.samples), rec.call_rate, rec.aaf)

	# Pandas DF to return:
	#	IF FILTER=PASS:
	#		CHROM REF ALT Nsamples CALL_RATE AAF (Allele Freq)

	sys.exit()
