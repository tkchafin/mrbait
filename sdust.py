class DustMasker:
	"""Python class for applying DUST algorithm to
	a nucleotide sequence, to find low-complexity
	regions"""

	def __init__ (seq, score):
		self.seq = seq
		self.score = score
		self.masked = None
		self.triplets = []

	def add_triplet(r,c,t):
		r = r + c[t]
		c[t] = c[t]+1


seq="AAATGTATTCGCGTCGCGCGCGGGCGCGCGGCGCGCGTATACTAGCTAGCTGACTGATCTAGCTATCGATCGACTG"
mask = DustMasker(seq, 10)
