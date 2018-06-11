.. include:: global.rst

.. _files:

Input files
===========

This section describes the input file types accepted by MrBait.

Assembled genomes
-----------------
MrBait only accepts genome assemblies formatted as FASTA. These can represent
contigs, scaffolds, or entire chromosomes. According to the FASTA specifications,
a sequence should begin with a header line, or short description (indicated by the
“>” symbol), followed by a second line containing sequence data. For MrBait, it does
not matter if the following lines are interleaved or on a single line, and any blank
lines in the file will be ignored, as will any leading or trailing whitespace.

An example FASTA-formatted sequence is given below.

.. code-block:: none
   :linenos:

   >chr1.scaffold1
   ATAGCTCGGCTACGTGATCGCGTGCTC-ATGCTAGCGCTNNNNNNNNATGATTGCTTTT
   TGTGTGTGCAAGCACTGCCGRGCTACGCGCTACTGCCRCCTAGTATGTGTGGCCGCTAC
   TAGTCCGCGCTAGCTtTtagatctcgtggcgccgcgcgcgtcgcacgatcgtacgcgcc
   >chr1.scaffold2
   ATCGTGCTGCGGCGCTGCCTCAGC…
   …
   …
   …

Annotating genomes with VCF
~~~~~~~~~~~~~~~~~~~~~~~~~~~

mrbait_ also supports supplementing genomic sequences with coordinate-reference SNP
data (e.g. obtained from population-level sequencing) using the `Variant Call Format <http://samtools.github.io/hts-specs/>`_:

.. code-block:: none
   :linenos:

   ##fileformat=VCFv4.2
   ##FORMAT=<ID=GT,Number=1,Type=Integer,Description="Genotype">
   ##FORMAT=<ID=GP,Number=G,Type=Float,Description="Genotype Probabilities">
   ##FORMAT=<ID=PL,Number=G,Type=Float,Description="Phred-scaled Genotype Likelihoods">
   #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMP001	SAMP002
   chr1.scaffold1	48	rs11449	G	A	.	PASS	.	GT	0/0	0/1
   chr1.scaffold1	47	rs11449	T	A	.	PASS	.	GT	0/0	0/1
   chr1.scaffold2	1 rs84825 A	T	.	PASS	.	GT:GP	0/1:.	0/1:0.03,0.97,0
   …
   …
