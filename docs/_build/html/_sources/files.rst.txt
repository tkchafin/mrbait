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

It is important to note that the VCF format can communicate much more information
than mrbait_ will utilize. The CHROM and POS columns will be parsed to locate the
reference position for each SNP, and the REF and ALT columns to write a new consensus
base at that position using IUPAC ambiguity codes (e.g. C/T = Y). More functionality
will be added in future versions of mrbait_.

It is highly recommended you add variant data if it is available, as it will be used
both for finding adequately conserved regions for bait design, as well as for filtering
target regions for those which capture flanking SNPs.

Annotating genomes with GFF
~~~~~~~~~~~~~~~~~~~~~~~~~~~

mrbait_ can also make use of genomic features provided using the Generic Feature Format (GFF),
independently or in addition to any variant data provided via VCF. mrbait assumes that input
GFF files follow the version 3 `GFF specification <https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md>`\_:

.. code-block:: none
   :linenos:

   ##gff-version 3
   chr1.scaffold1	.	gene	10	180	.	+	.	ID=gene0001;Alias=targets
   chr1.scaffold1	.	mRNA	20	180	.	+	.	ID=mrna0001;Parent=gene0001
   chr1.scaffold1	.	exon	10	128	.	+	.	ID=tfbs00001;Parent=gene0001
   …
   …

Columns should be separated by tabs and defined according to the GFF3 standard (e.g.
column 1 contains the sequence ID). mrbait will use the sequence ID (column 1) to map
coordinates in GFF columns 4 and 5 to the reference provided in your FASTA file, thus
these identifiers must be identical. mrbait will also categorize features internally by
the type (e.g. “exon”) given in column 3, and by any alias assigned in the attributes
column (column 9). All other columns are ignored. You can use either type or alias to
tell mrbait_ to target those features for bait design.

If you are not targeting all of a single type (e.g. CDS, or exon), you can either pre-filter
your GFF file prior to loading, or you can annotate features of interest using the Alias
attribute. 
