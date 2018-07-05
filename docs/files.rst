.. include:: global.rst

.. _files:

Input files
===========

This section describes the input file types accepted by MrBait.

Assembled genomes
-----------------
mrbait_ only accepts genome assemblies formatted as FASTA. These can represent
contigs, scaffolds, or entire chromosomes. According to the FASTA specifications,
a sequence should begin with a header line, or short description (indicated by the
“>” symbol), followed by a second line containing sequence data. It does
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

NOTE: When using VCF, the REF column is ignored. Instead, the reference allele will be
taken from the FASTA reference provided. For cases when the reference allele is an N or
gap (-), you can choose to either retain the N/gap allele, OR attempt to override it
using the ALT alleles provided in the VCF for that position (*--vcfALT*)

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

Multiple genome alignments
--------------------------

mrbait_ reads two different input file types for multiple genome alignments. These can
be provided using the Multiple Alignment Format (`MAF <https://genome.ucsc.edu/FAQ/FAQformat.html#format9.3>`_), or the eXtended Multi-FastA (`XMFA <https://asap.genetics.wisc.edu/software/mauve/mauve-user-guide/mauve-output-file-formats.php>`_) formats.

The MAF format is output by several multiple alignment programs, including `MAFFT <https://mafft.cbrc.jp/alignment/software/>`_
and `Mugsy <http://mugsy.sourceforge.net/>`_, and take the following general form:

.. code-block:: none
   :linenos:

   ##maf version=1 scoring=tba.v8
   # tba.v8 (((human chimp) baboon) (mouse rat))
   # multiz.v7
   # maf_project.v5 _tba_right.maf3 mouse _tba_C
   # single_cov2.v4 single_cov2 /dev/stdin

   a score=5062.0
   s hg16.chr7    27699739 6 + 158545518 RAAAGAGATGCTAAGCCAATGAGTTGATGTCTCTCAATGTGTG
   s panTro1.chr6 28862317 6 + 161576975 RAAAGAGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTGTG
   s baboon         241163 6 +   4622798 TAAAGAGATGCTAAGCCAATGAGTTGTTGTCTCTRAATGTGTG
   s mm4.chr6     53303881 6 + 151104725 TAAAGAGATGCTAAGCCAATGAGTTGTTGTCGCTCAATGTGTG
   s rn3.chr4     81444246 6 + 187371129 taaggaGATGCTAAGCCAATGAGTTGTTGTCGCTCAATGTGTG

   …
   …
   …

Comment lines (starting with “#”) are ignored by mrbait_. Alignment blocks (considered
by mrbait_ to each represent different loci) are started with “a”, followed by sequence
lines starting with “s”. Source, strand, and coordinate positions are not informative for
mrbait_, nor are lines starting with other letters (which can be used in the MAF
format to communicate additional information about the preceding sequence, such as
quality scores).

The eXtended Multi-FastA (XMFA) format output by the multiple-genome aligner
MAUVE (which outputs it as “.alignment”) is an extension of the standard FASTA format
to allow alignment blocks from many different loci, with header lines representing
identifiers for the aligned sequence, and start-end coordinates representing the alignment
block location within the genome, followed by the sequence:

.. code-block:: none
   :linenos:

   >1:1-230 +
   ATAGC-NAATC--GC…
   >2:210-440 -
   ATTGGCCAATCCCC…
   >3:3-230 +
   TTA-CCAAGC--GC…
   =
   …
   …

Alignment blocks are delimited by the “=” symbol. All alignment blocks are assumed
by mrbait_ to represent separate, discontinuous loci. Note that because
no individual 'alignment block' in the .xmfa file is guaranteed to contain the same
genome representatives, no reference coordinates are saved by mrbait. This means
that additional annotation via GFF or VCF cannot be added to whole-genome alignments
provided in .xmfa format.

Reduced representation data
---------------------------

Alignments from reduced-representation methods such as restriction-site associate
DNA sequencing methods (RADseq) can be input using the MAF or XMFA formats, or using
the “.loci” format output by the RADseq assembly pipeline `pyrad <https://github.com/dereneaton/pyrad>`_ or
its successor `ipyrad <https://github.com/dereneaton/ipyrad>`_. This format shows individual
loci delimited by a line starting with “//” which features additional annotation of
variants and parsimony-informative sites:

.. code-block:: none
   :linenos:

   >PopA001	GTGTGATAGTAGTGATGTATTTTATAATATATATTATCGGATAT……
   >PopA002 	GTGTGARAGTAGTGATGTATTTTATAATATATATTATCGGATAT……
   >PopB001 	GTGTGACAGTAGTGATGTATTTTATAATATATATTATCGGATAT……
   >PopB002 	GAGTGATAGTAGTGATGTATTTTATAATATATATTATCGGATAT……
   //            *      *                                        |1
   …
   …
   …

mrbait ignores annotation information (since it parses variants anyways to generate
a consensus sequence), and only uses the “//” delimiter to distinguish between alignment
blocks. Creating a .loci file from other formats can be accomplished relatively easily.
For example, a series of separate alignments (each as .fasta), could be converted to the
**.loci** format using the following **bash** command:

.. code-block:: bash

   for file in `ls example*.fasta`; do
     awk 'BEGIN{ORS=""}$1~/^\>/{print $01"\t";next}{print $0"\n"}' $file
     >> example.loci;
     echo "//" >> example.loci;
   done
