.. include:: global.rst

.. _usage:

Usage options
=============

mrbait_ reads all options and inputs using command-line arguments provided
after the program name. For a quick look at all options from the command line,
call the help menu by typing **mrbait -h** from the terminal.

Note that options requiring a floating point number (e.g. *-q*) allow inputs
from 0.0 to 1.0, and options requiring an integer (e.g. *-c*) allow inputs ranging
from 1 to infinity.

Main Parameters
---------------
General options
~~~~~~~~~~~~~~~

.. |br| raw:: html

    <br>

-r, --resume
   *Resume*: This flag is used to tell mrbait if you would like to resume work |br|
   following a particular step. Use this option in conjunction with the |br|
   *--db* flag to continue the pipeline if you would like to re-perform |br|
   filtering steps without needing to re-load and parse alignments|br||br|
   Usage: |br|
   -r 1: Continue pipeline after Step 1 (loading alignments) |br|
   -r 2: Continue fter Step 2 (target discovery) |br|
   -r 3: Continue after Step 3 (target filtering) |br|
   -r 4: Continue after Step 4 (bait discovery) |br|
   For example, -r 4 will tell mrbait to re-do bait filtering and output

--db   *Database*: Use this with the --resume flag to specify a .sqlite|br|
   database file from which to start the pipeline.
-T, --threads   *Threads*: Number of threads to use with processes that run |br|
   in parallel. This will also be passed to vsearch_ and/or blast_ if those are|br|
   being called. [default=1]
-h, --help   *Help*: Exit and display the help menu

Input Options
~~~~~~~~~~~~~

-M,--maf   *MAF input*: Use this to provide the path to the multiple alignment MAF file
-X,--xmfa   *XMFA input*: As an alternative to the MAF file, you can provide the |br| .xmfa file output by the aligner Mauve.
-L,--loci   *LOCI input*: Multiple alignments can also be provided using the |br|
.loci file output by the RADseq assembly pipeline pyRAD.
-A,--assembly   *FASTA input*: Genome assembly provided as FASTA
-V,--vcf   *VCF input*: For use with --assembly: VCF file containing variant data
-G,--gff   *GFF input*: For use with --assembly: GFF file containing feature data
