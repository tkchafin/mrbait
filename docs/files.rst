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
