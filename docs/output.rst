.. include:: global.rst

.. _output:

Output Files
============

Final output of baits will be formatted as FASTA and named $out_baits.fasta
(where $out is defined using the -o/--out flag). When the *-t/--print_tr* option is in
use, targets will also be output as $out_targets.fasta, with an additional field
in the header indicating if these targets passed or failed target selection and filtering.

By default, baits are reported with any ambiguity sequences included (e.g. as a
consensus sequence) like so:

.. code-block:: none
   :linenos:

   >Locus1_Target4_Bait1
   ATGTAATRAGGTATATG……
   >Locus1_Target4_Bait2
   TATGAATGTCGCGCGAT……
   …
   …
   …

 If using the *-x/--expand* option, ambiguities will be reported as all combinations, like so:

.. code-block:: none
   :linenos:

   >Locus2_Target4_Bait1.1
   ATGTAATAAGGTATATG……
   >Locus2_Target4_Bait1.1
   ATGTAATGAGGTATATG……
   >Locus1_Target4_Bait2.1
   TATGAATGTCGCGCGAT……
   …
   …
   …

 Baits can also be printed as reverse complement. For example, if the *--expand* option was
 specified, in addition to *--strand both*:

.. code-block:: none
   :linenos:

   >Locus2_Target4_Bait1.1
   ATGTAATAAGGTATATG……
   >Locus2_Target4_Bait1.1_revcomp
   TACATTATTCCATATAC……
   >Locus2_Target4_Bait1.1
   ATGTAATGAGGTATATG……
   >Locus2_Target4_Bait1.1_revcomp
   TACATTACTCCATATAC ……
   >Locus1_Target4_Bait2.1
   TATGAATGTCGCGCGAT……
   >Locus1_Target4_Bait2.1_revcomp
   ATACTTACAGCGCGCTA……
   …
   …

mrbait_ will also produce a .sqlite file (e.g. $out.sqlite) which can be used with
the *--resume* flag to restart the pipeline at different stages- for example to
re-perform bait filtering with different options. This stores the complete database,
including all consensus loci parsed from the alignment input files, all targets, and
all bait sequences (including those which failed filtering) and can be used independently
with your own SQLite queries.
