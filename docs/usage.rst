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

-r, --resume
   **Resume**: This flag is used to tell mrbait if you would like to resume work |br|
   following a particular step. Use this option in conjunction with the |br|
   **--db** flag to continue the pipeline if you would like to re-perform |br|
   filtering steps without needing to re-load and parse alignments

   Usage: |br|
   -r 1: Continue pipeline after Step 1 (loading alignments) |br|
   -r 2: Continue fter Step 2 (target discovery) |br|
   -r 3: Continue after Step 3 (target filtering) |br|
   -r 4: Continue after Step 4 (bait discovery) |br|
   For example, -r 4 will tell mrbait_ to re-do bait filtering and output

--db   **Database**: Use this with the --resume flag to specify a .sqlite |br|
   database file from which to start the pipeline.
-T, --threads   **Threads**: Number of threads to use with processes that run |br|
   in parallel. This will also be passed to vsearch_ and/or blast_ if those are |br|
   being called. [default=1]
-h, --help   **Help**: Exit and display the help menu

Input Options
~~~~~~~~~~~~~

-M, --maf
   **MAF input**: Use this to provide the path to the multiple alignment MAF file
-X, --xmfa
   **XMFA input**: As an alternative to the MAF file, you can provide the |br|
   .xmfa file output by the aligner Mauve.
-L, --loci
   **LOCI input**: Multiple alignments can also be provided using the |br|
   .loci file output by the RADseq assembly pipeline pyRAD.
-A, --assembly
   **FASTA input**: Genome assembly provided as FASTA
-V, --vcf
   **VCF input**: For use with --assembly: VCF file containing variant data
-G, --gff
   **GFF input**: For use with --assembly: GFF file containing feature data

Alignment filtering/ consensus options (use with -M, -X, -L)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-c, --cov
  **Coverage**: Minimum number of individuals/sequences per alignment, |br|
  for MAF, XMFA, or LOCI inputs [default=1]
-l, --len
  **Minimum length**: Minimum alignment length to attempt bait design [default=80]
-q, --tresh  **Bad base threshold**: Threshold proportion of gaps or N (ambiguous or |br|
  poor quality) characters to over-ride the consensus base. For example, *-q 0.2* |br|
  would be interpreted as 20% of bases at a nucleotide position must be an “N” |br|
  or gap character in order for that character to be represented as the consensus |br|
  base. [default=0.1]
-Q, --max_ambig  **Max bad bases**: Maximum allowable proportion of gap/N characters |br|
  allowed in a consensus sequence before it will be discarded. *-Q 0.5* means a |br|
  consensus sequence can be 50% N’s or gap characters (“-“) before being dropped |br|
  from consideration. [default=0.5]
-k, --mask  **Mask threshold**: Threshold proportion of masked characters per nucleotide |br|
  column to mask the consensus base call. For use when case represents masking |br|
  information (where lowercase = masked), as when using the *-xsmall* option in |br|
  `RepeatMasker <http://www.repeatmasker.org/>`_ to flag low-complexity or repetitive sequences. |br|
  Case will be retained in the consensus on a per-base basis according to this threshold. |br|
  [default=0.1]
-K, --max_mask  **Max masked bases**: Maximum allowable proportion of masked characters |br|
  allowed in a consensus sequence before it will be discarded. *-K 0.5* means a |br|
  consensus sequence can be 50% masked (lowercase) before being |br|
  dropped from consideration.[default=0.5].

  If lowercase bases do not contain masking information, set to *-K 1.0*


General Bait Design Options
~~~~~~~~~~~~~~~~~~~~~~~~~~~

-b, --bait  **Bait length**: This is the length of desired baits, and will be used for |br|
  bait design as well as the sliding window width for target region |br|
  discovery [default=80]
-w, --win_shift  **Sliding window shift distance**: Shift distance for sliding window used |br|
  to discover target regions. Generally, there should not be a reason to |br|
  alter this. If target discovery (step 2) is taking a very long time, adjusting |br|
  this may make it faster although it could result in more targets failing |br|
  filtering [default=1]
-v, --var_max  **Maximum SNPs per bait**: Maximum allowable variants allowed in a bait |br|
  sequence. These can be expanded in the final output as each possible |br|
  non-ambiguous bait sequence for synthesis. Use this when there are not enough |br|
  conserved regions to capture enough loci for your design. [default=0]
-n, --numN  **Maximum Ns per bait**: Maximum allowable ambiguous (N) bases allowed |br|
  per bait. This could be increased when there are too many poor quality bases |br|
  in your alignment to design a sufficient number of probes, although keep in mind |br|
  this will affect the specificity of your resulting probes. [default=0]
-g, --numG  **Maximum gaps per bait**: Maximum allowable gap characters allowed per bait. |br|
  If dealing with alignments containing many indels, it might be desirable to |br|
  allow a small number per bait sequence. These can be expanded in the final |br|
  output using the -x,--expand option, which will expand gap characters as |br|
  A, G, T, C, and absent. [default=0]

Target Region Options
~~~~~~~~~~~~~~~~~~~~~

-R, --mult_reg  **Multiple targets per locus**: By default, mrbait_ only chooses one target |br|
  region (e.g. conserved region for which baits could be designed) per locus/ |br|
  alignment. When multiple are discovered, they are ranked according to the |br|
  criterion selected with the *-S,--select_r* option. When *-R,--mult_reg* not |br|
  in use, only a single target region (and corresponding baits) is chosen per |br|
  alignment. [default=false]
-m, --min_mult  **Minimum length for multiple targets**: Specify this to set a minimum alignment |br|
  or locus length to allow multiple target regions to be selected. By default will |br|
  be set to the value of *-l,--len* (thus, when *-R,--mult_reg* is used, all loci |br|
  passing length filter will be allowed multiple targets).
-D, --dist_r
  **Distance between targets**: When *-R,--mult_reg* is in use, use this parameter |br|
  to specify the minimum distance between targets. When targets are in conflict |br|
  (e.g. they are less than *-D,--dist_r* bases apart), conflicts will be resolved |br|
  using the criterion set with *-S,--select_r*. [default=100]
-d, --flank_dist
  **Flanking distance for target filtering**: Distance from boundaries of target |br|
  region to parse for counting SNPs, ambiguities, gaps, etc when filtering |br|
  target regions (see *-S,--select_r* and *-F,--filter_r*) [default=500]. |br|
  Note that this value will tell mrbait_ to search ‘*d*’ bases to the left AND |br|
  right of each target region.

  Note that currently, the same *--flank_dist* value will be used for all filters.
-S, --select_r
  **Target selection criterion**: Method to resolve conflicts when targets are |br|
  too close together (e.g. when *-R* and *-D*), or when only choosing one target |br|
  per locus/alignment.

  Usage:
  *-S snp*: Select target with most SNPs within *d* bases |br|
  *-S bad*: Select target with least gaps/Ns within *d* bases |br|
  *-S cons*: Select target with least SNPs within *d* bases |br|
  *-S rand*: Randomly select a target [default] |br|

  Example: *-d 100 -S snp* to choose region with most SNPs within |br|
  100 flanking bases
-F, --filter_r
  **Target filtering criteria**: Method(s) used to filter all target regions. |br|
  Can be specified any number of times to use additional filtering criteria. |br|

  Usage: |br|
  *-F len=[x,y]*: Length between *x* (min) and *y* (max) |br|
  *-F gap=[x]*: Maximum of *x* indels in target region |br|
  *-F bad=[x]*: Maximum of *x* N characters in target region |br|
  *-F snp=[x,y]*: Between *x* (min) and *y* (max) SNPs w/in *d* |br|
  *-F mask=[x]*: Maximum of *x* N characters in target region |br|
  *-F gc=[x,y]*: G/C propotion between *x* (min) and *y* (max) |br|
  *-F rand=[x]*: Randomly retain *x* targets |br|
  *-F pw=[i,q]*: Pairwise alignment, removing when *i* percent identity over at least *q* proportion of the sequences |br|
  *-F blast_i=[i,q]*: Only retain BLAST hits with *i* percent identity over at least *q* query coverage |br|
  *-F blast_i=[i,q]*: Exclude BLAST hits with *i* percent identity over at least *q* query coverage |br|
  *-F gff=[type]*: Only retain targets within *d* bases of a GFF-annotated feature of type type. Only for use when *-A* and *-G* inputs provided. Use *-F gff=all* to target any type of annotated feature. |br|
  *-F gff_a=[alias]*: Only retain targets within *d* bases of a GFF-annotated feature of tagged with the Alias attribute matching alias. Only for use when *-A* and *-G* inputs provided. |br|

  Examples: |br|
  *-F snp=1,10 -d 100* to sample when 1-10 SNPs within 100 bases |br|
  *-F gc=0.2,0.8 -F rand=100* to keep 100 random targets w/ 20-80% GC |br|
  *-F mask=0.1* to remove targets with >10% masked bases |br|
  *-d 1000 -F gff=exon* to keep targets within 100 bases of an exon |br|

Bait Selection Options
~~~~~~~~~~~~~~~~~~~~~~

-s, --select_b  **Bait selection scheme**: Use this to specify the desired method |br|
  to design baits from passing target regions.

Usage: |br|
  *-s tile=[x]*: Tile baits over whole region, with *x* overlap |br|
  *-s center=[n,x]*: *n* centered baits with *x* overlap |br|
  *-s flank=[n,x]*: *n* terminal baits (each end) with *x* overlap
  Default behavior is to tile baits across all targets with 50% overlap

-f, --filter_b
  **Bait filtering criteria**: Method(s) used to filter baits. |br|
  Can be specified any number of times to use additional filtering criteria.

  Usage: |br|
  *-F mask=[x]*: Maximum of *x* N characters in target region
  *-F gc=[x,y]*: G/C propotion between *x* (min) and *y* (max)
  *-F rand=[x]*: Randomly retain *x* targets
  *-F pw=[i,q]*: Pairwise alignment, removing when *i* percent identity over at least *q* proportion of the sequences |br|
  *-F blast_i=[i,q]*: Only retain BLAST hits with *i* percent identity over at least *q* query coverage |br|
  *-F blast_i=[i,q]*: Exclude BLAST hits with *i* percent identity over at least *q* query coverage
