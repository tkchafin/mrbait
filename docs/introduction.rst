.. include:: global.rst

A variety of genome reduction methods have been implemented to reduce costs of
applying next-generation sequencing methods to non-model organisms, or projects
with large numbers of samples (e.g. those focusing on the population scale). These
can broadly be classified into those which use restriction enzymes and size selection
for subsampling genomic complexity [e.g. RADseq methods (Baird et al., 2008;
Peterson et al., 2012)], and those which enrich for fragments selected a priori
using biotinylated RNA ‘baits’ (Lemmon et al., 2012; McCormack et al., 2012). The
latter benefit from increased specificity, yet require some genomic information
for marker development. To mitigate, some take a hybrid approach by using baits to
enrich RAD loci which are most consistently recovered, or to maximize capture of
parsimony-informative variation (Ali et al., 2016; Hoffberg et al., 2016).

Applying these targeted-enrichment methods (either via RAD-capture methods,
ultra-conserved elements, or anchored-enrichment) requires first bioinformatic
processing to parse large alignments, identify candidate regions for bait-design,
and design of complementary oligonucleotide sequences for synthesis. Although some
developers are transparent in providing computational resources and workflows to
design such probe sets (e.g. see Faircloth 2017), a generalized and flexible pipeline
does not yet exist. The motivation behind MrBait was to provide such a resource,
which could be universally applied to differing bait enrichment strategies (e.g.
targeting ultra-conserved regions vs. functional elements), and facilitate diverse
quality control methods to mitigate non-target capture (contamination, etc),
target-target hybridization, ambiguous mapping, and enrichment of repetitive DNA.
