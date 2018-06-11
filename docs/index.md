# MrBait v. 1.1.2

Pipeline for universal design of target-enrichment capture probes from genomic data

Tyler Chafin

Department of Biological Sciences 
University of Arkansas


February 24, 2018
(Last updated: April 12, 2018)


MrBait code is open-source and freely available at3:
https://github.com/tkchafin/mrbait


Official releases can be found at3: 
https://github.com/tkchafin/mrbait/releases



Chafin, et al. 2018. MrBait: Universal identification and design of targeted-enrichment capture probes. Bioinformatics. Submitted. 
Questions or comments can be directed to tkchafin@uark.edu 
Software and documentation provided under the GNU Public License v3.0 and distributed “as is” without warranty of any kind. 


# Contents

1.	Introduction 	3

2.	Pipeline description 	3

3.	Getting started	5
3.1 Availability 	5
3.2 Dependencies 	5
3.3 Installation 	6
3.4 Running MrBait 	6

4.	Input files	7
4.1 Assembled genomes 	7
4.2 Multiple genome alignments 	8
4.3 Reduced-representation data 	9

5.	Usage options	10
5.1 Main parameters 	10
5.2 Filtering using VSEARCH 	15
5.3 Filtering using BLAST 	17

6.	Output files 	18

7.	Benchmarking and hardware requirements 	19

8.	Acknowledgements 	21

9.	References 	21


## Overview

MkDocs is a **fast**, **simple** and **downright gorgeous** static site
generator that's geared towards building project documentation. Documentation
source files are written in Markdown, and configured with a single YAML
configuration file.

### Host anywhere

MkDocs builds completely static HTML sites that you can host on GitHub pages,
Amazon S3, or [anywhere][deploy] else you choose.

### Great themes available

There's a stack of good looking themes available for MkDocs. Choose between
the built in themes: [mkdocs] and [readthedocs], select one of the 3rd
party themes in the [MkDocs wiki], or [build your own].

### Preview your site as you work

The built-in dev-server allows you to preview your documentation as you're
writing it. It will even auto-reload and refresh your browser whenever you save
