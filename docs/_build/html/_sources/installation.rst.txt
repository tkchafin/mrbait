.. include:: global.rst

.. _installation:

Getting Started
===============
MrBait has been tested on Mac and Linux operating systems and is primarily supported on
those platforms. However, Windows users can easily install using the built-in
Linux subsystem for Windows 10.

In-development code can be found on the Github page: https://github.com/tkchafin/mrbait

If you find any issues with the program, please email me at tkchafin@uark.edu or
submit as an ‘issue’ on Github: https://github.com/tkchafin/mrbait/issues, which can
also be used for submitting feature requests. When submitting bugs or issues, please
include input files, your command-line call, and any output MrBait produced to the screen
or output files.

Availability
------------
Functioning releases can be found at:
https://github.com/tkchafin/mrbait/releases

Dependencies
------------
MrBait is written for Python3, and requires Python version >= 3.6.0. The recommended
method of acquiring Python and all other dependencies is via the Anaconda distribution,
as outlined in Section 3.3. A full list of dependencies is given below.

Python >= 3.6: https://www.python.org/
SQLite3: https://www.sqlite.org
BioPython: http://biopython.org/
Pandas >=0.22: https://pandas.pydata.org/
NumPy: http://www.numpy.org/
PyVcf: https://pyvcf.readthedocs.io
NetworkX: https://networkx.github.io/

MrBait can optionally use the following programs during bait development:

NCBI BLAST+: https://blast.ncbi.nlm.nih.gov/Blast.cgi
VSEARCH: https://github.com/torognes/vsearch

For these utilities, please cite the following:
Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL. 2009.
BLAST+: architecture and applications. BMC Bioinformatics. 1-:410. Doi:10.1186/1471-2105-10.421
-and-
Rognes R, Flouri T, Nichols B, Quince C, Mahe F. 2016. VSEARCH: A versatile and open
source tool for metagenomics. PeerJ. 4:e2584. Doi: 10.7717/peerj.2584

Installation
------------
By far the easiest way to acquire and install MrBait is via _conda, a command line interface 
for managing and installing packages. Download and install Anaconda for Python 3.6 here:
https://www.anaconda.com/download/. If you are wanting a minimal environment, or a faster
install, you can also use the Miniconda distribution (https://conda.io/miniconda.html) with
the same commands. After installation, be sure to test that conda is installed by typing
‘conda info’, which will print information about your installation. Note, you may first need
to reload your bash environment by typing **source ~/.bashrc** or **source ~/.bash_profile** on Mac.
Assuming success, the installation process is then very straightforward:
