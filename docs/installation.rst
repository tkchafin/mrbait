.. include:: global.rst

.. _installation:

Getting Started
===============
mrbait_ has been tested on Mac and Linux operating systems and is primarily supported on
those platforms. However, Windows users can easily install using the built-in
Linux subsystem for Windows 10.

In-development code can be found on the Github page: https://github.com/tkchafin/mrbait

If you find any issues with the program, please email me at tkchafin@uark.edu or
submit as an issue on `Github <https://github.com/tkchafin/mrbait/issues>`_, which can
also be used for submitting feature requests. When submitting bugs or issues, please
include input files, your command-line call, and any output MrBait produced to the screen
or output files.

Availability
------------
Functioning releases can be found at:
https://github.com/tkchafin/mrbait/releases

Dependencies
------------
mrbait_ is written for Python3, and requires Python version >= 3.6.0. The recommended
method of acquiring Python and all other dependencies is via the Anaconda distribution,
as outlined in Section 3.3. A full list of dependencies is given below.

Python_ >= 3.6
SQLite3_
BioPython_
Pandas_ >=0.22
numpy_
pyVCF_
networkx_

mrbait can optionally use the following programs during bait development:

blast_
vsearch_

For these utilities, please cite the following:
Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL. 2009.
BLAST+: architecture and applications. BMC Bioinformatics. 1-:410. Doi:10.1186/1471-2105-10.421

Rognes R, Flouri T, Nichols B, Quince C, Mahe F. 2016. VSEARCH: A versatile and open
source tool for metagenomics. PeerJ. 4:e2584. Doi: 10.7717/peerj.2584

Installation
------------
By far the easiest way to acquire and install mrbait_ is via conda_, a command line interface
for managing and installing packages. Download and install anaconda_ for Python 3.6 here:
https://www.anaconda.com/download/. If you are wanting a minimal environment, or a faster
install, you can also use the Miniconda distribution (https://conda.io/miniconda.html) with
the same commands. After installation, be sure to test that conda is installed by typing
**conda info**, which will print information about your installation. Note, you may first need
to reload your bash environment by typing **source ~/.bashrc** or **source ~/.bash_profile** on Mac.
Assuming success, the installation process is then very straightforward:

.. code-block:: bash

   #This command tells conda that the code and dependencies for mrbait can
   #be found in ‘channels’ bioconda, conda-forge, and tylerkchafin.
   conda install mrbait -c tylerkchafin -c bioconda -c conda-forge

   #If you would like to instead install the latest development version, you can
   #clone the github repository and
   #install MrBait like so (assuming you have git installed):
   git clone https://github.com/tkchafin/mrbait.git
   cd mrbait
   python ./setup.py install


You will then need to manually install both vsearch_ and blast_, only if you install
directly from the GitHub source using the setup.py installation. These will be installed
for you if you used conda_.

**Windows users**: MrBait is installable using the built-in Linux subsystem for Windows 10.
I have only tested using the Ubuntu OS subsystem configuration but assume that other Linux
distros would work equally well. If you prefer, you can also use a Linux installation on a
virtual machine, or installed portably on a `USB-attached drive <https://tutorials.ubuntu.com/tutorial/tutorial-create-a-usb-stick-on-ubuntu#0>`_, although
this may impact performance. Contact me at tkchafin@uark.edu if you have any issues getting
mrbait_ installed, or feel free to launch an ‘Issue’ on the GitHub page.

**HPC users**: One of the reasons I recommend using conda to manage your Python environment,
is that it keeps your packages separate from the system environment, which you often will not
have permissions to modify. Anaconda will instead install your own local flavor of Python in
your home directory, where is will also install any additional packages you choose to add.

**BLAST and VSEARCH**: conda will also install both BLAST and VSEARCH and place them within
your conda environment. If you would like to manually manage versions of these programs, or
use an existing installation, you can provide the paths to those binaries using the *--vsearch*
and *--blastn* commands for mrbait_.

Running mrbait_
---------------
Assuming you have completed the recommended conda_ install, mrbait_ and it’s `Dependencies`_
should already be in your path and is now fully ready to go. You can verify successful
installation, and view the help menu, by typing:
mrbait -h

Instructions for bait design are provided as arguments (see Section 5 for thorough usage
instructions, and Section 8 for tutorials). For example, to generate baits of length
80, tiled across target regions with an overlap of 40 bases, from a Multiple Alignment
File (MAF) “example.maf”:

.. code-block:: bash

  mrbait -M example.maf -b 80 -s tile=40

Or, to also filter for only alignments including 5 or more individuals, and
of length >500:

.. code-block:: bash

  mrbait -M example.maf -b 80 -s tile=40 -l 500 -c 5
