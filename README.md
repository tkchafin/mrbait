# Mr.Bait 

[![version][version-badge]][CHANGELOG] [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Python program for universal design of targeted-enrichment probes from multiple input datatypes. Can be used to generate probe sets for ultraconserved elements, anchored enrichment, RAD-capture, etc. 

Citation: Chafin TK, Douglas MR, Douglas ME. 2018. MrBait: Universal identification and design of targeted-enrichment capture probes. Bioinformatics. https://doi.org/10.1093/bioinformatics/bty548

Documentation: https://mrbait.readthedocs.io/

Conda package: https://anaconda.org/tylerkchafin/mrbait

## Installation
I currently support native installation on Linux and Mac OS. MrBait is not natively supported on Windows, however will run using the Linux subsystem for Windows 10 (see: https://docs.microsoft.com/en-us/windows/wsl/install-win10). MrBait was tested using the Ubuntu Windows-10 subsystem. 


Regardless of OS, the easiest way to install is using the conda package manager (https://conda.io/docs/):
```conda install mrbait -c tylerkchafin -c bioconda -c conda-forge ```

If you are concerned about conflicting dependencies, you can also create a contained virtual conda environment for mrbait and associated dependencies: 
```
#conda call to create a blank python3.6 environment named "mrbait"
conda create -n mrbait python=3.6

#activate the mrbait conda environment.
source activate mrbait

#add some channels 
conda config --env --add channels bioconda --add channels conda-forge

#install mrbait
conda install -c tylerkchafin mrbait

#run mrbait:
mrbait -v example.VCF <etc>
...
...

#exit environment
source deactivate

#note that you can now reactivate your mrbait environment at any time by typing:
source activate mrbait 

#and then exit again just as before
source deactivate
```

[CHANGELOG]: ./CHANGELOG.md
[LICENSE]: ./LICENSE
[version-badge]: https://img.shields.io/badge/version-1.2.1-blue.svg
