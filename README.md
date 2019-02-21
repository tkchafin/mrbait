# Mr.Bait 

[![version][version-badge]][CHANGELOG] [![license][license-badge]][LICENSE]

Python program for universal design of targeted-enrichment probes from multiple input datatypes. Can be used to generate probe sets for ultraconserved elements, anchored enrichment, RAD-capture, etc. 

Citation: Chafin TK, Douglas MR, Douglas ME. 2018. MrBait: Universalidentification and design of targeted-enrichment capture probes. Bioinformatics. https://doi.org/10.1093/bioinformatics/bty548

Documentation: https://mrbait.readthedocs.io/

Conda package: https://anaconda.org/tylerkchafin/mrbait

## Installation
I currently supports native installation on Linux and Mac OS. MrBait is not natively supported on Windows, however will run using the Linux subsystem for Windows 10 (see: https://docs.microsoft.com/en-us/windows/wsl/install-win10). MrBait was tested using the Ubuntu Windows-10 subsystem. 


Regardless of OS, the easiest way to install is using the conda package manager (https://conda.io/docs/):
```conda install mrbait -c tylerkchafin -c bioconda -c conda-forge ```


[CHANGELOG]: ./CHANGELOG.md
[LICENSE]: ./LICENSE
[version-badge]: https://img.shields.io/badge/version-1.1.6-blue.svg
[license-badge]: 	https://img.shields.io/aur/license/yaourt.svg
