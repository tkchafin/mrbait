# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/).

## [Releases](https://github.com/tkchafin/mrbait/releases)

v1.1.6
- Bug fix in reporting for bait design
- Bug fix in parallel target selection
- Added --target_all option
- Added option to design x baits per target with n max overlap (-b calc)

v1.1.5
- Bug fix for parsing command line options (-v and -c not working properly)
- Bug fix for targets not being added when sliding window exceeds locus

v1.1.4 https://github.com/tkchafin/mrbait/releases/tag/v1.1.4
- Added --vcfALT option (changes behavior for building consensus from VCF)
- Fixed bug in VCF parsing when reference is N or gap
- Added --print_loc option
- Fixed issue where XMFA hidden 'chunk' files were not deleted

v1.1.3 https://github.com/tkchafin/mrbait/releases/tag/v1.1.3
- Moved docs to mrbait.readthedocs.io
- Bugfix with consensus calling when gap/N content high
- Bugfix for XMFA reading
- XMFA fie chunking and parallel parsing fixed and added
- Minor bug fix in argument parsing

v1.1.2 https://github.com/tkchafin/mrbait/releases/tag/v1.1.2
- Changed numbering of baits and targets in FASTA output
- Added now functional --print_tr option
- Other minor edits in response to reviewer criticisms

v1.1.0 https://github.com/tkchafin/mrbait/releases/tag/v1.1.0
- Now installable via conda
- Changes to directory structure
- Non-compatible with v1.0.0

v1.0.0 https://github.com/tkchafin/mrbait/releases/tag/v1.0.0
- Functional version, pre-packaging for conda and before benchmarking

Please see the [source code](https://github.com/tkchafin/mrbait) for the most up-to-date development version, and see the complete [development history](https://github.com/tkchafin/mrbait/commits/master) for all changes.
