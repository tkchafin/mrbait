.. mrbait documentation master file, created by
   sphinx-quickstart on Mon Jun 11 08:58:50 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. include:: global.rst

mrbait: Universal identification and design of targeted-enrichment capture probes
=================================================================================


mrbait_ is a software pipeline for identifying regions of interest in DNA
sequence data and designing probes to enrich them.

The motivation behind mrbait is ease and flexibility of use. As such, mrbait
allows a variety of input types and facilitates a diverse array of bait design
approaches, such as those targeting ultraconserved elements, RAD-capture methods,
or those targeting exons or other genomic elements. mrbait also enables fast and
efficient iterative design (e.g. to explore parameter settings) using native
Python parallelization and an SQL database back-end. In this documentation, you
can learn about the overall process employed by mrbait `Pipeline <pipeline>`_,
how to install mrbait for use on a personal desktop or remote workstation or HPC
`Getting Started <installation>`_, or

mrbait_ code is open-source and freely available at on `GitHub <https://github.com/tkchafin/mrbait>`_

Official releases can be found `here <https://github.com/tkchafin/mrbait/releases>`_



Having issues running or installing mrbait_? Contact me at tkchafin@uark.edu
or post an Issue on the `GitHub page <https://github.com/tkchafin/mrbait/issues>`_.

Citation: Chafin, et al. 2018. MrBait: Universal identification and design of
targeted-enrichment capture probes. Bioinformatics. Submitted

Software and documentation provided under the GNU Public License v3.0 and distributed
“as is” without warranty of any kind.

.. toctree::
   :maxdepth: 2

   introduction.rst
   pipeline.rst
   installation.rst



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
