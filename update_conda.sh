#!/bin/bash

conda build -c conda-forge -c bioconda -c tylerkchafin conda.recipe
conda install mrbait -c tylerkchafin -c bioconda -c conda-forge --force-reinstall
