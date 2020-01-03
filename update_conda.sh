#!/bin/bash

conda build conda.recipe
conda install mrbait -c tylerkchafin -c bioconda -c conda-forge --force-reinstall
