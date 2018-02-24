#!/usr/bin/python

from setuptools import setup, find_packages
import glob
import re

def requires():
    """ gets packages from requirements.txt """
    with open('requirements.txt') as infile:
        return infile.read().splitlines()

## Auto-update version from git repo tags
INITFILE = "ipyrad/__init__.py"
CUR_VERSION = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                    open(INITFILE, "r").read(),
                    re.M).group(1)

setup(
    name="mrbait",
    version=CUR_VERSION,
    url="https://github.com/tkchafin/mrbait",
    author="Tyler K. Chafin",
    author_email="tkchafin@uark.edu",
    description="Universal design of target-enrichment capture probes from genomic data",
    packages=find_packages(),
    install_requires=requires(),
    #dependencies=dependency_links(),
    entry_points={
            'console_scripts': [
                'ipyrad = mrbait.mrbait:main',
            ],
    },
    license='GPL',
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
    ],
)
