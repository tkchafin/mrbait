#!/usr/bin/python

from setuptools import setup, find_packages

def requires():
    """ gets packages from requirements.txt """
    with open('requirements.txt') as infile:
        return infile.read().splitlines()

setup(
    name="mrbait",
    install_requires=requires(),
    version="1.0.4",
    description="Universal design of target-enrichment capture probes from genomic data",
    author='Tyler K. Chafin',
    author_email='tkchafin@uark.edu',
    url='https://github.com/tkchafin/mrbait',
	packages=['mrbait'],
    include_package_data=True,
    entry_points={'console_scripts': ['mrbait = mrbait.mrbait:main']},
    zip_safe=False,
    classifiers=[
        'Programming Language :: Python :: 3.6',
    ]
)
