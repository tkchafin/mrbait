#!/usr/bin/python

from setuptools import setup, find_packages

def requires():
    """ gets packages from requirements.txt """
    with open('requirements.txt') as infile:
        return infile.read().splitlines()

setup(
    name="mrbait",
    packages=['mrbait'],
    version="1.0.2",
    description="Universal design of target-enrichment capture probes from genomic data",
    author='Tyler K. Chafin',
    author_email='tkchafin@uark.edu',
    url='https://github.com/tkchafin/mrbait',
    install_requires=requires(),
	packages=find_packages(),
    include_package_data=True,
    entry_points={'console_scripts': ['mrbait = mrbait.mrbait:main']},
    zip_safe=False,
    classifiers=[
        'Programming Language :: Python :: 3.6',
    ]
)
