#!/usr/bin/python

from setuptools import setup, find_packages

def requires():
    """ gets packages from requirements.txt """
    with open('requirements.txt') as infile:
        return infile.read().splitlines()

setup(
    name="mrbait",
    version="1.0.1",
    url="https://github.com/tkchafin/mrbait",
    author="Tyler K. Chafin",
    author_email="tkchafin@uark.edu",
    description="Interactive assembly and analysis of RADseq data sets",
    packages=find_packages(),
    install_requires=requires(),
    entry_points={
            'console_scripts': [
                'mrbait = mrbait.mrbait:main'
            ],
    },
    license='GPL',
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
    ],
)
