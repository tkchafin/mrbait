#!/usr/bin/python

from subprocess import Popen, PIPE, CalledProcessError
import os

"""Includes utilities for calling VSEARCH and parsing output of pairwise alignments"""

def allpairsGlobal(binary, threads, seqpath, qid, qcov, outpath):
    vsearch = [binary,
            "--allpairs_global", seqpath,
            "--threads", str(threads),
            "--id", str(qid),
            "--blast6out", ".temp.pw",
            "--rowlen", "0", "--self",
            "--target_cov", str(qcov),
            "--quiet"]
    command = " ".join(vsearch)

    #Vsearch subprocess
    proc = Popen(vsearch, stdout=PIPE, stdin=PIPE, env={'PATH': os.getenv('PATH')})

    #WRap to enable keyboard interrupe
    try:
        t = proc.communicate()[0]
    except KeyboardInterrupt:
        proc.kill()
        raise KeyboardInterrupt

    #Get return code from process
    if proc.returncode:
        raise CalledProcessError ("VSEARCH exited with non-zero status")

def sortByLength(binary, seqpath, outpath):
    vsearch = [binary,
            "--sortbylength", seqpath,
            "--output", outpath]
    command = " ".join(vsearch)

    #Vsearch subprocess
    proc = Popen(vsearch, stdout=PIPE, stdin=PIPE, env={'PATH': os.getenv('PATH')})

    #WRap to enable keyboard interrupe
    try:
        t = proc.communicate()[0]
    except KeyboardInterrupt:
        proc.kill()
        raise KeyboardInterrupt

    #Get return code from process
    if proc.returncode:
        raise CalledProcessError ("VSEARCH exited with non-zero status")
