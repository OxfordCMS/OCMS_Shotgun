'''
concatenate_fastq.py
====================

:Author: Sandi Yen
:Tags: Python

Overview
========

This script concatenates fastq files of paired-end reads into one fastq file. It's written as a pipeline so paired-end fastqs can be processed as a job.

Usage
=====

Pipeline that takes in all fastq.*.gz files in current directory and concatenates paired-end reads for each sample into one fastq file. Concatenated fastqs written to  concat_fastq.dir

Example::

    ocms concatenate_fastq make full


Configuration
-------------
No configuration required

Input files
-----------
Input files should be fastq.1.gz, fastq.2.gz


Requirements
------------


Pipeline output
===============
reads from fastq.1.gz and fastq.2.gz concatenated into a fastq.gz file


Glossary
========

..glossary::


Code
====

'''

import sys
import os
import re
import glob
from pathlib import Path
from ruffus import *
from cgatcore import pipeline as P 

# get all sequence files within directory to process
SEQUENCEFILES = ("*fastq.*gz")

SEQUENCEFILES_REGEX = regex(r"(\S+)\.(fastq.*gz)")

######################################################
######################################################
######################################################
# concatenate paired end fastq files

# produces a concat_fastq.dir which contains
# concatenated fastq files
######################################################
######################################################
######################################################

@follows(mkdir("concat_fastq.dir"))
@collate(SEQUENCEFILES, 
         SEQUENCEFILES_REGEX,
         r"concat_fastq.dir/\1.fastq.gz")

def concatenate_fastq(infiles, outfile):
    """Expects paired end fastq files to be in format fastq.1.gz, fastq.2.gz
    """

    infiles = ' '.join(infiles)
    tempfile = os.path.splitext(outfile)[0]
    # concatenate files together
    statement = '''zcat %(infiles)s >> %(tempfile)s && gzip %(tempfile)s'''
    P.run(statement)

@follows(concatenate_fastq)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
