'''
pipeline_concatfastq.py
====================

:Author: Sandi Yen
:Tags: Python

Overview
========

This pipeline concatenates fastq files of paired-end reads into one fastq file. It's written as a pipeline so paired-end fastqs can be processed as a job.

Usage
=====

Pipeline that takes in all fastq.*.gz files in current directory and concatenates paired-end reads for each sample into one fastq file. Concatenated fastqs written to  concat_fastq.dir

Example::

    ocms_shotgun concatfastq make full


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
from ruffus import *
from cgatcore import pipeline as P 
import ocmsshotgun.modules.Humann3 as H
import ocmstoolkit.modules.Utility as Utility
# no params in config needed but putting here in case that changes
PARAMS = P.get_parameters(["pipeline.yml"])

# set input directory
indir = PARAMS.get('general_input.dir', 'input.dir')

# get all sequence files within directory to process
SEQUENCEFILES = Utility.get_fastns(indir)

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
@transform(SEQUENCEFILES, 
         regex(fr"{indir}/(\S+)\.(fastq.1.gz)"),
         r"concat_fastq.dir/\1.fastq.gz")
def concatFastq(infile, outfile):
    """Expects paired end fastq files to be in format fastq.1.gz, fastq.2.gz
    """
    tool = H.Humann3(infile, outfile, **PARAMS)
    statement = tool.concat_fastq(Utility.MetaFastn(infile))
    P.run(statement)

@follows(concatFastq)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
