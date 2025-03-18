'''
pipeline_prokka.py
====================

:Author: Holly Roach and Marcin Pekalski
:Tags: Python

Overview
========

This script uses prokka to predict open reading frames and generate protein assemblies.

Usage
=====

Script takes in all contigs.fasta files in input.dir and runs prokka on these files to return a protein assembly

Example::

    ocms_shotgun pipeline_prokka make full


Configuration
-------------
ocms_shotgun pipeline_prokka config

Input files
-----------
Input files should be contigs.fasta files.


Requirements
------------
prokka/1.14.5-gompi-2022b

Pipeline output
===============
prokka_output.dir containing outputs generated prokka, including protein assemblies.


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
FASTAFILES = ("input.dir/*contigs.fasta")
FASTAFILES_REGEX = regex(r"input\.dir\/(\S+?)_.+\.fasta")

PARAMS = P.get_parameters(['pipeline.yml'])

######################################################
######################################################
######################################################
# subsamples fastq files to specified depth

# produces a prokka_output.dir which contains
# predicted protein assembilies
######################################################
######################################################
######################################################

@follows(mkdir("prokka_output.dir"))
@transform(FASTAFILES, 
         FASTAFILES_REGEX,
         r"prokka_output.dir/\1")

def run_prokka(infile, outfile):
    """Expects single end to be fastq.gz and paired end fastq files to be 
       in format fastq.1.gz, fastq.2.gz
    """

    print(f"infile: {infile}")      # infile: input.dir/samples_pooled_corrected.megahit.contigs.fasta
    print(f"outfile: {outfile}")    # outfile: prokka_output.dir/samples_pooled_corrected.megahit/

    # capture sample id from output dir name
    sample_id = re.sub("prokka_output.dir/", "", outfile)   # Sample ID: samples_pooled_corrected.megahit/
    print(f"Sample ID: {sample_id}")   

    # create statment for running prokka
    # prokka --cpus 0 --prefix {wildcards.sample} --locustag {wildcards.sample} --outdir {output.dir} --force {input}

    statement = f'''prokka --cpus 0 --prefix {sample_id} --locustag {sample_id} --outdir {outfile} --force {infile}'''

    # create script for slurm job submission
    P.run(statement,
          job_threads = PARAMS['job_threads'],
          job_memory = PARAMS['job_memory'])

@follows(run_prokka)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
