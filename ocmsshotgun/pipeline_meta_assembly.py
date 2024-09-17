"""
=============================
Metagenomic assembly pipeline
=============================

:Author: Jethro Johnson
:Release: $Id$
:Date: |today|
:Tags: Python

This pipeline imports cleaned, trimmed reads from OCMC_preprocess pipelines  and creates a metagenomic assembly. Pooling of samples is optional. 

Configuration
-------------

The pipeline requires a configured :file:`pipeline_mts.yml` file.
"""

###############################################################################
#load modules
from ruffus import *
from ruffus.combinatorics import *
import sys
import os
import re
import glob
import shutil
import pickle
import sqlite3
import collections
import numpy as np
import pandas as pd
import cgatcore.experiment as E
import cgat.GTF as GTF
import cgatcore.iotools as IOTools
import cgat.BamTools as BamTools
import cgat.FastaIterator as FastaIterator
import cgatcore.pipeline as P
import ocmsshotgun.modules.Utility as utility
#import PipelineMetaAssembly as PMA

###############################################################################

# Set up parameters
PARAMS = P.get_parameters(["pipeline_mts.yml"])

# Determine the location of the input files (reads)
try:
    PARAMS["input"]
except NameError:
    DATADIR = "."
else:
    if PARAMS["input"] == 0:
        DATADIR = "."
    elif PARAMS["input"] == 1:
        DATADIR = "data.dir"
    else:
        DATADIR = PARAMS["input"]

print(f"Using DATADIR: {DATADIR}")
print(f"Pipeline parameters: {PARAMS}")

###############################################################################

# Check that input files correspond
FASTQ1S = utility.check_input(DATADIR)

# Check if pooling is enabled in the config
if PARAMS['pool']['enable']:

    @follows(mkdir('input_pooled.dir'))
    @collate(os.path.join(DATADIR, '*.fastq.1.gz'),
             regex(PARAMS['preprocess']['pool_input_regex']),  # Using regex from YAML
             r"input_pooled.dir/" + PARAMS['preprocess']['pool_output_regex'])  # Using output pattern from YAML
    def poolSamples(infiles, out_fastq1):
        '''Collapse samples based on supplied regular expression and check for the presence of fastq3 files.'''

        print(f"Pooling samples from: {infiles} to {out_fastq1}")
        print("Checking directory existence...")

        if not os.path.exists('input_pooled.dir'):
            print("Error: Directory input_pooled.dir was not created.")
            return

        # Define the corresponding fastq.2.gz and fastq.3.gz output filenames
        out_fastq2 = P.snip(out_fastq1, '.1.gz') + '.2.gz'
        out_fastq3 = P.snip(out_fastq1, '.1.gz') + '.3.gz'

        in_fastqs1 = infiles
        in_fastqs2 = [P.snip(x, '.1.gz') + '.2.gz' for x in in_fastqs1]
        in_fastqs3 = [P.snip(x, '.1.gz') + '.3.gz' for x in in_fastqs1]

        # Ensure all paired files (fastq.2.gz) exist before proceeding
        for f in in_fastqs2:
            assert os.path.exists(f), f"Expected paired file {f} does not exist."

        # Debugging: Print paths being checked for .fastq.3.gz files
        for f in in_fastqs3:
            print(f"Checking existence of: {f}")

        # Check if any .fastq.3.gz files exist for this particular group
        include_fastq3 = [os.path.exists(f) for f in in_fastqs3]
        print(f"FASTQ3 files presence: {include_fastq3}")

        if len(in_fastqs1) == 1:
            # Create symlinks for a single input file
            os.symlink(os.path.join('..', in_fastqs1[0]), out_fastq1)
            os.symlink(os.path.join('..', in_fastqs2[0]), out_fastq2)
            if include_fastq3[0]:
                os.symlink(os.path.join('..', in_fastqs3[0]), out_fastq3)
        else:
            # Concatenate multiple input files
            in_fastqs1_str = ' '.join(in_fastqs1)
            in_fastqs2_str = ' '.join(in_fastqs2)

            statement = (
                f"cat {in_fastqs1_str} > {out_fastq1} && "
                f"cat {in_fastqs2_str} > {out_fastq2}"
            )

            if any(include_fastq3):
                in_fastqs3_str = ' '.join([f for i, f in enumerate(in_fastqs3) if include_fastq3[i]])
                if in_fastqs3_str:
                    statement += f" && cat {in_fastqs3_str} > {out_fastq3}"

            print(f"Running statement: {statement}")

            # Execute the concatenation command with specified resources
            P.run(statement,
                  job_threads=PARAMS['pool']['job_threads'],
                  job_memory=PARAMS['pool']['job_memory'],
                  job_options=PARAMS['pool'].get('options', ''))

            print(glob.glob('input_pooled.dir/*.fastq.1.gz'))

####################################################################################
####METAGENOME ASSEMBLY
###############################################################################
#This section of the pipeline runs different metagenome assemblers. MetaSPAdes includes read correction by default using the BayesHammer tool, which corrects errors in short reads before assembly. Megahit, on the other hand, incorporates an advanced algorithm that handles read errors during the assembly process itself, so a separate read error correction step is not required in this pipeline when using Megahit

@follows(poolSamples if PARAMS['pool']['enable'] else None)
@transform(os.path.join('input_pooled.dir/*.fastq.1.gz') if PARAMS['pool']['enable'] else os.path.join(DATADIR, '*.fastq.1.gz'),
           regex(r"(.+)_pooled.fastq.1.gz"),  # Adjusted to match 'pooled' files
           add_inputs(r"\1_pooled.fastq.2.gz"),
           r"assembly/\1_metaspades")
def run_metaSPAdes(infiles, outdir):
    '''Run MetaSPAdes on paired-end reads and handle fastq.3.gz as singletons if present.'''

    in_fastq1, in_fastq2 = infiles
    outdir = os.path.join(outdir, "metaspades")

    # Check if the fastq.3.gz file exists for this sample
    in_fastq3 = P.snip(in_fastq1, '.1.gz') + '.3.gz'
    include_fastq3 = os.path.exists(in_fastq3)

    # Base command for running MetaSPAdes with paired-end reads
    statement = (
        f"spades.py --meta -1 {in_fastq1} -2 {in_fastq2} -o {outdir} "
        f"-t {PARAMS['spades']['job_threads']} "
        f"-m {PARAMS['spades']['job_memory']} "
    ) 

    # If fastq.3.gz exists, treat it as singleton reads using the -s option
    if include_fastq3:
        statement += f"-s {in_fastq3} "

    # Add extra options from YAML if provided
    if PARAMS['pool'].get('options', ''):
        statement += PARAMS['pool']['options']

    print(f"Running MetaSPAdes: {statement}")

    # Execute the assembly command
    P.run(statement, 
          job_threads=PARAMS['pool']['job_threads'], 
          job_memory=PARAMS['pool']['job_memory'])

# Run the pipeline
pipeline_run()

