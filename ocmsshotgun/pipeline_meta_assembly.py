"""
=============================
Metagenomic Assembly Pipeline
=============================

:Author: Jethro Johnson
:Release: $Id$
:Date: |today|
:Tags: Python

This pipeline imports cleaned, trimmed reads from OCMC_preprocess pipelines and creates a metagenomic assembly. Pooling of samples is optional.

Configuration
-------------
The pipeline requires a configured :file:`pipeline_mts.yml` file.
"""

# Load necessary modules
from ruffus import *
from ruffus.combinatorics import *
import os
import glob
import shutil
import pandas as pd
import cgatcore.pipeline as P
import ocmsshotgun.modules.Utility as utility
import PipelineMetaAssembly as PMA

# Load pipeline parameters from the configuration file
PARAMS = P.get_parameters(["pipeline_mts.yml"])

# Determine the location of the input files (reads)
if PARAMS["input"] == 0:
    DATADIR = "."
elif PARAMS["input"] == 1:
    DATADIR = "data.dir"
else:
    DATADIR = PARAMS["input"]

print(f"Using DATADIR: {DATADIR}")
print(f"Pipeline parameters: {PARAMS}")

# Verify input FASTQ files
FASTQ1S = utility.check_input(DATADIR)

###############################################################################
# Pooling samples (optional)
###############################################################################

@active_if(PARAMS['pool']['enable'])  # Only run if pooling is enabled
@follows(mkdir('input_pooled.dir'))
@collate(os.path.join(DATADIR, '*.fastq.1.gz'),
         regex(PARAMS['preprocess']['pool_input_regex']),  # Use regex from YAML
         r"input_pooled.dir/" + PARAMS['preprocess']['pool_output_regex'])  # Output pattern from YAML
def poolSamples(infiles, out_fastq1):
    '''Pool samples based on the provided regular expression and handle paired reads.'''
    print(f"Pooling samples from: {infiles} to {out_fastq1}")

    out_fastq2 = P.snip(out_fastq1, '.1.gz') + '.2.gz'
    out_fastq3 = P.snip(out_fastq1, '.1.gz') + '.3.gz'

    in_fastqs2 = [P.snip(x, '.1.gz') + '.2.gz' for x in infiles]
    in_fastqs3 = [P.snip(x, '.1.gz') + '.3.gz' for x in infiles]

    # Ensure paired reads exist
    for f in in_fastqs2:
        assert os.path.exists(f), f"Expected paired file {f} does not exist."

    # Handle third paired reads if they exist
    include_fastq3 = [os.path.exists(f) for f in in_fastqs3]

    if len(infiles) == 1:
        # Create symlinks for single input
        os.symlink(os.path.join('..', infiles[0]), out_fastq1)
        os.symlink(os.path.join('..', in_fastqs2[0]), out_fastq2)
        if include_fastq3[0]:
            os.symlink(os.path.join('..', in_fastqs3[0]), out_fastq3)
    else:
        # Concatenate multiple input files
        in_fastqs1_str = ' '.join(infiles)
        in_fastqs2_str = ' '.join(in_fastqs2)

        statement = f"cat {in_fastqs1_str} > {out_fastq1} && cat {in_fastqs2_str} > {out_fastq2}"

        if any(include_fastq3):
            in_fastqs3_str = ' '.join([f for i, f in enumerate(in_fastqs3) if include_fastq3[i]])
            if in_fastqs3_str:
                statement += f" && cat {in_fastqs3_str} > {out_fastq3}"

        print(f"Running statement: {statement}")
        P.run(statement, job_threads=PARAMS['pool']['job_threads'], job_memory=PARAMS['pool']['job_memory'])

###############################################################################
# Read Correction with BayesHammer (SPAdes)
###############################################################################

@active_if(PARAMS['pool']['enable'])
@follows(poolSamples, mkdir('spades_read_correction.dir'))
@transform(
    os.path.join('input_pooled.dir', '*.fastq.1.gz'),
    regex(r'input_pooled.dir/(.+)_pooled.fastq.1.gz'),
    r'spades_read_correction.dir/\1.fastq.1.gz')
def runReadCorrection(infile, outfile):
    '''Run BayesHammer read correction on pooled samples using SPAdes.'''

    cluster_options = PARAMS['spades']['error_correction_run_options']
    assembler = PMA.SpadesReadCorrection()

    # Build and run the BayesHammer command
    statement = assembler.build(infile, outfile, **PARAMS)
    P.run(statement)

    # Fetch processed reads
    assembler = PMA.fetchSpadesProcessedReads()
    statement = assembler(infile, outfile)
    P.run(statement)

###############################################################################
# SPAdes Metagenome Assembly
###############################################################################

@follows(mkdir('spades_assembly.dir'))
@transform(runReadCorrection,
           regex(r'.+/(.+).fastq.1.gz'),
           r'spades_assembly.dir/\1.spades.contigs.fasta')
def assembleWithMetaSpades(infile, outfile):
    '''Run MetaSPAdes on corrected reads for metagenomic assembly.'''

    cluster_options = PARAMS['spades_run_options']
    assembler = PMA.runMetaSpades()

    # Build and run the SPAdes assembly command
    statement = assembler.build(infile, outfile, **PARAMS)
    P.run(statement)

###############################################################################
# Run the pipeline
###############################################################################
if __name__ == "__main__":
    pipeline_run()


