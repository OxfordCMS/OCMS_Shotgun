# Load necessary modules
from ruffus import *
from ruffus.combinatorics import *
import re
import os
import sys
import glob
import shutil
import sqlite3 as s3
import pandas as pd
import cgatcore.pipeline as P
import modules.Utility as utility
import modules.MetaAssembly as PMA

###############################################################################
# General setup
###############################################################################

# Load pipeline parameters from the configuration file
PARAMS = P.get_parameters(["pipeline.yml"])

# Location of the input files (reads)
if PARAMS['general']['input'] == 0:
    DATADIR = '.'
elif PARAMS['general']['input'] == 1:
    DATADIR = './data.dir'
else:
    DATADIR = PARAMS['general']['input']
assert os.path.isdir(DATADIR), f'Input directory does not exist: {DATADIR}'

# print(f"Using DATADIR: {DATADIR}")
# print(f"Pipeline parameters: {PARAMS}")

# Check the input files correspond
FASTQ1s = utility.check_input(DATADIR)

# utility function
def connect():
    '''Connect to sqlite db for storing pipeline datasets'''

    dbh = s3.connect(PARAMS['database']['name'])

    return dbh


###############################################################################
# Pooling samples (optional)
###############################################################################
@follows(mkdir('01_input_pooled.dir'))
@collate(FASTQ1s,
         regex(PARAMS['preprocess']['pool_input_regex']),  # Use regex from YAML
         r"01_input_pooled.dir/" + PARAMS['preprocess']['pool_output_regex'])  # Output pattern from YAML
def poolSamples(infiles, out_fastq1):
    '''Pool samples based on the provided regular expression and handle paired reads.'''
    # print(f"Pooling samples from: {infiles} to {out_fastq1}")

    samples = [utility.matchReference(x) for x in infiles]

    # A list of the fastq1 files
    in_fastqs1 = [fq.fastq1 for fq in samples]
    
    # Check for paired read files
    in_fastqs2 = [fq.fastq2 for fq in samples]
    if any(in_fastqs2):
        # Sanity check... all files need to be paired for pooling.
        for fq in samples:
            fq2 = P.snip(fq.fastq1, fq.fq1_suffix) + fq.fq2_suffix
            assert os.path.exists(fq2),  f"Expected paired file {fq2} does not exist."
        out_fastq2 = P.snip(out_fastq1, samples[0].fq1_suffix) + samples[0].fq2_suffix
            
    # Check for singleton read files (not essential for all samples to have singletons)
    in_fastqs3 = [fq.fastq3 for fq in samples]
    if any(in_fastqs3):
        out_fastq3 = P.snip(out_fastq1, samples[0].fq1_suffix) + samples[0].fq3_suffix

    if len(in_fastqs1) == 1:
        # Create symlinks for single input
        os.symlink(os.path.abspath(in_fastqs1[0]), out_fastq1)
        if any(in_fastqs2):
            os.symlink(os.path.abspath(in_fastqs2[0]), out_fastq2)
        if any(in_fastqs3):
            os.symlink(os.path.abspath(in_fastqs3[0]), out_fastq3)
        
    else:
        # Concatenate multiple input files
        in_fastqs1 = ' '.join(in_fastqs1)
        statement = "cat %(in_fastqs1)s > %(out_fastq1)s"
        if any(in_fastqs2):
            in_fastqs2 = ' '.join(in_fastqs2)
            statement += " && cat %(in_fastqs2)s > %(out_fastq2)s"
        if any(in_fastqs3):
            in_fastqs3 = ' '.join(in_fastqs3)
            statement += " && cat %(in_fastqs3)s > %(out_fastq3)s"

    P.run(statement)
    
    
###############################################################################
# Read Processing
###############################################################################
@follows(mkdir('02_processed_reads.dir'))
@transform(poolSamples,
           regex('.+/(.+).fastq.1.gz'),
           r'02_processed_reads.dir/\1.fastq.1.gz')
def runReadProcessing(infile, outfile):
    '''Run BayesHammer read correction on pooled or unpooled samples
    using SPAdes, based on YAML configuration.
    '''

    if PARAMS['preprocess']['error_correction']:

        cluster_options = PARAMS['spades']['error_correction_run_options']
        assembler = PMA.SpadesReadCorrection()

        # Build and run the BayesHammer command
        statement = assembler.build(infile, outfile, **PARAMS)
        P.run(statement,
              job_threads=PARAMS['spades']['ec_threads'],
              job_memory=PARAMS['spades']['ec_memory'],
              job_options=PARAMS.get('spades_ec_job_options', ''))

        # Fetch processed reads
        assembler = PMA.fetchSpadesProcessedReads()
        statement = assembler(infile, outfile)
        P.run(statement)

    else:
        os.symlink(infile, outfile)

    
###############################################################################
# Run Assembly Tools
###############################################################################
@follows(mkdir('03_spades_assembly.dir'))
@transform(runReadProcessing,
           regex('.+/(.+).fast(q|a).1.gz'),
           r'03_spades_assembly.dir/\1.spades.contigs.fasta')
def assembleWithMetaSpades(infile, outfile):
    '''Run SPAdes --meta without read error correction'''
    assembler = PMA.runMetaSpades()
    statement = assembler.build(infile, outfile, **PARAMS)
    
    P.run(statement,
          job_threads=PARAMS['spades']['meta_threads'],
          job_memory=PARAMS['spades']['meta_memory'],
          job_options=PARAMS.get('spades_meta_job_options', ''))

    
@follows(mkdir('megahit_assembly.dir'))
@transform(runReadProcessing,
           regex('.+/(.+).fast(q|a).1.gz'),
           r'megahit_assembly.dir/\1.megahit.contigs.fasta') 
def assembleWithMegaHit(infile, outfile):
    '''Run MEGAHIT'''
    cluster_options = PARAMS['megahit_cluster_options']
    assembler = PMA.runMegaHit()
    statement = assembler.build(infile, outfile, **PARAMS)
    
    P.run(statement, job_options=cluster_options)


ASSEMBLY_TARGETS = []
assembly_targets = {'metaspades': (assembleWithMetaSpades,),
                    'megahit': (assembleWithMegaHit,)
                    }

for tool in PARAMS['general']['assemblers'].split(','):
    ASSEMBLY_TARGETS.extend(assembly_targets[tool])

@follows(*ASSEMBLY_TARGETS)
def assembleMetaGenome():
    '''Perform assembly using tools specified in config file'''
    pass


###############################################################################
    
def main(argv=None):
    if argv is None:
        argv=sys.argv
    P.main(argv) 

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
