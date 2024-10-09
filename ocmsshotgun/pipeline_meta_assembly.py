# Load necessary modules
from ruffus import *
from ruffus.combinatorics import *
import re
import os
import sys
import glob
import shutil
import pandas as pd
import cgatcore.pipeline as P
import ocmsshotgun.modules.Utility as utility
import ocmsshotgun.modules.MetaAssembly as PMA

# Load pipeline parameters from the configuration file
PARAMS = P.get_parameters(["pipeline.yml"])

# Location of the input files (reads)
DATADIR = "data.dir"

print(f"Using DATADIR: {DATADIR}")
print(f"Pipeline parameters: {PARAMS}")

# Check the input files correspond
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
# Update DATADIR after pooling (if enabled)
###############################################################################

@follows(poolSamples)  
def updateDATADIR():
    global DATADIR
    if PARAMS['pool']['enable']:
        DATADIR = "input_pooled.dir"
        print(f"Pooling is enabled. Setting DATADIR to: {DATADIR}")
    else:
        print(f"Pooling is not enabled. Using original DATADIR: {DATADIR}")

# Add a task to use the updated DATADIR
@follows(updateDATADIR)
def downstreamTask():
    print(f"Running downstream tasks with DATADIR: {DATADIR}")

###############################################################################
# Read Processing
###############################################################################

@active_if(PARAMS['read_error_correction']['enable'])  # Only run if error correction is enabled
@follows(mkdir('processed_reads.dir'))
@transform(
    os.path.join(DATADIR, '*.fastq.1.gz'),
    regex(re.escape(DATADIR) + r'/(.+).fastq.1.gz'),
    r'processed_reads.dir/\1.fastq.1.gz')
def runReadProcessing(infile, outfile):
    '''Run BayesHammer read correction on pooled or unpooled samples using SPAdes, based on YAML configuration.'''

    cluster_options = PARAMS['spades']['error_correction_run_options']
    assembler = PMA.SpadesReadCorrection()

    # Build and run the BayesHammer command
    statement = assembler.build(infile, outfile, **PARAMS)
    print(f"BayesHammer statement: {statement}")
    P.run(statement)

    # Fetch processed reads
    assembler = PMA.fetchSpadesProcessedReads()
    statement = assembler(infile, outfile)
    P.run(statement)

def main(argv=None):
    if argv is None:
        argv=sys.argv
    P.main(argv) 

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
