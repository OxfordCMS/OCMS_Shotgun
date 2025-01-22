f
om ruffus import *
from cgatcore import pipeline as P
import sys
import ocmsshotgun.modules.MAG as mb
import ocmsshotgun.modules.Utility as utility
import os
import re
import itertools

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
i

# Check the input files correspond
FASTQ1s = utility.check_input(DATADIR)

# Check for FASTA files (pooled or unpooled)
def check_fasta_files(datadir=DATADIR, fasta_suffix='fasta'):
    """
    Check for FASTA files in the specified directory.

    Args:
        datadir (str): Directory to search for FASTA files.
        fasta_suffix (str): Suffix for FASTA files (default: 'fasta').

    Returns:
        dict: Contains information about FASTA files:
            - 'fasta_files': List of FASTA file paths.
            - 'is_pooled': Boolean indicating if the FASTA files are pooled or unpooled.
    """
    # Compile regex to match FASTA files
    fasta_regex = re.compile(f"(\S+).{re.escape(fasta_suffix)}$")

    # List all FASTA files in the directory
    fasta_files = [os.path.join(datadir, f)
                   for f in os.listdir(datadir) if fasta_regex.match(f)]

    if not fasta_files:
        raise ValueError("No FASTA files detected. Ensure files have the correct suffix.")

    # Check if FASTA files are pooled or unpooled
    fq1_stubs = [re.match(r"(\S+)\.fastq\.1\.gz", os.path.basename(f)).group(1) for f in FASTQ1s]
    fasta_stubs = [fasta_regex.match(os.path.basename(f)).group(1) for f in fasta_files]

    if len(fasta_files) == 1:
        is_pooled = True  # Single pooled FASTA file
    elif sorted(fq1_stubs) == sorted(fasta_stubs):
        is_pooled = False  # Unpooled FASTA files matching FASTQ samples
    else:
        raise ValueError("FASTA files do not correspond to FASTQ samples. Check input consistency.")

    return {
        'fasta_files': fasta_files,
        'is_pooled': is_pooled
    }

# Validate FASTA files
fasta_info = check_fasta_files(DATADIR)
print("FASTA files:", fasta_info['fasta_files'])
print("Is FASTA pooled?", fasta_info['is_pooled'])

# utility function
def connect():
    '''Connect to sqlite db for storing pipeline datasets'''

    dbh = s3.connect(PARAMS['database']['name'])

###############################################################################
# BAM Files Processing
###############################################################################
@follows(mkdir('01_processed_bam_files.dir'))
@transform(FASTQ1s,
           regex('.+/(.+).fastq.1.gz'),
           r"01_processed_bam_files.dir/\1.bam")
def processBAMfiles(fastq1, outfile):
    """
    Align paired-end FASTQ files to a reference FASTA and generate BAM files.
    """
    # Determine the corresponding FASTQ2 file
    fastq2 = fastq1.replace('.1.gz', '.2.gz')

    # Determine the appropriate FASTA file (pooled or unpooled)
    if fasta_info['is_pooled']:
        reference_fasta = fasta_info['fasta_files'][0]  # Single pooled FASTA
    else:
        # Match the FASTA file to the sample based on naming convention
        sample_name = re.match(r"(.+)\.fastq\.1\.gz", os.path.basename(fastq1)).group(1)
        reference_fasta = [f for f in fasta_info['fasta_files'] if sample_name in f]
        if not reference_fasta:
            raise ValueError(f"No corresponding FASTA file found for sample: {sample_name}")
        reference_fasta = reference_fasta[0]

    # Use the `samtools_align` class from `mb` module
    align_task = mb.samtools_align(
        fastq1=fastq1,
        fastq2=fastq2,
        outfile=outfile,
        reference_fasta=reference_fasta,
        **PARAMS
    )

    # Build the alignment statement
    statement = align_task.buildStatement()

    # Execute the alignment statement
    P.run(statement,
          job_threads=PARAMS['samtools_job_threads'],
          job_memory=PARAMS['samtools_job_memory'],
          job_options=PARAMS.get('samtools_job_options', ''))


##############################################################################

def main(argv=None):
    if argv is None:
        argv=sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))

