# module for utility functions useful in handling shotgun data

import os
import re
import itertools
import errno
from cgatcore import pipeline as P

# Check that the input files correspond
def check_input(datadir='.'):

    fq1_regex = re.compile('(\S+).(fastq.1.gz)')
    mask1 = list(map(lambda x: bool(fq1_regex.match(x)),
                     os.listdir(datadir)))
    fastq1s = [os.path.join(datadir, i) \
               for i in itertools.compress(os.listdir(datadir),
                                           mask1)]

    if sum(mask1):
        fq2_regex = re.compile('(\S+).(fastq.2.gz)')
        mask2 = list(map(lambda x: bool(fq2_regex.match(x)),
                         os.listdir(datadir)))
        fastq2s = [os.path.join(datadir, i) \
                   for i in itertools.compress(os.listdir(datadir), mask2)]
        if sum(mask2):
            assert sum(mask1) == sum(mask2), 'Not all input files have pairs'
            fq1_stubs = [fq1_regex.match(x).group(1) for x in fastq1s]
            fq2_stubs = [fq2_regex.match(x).group(1) for x in fastq2s]
            assert sorted(fq1_stubs) == sorted(fq2_stubs), \
                "First and second read pair files do not correspond"        
    else:
        raise ValueError("No input files detected in run directory."
                         " Check the file suffixes follow the notation"
                         " fastq.1.gz and fastq.2.gz.")

    return fastq1s

def symlnk(inf, outf):
    try:
        os.symlink(os.path.abspath(inf), outf)
    except OSError as e:
        if e.errno == errno.EEXIST:
            os.remove(outf)
            os.symlink(inf, outf)

class matchReference(object):
    """
    Base class for generating run statements to match mWGS reads to 
    reference sequences. Intended to work with single, paired, or
    paired + singleton fastq files. 

    Some options are  assumed to be passed via kwargs, as this and 
    inherited classes are written to work with a PARAMS dict 
    generated from a pipeline.yml config file.

    ** Options:
    fn_suffix - option to pass something other than .fastq.1.gz
    """

    def __init__(self, fastq1, outfile, **PARAMS):
        self.outdir = os.path.dirname(outfile)
        self.fastq1 = fastq1
        self.fastq2 = None
        self.fastq3 = None
        self.fq1_suffix = '.fastq.1.gz'
        self.fq2_suffix = '.fastq.2.gz'
        self.fq3_suffix = '.fastq.3.gz'
        self.outfile = outfile
        self.PARAMS = PARAMS

        # Assume that files are fastq and end .fastq.1.gz
        # Find mate pair file
        fastq2 = P.snip(self.fastq1, self.fq1_suffix) + self.fq2_suffix

        if os.path.exists(fastq2):
            self.fastq2 = fastq2        

        # Find singleton file
        fastq3 = P.snip(self.fastq1, self.fq1_suffix) + self.fq3_suffix

        if os.path.exists(fastq3):
            assert self.fastq2, "Can't have singletons without mate pairs"
            self.fastq3 = fastq3
