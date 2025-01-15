# module for utility functions useful in handling shotgun data

import os
import re
import glob
import itertools
import errno
import cgatcore.iotools as IOTools
from cgatcore import pipeline as P

# Check that the input files correspond
def check_input(datadir='.', paired=True):

    if paired:
        fq1_regex = re.compile('(\S+).(fastq.1.gz)')
    else:
        fq1_regex = re.compile('(\S+).(fastq.gz)')

    mask1 = list(map(lambda x: bool(fq1_regex.match(x)),
                     os.listdir(datadir)))
    fastq1s = [os.path.join(datadir, i) \
               for i in itertools.compress(os.listdir(datadir),
                                           mask1)]

    if paired:
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
    else:
        if not sum(mask1):
            raise ValueError("No input files detected in run directory."
                              "Check the file suffix is fastq.gz")
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

    Some elements pulled form CGATMetaSequencing by Matt Jackson.

    Some options are  assumed to be passed via kwargs, as this and 
    inherited classes are written to work with a PARAMS dict 
    generated from a pipeline.yml config file.
    """

    def __init__(self, fastq1, outfile, **PARAMS):
        self.outdir = os.path.dirname(outfile)
        self.indir = os.path.dirname(fastq1)
        self.fastq1 = fastq1
        self.fastq2 = None
        self.fastq3 = None
        self.fq1_suffix = None
        self.fq2_suffix = None
        self.fq3_suffix = None
        self.prefixstrip = None
        self.prefix = None
        self.outfile = outfile
        self.PARAMS = PARAMS
        self.head = []

        '''check file can be opened on init & capture header 
           (first 5 lines used for interleave and format checks)
        '''
        try:
            self.openfile = IOTools.open_file(self.fastq1)
        except FileNotFoundError as e:
            msg = "cannot open file {}".format(self.fastq1)
            raise Exception(msg) from e
        self.head = [self.openfile.readline().rstrip("\n") for x in range(5)]

        '''autocheck file format and pairedness, 
           read count must be specified seperately
        '''
        self.getSuffix()
        self.isPaired()
        if self.prefixstrip is None:
            self.prefix = os.path.basename(self.fastq1.rstrip(self.fq1_suffix))
        else:
            self.prefix = os.path.basename(self.fastq1.rstrip(self.prefixstrip))

    def isPaired(self):
        if self.fastq1.endswith(".1.gz"):
            paired_name = self.fastq1.replace(".1",".2")
            assert len(glob.glob(paired_name)) > 0, (
                f"cannot find read 2 file at location {paired_name}"
                f" associated with read 1 file {self.fastq1}")
            self.fastq2 = paired_name

    '''get fastq1 file suffix 
    '''
    # allow custom extension to be passed; default extension is ".fastsq.1.gz"
    # inf_suffix used to set fq1_suffix, fq2_suffix, fq3_suffix, and prefix
    def getSuffix(self):
        # set suffix
        # if self.fastq1.endswith(".fastq.1.gz"):
        if self.fastq2 is not None:
            self.fq1_suffix = ".fastq.1.gz"
            self.fq2_suffix = '.fastq.2.gz'
            self.fq3_suffix = '.fastq.3.gz'
        else:
            self.fq1_suffix = ".fastq.gz"
