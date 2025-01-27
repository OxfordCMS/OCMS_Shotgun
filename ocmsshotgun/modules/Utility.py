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

class metaFastn(object):
    """
    Base class for reading in fastn files. At the moment, only works for fastqs. 
    Intended to work with single, paired, or paired + singleton fastq files. 

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
            msg = f"cannot open file {self.fastq1}"
            raise Exception(msg) from e
        self.head = [self.openfile.readline().rstrip("\n") for x in range(5)]

        '''autocheck file format and pairedness, 
           read count must be specified seperately
        '''
        self.isPaired()
        self.hasSingleton()
        self.getFormat()
        self.getSuffix()
        if self.prefixstrip is None:
            self.prefix = os.path.basename(self.fastq1.rstrip(self.fq1_suffix))
        else:
            self.prefix = os.path.basename(self.fastq1.rstrip(self.prefixstrip))

    '''check if paired and if containts interleaved pairs or matching files'''
    def isPaired(self):
        if self.fastq1.endswith(".1.gz"):
            paired_name = self.fastq1.replace(".1",".2")
            assert len(glob.glob(paired_name)) > 0, (
                f"cannot find read 2 file at location {paired_name}"
                f" associated with read 1 file {self.fastq1}")
            self.fastq2 = paired_name

    '''check for singletons'''
    def hasSingleton(self):
        fq3_name = self.fastq1.replace(".1",".3")
        if os.path.exists(fq3_name):
            # check file is not empty
            if os.stat(fq3_name).st_size != 0:
                self.fastq3 = fq3_name

    '''check it is fasta or fastq and if compressed'''    
    def getFormat(self):
        extensions=("fasta","fastq")
        for i in extensions:    
            if self.fastq1.endswith((i+".1.gz",i+".gz")):
                if i == "fastq":
                    self.fileformat="fastq"
                else:
                    self.fileformat="fasta"
           
        msg = f"file {self.fastq1} is not of the correct format (fasta or fastq)."
        assert self.fileformat, msg
        if self.fileformat == "fasta":
            assert self.head[0][0] == ">", (
                "invalid header on first line for fasta format")
        else:
            assert self.head[0][0] == "@", (
                "invalid header on first line for fastq format")
            
    '''get fastq1 file suffix '''
    # allow custom extension to be passed; default extension is ".fastsq.1.gz"
    # inf_suffix used to set fq1_suffix, fq2_suffix, fq3_suffix, and prefix
    def getSuffix(self):
        # set suffix
        # if self.fastq1.endswith(".fastq.1.gz"):
        if self.fastq2 is not None:
            assert self.fastq1.endswith(f".{self.fileformat}.1.gz"), (
                f"Paired-end {self.fileformat} files must be in notation"
                f" '.{self.fileformat}.1.gz'")
            self.fq1_suffix = f".{self.fileformat}.1.gz"
            self.fq2_suffix = f'.{self.fileformat}.2.gz'
            self.fq3_suffix = f'.{self.fileformat}.3.gz'
        else:
            assert self.fastq1.endswith(f".{self.fileformat}.gz"), (
                "Single-end fastq files must be in notation "
                f"'{self.fileformat}.gz'"
            )
            self.fq1_suffix = f".{self.fileformat}.gz"
