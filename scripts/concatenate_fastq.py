'''
concatenate_fastq.py
====================

:Author: Sandi Yen
:Tags: Python

Purpose
-------
Concatenates fastq files of paired-end reads into one fastq file

Usage
-----
.. Example use case

Example::
    python cgat_script_template.py

Type::
    python cgat_script_template.py --help

for command line help.

Command line options
--------------------
'''

import sys
import os
import re
import glob
import cgatcore.experiment as E

def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    Expects paired end fastq files to be in format fastq.1.gz, fastq.2.gz
    """

    if argv is None:
        argv = sys.argv

    # set up command line parser
    parser = E.ArgumentParser(description = __doc__)

    
    parser.add_argument("-i", "--indir", dest="indir", type=str,
                        default=os.getcwd(),
                        help="supply input directory containing fastq files")
    parser.add_argument("-o", "--outdir", dest="outdir", type=str,
                        default="concat_fastq.dir",
                        help="supply output directory")

    # add common options (-h/--help, ,,,) and parse command line
    (args) = E.start(parser, argv=argv)

    # create output directory if doesn't exist
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    # get list of fastq files in indir
    file_list = glob.glob(os.path.join(args.indir, "*fastq*gz"))
    
    # extract sample names
    sample_name = [os.path.basename(x) for x in file_list]
    sample_name = [re.sub(r'\.fastq.*gz', '', x) for x in sample_name]
    sample_name = list(set(sample_name))

    # work with one sample at a time
    for sample in sample_name:
        curr_file = [x for x in file_list if sample in x]
        
        outfile = sample + '.fastq'
        outfile = os.path.join(args.outdir, outfile)
        infile = " ".join(curr_file)
        # concatenate files together
        statement = "zcat %s >> %s && gzip %s" % (infile, outfile, outfile)
        os.system(statement)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
