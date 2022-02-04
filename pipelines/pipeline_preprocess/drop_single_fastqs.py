'''
drop_fastqs.py
==============

:Author: Jethro Johnson

Purpose
-------

Provided with a fastq file fq1.

Provided with a text file (containing one fastq identifier per line), 
which specifiesd singletons to remove.

Will iterate over the fastq and drop any entries in the file. 
'''


import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Fastq as Fastq

import pandas as pd

import os,re,sys


def main(argv=None):

    if argv is None:
        argv = sys.argv


    parser = E.OptionParser(
        version="%prog version: $Id$",
        usage=globals()["__doc__"])

    parser.add_option("--fastq1", dest="fastq1")
    parser.add_option("--to-drop-single", dest='to_remove_singletons')
    parser.add_option("--fastq-out1", dest="fq_out1")
    parser.add_option("--fastq-drop1", dest="fq_dropped1")

    (options, args) = E.Start(parser)

    reads_to_remove = IOTools.openFile(options.to_remove_singletons).readlines()
    reads_to_remove = set([x.strip() for x in reads_to_remove])

    fastq_out = IOTools.openFile(options.fq_out1, 'w')
    fastq_host = IOTools.openFile(options.fq_dropped1, 'w')

    reads = 0
    dropped_reads = 0
    for read in Fastq.iterate(IOTools.openFile(fastq1)):
        reads += 1
        if read.identifier.split()[0] in reads_to_remove:
            fastq_host.write("@%s\n%s\n+\n%s\n" %
                             (read.identifier,
                              read.seq,
                              read.quals))
            dropped_reads += 1
        else:
            fastq_out.write("@%s\n%s\n+\n%s\n" %
                            (read.identifier,
                             read.seq,
                             read.quals))

    fastq_out.close()
    fastq_host.close()

    try:
        percent_dropped = dropped_reads/float(reads)*100
    except ZeroDivisionError:
        percent_dropped = 0.0
        
    E.info('Dropped %i of %i reads (%f percent)' \
           % (dropped_reads, reads, percent_dropped))

if __name__=="__main__":
    sys.exit(main(sys.argv))
            
