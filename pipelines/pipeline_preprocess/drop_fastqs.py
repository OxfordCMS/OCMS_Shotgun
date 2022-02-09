'''
drop_fastqs.py
==============

:Author: Jethro Johnson

Purpose
-------

Provided with three fastq files fq1, fq2, fq3, where 1 and 2 are paired, 
while 3 has singletons. 

Provided with two text files (containing one fastq identifier per line), 
which specify pairs and singletons to remove, respectively. 

Will iterate over the fastqs and drop any entries in the respective file. 
'''


import cgatcore.experiment as E
import cgatcore.iotools as IOTools
import cgat.Fastq as Fastq

import pandas as pd

import os,re,sys


def main(argv=None):

    if argv is None:
        argv = sys.argv


    parser = E.OptionParser(
        version="%prog version: $Id$",
        usage=globals()["__doc__"])

    parser.add_option("--fastq1", dest="fastq1")
    parser.add_option("--fastq2", dest="fastq2")
    parser.add_option("--fastq3", dest="fastq3")

    parser.add_option("--to-drop-paired", dest='to_remove_paired')
    parser.add_option("--to-drop-single", dest='to_remove_singletons')

    parser.add_option("--fastq-out1", dest="fq_out1")
    parser.add_option("--fastq-out2", dest="fq_out2")
    parser.add_option("--fastq-out3", dest="fq_out3")

    parser.add_option("--fastq-drop1", dest="fq_dropped1")
    parser.add_option("--fastq-drop2", dest="fq_dropped2")
    parser.add_option("--fastq-drop3", dest="fq_dropped3")

    (options, args) = E.Start(parser)

    # Fetch the reads to remove
    pairs_to_remove = IOTools.openFile(options.to_remove_paired).readlines()
    pairs_to_remove = set([x.strip() for x in pairs_to_remove])
    
    singles_to_remove = IOTools.openFile(options.to_remove_singletons).readlines()
    singles_to_remove = set([x.strip() for x in singles_to_remove])

    # open the outfiles
    fastq1_out = IOTools.openFile(options.fq_out1, 'w')
    fastq2_out = IOTools.openFile(options.fq_out2, 'w')
    fastq3_out = IOTools.openFile(options.fq_out3, 'w')

    fastq1_host = IOTools.openFile(options.fq_dropped1, 'w')
    fastq2_host = IOTools.openFile(options.fq_dropped2, 'w')
    fastq3_host = IOTools.openFile(options.fq_dropped3, 'w')

    dropped_pairs = 0
    pairs = 0
    # Drop the paired reads
    for read1, read2 in zip(Fastq.iterate(IOTools.openFile(options.fastq1)),
                            Fastq.iterate(IOTools.openFile(options.fastq2))):
        pairs +=1
        
        # bmtagger truncates fastq headers at space and won't accept
        # non-identical headers therefore, if one read matches, both
        # are chucked.
        r1_id = read1.identifier.split()[0]
        r2_id = read2.identifier.split()[0]

        if r1_id in pairs_to_remove:
            # Both are host
            fastq1_host.write("@%s\n%s\n+\n%s\n" %
                              (read1.identifier,
                               read1.seq,
                               read1.quals))
            fastq2_host.write("@%s\n%s\n+\n%s\n" %
                              (read2.identifier,
                               read2.seq,
                               read2.quals))
            dropped_pairs += 1
        else:
            # Neither are host
            fastq1_out.write("@%s\n%s\n+\n%s\n" %
                             (read1.identifier,
                              read1.seq,
                              read1.quals))
            fastq2_out.write("@%s\n%s\n+\n%s\n" %
                             (read2.identifier,
                              read2.seq,
                              read2.quals))            
    # Drop singletons
    singletons = 0
    dropped_singletons = 0
    for read in Fastq.iterate(IOTools.openFile(options.fastq3)):
        singletons += 1
        if read.identifier.split()[0] in singles_to_remove:
            fastq3_host.write("@%s\n%s\n+\n%s\n" %
                              (read.identifier,
                               read.seq,
                               read.quals))
            dropped_singletons +=1
        else:
            fastq3_out.write("@%s\n%s\n+\n%s\n" %
                             (read.identifier,
                              read.seq,
                              read.quals))

    fastq1_out.close()
    fastq2_out.close()
    fastq3_out.close()
    fastq1_host.close()
    fastq2_host.close()
    fastq3_host.close()

    try:
        percent_pairs = dropped_pairs/float(pairs)*100
    except ZeroDivisionError:
        percent_pairs = 0.0
    try:
        percent_singletons = dropped_singletons/float(singletons)*100
    except ZeroDivisionError:
        percent_singletons = 0.0
        
    E.info('Dropped %i of %i read pairs (%f percent)' \
           % (dropped_pairs, pairs, percent_pairs))
    E.info('Dropped %i of %i singletons (%f percent)' \
           % (dropped_singletons, singletons, percent_singletons))

if __name__=="__main__":
    sys.exit(main(sys.argv))
